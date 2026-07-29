#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_predictions.py — 用公开实验数据库交叉核验 admet-ai 的预测值。

背景：
  masld_screener 给出的候选是纯虚拟筛选结果（理化计算 + 相似度代理 + admet-ai 预测）。
  源库 SDF 本身没有实验值，因此本脚本从外部公开数据库拉取真实测定值，
  与 admet-ai 对候选的预测逐项比对，给出"预测 vs 实验"一致率，作为置信度依据。

比对端点：
  - hERG 心脏毒性（自动，快速）:
        预测 = admet-ai hERG 概率 >= 0.5 判为"有毒"；
        实验 = PubChem BioAssay 中 gene=KCNH2 且 Assay 名含 "hERG" 的 Active/Inactive。
  - AMES 致突变 / DILI 肝毒（--deep 模式，较慢）:
        需扫描完整 assaysummary（部分常见药摘要很大，可能超时）；
        预测同理；实验 = PubChem BioAssay 中 Ames/Salmonella 或 hepatotox/liver 的结局。
        默认不开启，以避免对常见药的大摘要造成长时间超时。

溶解度 / 口服生物利用度：
  PubChem BioAssay 一般不含这类连续实验值；脚本保留 admet-ai 预测列供参考，
  建议到 admetSAR 3.0 / ADMETNet / PKKB 做人工或批量核对（见 VALIDATION_PLAN.md）。

数据源：
  - admet-ai  (本地模型，与 masld_screener Phase 4 完全一致，复算逐端点概率)
  - PubChem PUG-REST  (assaysummary, 实验 BioAssay 结果，基因过滤)

用法：
  python validate_predictions.py -i ranked_candidates.csv -o validation_report.csv
  python validate_predictions.py -i ranked_candidates.csv -o validation_report.csv --deep --limit 20
"""

import argparse
import ast
import csv
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter

# ----------------------------------------------------------------------------
# 复用与主流程一致的端点映射逻辑（hepato -> 已修正为 dili，避免误匹配 Clearance_Hepatocyte_AZ）
# ----------------------------------------------------------------------------
_TOX_KEYWORDS = {
    "herg_inhibition": ["herg"],
    "ames_mutagenicity": ["ames", "mutagen"],
    "hepatotoxicity": ["dili", "hepatotox"],
    "oral_bioavailability": ["oral", "bioavailability"],
    "aqueous_solubility": ["solubility", "logs", "aqueous"],
}


def _find_column(columns, key):
    norm = [c.lower().replace("_", "").replace(" ", "").replace("-", "") for c in columns]
    kw = _TOX_KEYWORDS.get(key, [key])
    for i, c in enumerate(norm):
        if any(k.replace("_", "") in c for k in kw):
            return columns[i]
    return None


# ----------------------------------------------------------------------------
# 网络请求
# ----------------------------------------------------------------------------
_SSL_CTX = ssl.create_default_context()
_SSL_CTX.check_hostname = False
_SSL_CTX.verify_mode = ssl.CERT_NONE


def _get_json(url, timeout=12, retries=2):
    last = None
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MASLD-validator/1.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=_SSL_CTX) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
            time.sleep(0.4)
        except Exception as e:
            last = e
            time.sleep(0.4)
    return {"__error__": str(last)}


# ----------------------------------------------------------------------------
# PubChem 解析（兼容 classic Columns/Row/Cell 格式）
# ----------------------------------------------------------------------------
_OUTCOMES = {"active", "inactive", "unspecified", "probe", "inconclusive", "both"}


def _parse_classic(tbl):
    """解析 PubChem assaysummary 的 classic 格式，返回 [(outcome, name), ...]。"""
    if not isinstance(tbl, dict):
        return []
    cols = tbl.get("Columns", {}).get("Column", [])
    if isinstance(cols, dict):
        cols = cols.get("Column", [])
    cidx = {c: i for i, c in enumerate(cols)}
    oi = cidx.get("Activity Outcome")
    ai = cidx.get("Assay Name")
    ni = cidx.get("Activity Name")
    rows = []
    for row in tbl.get("Row", []):
        cell = row.get("Cell")
        if isinstance(cell, str):
            try:
                cell = ast.literal_eval(cell)
            except Exception:
                cell = [cell]
        if not isinstance(cell, (list, tuple)):
            continue
        outcome = str(cell[oi]) if oi is not None and oi < len(cell) else ""
        name = ""
        if ai is not None and ai < len(cell):
            name += str(cell[ai]) + " "
        if ni is not None and ni < len(cell):
            name += str(cell[ni])
        rows.append((outcome.strip(), name.strip()))
    return rows


def _pick_outcome(rows, keyword):
    """在 rows 中找名称含 keyword 的实验结局，优先 Active/Inactive，其次其它。"""
    hits = [(o, n) for o, n in rows if keyword in n.lower()]
    if not hits:
        return None
    for o, n in hits:
        if o.lower() in ("active", "inactive"):
            return o.lower()
    return hits[0][0].lower()


def pubchem_cid_from_smiles(smiles):
    enc = urllib.parse.quote(smiles)
    j = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{enc}/cids/JSON")
    if isinstance(j, dict) and j.get("IdentifierList"):
        return str(j["IdentifierList"]["CID"][0])
    j = _get_json(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{enc}/cids/JSON?identity_type=inchikey"
    )
    if isinstance(j, dict) and j.get("IdentifierList"):
        return str(j["IdentifierList"]["CID"][0])
    return None


def pubchem_experiment(cid, deep=False):
    """返回 {'herg':..., 'ames':..., 'dili':...} 实验结局（active/inactive/None）。"""
    res = {"herg": None, "ames": None, "dili": None}
    if not cid:
        return res
    # hERG：基因过滤（快速，~2s），再按 Assay 名含 hERG 取结局
    j = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON?gene=KCNH2")
    if isinstance(j, dict) and j.get("Table"):
        rows = _parse_classic(j["Table"])
        res["herg"] = _pick_outcome(rows, "herg")
    # AMES / DILI：仅在 --deep 下扫描完整摘要（可能超时）
    if deep:
        j = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON", timeout=20)
        if isinstance(j, dict) and j.get("Table"):
            rows = _parse_classic(j["Table"])
            res["ames"] = _pick_outcome(rows, "ames") or _pick_outcome(rows, "salmonella") or _pick_outcome(rows, "mutagen")
            res["dili"] = _pick_outcome(rows, "dili") or _pick_outcome(rows, "hepat") or _pick_outcome(rows, "liver")
    return res


# ----------------------------------------------------------------------------
# admet-ai 预测（复算）
# ----------------------------------------------------------------------------
def run_admet(smiles_list):
    from admet_ai import ADMETModel
    model = ADMETModel()
    return model.predict(smiles_list)


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="用 PubChem 实验数据交叉验证 admet-ai 预测")
    ap.add_argument("-i", "--input", required=True, help="ranked_candidates.csv 路径")
    ap.add_argument("-o", "--output", default="validation_report.csv", help="输出报告 CSV 路径")
    ap.add_argument("--limit", type=int, default=0, help="仅处理前 N 个候选（试点用）")
    ap.add_argument("--deep", action="store_true", help="额外扫描完整 assaysummary 以核对 AMES/DILI（较慢，可能超时）")
    args = ap.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        cands = [r for r in csv.DictReader(f)]
    if args.limit:
        cands = cands[: args.limit]
    print(f"候选数: {len(cands)}  (deep={args.deep})")

    smiles_list = [r["smiles"] for r in cands]

    print("正在运行 admet-ai 预测 ...")
    t0 = time.time()
    df = run_admet(smiles_list)
    print(f"  admet-ai 列名: {list(df.columns)}")

    col_map = {
        "herg_inhibition": _find_column(df.columns, "herg_inhibition"),
        "ames_mutagenicity": _find_column(df.columns, "ames_mutagenicity"),
        "hepatotoxicity": _find_column(df.columns, "hepatotoxicity"),
        "oral_bioavailability": _find_column(df.columns, "oral_bioavailability"),
        "aqueous_solubility": _find_column(df.columns, "aqueous_solubility"),
    }
    print(f"  端点列映射: {col_map}")
    row_map = {}
    for k, smi in enumerate(smiles_list):
        try:
            row_map[smi] = df.loc[smi]
        except KeyError:
            try:
                row_map[smi] = df.iloc[k]
            except Exception:
                row_map[smi] = None

    pred = {}
    for smi in smiles_list:
        row = row_map.get(smi)
        d = {}
        if row is not None:
            for key, col in col_map.items():
                if col and col in row.index:
                    try:
                        d[key] = float(row[col])
                    except Exception:
                        d[key] = None
        pred[smi] = d

    print("正在查询 PubChem 实验 BioAssay (并行, 4 线程) ...", flush=True)
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _work(smi):
        cid = pubchem_cid_from_smiles(smi)
        expd = pubchem_experiment(cid, deep=args.deep)
        return smi, cid, expd

    exp, cid_of = {}, {}
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(_work, smi): smi for smi in smiles_list}
        done = 0
        for fut in as_completed(futs):
            smi, cid, expd = fut.result()
            cid_of[smi] = cid
            exp[smi] = expd
            done += 1
            if done % 10 == 0 or done == len(smiles_list):
                print(f"  已完成 {done}/{len(smiles_list)}", flush=True)

    def tox_of(prob):
        return None if prob is None else ("toxic" if prob >= 0.5 else "safe")

    def exp_of(outcome):
        if outcome == "active":
            return "toxic"
        if outcome == "inactive":
            return "safe"
        return None

    out_rows = []
    stats = {k: {"n_exp": 0, "match": 0, "tp": 0, "fp": 0, "tn": 0, "fn": 0} for k in ("herg", "ames", "dili")}
    for r, smi in zip(cands, smiles_list):
        p, e, cid = pred[smi], exp[smi], cid_of[smi]
        row_out = {
            "name": r.get("name", ""),
            "cid": cid or "",
            "pred_hERG": ("%.3f" % p.get("herg_inhibition")) if p.get("herg_inhibition") is not None else "",
            "exp_hERG": e["herg"] or "",
            "pred_AMES": ("%.3f" % p.get("ames_mutagenicity")) if p.get("ames_mutagenicity") is not None else "",
            "exp_AMES": (e["ames"] or "") if args.deep else "(deep)",
            "pred_DILI": ("%.3f" % p.get("hepatotoxicity")) if p.get("hepatotoxicity") is not None else "",
            "exp_DILI": (e["dili"] or "") if args.deep else "(deep)",
            "pred_bioavailability": ("%.3f" % p.get("oral_bioavailability")) if p.get("oral_bioavailability") is not None else "",
            "pred_solubility": ("%.3f" % p.get("aqueous_solubility")) if p.get("aqueous_solubility") is not None else "",
        }
        for key, pk, ek in (("herg", "herg_inhibition", "herg"),
                            ("ames", "ames_mutagenicity", "ames"),
                            ("dili", "hepatotoxicity", "dili")):
            if not args.deep and key != "herg":
                continue
            pv = tox_of(p.get(pk))
            ev = exp_of(e[ek]) if e[ek] else None
            if ev is not None:
                stats[key]["n_exp"] += 1
                if pv is not None and pv == ev:
                    stats[key]["match"] += 1
                if pv == "toxic" and ev == "toxic":
                    stats[key]["tp"] += 1
                elif pv == "toxic" and ev == "safe":
                    stats[key]["fp"] += 1
                elif pv == "safe" and ev == "safe":
                    stats[key]["tn"] += 1
                elif pv == "safe" and ev == "toxic":
                    stats[key]["fn"] += 1
        out_rows.append(row_out)

    fields = ["name", "cid", "pred_hERG", "exp_hERG", "pred_AMES", "exp_AMES",
              "pred_DILI", "exp_DILI", "pred_bioavailability", "pred_solubility"]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in out_rows:
            w.writerow(row)
    print(f"  [输出] 已写入验证报告: {args.output} ({len(out_rows)} 行)")

    print(f"\n{'=' * 64}")
    print("  预测 vs 实验 一致性汇总")
    print(f"{'-' * 64}")
    for key in ("herg", "ames", "dili"):
        if not args.deep and key != "herg":
            print(f"  {key.upper():6}: 未运行（需 --deep 或人工核对 admetSAR3/ADMETNet）")
            continue
        s = stats[key]
        if s["n_exp"] == 0:
            print(f"  {key.upper():6}: 无实验数据可比对（PubChem 覆盖不足或超时）")
            continue
        acc = s["match"] / s["n_exp"]
        sens = s["tp"] / (s["tp"] + s["fn"]) if (s["tp"] + s["fn"]) else float("nan")
        spec = s["tn"] / (s["tn"] + s["fp"]) if (s["tn"] + s["fp"]) else float("nan")
        print(f"  {key.upper():6}: 实验覆盖 {s['n_exp']:2d} 个 | 一致率 {acc*100:5.1f}% | "
              f"灵敏度 {sens*100:5.1f}% | 特异度 {spec*100:5.1f}%")
        print(f"          TP={s['tp']} FP={s['fp']} TN={s['tn']} FN={s['fn']}")
    print(f"{'=' * 64}")
    print(f"  admet-ai 预测 + PubChem 查询耗时: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
