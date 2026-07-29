#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
masld_screener.py — 从 SDF 文件筛选具有降脂活性且低毒性的 MASLD 小分子。

五阶段筛选管线：
  Phase 1  SDF 读取 / 去盐 / 标准化 / 基于 canonical SMILES 去重
  Phase 2  理化性质计算与初筛（基于 config.PHYS_CHEM_FILTERS）
  Phase 3  与 9 个已知活性参考分子的 Morgan 指纹 Tanimoto 相似性搜索
  Phase 4  ADMET 毒性预测（admet-ai，缺失时跳过并以默认分填充）
  Phase 5  综合评分与排序，输出 ranked_candidates.csv / 可选 .sdf

Python 版本要求：>= 3.8
依赖（见 requirements.txt）：rdkit, admet-ai, pandas, numpy, tqdm

用法示例：
  python masld_screener.py input.sdf
  python masld_screener.py input.sdf -o result.csv --output_sdf
  python masld_screener.py input.sdf --sim_threshold 0.35
"""

import argparse
import csv
import os
import sys
import time
import warnings

from tqdm import tqdm

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, SaltRemover
try:
    from rdkit.Chem.rdMolStandardize import Normalize
except ImportError:  # 某些 RDKit 构建中将模块置于 rdkit.Chem.MolStandardize
    from rdkit.Chem.MolStandardize import rdMolStandardize as _rdMolStandardize
    Normalize = _rdMolStandardize.Normalize

# 屏蔽 RDKit 的冗余日志（如芳香性/价键警告），保持输出整洁
RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

import config
from reference_molecules import REFERENCE_MOLECULES

# ---------------------------------------------------------------------------
# 通用工具
# ---------------------------------------------------------------------------

def log_phase(name, n_in, n_out, elapsed):
    """打印单个阶段的过程日志。"""
    rate = (1.0 - n_out / n_in) * 100.0 if n_in else 0.0
    print(f"\n{'=' * 64}")
    print(f"  {name}")
    print(f"{'-' * 64}")
    print(f"  输入分子数 : {n_in}")
    print(f"  输出分子数 : {n_out}")
    print(f"  过滤率     : {rate:.1f}%")
    print(f"  耗时       : {elapsed:.2f}s")
    print(f"{'=' * 64}")


def get_name(mol, idx):
    """从分子对象中获取可读名称；优先 _Name，其次常见 ID 字段。"""
    for prop in ("_Name", "ID", "Name", "name", "CODE", "catalog"):
        try:
            if mol.HasProp(prop):
                val = mol.GetProp(prop).strip()
                if val:
                    return val
        except Exception:
            continue
    return f"Mol_{idx}"


def lipinski_violations(MW, LogP, HBD, HBA):
    """计算 Lipinski 五规则违规数（MW>500 / LogP>5 / HBD>5 / HBA>10）。"""
    v = 0
    if MW > 500:
        v += 1
    if LogP > 5:
        v += 1
    if HBD > 5:
        v += 1
    if HBA > 10:
        v += 1
    return v


# ---------------------------------------------------------------------------
# Phase 1: SDF 读取、去盐、标准化、去重
# ---------------------------------------------------------------------------

def _open_sdf_supplier(sdf_path):
    """打开 SDF 读取器。

    RDKit 的 C++ 实现在某些平台（如 Windows）对含非 ASCII 字符
    （中文、空格等）的文件名支持不佳，会抛出 OSError: Bad input file。
    此处先尝试直接打开；失败则将文件复制到临时 ASCII 名副本后重试，
    保证可移植性与鲁棒性。
    """
    try:
        sup = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=True)
        # 触发一次解析以确认文件确实可打开（惰性读取，需访问）
        _ = sup.GetItemText(0) if sup.GetNumEntries() > 0 else None
        return sup, None
    except Exception:
        pass
    # 回退：复制到临时 ASCII 文件名
    import tempfile, shutil
    tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False, prefix="masld_in_")
    tmp.close()
    try:
        shutil.copyfile(sdf_path, tmp.name)
        sup = Chem.SDMolSupplier(tmp.name, removeHs=False, sanitize=True)
        return sup, tmp.name
    except Exception as e:
        try:
            os.remove(tmp.name)
        except Exception:
            pass
        raise OSError(f"无法读取 SDF 文件: {sdf_path} ({e})")


def phase1_read_dedup(sdf_path):
    t0 = time.time()
    supplier, tmp_path = _open_sdf_supplier(sdf_path)
    remover = SaltRemover.SaltRemover()

    seen = set()
    records = []
    n_in = 0
    n_dup = 0
    n_bad = 0

    for i, mol in enumerate(tqdm(supplier, desc="Phase 1 读取SDF", unit="mol")):
        if mol is None:
            n_bad += 1
            continue
        n_in += 1

        name = get_name(mol, i)

        # 去盐：剥离常见盐/溶剂，但避免把整个分子删空
        try:
            mol = remover.StripMol(mol, dontRemoveEverything=True)
        except Exception:
            pass

        # 标准化：电荷/互变异构/价键等归一化
        try:
            mol = Normalize(mol)
        except Exception:
            pass

        # 基于 canonical SMILES 去重
        try:
            can = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            n_bad += 1
            continue
        if can is None:
            n_bad += 1
            continue
        if can in seen:
            n_dup += 1
            continue
        seen.add(can)

        records.append({"mol": mol, "canonical": can, "name": name, "idx": i})

    n_out = len(records)
    log_phase("Phase 1: SDF 读取 / 去盐 / 标准化 / 去重", n_in, n_out, time.time() - t0)
    print(f"  跳过无法解析的分子 : {n_bad}")
    print(f"  去除重复分子数     : {n_dup}")
    if tmp_path is not None:
        try:
            os.remove(tmp_path)
        except Exception:
            pass
    return records


# ---------------------------------------------------------------------------
# Phase 2: 理化性质计算与初筛
# ---------------------------------------------------------------------------

def phase2_physchem(records):
    t0 = time.time()
    n_in = len(records)

    filters = config.PHYS_CHEM_FILTERS
    total_items = len(filters)

    passed = []
    for rec in tqdm(records, desc="Phase 2 理化计算", unit="mol"):
        mol = rec["mol"]
        try:
            MW = Descriptors.MolWt(mol)
            LogP = Crippen.MolLogP(mol)
            TPSA = Descriptors.TPSA(mol)
            HBD = Lipinski.NumHDonors(mol)
            HBA = Lipinski.NumHAcceptors(mol)
            RotB = Lipinski.NumRotatableBonds(mol)
            viol = lipinski_violations(MW, LogP, HBD, HBA)
        except Exception:
            continue

        vals = {
            "MW": MW, "LogP": LogP, "TPSA": TPSA,
            "HBD": HBD, "HBA": HBA, "RotB": RotB,
            "Lipinski_violations": viol,
        }

        ok = True
        npass = 0
        for key, (lo, hi) in filters.items():
            v = vals[key]
            if lo <= v <= hi:
                npass += 1
            else:
                ok = False

        rec["descriptors"] = vals
        rec["physchem_pass"] = npass
        if ok:
            passed.append(rec)

    n_out = len(passed)
    log_phase("Phase 2: 理化性质计算与初筛", n_in, n_out, time.time() - t0)
    print(f"  理化筛选项数(每项一票) : {total_items}")
    return passed


# ---------------------------------------------------------------------------
# Phase 3: 已知活性分子相似性搜索
# ---------------------------------------------------------------------------

def phase3_similarity(records, sim_threshold=None, radius=None, n_bits=None):
    t0 = time.time()
    n_in = len(records)
    if sim_threshold is None:
        sim_threshold = config.SIMILARITY_THRESHOLD
    if radius is None:
        radius = config.FINGERPRINT_RADIUS
    if n_bits is None:
        n_bits = config.FINGERPRINT_BITS

    # 构建 Morgan 指纹生成器（现代 API，缺失时回退到旧接口）
    try:
        from rdkit.Chem import rdFingerprintGenerator
        _gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        def _fp(m):
            return _gen.GetFingerprint(m)
    except Exception:
        def _fp(m):
            return AllChem.GetMorganFingerprintAsBitVst(m, radius, nBits=n_bits)

    # 参考分子指纹
    ref_fps = []
    for r in REFERENCE_MOLECULES:
        rmol = Chem.MolFromSmiles(r["smiles"])
        if rmol is None:
            print(f"  [警告] 参考分子 {r['name']} SMILES 无法解析，已跳过")
            continue
        fp = _fp(rmol)
        ref_fps.append((r["name"], r["target"], fp))

    passed = []
    for rec in tqdm(records, desc="Phase 3 相似性搜索", unit="mol"):
        mol = rec["mol"]
        try:
            fp = _fp(mol)
        except Exception:
            continue

        best_tc = -1.0
        best_name = ""
        best_target = ""
        for nm, tg, rfp in ref_fps:
            tc = DataStructs.TanimotoSimilarity(fp, rfp)
            if tc > best_tc:
                best_tc = tc
                best_name = nm
                best_target = tg

        rec["max_tc"] = best_tc
        rec["best_match"] = best_name
        rec["best_target"] = best_target

        if best_tc >= sim_threshold:
            passed.append(rec)

    n_out = len(passed)
    log_phase(f"Phase 3: 相似性搜索 (阈值 Tanimoto >= {sim_threshold})", n_in, n_out, time.time() - t0)
    print(f"  参考分子数 : {len(ref_fps)}")
    return passed


# ---------------------------------------------------------------------------
# Phase 4: ADMET 毒性预测
# ---------------------------------------------------------------------------

# 毒性筛选关键字 -> 候选列名匹配词
_TOX_KEYWORDS = {
    "herg_inhibition": ["herg"],
    "ames_mutagenicity": ["ames", "mutagen"],
    "hepatotoxicity": ["hepato"],
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


def _check_cond(val, cond):
    """根据筛选条件判断单条端点是否通过。"""
    # 布尔型条件：config 中 False 表示期望“安全/阴性”
    if isinstance(cond, bool):
        if isinstance(val, bool):
            return (val == False)  # 期望阴性
        if isinstance(val, str):
            low = val.lower()
            if any(w in low for w in ("neg", "non", "inactive", "false", "no")):
                return True
            if any(w in low for w in ("pos", "tox", "active", "true")):
                return False
            return False
        try:
            return float(val) < 0.5
        except Exception:
            return False
    # 字符串阈值条件：形如 ">0.3" / ">-5"
    if isinstance(cond, str) and cond.startswith(">"):
        try:
            thr = float(cond[1:])
            return float(val) > thr
        except Exception:
            return False
    return False


def phase4_admet(records):
    t0 = time.time()
    n_in = len(records)

    try:
        from admet_ai import ADMETModel
    except Exception:
        print(f"\n{'=' * 64}")
        print("  Phase 4: ADMET 毒性预测（admet-ai 未安装 -> 跳过）")
        print(f"{'-' * 64}")
        print("  [警告] 未检测到 admet-ai，跳过毒性预测阶段。")
        print("          toxic_score 默认设为 0.5（中性分）。")
        print("          安装方法: pip install admet-ai")
        print(f"  输入分子数 : {n_in}   输出分子数 : {n_in}")
        print(f"  耗时 : {time.time() - t0:.2f}s")
        print(f"{'=' * 64}")
        for rec in records:
            rec["toxic_score"] = 0.5
            rec["toxic_mappable"] = 0
        return records

    print(f"\n{'=' * 64}")
    print("  Phase 4: ADMET 毒性预测 (admet-ai)")
    print(f"{'-' * 64}")
    try:
        model = ADMETModel()
        smiles_list = [r["canonical"] for r in records]
        df = model.predict(smiles_list)
    except Exception as e:
        print(f"  [警告] admet-ai 预测失败，跳过 Phase 4: {e}")
        print("          toxic_score 默认设为 0.5。")
        for rec in records:
            rec["toxic_score"] = 0.5
            rec["toxic_mappable"] = 0
        return records

    # 打印列名供调试（符合需求）
    print("  admet-ai 输出列名:", list(df.columns))

    # 建立 SMILES -> 行 的映射（兼容 index=smiles 与按位置两种返回形式）
    row_map = {}
    for k, smi in enumerate(smiles_list):
        try:
            row_map[smi] = df.loc[smi]
        except KeyError:
            try:
                row_map[smi] = df.iloc[k]
            except Exception:
                row_map[smi] = None

    filters = config.TOXICITY_FILTERS
    col_for = {key: _find_column(df.columns, key) for key in filters}
    mappable = [k for k, c in col_for.items() if c is not None]
    print(f"  可映射的毒性端点 : {mappable if mappable else '无'}")

    n_full_pass = 0
    for rec in tqdm(records, desc="Phase 4 毒性评估", unit="mol"):
        row = row_map.get(rec["canonical"])
        if row is None:
            rec["toxic_score"] = 0.5
            rec["toxic_mappable"] = 0
            continue  # 无法评估时保留，交由综合评分中性处理
        npass = 0
        nmap = 0
        for key, cond in filters.items():
            col = col_for[key]
            if col is None:
                continue
            nmap += 1
            val = row[col]
            if _check_cond(val, cond):
                npass += 1
        rec["toxic_mappable"] = nmap
        rec["toxic_score"] = (npass / nmap) if nmap else 0.5
        if nmap > 0 and npass == nmap:
            n_full_pass += 1

    # 软筛选：保留全部候选，toxic_score 参与 Phase 5 综合评分；
    # 同时统计“完全通过全部可映射毒性端点”的分子数（供用户参考）。
    n_out = len(records)
    log_phase("Phase 4: ADMET 毒性预测与评分（软筛选）", n_in, n_out, time.time() - t0)
    print(f"  完全通过全部毒性端点(可映射)的分子数 : {n_full_pass}")
    return records


# ---------------------------------------------------------------------------
# Phase 5: 综合评分与排序
# ---------------------------------------------------------------------------

def phase5_score_and_rank(records):
    t0 = time.time()
    n_in = len(records)
    w = config.SCORE_WEIGHTS
    total_phys = len(config.PHYS_CHEM_FILTERS)

    rows = []
    for rec in records:
        d = rec.get("descriptors", {})
        max_tc = rec.get("max_tc", 0.0)
        phys_pass = rec.get("physchem_pass", 0)
        toxic = rec.get("toxic_score", 0.5)

        sim_score = min(max_tc / 0.7, 1.0)
        phys_score = (phys_pass / total_phys) if total_phys else 0.0
        final = w["similarity"] * sim_score + w["physchem"] * phys_score + w["toxicity"] * toxic

        rows.append({
            "name": rec["name"],
            "smiles": rec["canonical"],
            "MW": round(d.get("MW", float("nan")), 2),
            "LogP": round(d.get("LogP", float("nan")), 2),
            "TPSA": round(d.get("TPSA", float("nan")), 2),
            "HBD": d.get("HBD", ""),
            "HBA": d.get("HBA", ""),
            "RotB": d.get("RotB", ""),
            "Lipinski_violations": d.get("Lipinski_violations", ""),
            "max_tc": round(max_tc, 4),
            "best_match": rec.get("best_match", ""),
            "best_target": rec.get("best_target", ""),
            "physchem_compliance": round(phys_score, 4),
            "toxic_score": round(toxic, 4),
            "final_score": round(final, 4),
        })

    # 综合评分降序
    rows.sort(key=lambda x: x["final_score"], reverse=True)
    for i, r in enumerate(rows, start=1):
        r["rank"] = i

    n_out = len(rows)
    log_phase("Phase 5: 综合评分与排序", n_in, n_out, time.time() - t0)
    return rows


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------

CSV_COLUMNS = [
    "rank", "name", "smiles", "MW", "LogP", "TPSA", "HBD", "HBA", "RotB",
    "Lipinski_violations", "max_tc", "best_match", "best_target",
    "physchem_compliance", "toxic_score", "final_score",
]


def write_csv(rows, out_csv):
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
    print(f"  [输出] 已写入 CSV : {out_csv}  ({len(rows)} 行)")


def write_sdf(records_ranked, out_sdf):
    """将排序后的候选分子写出为 SDF，属性字段写入评分信息。"""
    # records_ranked: 与 rows 对应的原始 rec 列表（按排名顺序）
    writer = Chem.SDWriter(out_sdf)
    for rank, rec in enumerate(records_ranked, start=1):
        mol = rec["mol"]
        try:
            mol.SetProp("rank", str(rank))
            mol.SetProp("final_score", f"{rec.get('final_score', 0):.4f}")
            mol.SetProp("max_tc", f"{rec.get('max_tc', 0):.4f}")
            mol.SetProp("best_match", str(rec.get("best_match", "")))
            mol.SetProp("best_target", str(rec.get("best_target", "")))
            mol.SetProp("physchem_compliance", f"{rec.get('physchem_compliance', 0):.4f}")
            mol.SetProp("toxic_score", f"{rec.get('toxic_score', 0.5):.4f}")
            writer.write(mol)
        except Exception:
            continue
    writer.close()
    print(f"  [输出] 已写入 SDF : {out_sdf}")


def print_top10(rows):
    print(f"\n{'=' * 80}")
    print("  Top 10 候选分子摘要")
    print(f"{'-' * 80}")
    header = f"{'rank':>4}  {'name':24} {'final':>6} {'maxTC':>7}  {'target':20} {'MW':>7}"
    print(header)
    print("-" * 80)
    for r in rows[:10]:
        nm = r["name"][:24]
        tg = str(r["best_target"])[:20]
        print(f"{r['rank']:>4}  {nm:24} {r['final_score']:>6.3f} {r['max_tc']:>7.3f}  {tg:20} {r['MW']:>7}")
    print(f"{'=' * 80}")


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="从 SDF 文件筛选 MASLD 降脂活性小分子（去重标准化→理化初筛→相似性搜索→ADMET毒性预测→综合排序）"
    )
    parser.add_argument("input_sdf", help="输入 SDF 文件路径")
    parser.add_argument("-o", "--output", default="ranked_candidates.csv", help="输出 CSV 路径 (默认 ranked_candidates.csv)")
    parser.add_argument("--output_sdf", action="store_true", help="同时输出排序后的 SDF 文件")
    parser.add_argument("--sim_threshold", type=float, default=None, help="相似性保留阈值 (默认使用 config.SIMILARITY_THRESHOLD)")
    args = parser.parse_args()

    if not os.path.isfile(args.input_sdf):
        print(f"[错误] 找不到输入文件: {args.input_sdf}")
        sys.exit(1)

    print(f"\n########## MASLD Screener 开始 ##########")
    print(f"输入文件 : {args.input_sdf}")
    print(f"输出文件 : {args.output}")

    t_total = time.time()

    # Phase 1
    recs = phase1_read_dedup(args.input_sdf)
    if not recs:
        print("[提示] Phase 1 后无分子，流程结束。")
        write_csv([], args.output)
        return

    # Phase 2
    recs = phase2_physchem(recs)
    if not recs:
        print("[提示] Phase 2 理化初筛后无分子，流程结束。")
        write_csv([], args.output)
        return

    # Phase 3
    recs = phase3_similarity(recs, sim_threshold=args.sim_threshold)
    if not recs:
        print("[提示] Phase 3 相似性搜索后无分子，流程结束。")
        write_csv([], args.output)
        return

    # Phase 4
    recs = phase4_admet(recs)
    if not recs:
        print("[提示] Phase 4 毒性筛选后无分子，流程结束。")
        write_csv([], args.output)
        return

    # Phase 5
    rows = phase5_score_and_rank(recs)

    # 输出
    write_csv(rows, args.output)
    if args.output_sdf:
        sdf_path = os.path.splitext(args.output)[0] + ".sdf"
        write_sdf(recs, sdf_path)

    print_top10(rows)

    print(f"\n总耗时 : {time.time() - t_total:.2f}s")
    print(f"########## MASLD Screener 完成 ##########")


if __name__ == "__main__":
    main()
