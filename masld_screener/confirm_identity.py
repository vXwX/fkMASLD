#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
逐一核验 T001_ranked.csv 中每个候选的身份：
- 按 InChIKey 精确（27 位）查询 PubChem，拿 CID + 收录名（Title / 第一条通用名同义词）
- 精确不中时，用连接性首块（14 位）查询，判断"骨架已知、立体异构未收录"
- 输出：known_flag / known_as / pubchem_cid / pubchem_name 四列

判定口径（known_flag）：
  参考药本体            max_tc≈1.0 且是我们自己的锚点药
  已知化合物(PubChem)   InChIKey 精确命中 PubChem，且有通用名
  已收录(仅系统名)      InChIKey 精确命中，但无通用名（仅 IUPAC/CID 名）
  骨架已知(立体异构未收录) 精确不中，但连接性层命中 → 该平面结构已知，此立体异构未单独收录
  疑似新颖(未收录)      精确 + 连接性均未命中 PubChem
"""
import sys, csv, json, time, ssl, argparse, urllib.request, urllib.error, urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import inchi

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# 手工确认的锚点/参考药（TC=1.0 或已核验），核验时保留并优先展示
CURATED = {
    "T7395": ("参考药本体", "Lanifibranor (PPARα/γ/δ, 临床III期)"),
    "T1789": ("参考药本体", "Obeticholic acid OCA (FXR, III期)"),
    "T3595": ("参考药本体", "Resmetirom (THR-β, 2024获批 Rezdiffra)"),
    "T7184": ("参考药本体", "Firsocostat GS-0976 (ACC, Phase 2)"),
    "TQ0243": ("已知类似物", "Firsocostat 类似物 (噻吩并嘧啶酰胺)"),
}


def _get(url, timeout=15, retries=2):
    last = None
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "identity-check/1.0"})
            r = urllib.request.urlopen(req, timeout=timeout, context=CTX)
            return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, b""
            last = e
        except Exception as e:
            last = e
        time.sleep(0.4)
    return None, (repr(last).encode() if last else b"")


def cids_by_inchikey(ik):
    st, data = _get(f"{BASE}/compound/inchikey/{ik}/cids/JSON")
    if st == 200:
        try:
            return json.loads(data).get("IdentifierList", {}).get("CID", [])
        except Exception:
            return []
    return []


def title_of(cid):
    st, data = _get(f"{BASE}/compound/cid/{cid}/property/Title/JSON")
    if st == 200:
        try:
            props = json.loads(data).get("PropertyTable", {}).get("Properties", [])
            if props:
                return props[0].get("Title", "")
        except Exception:
            pass
    return ""


def first_synonym(cid):
    st, data = _get(f"{BASE}/compound/cid/{cid}/synonyms/JSON")
    if st == 200:
        try:
            infos = json.loads(data).get("InformationList", {}).get("Information", [])
            if infos:
                syns = infos[0].get("Synonym", [])
                # 优先选一个像"通用名"的（含字母、非纯数字/CAS）
                for s in syns[:8]:
                    ss = s.strip()
                    if any(c.isalpha() for c in ss) and not ss.replace("-", "").isdigit():
                        return ss
                return syns[0] if syns else ""
        except Exception:
            pass
    return ""


def is_systematic(name):
    """粗略判断是否只是系统名（IUPAC）而非通用名。"""
    if not name:
        return True
    n = name.lower()
    markers = ["[", "]", "(1", "(2", "(3", "(4", "yl)", "oxy", "amino", "hydroxy",
               "carboxylic", "methyl", "ethyl", "-di", "-tri", "acid,", "beta-",
               "alpha-", "1,", "2,", "3,", "4,", "5,"]
    hits = sum(1 for m in markers if m in n)
    return hits >= 2 and len(name) > 25


def _post(url, data, timeout=25, retries=2):
    last = None
    for _ in range(retries):
        try:
            req = urllib.request.Request(
                url, data=urllib.parse.urlencode(data).encode(),
                headers={"User-Agent": "identity-check/1.0",
                         "Content-Type": "application/x-www-form-urlencoded"})
            r = urllib.request.urlopen(req, timeout=timeout, context=CTX)
            return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, b""
            last = e
        except Exception as e:
            last = e
        time.sleep(0.4)
    return None, b""


def cids_same_connectivity(smi):
    """用 fastidentity (SMILES POST) 按 same_connectivity 检索同平面骨架。"""
    url = (f"{BASE}/compound/fastidentity/smiles/cids/JSON"
           f"?identity_type=same_connectivity")
    st, data = _post(url, {"smiles": smi})
    if st == 200:
        try:
            return json.loads(data).get("IdentifierList", {}).get("CID", [])
        except Exception:
            return []
    return []


def _name_for(cid):
    title = title_of(cid)
    syn = first_synonym(cid)
    name = title or syn
    # 若 Title 像系统名而同义词像通用名，优先通用名
    if title and syn and (not is_systematic(syn)) and is_systematic(title):
        name = syn
    return name


def classify(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return dict(flag="解析失败", known_as="", cid="", pname="", ik="")
    ik = inchi.MolToInchiKey(m)
    # 1) InChIKey 精确匹配（含立体）
    cids = cids_by_inchikey(ik)
    if cids:
        cid = cids[0]
        name = _name_for(cid)
        if is_systematic(name):
            return dict(flag="已收录(仅系统名)", known_as=(name or "")[:90], cid=str(cid),
                        pname=name, ik=ik)
        return dict(flag="已知化合物(PubChem)", known_as=name, cid=str(cid),
                    pname=name, ik=ik)
    # 2) 连接性层匹配（同骨架、立体异构可能未单独收录）
    ccids = cids_same_connectivity(smi)
    if ccids:
        cid = ccids[0]
        name = _name_for(cid)
        tag = (name or "").strip()
        return dict(flag="骨架已知(立体异构未收录)",
                    known_as=(f"同骨架: {tag}" if tag else "同平面结构已收录"),
                    cid=str(cid), pname=name, ik=ik)
    # 3) 精确 + 连接性均未命中
    return dict(flag="疑似新颖(未收录)", known_as="", cid="", pname="", ik=ik)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="masld_screener/T001_ranked.csv")
    ap.add_argument("-o", "--output", default="masld_screener/T001_ranked.csv")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.input, encoding="utf-8")))
    print(f"读取 {len(rows)} 个候选，开始 PubChem 逐一核验 (InChIKey 精确)...", flush=True)

    def work(r):
        res = classify(r["smiles"])
        return r["name"], res

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(work, r): r for r in rows}
        done = 0
        for fut in as_completed(futs):
            nm, res = fut.result()
            results[nm] = res
            done += 1
            if done % 10 == 0 or done == len(rows):
                print(f"  已核验 {done}/{len(rows)}", flush=True)

    # 写回
    for r in rows:
        nm = r["name"]
        res = results.get(nm, {})
        cur = CURATED.get(nm)
        if cur:
            r["known_flag"] = cur[0]
            # 参考药也附上 PubChem 名/CID
            r["known_as"] = cur[1]
        else:
            r["known_flag"] = res.get("flag", "疑似新颖(未收录)")
            r["known_as"] = res.get("known_as", "")
        r["pubchem_cid"] = res.get("cid", "")
        r["pubchem_name"] = res.get("pname", "")

    fields = [c for c in rows[0].keys() if c not in ("pubchem_cid", "pubchem_name")]
    # 确保新列在末尾
    for c in ("known_flag", "known_as"):
        if c not in fields:
            fields.append(c)
    fields = list(dict.fromkeys(fields)) + ["pubchem_cid", "pubchem_name"]

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    # 汇总
    from collections import Counter
    c = Counter(r["known_flag"] for r in rows)
    print("\n=== 核验汇总 ===")
    for k, v in c.most_common():
        print(f"  {k}: {v}")
    n_pub = sum(1 for r in rows if r.get("pubchem_cid"))
    print(f"  PubChem 精确命中 CID: {n_pub}/{len(rows)}")
    print(f"\n已写回 {args.output}")


if __name__ == "__main__":
    main()
