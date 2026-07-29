#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对 81 候选的 physchem 描述符做独立对账：
RDKit 计算值 (CSV 内)  vs  PubChem PUG-REST 报告值 (MW / XLogP / TPSA / HBD / HBA / RotB)
注意：PubChem 的 MW/XLogP/TPSA 也是 in-silico 计算值，非实验实测；
      本脚本用于"两套独立算法一致性校验"，严重不符才提示结构/SMILES 问题。
网络需放行（sandbox 外执行）。
"""
import csv, json, sys, time, urllib.request, urllib.error

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
SRC = "T001_ranked.csv"
OUT = "physchem_crosscheck.csv"

# 偏差阈值
TH_MW = 1.0        # Da 绝对差 -> 结构存疑
TH_TPSA = 5.0      # 绝对差 -> 存疑
TH_LOGP = 1.5      # 绝对差 -> 显著(算法不同, 容许更大)
TH_COUNT = 0       # HBD/HBA/RotB 不等 -> 计数差异

def fetch_props(cids):
    """批量拉 PubChem 属性, 返回 {cid: dict}"""
    url = f"{BASE}/compound/cid/{','.join(map(str,cids))}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MASLD-screen/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            out = {}
            for p in data.get("PropertyTable", {}).get("Properties", []):
                out[p["CID"]] = p
            return out
        except Exception as e:
            print(f"  batch fetch attempt {attempt+1} failed: {e}", file=sys.stderr)
            time.sleep(2)
    # 回退: 逐 CID
    out = {}
    for cid in cids:
        try:
            u = f"{BASE}/compound/cid/{cid}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
            req = urllib.request.Request(u, headers={"User-Agent": "MASLD-screen/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read().decode())
            ps = d.get("PropertyTable", {}).get("Properties", [])
            if ps:
                out[cid] = ps[0]
        except Exception as e:
            print(f"  cid {cid} fetch failed: {e}", file=sys.stderr)
        time.sleep(0.1)
    return out

def f(x):
    try:
        return float(x)
    except:
        return None

def main():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    cids = [int(r["pubchem_cid"]) for r in rows if r["pubchem_cid"].strip()]
    print(f"fetching PubChem properties for {len(cids)} cids ...")
    pc = fetch_props(cids)
    print(f"  got {len(pc)}/{len(cids)}")

    out_rows = []
    n_flag = 0
    for r in rows:
        cid = int(r["pubchem_cid"]) if r["pubchem_cid"].strip() else None
        p = pc.get(cid, {})
        pc_mw = f(p.get("MolecularWeight"))
        pc_xlogp = f(p.get("XLogP"))
        pc_tpsa = f(p.get("TPSA"))
        pc_hbd = f(p.get("HBondDonorCount"))
        pc_hba = f(p.get("HBondAcceptorCount"))
        pc_rot = f(p.get("RotatableBondCount"))

        rk_mw = f(r["MW"]); rk_logp = f(r["LogP"]); rk_tpsa = f(r["TPSA"])
        rk_hbd = f(r["HBD"]); rk_hba = f(r["HBA"]); rk_rot = f(r["RotB"])

        d_mw = (rk_mw - pc_mw) if (rk_mw is not None and pc_mw is not None) else None
        d_tpsa = (rk_tpsa - pc_tpsa) if (rk_tpsa is not None and pc_tpsa is not None) else None
        d_logp = (rk_logp - pc_xlogp) if (rk_logp is not None and pc_xlogp is not None) else None

        flags = []
        if d_mw is not None and abs(d_mw) > TH_MW:
            flags.append(f"MW差{abs(d_mw):.1f}")
        if d_tpsa is not None and abs(d_tpsa) > TH_TPSA:
            flags.append(f"TPSA差{abs(d_tpsa):.1f}")
        if d_logp is not None and abs(d_logp) > TH_LOGP:
            flags.append(f"LogP差{abs(d_logp):.1f}")
        for nm, a, b in [("HBD", rk_hbd, pc_hbd), ("HBA", rk_hba, pc_hba), ("RotB", rk_rot, pc_rot)]:
            if a is not None and b is not None and abs(a - b) > TH_COUNT:
                flags.append(f"{nm}:{int(a)}vs{int(b)}")
        flag = ";".join(flags)
        if flag:
            n_flag += 1

        out_rows.append({
            "rank": r["rank"], "name": r["name"], "cid": cid or "",
            "rdkit_MW": r["MW"], "pc_MW": p.get("MolecularWeight", ""),
            "d_MW": (f"{d_mw:.2f}" if d_mw is not None else ""),
            "rdkit_LogP": r["LogP"], "pc_XLogP": p.get("XLogP", ""),
            "d_LogP": (f"{d_logp:.2f}" if d_logp is not None else ""),
            "rdkit_TPSA": r["TPSA"], "pc_TPSA": p.get("TPSA", ""),
            "d_TPSA": (f"{d_tpsa:.2f}" if d_tpsa is not None else ""),
            "rdkit_HBD": r["HBD"], "pc_HBD": p.get("HBondDonorCount", ""),
            "rdkit_HBA": r["HBA"], "pc_HBA": p.get("HBondAcceptorCount", ""),
            "rdkit_RotB": r["RotB"], "pc_RotB": p.get("RotatableBondCount", ""),
            "flag": flag,
        })

    fields = ["rank","name","cid","rdkit_MW","pc_MW","d_MW","rdkit_LogP","pc_XLogP","d_LogP",
              "rdkit_TPSA","pc_TPSA","d_TPSA","rdkit_HBD","pc_HBD","rdkit_HBA","pc_HBA",
              "rdkit_RotB","pc_RotB","flag"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    print(f"\n=== 对账完成 ===")
    print(f"输出: {OUT}")
    print(f"被标记(任一偏差超阈值)的行数: {n_flag}/{len(rows)}")
    if n_flag:
        print("标记明细:")
        for o in out_rows:
            if o["flag"]:
                print(f"  #{o['rank']:>3} {o['name']:<10} -> {o['flag']}")
    else:
        print("两套算法在所有化合物上一致, 无结构/SMILES 存疑。")

if __name__ == "__main__":
    main()
