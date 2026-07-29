#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""给 physchem_crosscheck.csv 增加 severity / interpretation 列, 避免'74/81 被标记'被误读为数据错误。"""
import csv

SRC = "physchem_crosscheck.csv"

def classify(flag):
    if not flag:
        return "OK", "RDKit 与 PubChem 完全一致"
    parts = flag.split(";")
    has_mw = any(p.startswith("MW差") for p in parts)
    has_tpsa = any(p.startswith("TPSA差") for p in parts)
    has_logp = any(p.startswith("LogP差") for p in parts)
    has_count = any(p.startswith(("HBD:", "HBA:", "RotB:")) for p in parts)
    sev = []
    notes = []
    if has_mw:
        sev.append("结构存疑"); notes.append("MW 不符 -> 需核对 SMILES/分子式")
    if has_tpsa:
        sev.append("TPSA方法学"); notes.append("TPSA 差异多因盐/电荷/互变异构的取质子态不同, 非错误")
    if has_logp:
        sev.append("LogP算法差"); notes.append("RDKit MolLogP 与 PubChem XLogP3 为不同算法, 差异属正常")
    if has_count and not (has_mw or has_tpsa or has_logp):
        sev.append("计数惯例"); notes.append("HBD/HBA/RotB 计数口径不同(RDKit vs PubChem), 不影响筛选")
    if not sev:
        sev.append("其它")
    return "/".join(sev), "; ".join(notes)

rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
missing = [r for r in rows if not r["pc_MW"].strip()]
for r in rows:
    sev, note = classify(r["flag"])
    r["severity"] = sev
    r["interpretation"] = note

fields = list(rows[0].keys())
# 调整列序: severity/interp 放最后(已是最后两列)
with open(SRC, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

from collections import Counter
c = Counter(r["severity"] for r in rows)
print("=== 对账分类汇总 ===")
for k, v in c.most_common():
    print(f"  {k}: {v}")
print(f"\nMW 不符(结构存疑)数量: {sum(1 for r in rows if '结构存疑' in r['severity'])}")
print(f"缺失 PubChem 属性的 CID: {[ (r['rank'], r['name'], r['cid']) for r in missing ]}")
print("\n需注意的 TPSA/LogP 显著项(非错误, 方法学差异):")
for r in rows:
    if "TPSA方法学" in r["severity"] or "LogP算法差" in r["severity"]:
        print(f"  #{r['rank']:>3} {r['name']:<10} flag={r['flag']}")
