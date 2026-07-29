#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""对 confirm_identity 的结果做人工精修：补全空名、修正真实身份、标注结构假阳性。"""
import csv

SRC = "masld_screener/T001_ranked.csv"

# 1) 空 known_as 补全（有 CID 但同义词库无通用名，多为钠盐/系统名）
FILL = {
    "TN2215": "牛磺酸结合胆酸钠盐类 (PubChem 已收录, 无通用名)",
    "TN2349": "牛磺酸结合胆酸钠盐类 (PubChem 已收录, 无通用名)",
}

# 2) 修正/明确真实身份（PubChem 名比手工标签更权威）
OVERRIDE_AS = {
    "TQ0243": "ND-646 (已知 ACC 抑制剂; 原标为 Firsocostat 类似物)",
}
OVERRIDE_FLAG = {
    "TQ0243": "已知药/工具化合物",
}

# 3) 结构相似但适应症无关的已知药 —— 结构假阳性（务必点名，避免误判为 MASLD 候选）
IRRELEVANT = {
    "T1172":  "⚠️结构假阳性: Diclazuril 是抗球虫兽药, 仅与 Resmetirom 共享三嗪二酮子结构, 与 MASLD/THR-β 无关",
    "T19526": "⚠️结构假阳性: Pradefovir 是抗乙肝前药(adefovir 前药), 仅与 VK2809 共享环磷酸酯前药子结构, 与 THR-β 无关",
    "T10953": "⚠️结构假阳性: Dafadine-A 是 CYP/DAF-9 工具抑制剂, 仅与 Denifanstat 结构相似, 与 FASN/MASLD 无关",
}

rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
fields = list(rows[0].keys())
if "relevance_note" not in fields:
    fields.append("relevance_note")

for r in rows:
    nm = r["name"]
    if nm in FILL and not r.get("known_as"):
        r["known_as"] = FILL[nm]
    if nm in OVERRIDE_AS:
        r["known_as"] = OVERRIDE_AS[nm]
    if nm in OVERRIDE_FLAG:
        r["known_flag"] = OVERRIDE_FLAG[nm]
    r["relevance_note"] = IRRELEVANT.get(nm, "")

with open(SRC, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print("补丁完成。relevance_note 已标注结构假阳性:",
      [r["name"] for r in rows if r["relevance_note"]])
print("空 known_as 剩余:", [r["name"] for r in rows if not r["known_as"]])
