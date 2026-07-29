#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""给 T001_ranked.csv 增补一列 toxic_score_basis, 明确 toxic_score / final_score 为
ADMET-AI 机器学习预测, 非实验实测毒理数据。不改任何现有列与数值。"""
import csv, shutil, os

SRC = "T001_ranked.csv"
BACKUP_DIR = "../_backups"
COL = "toxic_score_basis"
VAL = "ADMET-AI预测得分·非实验实测毒理(hERG/AMES/DILI实验库零命中;final_score由其推导)"

os.makedirs(BACKUP_DIR, exist_ok=True)
backup = os.path.join(BACKUP_DIR, "T001_ranked.backup_before_scorebasis.csv")
shutil.copy(SRC, backup)
print(f"已备份原 CSV -> {backup}")

rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
fields = list(rows[0].keys())
if COL not in fields:
    fields.append(COL)
for r in rows:
    r[COL] = VAL

with open(SRC, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"已新增列 {COL}, 共 {len(rows)} 行, 现有列数 = {len(fields)}")
