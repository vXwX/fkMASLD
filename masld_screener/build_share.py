#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打包共享 zip：masld_screener/ 全量(排除 _backups / backup csv / __pycache__) + SHARE_NOTES.md。"""
import os, zipfile, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # D:\AI4S
SRC = os.path.join(ROOT, "masld_screener")
NOTES = os.path.join(ROOT, "SHARE_NOTES.md")
OUT = os.path.join(ROOT, "masld_screener_share_20260729.zip")

EXCLUDE_DIRS = {"_backups", "__pycache__", ".git"}
EXCLUDE_SUFFIX = (".backup", ".backup_before_identity.csv", ".backup_before_scorebasis.csv")

def keep(path):
    bn = os.path.basename(path)
    if any(bn.endswith(s) for s in EXCLUDE_SUFFIX):
        return False
    return True

def main():
    if os.path.exists(OUT):
        os.remove(OUT)
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
        # masld_screener
        for dp, dns, fns in os.walk(SRC):
            dns[:] = [d for d in dns if d not in EXCLUDE_DIRS]
            for fn in fns:
                fp = os.path.join(dp, fn)
                if not keep(fp):
                    continue
                arc = os.path.relpath(fp, ROOT)
                z.write(fp, arc)
        # notes
        if os.path.exists(NOTES):
            z.write(NOTES, "SHARE_NOTES.md")
    size = os.path.getsize(OUT)
    # list contents
    with zipfile.ZipFile(OUT) as z:
        names = z.namelist()
    print(f"打包完成: {OUT} ({size} bytes, {len(names)} 个条目)")
    # 确认无 backup
    backs = [n for n in names if "backup" in n.lower() or "_backups" in n.lower()]
    print("包含 backup 文件:", backs if backs else "无")
    # 确认关键文件在场
    for key in ["T001_ranked.csv", "physchem_crosscheck.csv", "crosscheck_physchem.py",
                "confirm_identity.py", "SHARE_NOTES.md"]:
        hit = [n for n in names if n.endswith(key)]
        print(f"  {key}: {'OK' if hit else '缺失!!!'}")

if __name__ == "__main__":
    main()
