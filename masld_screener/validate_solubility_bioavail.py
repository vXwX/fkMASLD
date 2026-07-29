# -*- coding: utf-8 -*-
"""
validate_solubility_bioavail.py
用公开实验数据集交叉核验 admet-ai 的"溶解度 / 口服生物利用度"预测。

数据源（均为 admet-ai 对应端点的原始训练/基准实验集）：
  - AqSolDB        (9982 化合物, 实验水溶解度 LogS, 连续值)  -> admet-ai 端点 Solubility_AqSolDB
  - Bioavailability_Ma (640 化合物, 口服生物利用度 0/1 二分类) -> admet-ai 端点 Bioavailability_Ma
  经 Harvard Dataverse (TDC) 下载: datafile 4259610 / 4259567

匹配方式：InChIKey 精确匹配（27 位全串）。同时记录 InChIKey 首块（14 位连接性）匹配，
供人工判断立体异构差异。

⚠️ 关键局限（务必写进交付说明）：
  AqSolDB / Bioavailability_Ma 正是 admet-ai 这两个端点的**训练数据集**。
  若候选命中，说明该分子极可能在训练集中——此时"预测≈实验"反映的是**模型记忆**，
  并非独立泛化能力的证明。因此本比对结论应表述为：
  "命中候选的 admet-ai 报告值与公开实验值是否自洽"，而非"模型预测多准"。

用法:
  python validate_solubility_bioavail.py -i T001_ranked.csv -o solubility_bioavail_report.csv
"""
import argparse, csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DS_DIR = os.path.join(ROOT, "exp_datasets")


def _find_column(columns, key):
    norm = [c.lower().replace("_", "").replace(" ", "").replace("-", "") for c in columns]
    k = key.lower().replace("_", "").replace(" ", "").replace("-", "")
    for i, c in enumerate(norm):
        if k in c:
            return columns[i]
    return None


def inchikey(smiles):
    from rdkit import Chem
    from rdkit.Chem import inchi
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        return inchi.MolToInchiKey(m)
    except Exception:
        return None


def load_dataset(path):
    """返回 {full_inchikey: (Y, drug_id, smiles)} 与 {first_block: [(Y, drug_id, smiles)...]}"""
    full, block = {}, {}
    with open(path, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            smi = (row.get("Drug") or "").strip().strip('"')
            y = row.get("Y")
            did = (row.get("Drug_ID") or "").strip().strip('"')
            if not smi:
                continue
            ik = inchikey(smi)
            if not ik:
                continue
            full[ik] = (y, did, smi)
            block.setdefault(ik.split("-")[0], []).append((y, did, smi))
    return full, block


def run_admet(smiles_list):
    from admet_ai import ADMETModel
    model = ADMETModel()
    return model.predict(smiles_list)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", default="solubility_bioavail_report.csv")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.input, encoding="utf-8")))
    scol = _find_column(list(rows[0].keys()), "smiles")
    ncol = _find_column(list(rows[0].keys()), "name") or scol
    smiles_list = [r[scol] for r in rows]
    names = [r[ncol] for r in rows]
    print(f"候选数: {len(smiles_list)}")

    print("加载实验数据集 (AqSolDB / Bioavailability_Ma) ...", flush=True)
    sol_full, sol_block = load_dataset(os.path.join(DS_DIR, "solubility_aqsoldb.raw"))
    bio_full, bio_block = load_dataset(os.path.join(DS_DIR, "bioavailability_ma.raw"))
    print(f"  AqSolDB: {len(sol_full)} 条 | Bioavailability_Ma: {len(bio_full)} 条", flush=True)

    print("运行 admet-ai 预测 ...", flush=True)
    df = run_admet(smiles_list)
    sol_col = _find_column(df.columns, "solubility_aqsoldb") or _find_column(df.columns, "aqueoussolubility") or _find_column(df.columns, "solubility")
    bio_col = _find_column(df.columns, "bioavailability_ma") or _find_column(df.columns, "bioavailability")
    print(f"  admet-ai 溶解度列: {sol_col} | 生物利用度列: {bio_col}", flush=True)
    pred = df.reset_index(drop=True)

    out = []
    for i, (name, smi) in enumerate(zip(names, smiles_list)):
        ik = inchikey(smi)
        block = ik.split("-")[0] if ik else ""
        pred_sol = pred.iloc[i][sol_col] if sol_col else ""
        pred_bio = pred.iloc[i][bio_col] if bio_col else ""

        # solubility exact / block match
        s_exp, s_kind, s_src = "", "", ""
        if ik in sol_full:
            s_exp, s_src, _ = sol_full[ik]; s_kind = "exact"
        elif block in sol_block:
            s_exp, s_src, _ = sol_block[block][0]; s_kind = "block"

        # bioavailability exact / block match
        b_exp, b_kind, b_src = "", "", ""
        if ik in bio_full:
            b_exp, b_src, _ = bio_full[ik]; b_kind = "exact"
        elif block in bio_block:
            b_exp, b_src, _ = bio_block[block][0]; b_kind = "block"

        # solubility abs error
        s_err = ""
        if s_exp != "" and pred_sol != "":
            try:
                s_err = round(abs(float(pred_sol) - float(s_exp)), 3)
            except Exception:
                s_err = ""
        # bioavailability agreement (pred>=0.5 -> 1)
        b_agree = ""
        if b_exp != "" and pred_bio != "":
            try:
                b_pred_cls = 1 if float(pred_bio) >= 0.5 else 0
                b_agree = "MATCH" if b_pred_cls == int(float(b_exp)) else "MISMATCH"
            except Exception:
                b_agree = ""

        out.append({
            "name": name, "inchikey": ik or "",
            "pred_LogS": round(float(pred_sol), 3) if pred_sol != "" else "",
            "exp_LogS": s_exp, "sol_match": s_kind, "sol_abs_err": s_err, "sol_src": s_src,
            "pred_bioavail_prob": round(float(pred_bio), 3) if pred_bio != "" else "",
            "exp_bioavail_cls": b_exp, "bio_match": b_kind, "bio_agree": b_agree, "bio_src": b_src,
        })

    fields = ["name", "inchikey", "pred_LogS", "exp_LogS", "sol_match", "sol_abs_err", "sol_src",
              "pred_bioavail_prob", "exp_bioavail_cls", "bio_match", "bio_agree", "bio_src"]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(out)

    # summary
    sol_hits = [r for r in out if r["sol_match"]]
    bio_hits = [r for r in out if r["bio_match"]]
    sol_exact = [r for r in sol_hits if r["sol_match"] == "exact"]
    bio_exact = [r for r in bio_hits if r["bio_match"] == "exact"]
    errs = [r["sol_abs_err"] for r in sol_exact if r["sol_abs_err"] != ""]
    print("\n========== 溶解度 / 生物利用度 实验比对汇总 ==========")
    print(f"溶解度 AqSolDB 命中: 精确 {len(sol_exact)} / 连接性 {len(sol_hits)-len(sol_exact)} (共 {len(sol_hits)}/{len(out)})")
    if errs:
        errs_sorted = sorted(errs)
        mae = sum(errs) / len(errs)
        med = errs_sorted[len(errs_sorted)//2]
        within1 = sum(1 for e in errs if e <= 1.0)
        print(f"  精确命中 LogS 绝对误差: MAE={mae:.2f}  中位={med:.2f}  |误差<=1 log 单位|: {within1}/{len(errs)}")
    print(f"生物利用度 Ma 命中: 精确 {len(bio_exact)} / 连接性 {len(bio_hits)-len(bio_exact)} (共 {len(bio_hits)}/{len(out)})")
    b_ok = [r for r in bio_exact if r["bio_agree"] == "MATCH"]
    if bio_exact:
        print(f"  精确命中分类一致: {len(b_ok)}/{len(bio_exact)}")
    print("\n[数据泄漏提醒] AqSolDB / Bioavailability_Ma 即 admet-ai 训练集，命中=模型很可能已见过该分子；")
    print("  故上述一致性反映'报告值是否与公开实验自洽'，非独立泛化准确度。")
    print(f"\n报告已写入: {args.output}")


if __name__ == "__main__":
    main()
