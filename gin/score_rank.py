import os

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(SCRIPT_DIR, "results")
OUT_MOL = os.path.join(RESULT_DIR, "mols.csv")
OUT_SCORES = os.path.join(RESULT_DIR, "all_scores.csv")
OUT_RANK = os.path.join(RESULT_DIR, "predicted_inactive_ranked.csv")

ACT_LOW, ACT_HIGH = 0.05, 0.30  # activity thresholds (same as ranking pipeline)


def act_level(s):
    if s <= ACT_LOW:
        return "低"
    if s <= ACT_HIGH:
        return "中"
    return "高"


def main():
    scores = pd.read_csv(OUT_SCORES)
    mols = pd.read_csv(OUT_MOL, usecols=["ID", "has_active"])

    # keep only molecules without any activity data (prediction targets)
    targets = mols[mols["has_active"].isna()]["ID"].tolist()
    df = scores[scores["ID"].isin(targets)].copy()
    print(f"target molecules (no activity data): {len(df)}")

    # activity score = P_act * strength
    df["active_score"] = df["P_act"] * df["STRENGTH"]
    # toxicity signals
    df["P_tox_any"] = df["P_tox_1"] + df["P_tox_2"]   # any toxicity probability
    df["P_tox_high"] = df["P_tox_2"]                  # high-risk class probability
    df["toxicity_score"] = df["P_tox_any"]
    # recommendation score
    df["recommend_score"] = df["active_score"] * (1.0 - df["toxicity_score"])

    # activity level from thresholds
    df["activity_level"] = df["active_score"].apply(act_level)

    # toxicity levels: quantile-based within target population (model is biased
    # toward high-risk class, so relative thresholds give a meaningful spread)
    t_high = df["P_tox_high"]
    t_any = df["toxicity_score"]
    hi = np.percentile(t_high, 66.7)
    lo = np.percentile(t_any, 33.3)

    def tox_level(row):
        if row["P_tox_high"] > hi:
            return "高毒"
        if row["toxicity_score"] > lo:
            return "毒"
        return "低毒"

    df["toxicity_level"] = df.apply(tox_level, axis=1)

    df = df.sort_values("recommend_score", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)

    cols = [
        "rank", "ID", "SMILES", "active_score", "activity_level",
        "toxicity_score", "P_tox_high", "toxicity_level", "recommend_score",
    ]
    df[cols].to_csv(OUT_RANK, index=False)

    print(f"Saved {OUT_RANK} ({len(df)} rows)")
    print("\nActivity level distribution:")
    print(df["activity_level"].value_counts().to_string())
    print("\nToxicity level distribution:")
    print(df["toxicity_level"].value_counts().to_string())
    print("\nTop 20 by recommend score:")
    print(df[cols[:9]].head(20).to_string(index=False))
    print(f"\nthresholds: high_risk>{hi:.3f}, tox>{lo:.3f}")


if __name__ == "__main__":
    main()