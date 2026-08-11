import argparse
import os
import pickle

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from gin_dataset import MolGraphDataset, collate_graphs
from gin_model import MultiTaskGIN
from config import ATOM_DIM, BOND_DIM, HIDDEN, N_LAYERS

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(SCRIPT_DIR, "results")
OUT_MOL = os.path.join(RESULT_DIR, "mols.csv")
OUT_GRAPH = os.path.join(RESULT_DIR, "graphs.pkl")
MODEL_PATH = os.path.join(RESULT_DIR, "multitask_gin.pt")
OUT_SCORES = os.path.join(RESULT_DIR, "all_scores.csv")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--batch", type=int, default=512)
    parser.add_argument("--model", default=MODEL_PATH,
                        help="Path to model state_dict (default: CV best model).")
    parser.add_argument("--only-predict-targets", action="store_true",
                        help="Only score molecules without activity labels (the 8464).")
    args = parser.parse_args()

    device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu")
    print("device:", device)

    mol_df = pd.read_csv(OUT_MOL)
    with open(OUT_GRAPH, "rb") as f:
        graph = pickle.load(f)

    if args.only_predict_targets:
        mol_df = mol_df[mol_df["has_active"].isna()].copy()
        print("predicting only target molecules:", len(mol_df))
    else:
        print("predicting all molecules:", len(mol_df))

    # build records (graph + placeholder labels)
    records = []
    for _, row in mol_df.iterrows():
        records.append({
            "graph": graph[row["ID"]],
            "has_active": None,
            "act_ratio": None,
            "has_tox": None,
            "tox_level": None,
        })
    ds = MolGraphDataset(records)
    loader = DataLoader(ds, batch_size=args.batch, shuffle=False, num_workers=0, collate_fn=collate_graphs)

    # infer input dim from first record and verify it matches the shared config
    in_dim = records[0]["graph"]["x"].shape[1]
    assert in_dim == ATOM_DIM, f"graph x dim {in_dim} != config ATOM_DIM {ATOM_DIM}"
    model = MultiTaskGIN(in_dim, BOND_DIM, hidden=HIDDEN, n_layers=N_LAYERS, n_tox=3).to(device)
    model.load_state_dict(torch.load(args.model, map_location=device))
    model.eval()

    p_act_all, strength_all, tox_probs_all = [], [], []
    with torch.no_grad():
        for batch in loader:
            x = batch["x"].to(device)
            ei = batch["edge_index"].to(device)
            ea = batch["edge_attr"].to(device)
            bj = batch["batch_idx"].to(device)
            ng = batch["num_graphs"]
            p_act, strength, tox_logits = model(x, ei, ea, bj, ng)
            tox_probs = torch.softmax(tox_logits, dim=-1)
            p_act_all.append(p_act.detach().cpu().numpy())
            strength_all.append(strength.detach().cpu().numpy())
            tox_probs_all.append(tox_probs.detach().cpu().numpy())

    p_act = np.concatenate(p_act_all)
    strength = np.concatenate(strength_all)
    tox_probs = np.concatenate(tox_probs_all)

    out = mol_df[["ID", "SMILES"]].copy()
    out["P_act"] = p_act
    out["STRENGTH"] = strength
    out["P_tox_0"] = tox_probs[:, 0]
    out["P_tox_1"] = tox_probs[:, 1]
    out["P_tox_2"] = tox_probs[:, 2]
    out.to_csv(OUT_SCORES, index=False)
    print(f"Saved {OUT_SCORES} ({len(out)} rows)")


if __name__ == "__main__":
    main()