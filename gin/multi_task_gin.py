import argparse
import json
import os
import pickle

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import scipy.stats as ss
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from torch.utils.data import DataLoader

from gin_dataset import MolGraphDataset, collate_graphs
from gin_model import MultiTaskGIN
from config import ATOM_DIM, BOND_DIM, HIDDEN, LR, N_LAYERS, PATIENCE

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(SCRIPT_DIR, "results")
OUT_MOL = os.path.join(RESULT_DIR, "mols.csv")
OUT_GRAPH = os.path.join(RESULT_DIR, "graphs.pkl")

MODEL_PATH = os.path.join(RESULT_DIR, "multitask_gin.pt")
FULL_MODEL_PATH = os.path.join(RESULT_DIR, "multitask_gin_full.pt")
REPORT_PATH = os.path.join(RESULT_DIR, "validation_report.json")


def build_records(mol_df, graph):
    records = []
    for _, row in mol_df.iterrows():
        m_id = row["ID"]
        records.append({
            "graph": graph[m_id],
            "has_active": None if pd.isna(row["has_active"]) else int(row["has_active"]),
            "act_ratio": float(row["act_ratio"]) if not pd.isna(row["act_ratio"]) else None,
            "has_tox": None if pd.isna(row["has_tox"]) else int(row["has_tox"]),
            "tox_level": None if pd.isna(row["tox_level"]) else int(row["tox_level"]),
            "group": row["active_targets"] if pd.notna(row["active_targets"]) and str(row["active_targets"]) else "no_target",
        })
    return records


def loss_fn(batch, model, device):
    x = batch["x"].to(device)
    ei = batch["edge_index"].to(device)
    ea = batch["edge_attr"].to(device)
    bj = batch["batch_idx"].to(device)
    ng = batch["num_graphs"]
    has_active = batch["has_active"].to(device)
    act_ratio = batch["act_ratio"].to(device)
    tox_level = batch["tox_level"].to(device)

    p_act, strength, tox_logits = model(x, ei, ea, bj, ng)

    # activity classification (BCE on labeled)
    act_mask = ~torch.isnan(has_active)
    L_clf = nn.functional.binary_cross_entropy(p_act[act_mask], has_active[act_mask])

    # activity regression (on positive samples only), target = act_ratio in [0,1]
    pos_mask = act_mask & (has_active > 0)
    L_reg = torch.tensor(0.0, device=device)
    if pos_mask.any():
        L_reg = nn.functional.mse_loss(strength[pos_mask], act_ratio[pos_mask])

    # toxicity multiclass CE (only labeled)
    tox_mask = tox_level >= 0
    L_tox = torch.tensor(0.0, device=device)
    if tox_mask.sum() > 0:
        L_tox = nn.functional.cross_entropy(tox_logits[tox_mask], tox_level[tox_mask])

    return L_clf + L_reg + L_tox, L_clf, L_reg, L_tox


def evaluate(model, loader, device):
    model.eval()
    ys, preds = [], []
    reg_t, reg_p = [], []
    tox_t, tox_p = [], []
    with torch.no_grad():
        for batch in loader:
            x = batch["x"].to(device)
            ei = batch["edge_index"].to(device)
            ea = batch["edge_attr"].to(device)
            bj = batch["batch_idx"].to(device)
            ng = batch["num_graphs"]
            p_act, strength, tox_logits = model(x, ei, ea, bj, ng)
            act_mask = ~torch.isnan(batch["has_active"])
            if act_mask.sum() > 0:
                ys.append(batch["has_active"][act_mask].numpy())
                preds.append(p_act[act_mask].detach().cpu().numpy())
            pos_mask = act_mask & (batch["has_active"] > 0)
            if pos_mask.sum() > 0:
                reg_t.append(batch["act_ratio"][pos_mask].numpy())
                reg_p.append(strength[pos_mask].detach().cpu().numpy())
            tox_mask = batch["tox_level"] >= 0
            if tox_mask.sum() > 0:
                tox_t.append(batch["tox_level"][tox_mask].numpy())
                tox_p.append(torch.softmax(tox_logits[tox_mask], dim=-1).detach().cpu().numpy())
    return (np.concatenate(ys) if ys else None, np.concatenate(preds) if preds else None,
            np.concatenate(reg_t) if reg_t else None, np.concatenate(reg_p) if reg_p else None,
            np.concatenate(tox_t) if tox_t else None, np.concatenate(tox_p, axis=0) if tox_p else None)


def metrics(y, p, rt, rp, tt, tp):
    m = {}
    if y is not None and len(np.unique(y)) > 1:
        m["auc_roc_act"] = float(roc_auc_score(y, p))
        m["auc_pr_act"] = float(average_precision_score(y, p))
    if rt is not None and len(rt) > 1 and np.std(rt) > 0 and np.std(rp) > 0:
        m["pearson_act"] = float(ss.pearsonr(rt, rp)[0])
        m["spearman_act"] = float(ss.spearmanr(rt, rp)[0])
        m["rmse_reg"] = float(np.sqrt(np.mean((rt - rp) ** 2)))
    if tt is not None and len(np.unique(tt)) > 1:
        from sklearn.metrics import accuracy_score
        m["tox_acc"] = float(accuracy_score(tt, np.argmax(tp, axis=1)))
        # macro AUC one-vs-rest
        from sklearn.metrics import roc_auc_score as ras
        macro = []
        for c in np.unique(tt):
            yb = (tt == c).astype(int)
            if len(np.unique(yb)) > 1:
                macro.append(ras(yb, tp[:, int(c)]))
        if macro:
            m["tox_auc_macro"] = float(np.mean(macro))
    return m


def train_full(records, args, device):
    """Train one final model on all labeled data (no CV, no early stopping)."""
    print(f"\n=== Full retrain on all {len(records)} labeled records ===")
    ds = MolGraphDataset(records)
    loader = DataLoader(ds, batch_size=args.batch, shuffle=True, num_workers=0, collate_fn=collate_graphs)

    model = MultiTaskGIN(ATOM_DIM, BOND_DIM, hidden=HIDDEN, n_layers=N_LAYERS, n_tox=3).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)

    for ep in range(args.epochs):
        model.train()
        tot = tot_c = tot_r = tot_t = 0.0
        nb = 0
        for batch in loader:
            opt.zero_grad()
            L, Lc, Lr, Lt = loss_fn(batch, model, device)
            L.backward()
            opt.step()
            tot += L.item(); tot_c += Lc.item(); tot_r += Lr.item(); tot_t += Lt.item()
            nb += 1
        if ep == 0 or (ep + 1) % 5 == 0:
            print(f"  ep{ep}: L={tot/nb:.3f} clf={tot_c/nb:.3f} reg={tot_r/nb:.4f} tox={tot_t/nb:.3f}")

    state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    torch.save(state, FULL_MODEL_PATH)
    print(f"Saved full-retrain model: {FULL_MODEL_PATH}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--full", action="store_true",
                        help="Retrain one final model on all labeled data (no CV).")
    args = parser.parse_args()

    device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu")
    print("device:", device, "| cuda:", torch.cuda.is_available())

    mol_df = pd.read_csv(OUT_MOL)
    with open(OUT_GRAPH, "rb") as f:
        graph = pickle.load(f)

    records = [r for r in build_records(mol_df, graph) if r["has_active"] is not None]
    print("labeled records:", len(records))
    if not records:
        print("ERROR: no labeled records")
        return

    if args.full:
        train_full(records, args, device)
        return

    groups = [r["group"] for r in records]
    gkf = GroupKFold(n_splits=args.folds)
    split_data = list(gkf.split(records, groups=groups))

    fold_reports = []
    best_overall_auc = -1
    best_state = None

    for fold, (tr_idx, va_idx) in enumerate(split_data):
        print(f"\n=== Fold {fold+1}/{args.folds} (train {len(tr_idx)}, val {len(va_idx)}) ===")
        tr_ds = MolGraphDataset([records[i] for i in tr_idx])
        va_ds = MolGraphDataset([records[i] for i in va_idx])
        tr_loader = DataLoader(tr_ds, batch_size=args.batch, shuffle=True, num_workers=0, collate_fn=collate_graphs)
        va_loader = DataLoader(va_ds, batch_size=512, shuffle=False, num_workers=0, collate_fn=collate_graphs)

        model = MultiTaskGIN(ATOM_DIM, BOND_DIM, hidden=HIDDEN, n_layers=N_LAYERS, n_tox=3).to(device)
        opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-4)
        best_auc = -1
        best_fold_state = None
        patience = 0

        for ep in range(args.epochs):
            model.train()
            tot = tot_c = tot_r = tot_t = 0.0
            nb = 0
            for batch in tr_loader:
                opt.zero_grad()
                L, Lc, Lr, Lt = loss_fn(batch, model, device)
                L.backward()
                opt.step()
                tot += L.item(); tot_c += Lc.item(); tot_r += Lr.item(); tot_t += Lt.item()
                nb += 1
            y, p, rt, rp, tt, tp = evaluate(model, va_loader, device)
            m = metrics(y, p, rt, rp, tt, tp)
            auc = m.get("auc_roc_act", -1)
            if auc > best_auc:
                best_auc = auc
                best_fold_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                patience = 0
            else:
                patience += 1
            if fold == 0 and (ep == 0 or (ep + 1) % 5 == 0):
                print(f"  ep{ep}: L={tot/nb:.3f} clf={tot_c/nb:.3f} reg={tot_r/nb:.4f} tox={tot_t/nb:.3f} | va_act_auc={auc:.4f}")
            if patience >= PATIENCE:
                print(f"  early stop at ep{ep}, best_auc={best_auc:.4f}")
                break

        model.load_state_dict(best_fold_state)
        y, p, rt, rp, tt, tp = evaluate(model, va_loader, device)
        m = metrics(y, p, rt, rp, tt, tp)
        m["fold"] = fold
        fold_reports.append(m)
        print("  val: ", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()})

        if best_auc > best_overall_auc:
            best_overall_auc = best_auc
            best_state = best_fold_state

    torch.save(best_state, MODEL_PATH)
    print(f"\nSaved best model: {MODEL_PATH} (val_auc={best_overall_auc:.4f})")

    # aggregate
    report = {
        "folds": fold_reports,
        "best_val_auc": round(best_overall_auc, 4),
    }
    for key in ["auc_roc_act", "auc_pr_act", "pearson_act", "spearman_act", "rmse_reg", "tox_acc", "tox_auc_macro"]:
        vals = [f[key] for f in fold_reports if key in f]
        if vals:
            report["mean_" + key] = round(float(np.mean(vals)), 4)
            report["std_" + key] = round(float(np.std(vals)), 4)
    with open(REPORT_PATH, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print("Saved", REPORT_PATH)


if __name__ == "__main__":
    main()