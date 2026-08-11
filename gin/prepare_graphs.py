import os
import pickle

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from config import ATOM_DIM, BOND_DIM, MAX_ATOM_NUM, MAX_DEG

RDLogger.DisableLog("rdApp.*")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
RESULT_DIR = os.path.join(SCRIPT_DIR, "results")
os.makedirs(RESULT_DIR, exist_ok=True)

MOLECULE_CSV = os.path.join(PROJECT_DIR, "molecule", "T001_TargetMol_SMILES.csv")
ACTIVITY_CSV = os.path.join(PROJECT_DIR, "data", "activity.csv")
TOXICITY_CSV = os.path.join(PROJECT_DIR, "data", "toxicity.csv")

OUT_MOL = os.path.join(RESULT_DIR, "mols.csv")
OUT_GRAPH = os.path.join(RESULT_DIR, "graphs.pkl")
OUT_LABELS = os.path.join(RESULT_DIR, "labels.pkl")

HIGH_TOX_KEYWORDS = {
    "herg", "qt prolongation", "cardiotox", "hepatotox", "nephrotox",
    "neurotox", "genotox", "death", "lethal", "carcinogen", "mutagen",
    "teratogen",
}


def build_labels():
    """Aggregate activity & toxicity CSVs into per-molecule labels."""
    import pandas as pd

    act_count = {}
    act_total = {}
    act_targets = {}
    for chunk in pd.read_csv(ACTIVITY_CSV, usecols=["ID", "Activity_Outcome", "Target_GeneID"],
                             chunksize=400000):
        for mol_id, g in chunk.groupby("ID"):
            act_total[mol_id] = act_total.get(mol_id, 0) + len(g)
            act_count[mol_id] = act_count.get(mol_id, 0) + int(
                (g["Activity_Outcome"].str.lower() == "active").sum()
            )
            act_targets.setdefault(mol_id, set())
            tgt = g["Target_GeneID"].dropna().astype(str).unique().tolist()
            act_targets[mol_id].update(t for t in tgt if t.strip())

    tox_total = {}
    tox_active = {}
    tox_highrisk = {}
    for chunk in pd.read_csv(TOXICITY_CSV, usecols=["ID", "Activity_Outcome", "Assay_Name"],
                             chunksize=400000):
        for mol_id, g in chunk.groupby("ID"):
            tox_total[mol_id] = tox_total.get(mol_id, 0) + len(g)
            tox_active[mol_id] = tox_active.get(mol_id, 0) + int(
                (g["Activity_Outcome"].str.lower() == "active").sum()
            )
            names = " ".join(g["Assay_Name"].fillna("").astype(str).str.lower().tolist())
            if any(kw in names for kw in HIGH_TOX_KEYWORDS):
                tox_highrisk[mol_id] = True

    labels = {}
    for mol_id in set(act_count) | set(tox_total):
        ac = act_count[mol_id]
        at = act_total[mol_id]
        tt = tox_total.get(mol_id, 0)
        tva = tox_active.get(mol_id, 0)
        thr = tox_highrisk.get(mol_id, False)
        rat = ac / at if at else float("nan")
        labels[mol_id] = {
            "has_active": int(ac > 0),
            "act_ratio": rat,
            "active_targets": ";".join(sorted(act_targets.get(mol_id, ()))),
            "has_tox": int(tva > 0) if tt else None,
            "tox_level": int(2 if thr else (1 if tva > 0 else 0)) if tt else None,
        }
    return labels


def atom_to_features(atom):
    anum = atom.GetAtomicNum()
    if anum >= MAX_ATOM_NUM:
        anum = 0
    v = [0.0] * ATOM_DIM
    v[anum] = 1.0
    deg = min(atom.GetTotalDegree(), MAX_DEG - 1)
    v[MAX_ATOM_NUM + deg] = 1.0
    base = MAX_ATOM_NUM + MAX_DEG
    v[base + 0] = float(atom.GetFormalCharge())
    v[base + 1] = float(atom.GetTotalNumHs()) / 4.0
    v[base + 2] = float(atom.GetIsAromatic())
    v[base + 3] = float(atom.IsInRing())
    v[base + 4] = float(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
    return (v, base)


def bond_to_features(bond):
    v = [0.0] * BOND_DIM
    bt = bond.GetBondType()
    if bt == Chem.BondType.SINGLE:
        v[0] = 1.0
    elif bt == Chem.BondType.DOUBLE:
        v[1] = 1.0
    elif bt == Chem.BondType.TRIPLE:
        v[2] = 1.0
    elif bt == Chem.BondType.AROMATIC:
        v[3] = 1.0
    v[4] = float(bond.GetIsConjugated())
    v[5] = float(bond.GetIsAromatic())
    return v


def mol_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    xs = []
    for a in mol.GetAtoms():
        vec, base = atom_to_features(a)
        xs.append(vec)
    x = np.array(xs, dtype=np.float32)

    rows, cols, es = [], [], []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        rows.append(i)
        cols.append(j)
        rows.append(j)
        cols.append(i)
        es.append(bond_to_features(bond))
        es.append(bond_to_features(bond))
    edge_index = np.stack([rows, cols], axis=0) if rows else np.empty((2, 0), dtype=np.int64)
    edge_attr = np.array(es, dtype=np.float32) if es else np.empty((0, BOND_DIM), dtype=np.float32)
    return {
        "x": x,
        "edge_index": edge_index,
        "edge_attr": edge_attr,
    }


def main():
    labels = build_labels()
    print(f"[prepare] labels built for {len(labels)} molecules")

    mol_df = pd.read_csv(MOLECULE_CSV)
    graphs = {}
    parsed, failed = 0, 0
    missing_rows = []

    records = []
    for idx, row in mol_df.iterrows():
        m_id = row["ID"]
        lb = labels.get(m_id, {})
        g = mol_to_graph(row["SMILES"])
        if g is None:
            failed += 1
            missing_rows.append(m_id)
            continue
        parsed += 1
        graphs[m_id] = g
        records.append({
            "ID": m_id,
            "SMILES": row["SMILES"],
            "has_active": lb.get("has_active"),
            "act_ratio": lb.get("act_ratio", float("nan")),
            "has_tox": lb.get("has_tox"),
            "tox_level": lb.get("tox_level"),
            "active_targets": lb.get("active_targets", ""),
        })

    mol_out = pd.DataFrame(records)
    mol_out.to_csv(OUT_MOL, index=False)
    with open(OUT_GRAPH, "wb") as f:
        pickle.dump(graphs, f)
    with open(OUT_LABELS, "wb") as f:
        pickle.dump(labels, f)

    print(f"[prepare] parsed {parsed}, failed {failed}")
    if missing_rows:
        print("failed sample:", missing_rows[:10])
    print(f"[prepare] saved {OUT_MOL}")
    print(f"[prepare] saved graphs for {len(graphs)} molecules")
    if not graphs:
        print("WARNING: no graphs parsed!")
        return
    sample_key = next(iter(graphs))
    print("[prepare] sample keys:", list(graphs[sample_key].keys()),
          "x shape:", graphs[sample_key]["x"].shape,
          "edges:", graphs[sample_key]["edge_index"].shape)


if __name__ == "__main__":
    main()