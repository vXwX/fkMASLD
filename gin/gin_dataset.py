import torch
from torch.utils.data import Dataset


class MolGraphDataset(Dataset):
    """Dataset over a list of molecule records. Each record: graph dict + labels."""

    def __init__(self, records):
        self.records = records

    def __len__(self):
        return len(self.records)

    def __getitem__(self, idx):
        rec = self.records[idx]
        g = rec["graph"]
        # per-sample labels
        return {
            "x": torch.as_tensor(g["x"], dtype=torch.float32),
            "edge_index": torch.as_tensor(g["edge_index"], dtype=torch.long),
            "edge_attr": torch.as_tensor(g["edge_attr"], dtype=torch.float32),
            "has_active": torch.as_tensor(rec["has_active"], dtype=torch.float32) if rec["has_active"] is not None else torch.tensor(float("nan")),
            "act_ratio": torch.as_tensor(rec.get("act_ratio", float("nan")), dtype=torch.float32) if rec.get("act_ratio") is not None else torch.tensor(float("nan")),
            "has_tox": torch.as_tensor(rec["has_tox"], dtype=torch.float32) if rec["has_tox"] is not None else torch.tensor(float("nan")),
            "tox_level": torch.as_tensor(rec["tox_level"], dtype=torch.long) if rec["tox_level"] is not None else torch.tensor(-1),
        }


def collate_graphs(batch):
    """Collate variable-size graphs into a single batched graph with offsets."""
    n_graphs = len(batch)
    xs, rows, cols, es = [], [], [], []
    batch_idx = []
    has_act, act_ratio, has_tox, tox_level = [], [], [], []
    offset = 0
    for i, item in enumerate(batch):
        x = item["x"]
        ei = item["edge_index"]
        xs.append(x)
        rows.append(ei[0] + offset)
        cols.append(ei[1] + offset)
        es.append(item["edge_attr"])
        batch_idx.append(torch.full((x.size(0),), i, dtype=torch.long))
        has_act.append(item["has_active"].unsqueeze(0))
        act_ratio.append(item["act_ratio"].unsqueeze(0))
        has_tox.append(item["has_tox"].unsqueeze(0))
        tox_level.append(item["tox_level"].unsqueeze(0))
        offset += x.size(0)

    x = torch.cat(xs, dim=0)
    edge_index = torch.stack([torch.cat(rows), torch.cat(cols)], dim=0)
    edge_attr = torch.cat(es, dim=0)
    batch_idx = torch.cat(batch_idx, dim=0)
    has_act = torch.cat(has_act, dim=0)
    act_ratio = torch.cat(act_ratio, dim=0)
    has_tox = torch.cat(has_tox, dim=0)
    tox_level = torch.cat(tox_level, dim=0)
    return {
        "x": x,
        "edge_index": edge_index,
        "edge_attr": edge_attr,
        "batch_idx": batch_idx,
        "num_graphs": n_graphs,
        "has_active": has_act,
        "act_ratio": act_ratio,
        "has_tox": has_tox,
        "tox_level": tox_level,
    }