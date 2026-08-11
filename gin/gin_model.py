import torch
import torch.nn as nn


def graph_pool(node_emb, batch_idx, num_graphs):
    """Mean + max pooling per graph. batch_idx: graph index per node."""
    dev = node_emb.device
    mean = torch.zeros(num_graphs, node_emb.size(1), device=dev)
    cnt = torch.zeros(num_graphs, 1, device=dev)
    ones = torch.ones_like(batch_idx, dtype=node_emb.dtype)
    mean.index_add_(0, batch_idx, node_emb)
    cnt.index_add_(0, batch_idx, ones.unsqueeze(1))
    mean = mean / cnt.clamp(min=1)

    mx = torch.full((num_graphs, node_emb.size(1)), float("-inf"), device=dev)
    mx = mx.scatter_reduce_(0, batch_idx.unsqueeze(1).expand_as(node_emb), node_emb, reduce="amax")
    mx = torch.nan_to_num(mx, nan=0.0, posinf=0.0, neginf=0.0)
    return torch.cat([mean, mx], dim=1)


class GINConv(nn.Module):
    """Graph Isomorphism Network convolution (sum aggregation + MLP).

    Node features are projected to `out_dim`, edge features are projected and
    added to the source node's message, then sum-aggregated and passed through
    (BN -> ReLU -> Linear). Structure:
        x' = mlp2( relu( bn( (1+eps)*proj(x) + sum_{j->i} (proj(x_j) + edge(x_ij)) ) ) )
    """

    def __init__(self, in_dim, out_dim, edge_dim):
        super().__init__()
        self.eps = nn.Parameter(torch.zeros(1))
        self.node_proj = nn.Linear(in_dim, out_dim)
        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_dim, out_dim),
            nn.ReLU(),
        )
        self.bn = nn.BatchNorm1d(out_dim)
        self.act = nn.ReLU()
        self.mlp2 = nn.Linear(out_dim, out_dim)

    def forward(self, x, edge_index, edge_attr):
        row, col = edge_index
        xp = self.node_proj(x)
        edge_msgs = self.edge_mlp(edge_attr)
        agg = torch.zeros_like(xp)
        agg.index_add_(0, row, xp[col] + edge_msgs)
        h = (1.0 + self.eps) * xp + agg
        h = self.bn(h)
        return self.mlp2(self.act(h))


class MultiTaskGIN(nn.Module):
    def __init__(self, in_dim, bond_dim, hidden=128, n_layers=3, n_tox=3):
        super().__init__()
        self.gin_layers = nn.ModuleList()
        prev = in_dim
        for _ in range(n_layers):
            self.gin_layers.append(GINConv(prev, hidden, bond_dim))
            prev = hidden
        self.dropout = nn.Dropout(p=0.2)
        self.pool_lin = nn.Linear(hidden * 2, hidden)
        # heads
        self.act_clf = nn.Sequential(nn.Linear(hidden, 64), nn.ReLU(), nn.Dropout(0.1), nn.Linear(64, 1))
        self.act_reg = nn.Sequential(nn.Linear(hidden, 64), nn.ReLU(), nn.Dropout(0.1), nn.Linear(64, 1))
        self.tox_head = nn.Sequential(nn.Linear(hidden, 64), nn.ReLU(), nn.Dropout(0.1), nn.Linear(64, n_tox))

    def embed(self, x, edge_index, edge_attr, batch_idx, num_graphs):
        h = x
        for layer in self.gin_layers:
            h = layer(h, edge_index, edge_attr)
            h = self.dropout(h)
        h = graph_pool(h, batch_idx, num_graphs)
        return self.pool_lin(h)

    def forward(self, x, edge_index, edge_attr, batch_idx, num_graphs):
        h = self.embed(x, edge_index, edge_attr, batch_idx, num_graphs)
        p_act = torch.sigmoid(self.act_clf(h)).squeeze(-1)
        strength = torch.sigmoid(self.act_reg(h)).squeeze(-1)
        tox_logits = self.tox_head(h)
        return p_act, strength, tox_logits

    def predict(self, x, edge_index, edge_attr, batch_idx, num_graphs):
        """Inference: return (p_act, strength, tox_probs)."""
        self.eval()
        with torch.no_grad():
            h = self.embed(x, edge_index, edge_attr, batch_idx, num_graphs)
            p_act = torch.sigmoid(self.act_clf(h)).squeeze(-1)
            strength = torch.sigmoid(self.act_reg(h)).squeeze(-1)
            tox_probs = torch.softmax(self.tox_head(h), dim=-1)
        return p_act, strength, tox_probs
