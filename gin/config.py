"""Shared constants for the GIN multi-task pipeline.

These must stay in sync with gin/prepare_graphs.py feature construction:
  ATOM_DIM = MAX_ATOM_NUM + MAX_DEG + 7  (one-hot atom num + one-hot degree + 7 scalar feats)
  BOND_DIM = 6                            (single/double/triple/aromatic/conjugated/aromatic-bit)
"""

MAX_ATOM_NUM = 100
MAX_DEG = 6
ATOM_DIM = MAX_ATOM_NUM + MAX_DEG + 7
BOND_DIM = 6

HIDDEN = 128
N_LAYERS = 3
LR = 1e-3
PATIENCE = 8