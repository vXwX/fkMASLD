# MASLD Drug Discovery Screening Pipeline

## Usage
在 opencode 会话中直接说"按 MASLD 药物筛选流程处理"即可自动加载本 skill。

## Description
Complete virtual screening pipeline for MASLD (Metabolic dysfunction-Associated Steatotic Liver Disease) drug discovery. Takes a compound library (SDF), extracts SMILES, fetches PubChem bioactivity/toxicity data, computes physicochemical properties, ranks molecules by activity-toxicity profile, identifies related target proteins, detects binding pockets with D3Pockets, scores protein-ligand interactions with DrugCLIP, and uses a GIN multi-task model to predict activity × toxicity for the 8464 molecules lacking PubChem data, producing a final high-activity low-toxicity top10.

## Directory Structure
```
project/
├── molecule/              # Compound library
│   ├── extract_smiles.py  # SDF → SMILES CSV with dedup
│   └── T001_TargetMol_SMILES.csv  # 22868 molecules
├── data/
│   ├── bioassay/          # PubChem bioactivity/toxicity data
│   ├── property/          # Physicochemical properties
│   ├── supp/              # Missing molecule records
│   ├── fetch_bioassay.py  # Fetch activity/toxicity from PubChem
│   └── fetch_properties.py # Fetch properties from PubChem
├── rank/
│   ├── rank_molecules.py  # Activity × Toxicity 9-grid classification
│   ├── ranked_molecules.csv  # 14404 ranked molecules
│   ├── export_targets.py  # Target frequency summary
│   └── map_targets.py     # RefSeq → Gene Symbol mapping
├── protein/
│   ├── fetch_pdb_structures.py  # Download PDB + generate apo structures
│   ├── pdb_apo/           # 18 apo protein structures
│   ├── enrichment/        # GO/KEGG/Reactome enrichment
│   └── relpro/            # PPI network expansion via STRING
├── pocket/
│   ├── run_d3pockets.py   # D3Pockets cavity detection
│   ├── pocket_output/     # Detected pockets per protein
│   └── score/
│       ├── prepare_mols.py    # SMILES → 3D conformations → LMDB
│       ├── run_scoring.py     # DrugCLIP scoring (GPU)
│       ├── merge_scores.py    # Merge all gene scores
│       └── results/           # Per-gene and merged scores
├── masld_screener/        # Alternative MASLD screening pipeline
├── gin/                   # GIN multi-task model @ DrugClip
│   ├── config.py              # Shared constants (ATOM_DIM/BOND_DIM/HIDDEN/…)
│   ├── prepare_graphs.py      # SMILES → molecule graph + labels (activity/toxicity)
│   ├── gin_model.py           # MultiTaskGIN model
│   ├── gin_dataset.py         # Graph dataset
│   ├── multi_task_gin.py      # Train 4-fold GroupKFold (grouped by activity target)
│   ├── predict.py             # Predict P_act/STRENGTH/P_tox_0-2 for all molecules
│   ├── score_rank.py          # Rank no-data molecules by active×(1−tox)
│   └── results/
│       ├── predicted_inactive_ranked.csv  # 8464 no-data molecules ranked
│       ├── validation_report.json         # Activity AUC/Spearman, toxicity AUC
│       └── top10.csv                     # Final high-act low-tox Top10
├── cross_analysis.py      # Cross-analysis: rank × DrugCLIP scores
├── top_candidates.csv     # All scored pairs (97056)
├── top_candidates_unique.csv  # Best score per molecule (5392)
└── top_candidates_confirmed_targets.csv  # Molecules with validated targets
```

## Pipeline Steps

### 0. Environment Setup
```bash
# env: general Python tasks (RDKit, tqdm, etc.)
conda activate env
# D3Pockets: GPU pocket detection
conda activate D3Pockets
# DrugClip: DrugCLIP scoring + GIN multi-task training (PyTorch + CUDA + RDKit + sklearn)
conda activate DrugClip
```

### 1. Extract SMILES from SDF
```bash
conda run -n env python molecule/extract_smiles.py
```

### 2. Fetch PubChem Data
```bash
conda run -n env python data/fetch_properties.py
conda run -n env python data/fetch_bioassay.py
```

### 3. Activity-Toxicity Ranking
```bash
conda run -n env python rank/rank_molecules.py
```

### 4. Target Enrichment & PPI
```bash
conda run -n env python protein/enrichment/enrich_targets.py
conda run -n env python protein/relpro/find_related_proteins.py
```

### 5. Download Protein Structures
```bash
conda run -n env python protein/fetch_pdb_structures.py
```

### 6. Pocket Detection (D3Pockets)
```bash
conda run -n D3Pockets python pocket/run_d3pockets.py
```

### 7. DrugCLIP Scoring (GPU)
```bash
conda run -n env python pocket/score/prepare_mols.py
CUDA_VISIBLE_DEVICES=0 conda run -n DrugClip python pocket/score/run_scoring.py
conda run -n env python pocket/score/merge_scores.py
```

### 8. Cross-Analysis
```bash
conda run -n env python cross_analysis.py
```

### 9. GIN Prediction for No-Data Molecules (DrugClip)
For the 8,464 molecules without PubChem activity data, train a multi-task GIN
(activity: P_act + STRENGTH; toxicity: P_tox_0/1/2) with 4-fold GroupKFold grouped by
activity target to prevent leakage, then rank them.

```bash
conda run -n env python gin/prepare_graphs.py
conda run -n DrugClip python gin/multi_task_gin.py
conda run -n DrugClip python gin/predict.py
conda run -n DrugClip python gin/score_rank.py
```

- `gin/results/predicted_inactive_ranked.csv` — 8,464 ranked (score = P_act × STRENGTH × (1 − tox))
- `gin/results/validation_report.json` — 4-fold activity AUC/Spearman + toxicity AUC
- `gin/results/top10.csv` — final Top10 (active_score ≥ 90分位 AND toxicity_score ≤ 10分位), sorted by recommend_score

## Key Classification Rules

### Activity
| Level | Condition |
|-------|-----------|
| 高活性 | Active ratio > 30% AND ≥ 2 targets |
| 活性 | Active ratio > 5% OR ≥ 1 target |
| 低活性 | Rest |

### Toxicity
| Level | Condition |
|-------|-----------|
| 高毒性 | Contains high-risk keywords (hERG, hepatotox, etc.) |
| 毒性 | Any Active record in toxicity assays |
| 低毒性 | All Inactive or no data |

## Key Scripts Location
- `cross_analysis.py` — project root
- `export_confirmed_targets.py` — project root
- `rank/rank_molecules.py` — rank directory
- `rank/export_targets.py` — rank directory
- `rank/map_targets.py` — rank directory
- `protein/fetch_pdb_structures.py` — protein directory
- `pocket/run_d3pockets.py` — pocket directory
- `pocket/score/prepare_mols.py` — pocket/score directory
- `pocket/score/run_scoring.py` — pocket/score directory
- `pocket/score/merge_scores.py` — pocket/score directory
- `gin/config.py` — gin directory
- `gin/prepare_graphs.py` — gin directory
- `gin/multi_task_gin.py` — gin directory
- `gin/predict.py` — gin directory
- `gin/score_rank.py` — gin directory