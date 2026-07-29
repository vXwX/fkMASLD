# MASLD Drug Discovery Screening Pipeline

## Description
Complete virtual screening pipeline for MASLD (Metabolic dysfunction-Associated Steatotic Liver Disease) drug discovery. Takes a compound library (SDF), extracts SMILES, fetches PubChem bioactivity/toxicity data, computes physicochemical properties, ranks molecules by activity-toxicity profile, identifies related target proteins, detects binding pockets with D3Pockets, and scores protein-ligand interactions with DrugCLIP.

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
├── cross_analysis.py      # Cross-analysis: rank × DrugCLIP scores
├── top_candidates.csv     # All scored pairs (97056)
├── top_candidates_unique.csv  # Best score per molecule (5392)
└── top_candidates_confirmed_targets.csv  # Molecules with validated targets
```

## Pipeline Steps

### 0. Environment Setup
```bash
# tmp: general Python tasks (RDKit, tqdm, etc.)
conda activate tmp
# D3Pockets2J: GPU pocket detection
conda activate D3Pockets2J
# LigUnity: DrugCLIP scoring
conda activate LigUnity
```

### 1. Extract SMILES from SDF
```bash
conda run -n tmp python molecule/extract_smiles.py
# Output: molecule/T001_TargetMol_SMILES.csv (22688 molecules, deduplicated)
```

### 2. Fetch PubChem Data
```bash
conda run -n tmp python data/fetch_properties.py   # Physicochemical properties
conda run -n tmp python data/fetch_bioassay.py     # Bioactivity + toxicity
```

### 3. Activity-Toxicity Ranking
```bash
conda run -n tmp python rank/rank_molecules.py
# Output: rank/ranked_molecules.csv with 9-grid classification
```

### 4. Target Enrichment & PPI
```bash
conda run -n tmp python protein/enrichment/enrich_targets.py  # GO/KEGG/Reactome
conda run -n tmp python protein/relpro/find_related_proteins.py  # STRING PPI
```

### 5. Download Protein Structures
```bash
conda run -n tmp python protein/fetch_pdb_structures.py
# Output: protein/pdb_apo/{gene}_apo.pdb (18 proteins)
```

### 6. Pocket Detection (D3Pockets)
```bash
conda run -n D3Pockets2J python pocket/run_d3pockets.py
# Output: pocket/pocket_output/{gene}/pocket_detected/
```

### 7. DrugCLIP Scoring (GPU)
```bash
# Step 7a: Generate molecule conformations
conda run -n tmp python pocket/score/prepare_mols.py

# Step 7b: Run scoring per gene
CUDA_VISIBLE_DEVICES=0 conda run -n LigUnity python pocket/score/run_scoring.py

# Step 7c: Merge results
conda run -n tmp python pocket/score/merge_scores.py
```

### 8. Cross-Analysis
```bash
conda run -n tmp python cross_analysis.py
# Output: top_candidates.csv, top_candidates_unique.csv
```

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