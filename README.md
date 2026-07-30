# MASLD Drug Discovery Screening Pipeline

MASLD (Metabolic dysfunction-Associated Steatotic Liver Disease) 药物发现虚拟筛选流程。从 22966 个化合物的 SDF 库出发，经过 SMILES 提取、PubChem 活性/毒性数据获取、理化性质计算、活性-毒性分层排序、靶点分析、蛋白结构准备、口袋检测、DrugCLIP 分子对接打分，最终筛选出 Top 10 候选分子。

## 目录结构

```
├── molecule/                    # 化合物库
│   ├── extract_smiles.py        # SDF → SMILES 提取（去重）
│   └── T001_TargetMol_SMILES.csv  # 22868 个分子
├── data/
│   ├── bioassay/                # PubChem 活性/毒性数据
│   │   ├── activity.csv         # 523 万条活性数据
│   │   └── toxicity.csv         # 78 万条毒性数据
│   ├── property/                # 理化性质（22868 分子）
│   ├── fetch_bioassay.py        # PubChem 活性/毒性获取
│   └── fetch_properties.py      # PubChem 理化性质获取
├── rank/
│   ├── rank_molecules.py        # 活性×毒性 9宫格分层排序
│   ├── ranked_molecules.csv     # 14404 分子排序结果
│   ├── export_targets.py        # 靶点频次统计
│   └── map_targets.py           # RefSeq → 基因名映射
├── protein/
│   ├── pdb_apo/                 # 18 个空蛋白 PDB 结构
│   ├── enrichment/              # GO/KEGG/Reactome 富集分析
│   ├── relpro/                  # STRING PPI 网络扩展
│   └── fetch_pdb_structures.py  # PDB 下载+空蛋白生成
├── pocket/
│   ├── run_d3pockets.py         # D3Pockets 空腔检测
│   ├── pocket_output/           # 各蛋白检测到的口袋
│   └── score/
│       ├── prepare_mols.py      # SMILES → 3D 构象 → LMDB
│       ├── run_scoring.py       # DrugCLIP 打分 (GPU)
│       ├── merge_scores.py      # 合并所有基因打分
│       └── results/             # 409608 条打分记录
├── masld_screener/              # 备选 MASLD 筛选管道
│   └── T001_ranked.csv          # 81 个候选分子
├── top_candidates.csv           # 97056 分子-蛋白对
├── top_candidates_unique.csv    # 5392 个唯一分子
├── top_candidates_confirmed_targets.csv  # 64 个实验验证靶点分子
├── top_candidates_vs_masld.csv  # 两条管道交集（16 个）
├── top_10_final.csv             # 最终 Top 10 候选
├── classification_criteria.csv  # 活性/毒性分层标准
├── cross_analysis.py            # 交叉分析
└── pipeline.tex                 # 流程图（Overleaf 可用）
```

## 分层标准

| 层级 | 活性条件 | 毒性条件 |
|------|----------|----------|
| 高 | Positive 占比 > 30% 且 ≥ 2 个靶点 | 含高危关键词 (hERG/肝毒性/致癌等) |
| 中 | Positive 占比 > 5% 或 ≥ 1 个靶点 | 任意 Positive 记录 |
| 低 | 全部 Negative 或无数据 | 全部 Negative 或无数据 |

## 运行流程

1. **分子准备**: `conda run -n env python molecule/extract_smiles.py`
2. **数据获取**: `conda run -n env python data/fetch_properties.py` + `fetch_bioassay.py`
3. **分层排序**: `conda run -n env python rank/rank_molecules.py`
4. **蛋白结构**: `conda run -n env python protein/fetch_pdb_structures.py`
5. **口袋检测**: `conda run -n D3Pockets python pocket/run_d3pockets.py`
6. **DrugCLIP**: `conda run -n env python pocket/score/prepare_mols.py` → `CUDA_VISIBLE_DEVICES=0 conda run -n DrugClip python pocket/score/run_scoring.py`
7. **交叉分析**: `conda run -n env python cross_analysis.py`

## 环境要求

- **env**: RDKit, tqdm, lmdb, mygene, numpy, pandas, biopandas
- **D3Pockets**: D3Pockets
- **DrugClip**: DrugCLIP (PyTorch, CUDA)

## 最终结果

`top_10_final.csv` 包含 10 个候选分子，在 DrugCLIP 打分和 MASLD 筛选中综合排名最高，且排除毒性分子。详见 [PIPELINE.md](PIPELINE.md) 和 [SKILL.md](SKILL.md)。
