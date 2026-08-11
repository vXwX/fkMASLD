# MASLD Drug Discovery Screening Pipeline

MASLD (Metabolic dysfunction-Associated Steatotic Liver Disease) 药物发现虚拟筛选流程。从 22966 个化合物的 SDF 库出发，经过 SMILES 提取、PubChem 活性/毒性数据获取、理化性质计算、活性-毒性分层排序、靶点分析、蛋白结构准备、口袋检测、DrugCLIP 分子对接打分，最终筛选出 Top 10 候选分子。

## 目录结构

```
├── molecule/           # 化合物库 (22868 SMILES)
├── data/               # PubChem 数据 (活性/毒性/理化性质)
├── rank/               # 活性×毒性 分层排序
├── protein/            # 靶点分析 + 18 空蛋白 PDB
├── pocket/             # D3Pockets 口袋检测 + DrugCLIP 打分
├── gin/                # GIN 多任务模型: 对无数据化合物预测活性×毒性并排序
│   ├── results/predicted_inactive_ranked.csv  # 8464 无数据分子排序结果
├── Top10.csv    # 最终 Top 10 候选分子
└── SKILL.md            # opencode Skill 定义
```

## 分层标准

| 层级 | 活性条件 | 毒性条件 |
|------|----------|----------|
| 高 | Positive 占比 > 30% 且 ≥ 2 个靶点 | 含高危关键词 (hERG/肝毒性/致癌等) |
| 中 | Positive 占比 > 5% 或 ≥ 1 个靶点 | 任意 Positive 记录 |
| 低 | 全部 Negative 或无数据 | 全部 Negative 或无数据 |

## 使用 Agent

在 AI agent（如 opencode）会话中直接说"按 MASLD 药物筛选流程处理"即可自动加载本流程。
`SKILL.md` 定义了完整的 pipeline 描述，需将 `SKILL.md` 放置在 agent 可识别的 skills 目录下
（如 `.opencode/skills/` 或 `~/.config/opencode/skill/`）。

## 手动运行流程

1. **分子准备**: `conda run -n env python molecule/extract_smiles.py`
2. **数据获取**: `conda run -n env python data/fetch_properties.py` + `fetch_bioassay.py`
3. **分层排序**: `conda run -n env python rank/rank_molecules.py`
4. **蛋白结构**: `conda run -n env python protein/fetch_pdb_structures.py`
5. **口袋检测**: `conda run -n D3Pockets python pocket/run_d3pockets.py`
6. **DrugCLIP**: `conda run -n env python pocket/score/prepare_mols.py` → `CUDA_VISIBLE_DEVICES=0 conda run -n DrugClip python pocket/score/run_scoring.py`
7. **交叉分析**: `conda run -n env python cross_analysis.py`

## GIN 无数据化合物排序流程 (LigUnity 环境)

对 8,464 个无活性/毒性数据的化合物, 用 GIN 多任务模型预测活性分与毒性分并排序。

```bash
# 1. 构建分子图 + 标签
conda run -n LigUnity python gin/prepare_graphs.py
# 2. 训练多任务 GIN (4折GroupKFold, 按活性靶点分组防泄漏)
conda run -n LigUnity python gin/multi_task_gin.py
# 3. 对所有分子打分
conda run -n LigUnity python gin/predict.py
# 4. 生成 8,464 无数据分子排序
conda run -n LigUnity python gin/score_rank.py
```

输出: `gin/results/predicted_inactive_ranked.csv` (排序分 = 活性分 × (1 − 毒性分))
验证报告: `gin/results/validation_report.json` (活性 AUC/Spearman、毒性 AUC)

## 环境要求

- **env**: RDKit, tqdm, lmdb, mygene, numpy, pandas, biopandas
- **D3Pockets**: D3Pockets
- **DrugClip**: DrugCLIP (PyTorch, CUDA)
- **LigUnity**: GIN 模型训练 (PyTorch + CUDA + RDKit + sklearn)

## 最终结果

`top_10_final.csv` 包含 10 个候选分子，在 DrugCLIP 打分和 MASLD 筛选中综合排名最高，且排除毒性分子。详见 [PIPELINE.md](PIPELINE.md) 和 [SKILL.md](SKILL.md)。
