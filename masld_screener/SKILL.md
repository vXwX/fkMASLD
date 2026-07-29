---
name: masld-screener
description: >-
  从 SDF 文件筛选 MASLD（代谢功能障碍相关脂肪性肝病）降脂活性小分子。
  执行五阶段管线：去重标准化 → 理化初筛 → 相似性搜索 → ADMET 毒性预测 → 综合评分排序。
  当用户要求筛选降脂活性小分子、进行 MASLD 药物筛选、或对 SDF 分子库做降脂/低毒过滤时自动触发。
---

# masld-screener

从大规模 SDF 分子库中，自动筛选具有 **降脂活性潜力** 且 **低毒性** 的口服可及小分子，用于
MASLD（代谢功能障碍相关脂肪性肝病，旧称 NAFLD/NASH）先导化合物的优先排序。

## 触发条件

- 用户要求「筛选降脂活性小分子 / 降脂候选药物」
- 用户要求进行「MASLD 药物筛选 / NASH 先导化合物筛选」
- 用户希望对某个 `.sdf` 分子库做「降脂 + 低毒」过滤或优先排序
- 用户提供 SDF 文件并希望按已知临床阶段活性分子做相似性检索

## 依赖

- `rdkit`：分子读取、去盐、标准化、指纹与理化性质
- `admet-ai`：ADMET/毒性端点预测（**可选**；未安装时自动跳过 Phase 4，toxic 默认 0.5）
- `pandas`、`numpy`：数据处理
- `tqdm`：进度条

安装：`pip install -r requirements.txt`

## 输入

- 工作空间中的 `.sdf` 文件路径（可为厂商化合物库，如 TargetMol、Enamine、ChemDiv 等）

## 输出

- `ranked_candidates.csv`：按综合评分降序排列的候选分子，字段含
  `rank, name, smiles, MW, LogP, TPSA, HBD, HBA, RotB, Lipinski_violations,
   max_tc, best_match, best_target, physchem_compliance, toxic_score, final_score`
- `ranked_candidates.sdf`（可选，`--output_sdf`）：排序后的候选分子，属性字段写入评分

## 五阶段管线

1. **去重与标准化**：`SDMolSupplier` 读取 → 去盐（`SaltRemover`）→ 标准化（`rdMolStandardize.Normalize`）
   → 基于 canonical SMILES 去重。
2. **理化初筛**：计算 MW / LogP / TPSA / HBD / HBA / RotB / Lipinski 违规数，按 `config.PHYS_CHEM_FILTERS` 区间过滤。
3. **相似性搜索**：对 9 个临床阶段 MASLD 参考分子生成 Morgan 指纹（ECFP4, r=2, 2048 bit），
   计算每个候选与参考的 Tanimoto 相似度，保留 `max_tc >= SIMILARITY_THRESHOLD` 的分子并记录最佳匹配。
4. **ADMET 毒性预测**：用 `admet-ai` 预测 hERG / Ames / 肝毒性 / 口服生物利用度 / 水溶性，动态检测列名后按
   `config.TOXICITY_FILTERS` 过滤（未安装则跳过）。
5. **综合评分排序**：`final = 0.4×相似性 + 0.3×理化合规度 + 0.3×毒性安全度`，降序输出。

## 参考分子（见 reference_molecules.py）

Resmetirom (THR-β, 已批准)、VK2809 (THR-β, Ph2)、Obeticholic acid (FXR, 已批准)、
Cilofexor (FXR, Ph3)、Lanifibranor (PPARα/γ/δ, Ph3)、Firsocostat (ACC, Ph2)、
Denifanstat (FASN, Ph3)、BI-1467335 (AOC3/VAP-1, Ph2)、Aramchol (SCD1, Ph2)。

## 快速使用

```bash
cd masld_screener
python masld_screener.py /path/to/library.sdf -o ranked_candidates.csv --output_sdf
```

## 阈值调整

全部阈值集中于 `config.py`，修改即可调整严格程度，无需改动主脚本。
