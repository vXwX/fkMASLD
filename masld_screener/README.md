# masld_screener — MASLD 降脂活性小分子筛选 Skill

从大规模 SDF 分子库中，自动筛选具有 **降脂活性潜力** 且 **低毒性** 的口服可及小分子，
用于 MASLD（代谢功能障碍相关脂肪性肝病，旧称 NAFLD/NASH）先导化合物的优先排序。

---

## 1. 功能简介

`masld_screener` 实现了一条完整、可复现的虚拟筛选管线，包含五个阶段：

| 阶段 | 功能 | 关键阈值来源 |
|------|------|--------------|
| 1 | SDF 读取、去盐、标准化、基于 canonical SMILES 去重 | `config.DEDUP_METHOD` |
| 2 | 理化性质计算（MW/LogP/TPSA/HBD/HBA/RotB/Lipinski）与区间初筛 | `config.PHYS_CHEM_FILTERS` |
| 3 | 与 9 个临床阶段 MASLD 参考分子的 Morgan 指纹 Tanimoto 相似性检索 | `config.SIMILARITY_THRESHOLD` 等 |
| 4 | ADMET 毒性预测（hERG / Ames / 肝毒 / 口服利用度 / 水溶性） | `config.TOXICITY_FILTERS` |
| 5 | 综合评分（相似性 0.4 + 理化 0.3 + 毒性 0.3）并降序排序 | `config.SCORE_WEIGHTS` |

最终输出按综合评分降序排列的候选分子清单，便于优先进行实验验证。

---

## 2. 安装依赖

要求 **Python >= 3.8**。建议在独立虚拟环境中安装：

```bash
# 建议新建虚拟环境（任选其一）
python -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate

# 安装依赖
pip install -r requirements.txt
```

> `admet-ai` 为可选项：若未安装，脚本会在 Phase 4 打印警告并自动跳过毒性预测，
> 将 `toxic_score` 设为中性值 0.5，其余阶段不受影响。

`requirements.txt` 内容：

```
rdkit
admet-ai
pandas
numpy
tqdm
```

---

## 3. 使用方法（命令行示例）

```bash
# 基本用法：输入 SDF，输出 ranked_candidates.csv
python masld_screener.py library.sdf

# 指定输出文件名，并同时导出排序后的 SDF
python masld_screener.py library.sdf -o result.csv --output_sdf

# 提高相似性门槛（保留与参考分子更相似的候选）
python masld_screener.py library.sdf --sim_threshold 0.35
```

参数说明：

- `input_sdf`：必填，输入 SDF 文件路径。
- `-o / --output`：输出 CSV 路径，默认 `ranked_candidates.csv`。
- `--output_sdf`：开关，导出排序后的候选分子 SDF（与 CSV 同名、扩展名 `.sdf`）。
- `--sim_threshold`：覆盖 `config.SIMILARITY_THRESHOLD`，用于临时调整相似性门槛。

---

## 4. 各阶段筛选逻辑说明

### Phase 1 — 去重与标准化
- 使用 RDKit `SDMolSupplier` 读取 SDF（保留原始名称/ID 属性）。
- `SaltRemover` 去除常见盐与溶剂片段（`dontRemoveEverything=True` 防止删空）。
- `rdMolStandardize.Normalize` 统一电荷/互变异构/价键表示。
- 以 `Chem.MolToSmiles(mol, canonical=True)` 生成规范 SMILES，用 `set` 去重。

### Phase 2 — 理化初筛
- 计算：MW（`Descriptors.MolWt`）、LogP（`Crippen.MolLogP`）、TPSA（`Descriptors.TPSA`）、
  HBD/HBA（`Lipinski.NumHDonors/NumHAcceptors`）、RotB（`Lipinski.NumRotatableBonds`）、
  Lipinski 违规数（MW>500 / LogP>5 / HBD>5 / HBA>10 的计数）。
- 逐项检查是否落入 `PHYS_CHEM_FILTERS` 的 `(min, max)` 区间，全部通过才保留。
- 同时记录「通过项数 / 总项数」用于 Phase 5 的理化合规度评分。

### Phase 3 — 相似性搜索
- 对 9 个参考分子生成 Morgan 指纹（`radius=2`, `nBits=2048`，即 ECFP4）。
- 对每个候选计算与全部参考的 Tanimoto 相似度，取最大值 `max_tc`。
- 保留 `max_tc >= SIMILARITY_THRESHOLD` 的分子，并记录最佳匹配分子名与靶点。

### Phase 4 — ADMET 毒性预测
- 使用 `admet-ai` 的 `ADMETModel.predict` 批量预测。
- **动态检测列名**：打印 `df.columns` 供调试，按关键字（herg / ames / hepato / oral / solubility 等）
  自动匹配端点列。
- 依据 `TOXICITY_FILTERS` 判定：布尔型（`False`=期望安全）、阈值型（`">0.3"` / `">-5"`）。
- 无法映射的端点不据此淘汰分子；毒性安全度 = 通过项数 / 可映射项数。
- `admet-ai` 缺失或预测失败时，自动跳过，`toxic_score` 默认 0.5。

### Phase 5 — 综合评分与排序
- 相似性得分 = `min(max_tc / 0.7, 1.0)`
- 理化合规度 = `通过项数 / 总项数`
- 毒性安全度 = `通过毒性检查项数 / 总检查项数`
- 综合评分 = `0.4×相似性 + 0.3×理化合规度 + 0.3×毒性安全度`
- 按综合评分降序排列，输出 Top 10 摘要。

---

## 5. 输出文件说明

### ranked_candidates.csv
列：`rank, name, smiles, MW, LogP, TPSA, HBD, HBA, RotB, Lipinski_violations,
max_tc, best_match, best_target, physchem_compliance, toxic_score, final_score`

- `max_tc`：与参考分子的最大 Tanimoto 相似度
- `best_match` / `best_target`：最相似参考分子名 / 靶点
- `physchem_compliance`：理化合规度（0–1）
- `toxic_score`：毒性安全度（0–1；未跑 Phase 4 时为 0.5）
- `final_score`：综合评分（0–1，越高越优先）

### ranked_candidates.sdf（可选）
排序后的候选分子，SDF 属性字段写入 `rank / final_score / max_tc / best_match /
best_target / physchem_compliance / toxic_score`，便于直接在分子可视化软件中查看。

---

## 6. 如何调整阈值

所有阈值集中在 **`config.py`**，修改后无需改动 `masld_screener.py`：

- 放宽/收紧理化范围：编辑 `PHYS_CHEM_FILTERS`（如将 `MW` 上限从 600 调到 700）。
- 调整相似性门槛：编辑 `SIMILARITY_THRESHOLD`（值越小保留越多）。
- 修改指纹参数：`FINGERPRINT_RADIUS` / `FINGERPRINT_BITS`。
- 调整毒性判定：`TOXICITY_FILTERS`（端点与条件）。
- 调整评分权重：`SCORE_WEIGHTS`（三者之和建议为 1.0）。

如需更换/增补参考分子，编辑 `reference_molecules.py`（每个分子为
`{"smiles": "...", "target": "...", "stage": "..."}` 字典）。

---

## 7. 可移植性说明

- 整个 `masld_screener/` 文件夹为**自包含单元**，可直接复制到任何装有 Python (>=3.8)
  的机器使用。
- 所有路径均使用相对路径或 `os.path` 处理，**不依赖任何工作空间绝对路径**。
- 依赖完整列于 `requirements.txt`，脚本头部注明 Python 版本要求。
- 参考分子 SMILES 已内嵌于 `reference_molecules.py`，无需联网即可运行（Phase 4 的
  `admet-ai` 在联网安装后本地推理，亦无需实时联网）。
- `admet-ai` 为可选依赖，缺失时管线其余部分完全可运行。

### 已知注意事项
- **非 ASCII 文件名兼容**：RDKit 的 C++ `SDMolSupplier` 在部分平台（如 Windows）对含
  中文、空格等非 ASCII 字符的文件名支持不佳，会报 `OSError: Bad input file`。
  本工具已内置自动回退——当直接打开失败时，会临时复制为 ASCII 命名副本再读取，
  运行结束后自动清理，**无需用户干预**。如希望进一步提速，可将输入 SDF 预先命名为
  纯英文/数字文件名（如 `library.sdf`）。
- **大文件**：对百万级分子的超大型 SDF，建议在内存充足的环境下运行；Phase 4 的
  `admet-ai` 需加载本地模型，耗时随候选数量线性增长。
