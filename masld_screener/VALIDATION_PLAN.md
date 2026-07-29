# 实验值交叉验证方案（VALIDATION_PLAN）

本 Skill 的筛选结果（理化计算 + 相似度代理 + admet-ai 预测）属于**纯虚拟筛选**，
源库 SDF 本身不含任何实验活性/毒性字段。为回答"预测到底可不可信"，
下面给出可行的实验值来源与验证路径。

## 1. 关键认知

- **admet-ai 输出的是 ML 预测值**，不是实验测定值。
- **理化性质 / 结构相似度是 RDKit 的确定性计算**，不是模型猜测。
- **"降脂活性"是相似度代理（类比假设）**，并非对候选自身靶点的直接活性预测。
- 整条流程无任何湿实验数据；结果仅用于**优先级排序**，最终需实验确认。

## 2. 公开实验数据库（可用于交叉验证）

### 2.1 实验 ADMET（毒性/药代）—— 直接对标 admet-ai 的预测端点

| 数据库 | 网址 | 实验内容 | 用途 |
|--------|------|---------|------|
| admetSAR 3.0 | lmmd.ecust.edu.cn/admetsar3 | 37 万条实验 ADMET，覆盖 10.4 万化合物；含 hERG / AMES / DILI / 溶解度 / 口服生物利用度 | 最对口；可下载数据集批量比对 |
| ADMETNet | bioinf.xmu.edu.cn/ADMETNet | 逐化合物实验页（溶解度、生物利用度、AMES、肝毒等） | 单分子点查最快 |
| PKKB | cadd.suda.edu.cn/admet | 1685 个药物的实验 ADMET，11 个子集可免费下载 | 模型基准（AUROC/PR）干净集 |
| ADMETlab 3.0 | rdnase.scbdd.com | 40 万+ 条目，支持 ADME/T 相似度 read-across | 对无实验值者用近邻外推 |

### 2.2 实验生物活性（IC50/EC50/Kd/Ki）—— 验证"相似度→降脂活性"代理

| 数据库 | 网址 | 内容 |
|--------|------|------|
| ChEMBL | ebi.ac.uk/chembl | 已发表文献的 IC50/EC50/Kd/Ki，按靶点（PPARα/γ/δ、FXR、THR-β 等）组织 |
| BindingDB | bindingdb.org | 蛋白-配体实验结合亲和力 |
| PubChem BioAssay | pubchem.ncbi.nlm.nih.gov | 上百万条实验 Assay 结果（PUG-REST 可程序化访问） |
| Guide to Pharmacology | guidetopharmacology.org | 策展级配体活性表（参考药实测活性确认） |

## 3. 自动验证脚本（已提供：`validate_predictions.py`）

脚本对候选复算 admet-ai 逐端点概率，并调用 **PubChem PUG-REST** 拉取实验 BioAssay 结果做比对。

- **hERG（自动、快速）**：`assaysummary?gene=KCNH2` 基因过滤 + Assay 名含 "hERG" 的
  Active/Inactive 结局。`预测=概率≥0.5 判有毒`，`实验=Active 判有毒`。
- **AMES / DILI（可选 `--deep`，较慢）**：扫描完整 assaysummary（部分常见药摘要很大可能超时），
  按 Ames/Salmonella/mutagen 与 hepatotox/liver/dili 关键词取结局。
- **溶解度 / 口服生物利用度**：PubChem BioAssay 一般不含连续实验值，脚本保留预测列，
  建议到 admetSAR3 / ADMETNet / PKKB 做人工或批量核对。

用法：
```bash
python validate_predictions.py -i T001_ranked.csv -o validation_report.csv
python validate_predictions.py -i T001_ranked.csv -o validation_report.csv --deep --limit 20
```

## 4. 已知局限

- PubChem 实验覆盖对**常见药/内源物高**、对全新骨架化合物可能很低。
- 一个化合物在 PubChem 中可能有多个 hERG Assay，脚本取首个 Active/Inactive 作为代表，
  不同 Assay 浓度/体系可能给出不同结论（以 "Unspecified/Inconclusive" 居多的不计入比对）。
- `--deep` 模式对大摘要分子会超时，相应端点标记为无数据；这不是预测错误，而是覆盖不足。
- 完整、严格的置信度（AUROC/PR）建议用 admetSAR3 下载集或内部 Assay 标注数据重算。

## 5. 推荐落地步骤（给对接方）

1. 用 `validate_predictions.py` 得到 hERG 自动一致率（已包含在分享包 `validation_report.csv`）。
2. 对 Top 候选的 CAS/SMILES，到 **admetSAR3 / ADMETNet** 拉实验 hERG/AMES/DILI/溶解度/生物利用度，
   与 `validation_report.csv` 预测列人工或批量比对，算一致率。
3. 用 **ChEMBL** 查参考药（Resmetirom/OCA/Lanifibranor…）对其靶点的实测 IC50，确认锚点活性真实；
   并查候选是否已有对应靶点实测活性。
4. 最终以**湿实验（体外 Assay / 体内）**确认 Top 候选。

## 6. hERG 自动验证结果（已执行，2026-07-28）

命令：`validate_predictions.py -i T001_ranked.csv -o validation_report.csv`（默认 hERG 模式，非 `--deep`）

| 端点 | 实验覆盖 | 一致率 | 特异度 | 灵敏度 | TP/FP/TN/FN |
|------|---------|--------|--------|--------|-------------|
| hERG | 10/81 | 90.0% | 90.0% | 不可评估(无 active) | 0/1/9/0 |
| AMES | 13/81 | 100.0% | 100.0% | 不可评估(无 active) | 0/0/13/0 |
| DILI | 8/81 | 75.0% | 83.3% | **50.0%** | 1/1/5/1 |

命令：`validate_predictions.py -i T001_ranked.csv -o validation_report.csv --deep`（hERG 基因过滤 + AMES/DILI 完整摘要；2026-07-29 执行）

明细与关键发现：
- **hERG**：10 个有数据（均 inactive）。9 一致 + 1 假阳性 **T10953**（pred 0.750 有毒 vs 实验 inactive）。1 个 "unspecified" 不计（T7184）。
- **AMES**：13 个有数据（均 inactive）。13/13 一致，特异度 100%。
  ⚠️ 但 admet-ai 把 **T19526(0.956)、TN2215(0.664)、TN2349(0.664)** 判为 AMES 致突变"有毒"，而 PubChem **无这 3 个的实验 AMES 数据** → 属"高预测风险、未实验确认"，需优先人工核对。
- **DILI（最关键的肝毒端点，且首次可评估灵敏度）**：
  - **T3595(Resmetirom)**：pred 0.999 有毒 ↔ 实验 active → **真阳性(TP)**，模型抓对了。
  - **T11954**：pred 0.402 判"安全" ↔ 实验 active（肝毒）→ **假阴性(FN)**！admet-ai 险些漏掉一个真实肝毒化合物（0.402 贴近 0.5 阈值）。
  - **T10953**：pred 0.887 有毒 ↔ 实验 inactive → **假阳性(FP)**。
  - 灵敏度仅 50%（2 个真实肝毒里只抓到 1 个）→ 对 MASLD（肝病）筛选而言，肝毒漏检风险必须正视。
  - 另有 16 个候选的 DILI 实验结局为 "unspecified" 被排除，实际覆盖更薄。

**核心局限（务必同步给对接方）：**
- 各端点 PubChem 覆盖极低（hERG 10/81、AMES 13/81、DILI 8/81），且 hERG/AMES 覆盖集**全为 inactive**，
  故 hERG/AMES 只能证明"特异度（少误报）"，灵敏度仍不可评估；**仅 DILI 因恰好有 active 样本，灵敏度可估为 50%**。
- **溶解度 / 口服生物利用度**：PubChem BioAssay 一般不含连续实验值。**已改用公开实验数据集
  （AqSolDB 9980 条实验水溶解度 + Bioavailability_Ma 640 条实验口服生物利用度）批量核对**——
  详见第 8 节。溶解度获得 26/81 命中；生物利用度 0 命中。
- 候选集中存在重复结构：T5072 / T2831 共享 CID 10140；T11557 / T17180 共享 CID 72947731
  （同 SMILES 自然映射到同一 PubChem CID，属正常去重现象，但提示源库有重复条目）。

## 7. hepatotoxicity 端点映射修正与排名重算

- **问题**：早期版本 `masld_screener.py` 与 `validate_predictions.py` 用关键字 `"hepato"` 匹配 admet-ai 输出列，
  误匹配到 `Clearance_Hepatocyte_AZ`（清除率，非 0–1 概率），导致 `toxic_score` 中肝毒贡献系统性偏低。
- **修正**：改为 `["dili", "hepatotox"]`（见 `masld_screener.py` 第 308 行注释）。该列现为 admet-ai 的 DILI 概率，
  更贴近"肝毒"语义；`Clearance` 类列不再参与肝毒评分。
- **重算**：已用修正后脚本对全量 22966 个分子重跑筛选，生成新版 `T001_ranked.csv`
  （任务 inFVzC，2026-07-28）；分享包内纳入的是此**修正版**排名。
- 注意：`validation_report.csv` 中的 `pred_DILI` 已是修正映射下的 DILI 概率，与排名一致。

## 8. 溶解度 / 口服生物利用度 实验比对（已执行，2026-07-29）

命令：`validate_solubility_bioavail.py -i T001_ranked.csv -o solubility_bioavail_report.csv`

**数据源**（均从 Harvard Dataverse / TDC 下载，见 `exp_datasets/`）：
- **AqSolDB**（datafile 4259610，9980 条实验水溶解度 LogS，连续值）→ 对应 admet-ai 端点 `Solubility_AqSolDB`
- **Bioavailability_Ma**（datafile 4259567，640 条实验口服生物利用度 0/1）→ 对应 admet-ai 端点 `Bioavailability_Ma`

匹配方式：InChIKey 精确（27 位）+ 连接性首块（14 位，供立体异构参考）。

| 端点 | 命中(精确/连接性/合计) | 一致性指标 |
|------|----------------------|-----------|
| 溶解度 AqSolDB | 8 / 18 / **26 / 81** | 精确命中 LogS 绝对误差 MAE=**0.92**，中位 0.97，误差≤1 log 单位 **5/8** |
| 生物利用度 Ma | 0 / 0 / **0 / 81** | **无候选命中实验集，无法用实验值验证** |

**8 个精确命中溶解度明细（全部为胆汁酸类已知化合物）：**

| 候选 | 身份 | pred LogS | exp LogS | 绝对误差 |
|------|------|-----------|----------|---------|
| T2202 | 石胆酸 lithocholic acid | -5.74 | -6.00 | **0.26** ✅最准 |
| T0700 | 熊脱氧胆酸 UDCA | -4.91 | -4.29 | 0.62 |
| T2963 | 胆酸 cholic acid | -4.06 | -3.37 | 0.69 |
| T4588 | 甘氨鹅脱氧胆酸 | -4.45 | -5.15 | 0.70 |
| T2965 | 去氧胆酸 deoxycholic acid | -4.92 | -3.95 | 0.97 |
| T0847 | 鹅脱氧胆酸 CDCA | -4.89 | -3.64 | 1.25 |
| T2831 | 甘氨胆酸 glycocholic acid | -3.74 | -5.15 | 1.41 |
| T5072 | 甘氨胆酸 glycocholic acid | -3.72 | -5.15 | **1.43** ⚠️偏差最大 |

**核心局限（务必同步给对接方）：**
- **⚠️ 数据泄漏警示**：AqSolDB / Bioavailability_Ma **正是 admet-ai 这两个端点的训练数据集**。
  命中候选（8 个胆汁酸）极可能在训练集内，故"预测≈实验"反映的是**模型记忆**，
  **不是**对新颖候选的独立泛化准确度。结论应表述为"报告的溶解度值与公开实验自洽（误差~1 log 单位）"，
  而非"模型溶解度预测很准"。
- **溶解度命中全部落在已知胆汁酸上**（与 `known_flag` 标记一致），排名靠前的**新颖候选无一命中实验溶解度**，
  故新颖候选的溶解度仍只有预测值、无实验佐证。
- **口服生物利用度 0 命中**：81 候选无一在 Ma 实验集内 → 生物利用度端点**完全无实验佐证**，
  仍为纯预测值，需湿实验（Caco-2 / 大鼠 PK）确认。
- 误差量级参考：admet-ai 官方报告 AqSolDB 测试集 MAE≈0.9–1.0 log 单位；本次命中子集 MAE 0.92 与之一致。
