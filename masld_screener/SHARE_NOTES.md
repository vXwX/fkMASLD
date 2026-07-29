# MASLD 虚拟筛选 — 共享包说明（2026-07-29 · v6）

本压缩包是 `masld_screener` Skill 的完整交付物，供对接方直接复用或继续开发。

## 包内内容

| 文件/目录 | 说明 |
|-----------|------|
| `masld_screener/` | Skill 全部源码与产物（可直接 `python masld_screener.py <sdf> -o out.csv` 运行） |
| `masld_screener/masld_screener.py` | 主筛选流程：理化计算(RDKit) + 结构相似度(分子骨架/靶点锚点) + admet-ai 毒性预测 |
| `masld_screener/validate_predictions.py` | 毒性端点交叉验证脚本（PubChem BioAssay：hERG/AMES/DILI） |
| `masld_screener/validate_solubility_bioavail.py` | **本次新增**：溶解度/生物利用度实验比对脚本（AqSolDB + Bioavailability_Ma） |
| `masld_screener/T001_ranked.csv` | **筛选结果（修正版）**：22966 个分子的 Top 81 候选，含理化/相似度/毒性/综合分；身份列 `known_flag`/`known_as`/`pubchem_cid`/`pubchem_name`/`relevance_note`，并新增 `toxic_score_basis`（标注毒性分为预测值非实测）；**81 个已逐一 PubChem 核验身份** |
| `masld_screener/confirm_identity.py` | **本次新增**：逐一 PubChem 身份核验脚本（InChIKey 精确 + 连接性层） |
| `masld_screener/crosscheck_physchem.py` | **v6 新增**：physchem 描述符独立对账脚本（RDKit 计算值 vs PubChem 报告值 MW/XLogP/TPSA/HBD/HBA/RotB） |
| `masld_screener/finalize_crosscheck.py` | **v6 新增**：对账结果分级（OK/计数惯例/TPSA方法学/LogP算法差/结构存疑） |
| `masld_screener/physchem_crosscheck.csv` | **v6 新增**：81 候选双源 physchem 对账表，含 `severity`/`interpretation` 列 |
| `masld_screener/T001_run.log` | 全量筛选运行日志（两万多个分子的筛选过程可追溯） |
| `masld_screener/T001_ranked.sdf` | 同批候选的结构 SDF |
| `masld_screener/validation_report.csv` | hERG/AMES/DILI 预测 vs 实验 一致性报告（PubChem） |
| `masld_screener/solubility_bioavail_report.csv` | **本次新增**：溶解度/生物利用度 预测 vs 实验 比对（AqSolDB/Ma） |
| `masld_screener/VALIDATION_PLAN.md` | **本次更新**：实验数据库清单、五端点验证结果与局限（第 6/8 节） |
| `masld_screener/DB_CONFIRMATION.md` | 9 个参考分子及 Top 候选的身份/结构/标注联网核验 |
| `masld_screener/SKILL.md` / `README.md` / `config.py` / `reference_molecules.py` | Skill 元信息与配置 |
| `masld_screener/test_*.sdf` / `test_ranked.csv` | 10 分子小样例，便于快速验证流程 |

## 本次（2026-07-29）关键变更

0. **🔴 81 个候选逐一 PubChem 身份核验（本次重点）**：用 `confirm_identity.py` 按 InChIKey 精确匹配，
   确认 **81/81 全部为 PubChem 已收录已知化合物，无一新颖**；CSV 补 `pubchem_cid`/`pubchem_name`/`relevance_note` 列；
   标注 3 个"结构假阳性"（Diclazuril/Pradefovir/Dafadine-A，其他适应症已知药）。详见下方"已知化合物逐一核验结论"。
1. **实验值交叉验证已扩展为三端点（hERG / AMES / DILI）**：用 PubChem PUG-REST
   （hERG 基因过滤 KCNH2 + AMES/DILI 完整摘要 `--deep`）拉取实验 BioAssay，与 admet-ai 预测逐项比对。
   - **hERG**：覆盖 10/81，一致率 90%，特异度 90%（灵敏度不可评估，覆盖集无 active）。
   - **AMES**：覆盖 13/81，一致率 100%，特异度 100%（灵敏度不可评估）；T19526/TN2215/TN2349 被预测为 AMES 有毒但无实验数据。
   - **DILI**：覆盖 8/81，一致率 75%，特异度 83%、**灵敏度 50%**（首次可评估）；含 1 假阴性 T11954（漏检真实肝毒）、1 假阳性 T10953。
   - 详见 `VALIDATION_PLAN.md` 第 6 节。
2. **溶解度 / 口服生物利用度 也补上了实验比对**（用公开实验数据集，非 PubChem）：
   下载 **AqSolDB**（9980 条实验水溶解度）与 **Bioavailability_Ma**（640 条实验口服生物利用度），按 InChIKey 匹配 81 候选。
   - **溶解度**：命中 **26/81**（8 精确 + 18 连接性）；精确命中 LogS 绝对误差 MAE=**0.92** log 单位（与 admet-ai 官方 ~0.9 一致）。
     ⚠️ 但 8 个精确命中**全是胆汁酸类已知化合物**，且 AqSolDB 正是 admet-ai 训练集 → 属**模型记忆**，非新颖候选的泛化验证。
   - **生物利用度**：**0/81 命中**，81 候选无一在 Ma 实验集内 → 该端点**完全无实验佐证**，仍是纯预测值。
   - 详见 `VALIDATION_PLAN.md` 第 8 节。
3. **hepatotoxicity 端点映射修正**：早期用 `"hepato"` 误匹配到 admet-ai 的 `Clearance_Hepatocyte_AZ` 列，
   导致 `toxic_score` 系统性偏低；已改为 `["dili","hepatotox"]`。`T001_ranked.csv` 已用修正版重算。
4. **🔵 physchem 描述符独立对账（v6 新增，回应"每列是否都核对过数据库"）**：用 `crosscheck_physchem.py`
   把 CSV 内 RDKit 计算的 `MW/LogP/TPSA/HBD/HBA/RotB` 与 PubChem PUG-REST 报告的同字段逐一比对（81/81 全部取到 MW）。
   - **MW 不符（结构存疑）= 0**：81 个分子量全部与 PubChem 一致 → 结构/SMILES 正确，physchem 列非杜撰。
   - 其余"标记"均为**可解释的方法学差异**，非数据错误：计数惯例差异 49（RDKit vs PubChem 受体/供体计数口径，
     如羧酸双氧计 1 vs 2）、TPSA 方法学差异 17（盐/电荷/互变异构取质子态不同，如牛磺胆酸钠盐差≈1 个磺酸酯基团）、
     LogP 算法差 8（RDKit MolLogP vs PubChem XLogP3 本就是不同模型）。
   - ⚠️ 注意：PubChem 的 MW/XLogP/TPSA 本身也是**计算值非实验实测**，本对账是"两套独立算法一致性校验"，
     严重不符才提示结构问题——本次无严重不符。详见下方"physchem 独立对账结论"。
5. **🔵 `toxic_score_basis` 列新增（v6）**：明确 `toxic_score` 与 `final_score` 来自 **ADMET-AI 机器学习预测**，
   **非实验实测毒理数据**（hERG/AMES/DILI 实验库零命中）；发送/对接时不得把这两列写成"实测值"。

## 🔴 已知化合物逐一核验结论（最重要，发送前必看）

**已用 `confirm_identity.py` 对全部 81 个候选按 InChIKey 精确匹配 PubChem 逐一核验身份。结论：**

> **81 / 81（100%）候选在 PubChem 均有精确匹配（含立体化学），全部是已收录/已知化合物，无一新颖骨架。**

`T001_ranked.csv` 现含 5 列身份信息：`known_flag` / `known_as` / `pubchem_cid` / `pubchem_name` / `relevance_note`。分布：

- **参考药本体（4 个）**：T7395=Lanifibranor、T1789=OCA、T3595=Resmetirom、T7184=Firsocostat（TC=1.0，即锚点药自身）。
- **已知药/工具化合物（多个）**：Denifanstat(T15271)、Sobetirome(T5313)、MB-07344/VK2809 活性体(T11954)、INT-777(T11662L, TGR5)、BAR501/BAR502(T4083/TQ0252)、TVB-3664(T17181)、ND-646(TQ0243, ACC)、Terbufibrol(T13919)、Bsh-IN-1(T10623) 等。
- **已知胆汁酸及其衍生/结合物（绝大多数）**：UDCA、CDCA、胆酸、去氧胆酸、石胆酸、各类甘氨/牛磺结合物、酮基/异构体等——全部是内源或已收录胆汁酸。
- **已收录(仅系统名)（12 个）**：PubChem 有 CID 但只挂 IUPAC 系统名（如钠盐、立体异构体），**仍属已知收录化合物**，非新颖。

### ⚠️ 3 个"结构假阳性"必须点名（`relevance_note` 列已标注）

以下 3 个虽是已知药，但属**其他适应症**，仅因子结构相似被 Tanimoto 拉进来，**与 MASLD 无关**，不应视作有效候选：

| 候选 | 真实身份 | 为何误入 |
|------|---------|---------|
| **T1172** | Diclazuril（抗球虫兽药） | 与 Resmetirom 共享三嗪二酮子结构 |
| **T19526** | Pradefovir（抗乙肝前药） | 与 VK2809 共享环磷酸酯前药子结构 |
| **T10953** | Dafadine-A（CYP/DAF-9 工具抑制剂） | 与 Denifanstat 结构相似 |

**总体含义**：本次筛选命中的是"已知化学空间"（现货库本就是已合成/已收录化合物），**没有新颖分子**——这对"药物重定位/已知活性验证"有价值，但若目标是"发现全新骨架"，则需换新颖化合物库或做骨架跃迁。IP 与新颖性风险务必与对接方讲清。

## 🔵 physchem 独立对账结论（v6 新增，回应"每列是否都核对过数据库"）

**否——只有身份列（已知_flag/known_as/pubchem_cid/pubchem_name/relevance_note）做了数据库核验；physchem 与评分列是计算/预测值，未逐库核对。但 physchem 列已做"双源一致性对账"：**

用 `crosscheck_physchem.py` 对每个候选取 PubChem PUG-REST 报告的 `MolecularWeight / XLogP / TPSA / HBondDonorCount / HBondAcceptorCount / RotatableBondCount`，与 CSV 内 RDKit 计算值逐一比对：

| 结果 | 数量 | 含义 |
|------|------|------|
| MW 不符（结构存疑） | **0 / 81** | 分子量全与 PubChem 一致 → 结构/SMILES 正确，physchem 列非杜撰 |
| 计数惯例差异（良性） | 49 | RDKit vs PubChem 受体/供体/可旋转键计数口径不同，不影响筛选 |
| TPSA 方法学差异 | 17 | 盐/电荷/互变异构取质子态不同（如牛磺胆酸钠盐差≈1 个磺酸酯基团 ~22；Lanifibranor 差 36.7 因 PubChem 规范式取质子态不同） |
| LogP 算法差 | 8 | RDKit MolLogP 与 PubChem XLogP3 本就是不同预测模型，差异正常 |
| 完全一致（OK） | 7 | 两套算法零差异 |

> ⚠️ 关键澄清：PubChem 的 MW / XLogP / TPSA 本身也是**计算值（in-silico），非实验实测**。本对账是"两套独立算法一致性校验"——只有严重不符才提示结构/SMILES 问题；本次**无严重不符**，故 physchem 列可放心用于排序，但仍是预测性质，最终以湿实验为准。

完整逐行数据见 `masld_screener/physchem_crosscheck.csv`（含 `severity` / `interpretation` 列）。

## 实验数据库检索完成度（务必如实转告）

- **已实际检索并比对（数据源 PubChem BioAssay）**：hERG（基因过滤 KCNH2）+ AMES + DILI（`--deep` 扫描完整摘要）三个毒性端点：
  - **hERG**：覆盖 10/81，一致率 90%，特异度 90%（灵敏度因覆盖集无 active 样本不可评估）。
  - **AMES**：覆盖 13/81，一致率 100%，特异度 100%（灵敏度不可评估）；但 **T19526 / TN2215 / TN2349** 被 admet-ai 判为 AMES 致突变"有毒"却无实验数据 → 属"未确认高风险"，需优先人工核对。
  - **DILI（肝毒，最关键）**：覆盖 8/81，一致率 75%，特异度 83.3%、**灵敏度 50%**（首次可评估）；含 1 假阴性 **T11954**（admet-ai 判安全却实验 active，险些漏检真实肝毒）、1 假阳性 **T10953**。
- **已实际检索并比对（数据源 AqSolDB / Bioavailability_Ma 公开实验集，经 Harvard Dataverse/TDC 下载）**：
  - **溶解度**：命中 26/81（8 精确 + 18 连接性），精确命中 LogS 绝对误差 MAE=0.92 log 单位；但命中全为胆汁酸已知化合物，且属 admet-ai 训练集（数据泄漏，见下）。
  - **口服生物利用度**：0/81 命中，无实验佐证。
- **⚠️ 数据泄漏警示**：AqSolDB / Bioavailability_Ma 即 admet-ai 对应端点的**训练集**，命中候选 = 模型很可能已见过 →
  溶解度"预测≈实验"反映的是**记忆而非泛化**，不能当作新颖候选溶解度准确度的证据。
- **仍未检索的源**：admetSAR 3.0 网站可访问但检索需前端 JS 交互、resource 页无可下载数据集；ChEMBL REST 在本环境 404；
  **ChEMBL / BindingDB 实验活性**（验证"降脂活性"代理）仍未查。
- 结论：**五个 ADMET 端点（hERG/AMES/DILI/溶解度/生物利用度）均已尽本环境所能做实验比对**；
  尚缺的是：① 含 active 样本的毒性实验集（补 hERG/AMES 灵敏度）；② 新颖候选的独立（非训练集）溶解度实测；
  ③ 生物利用度实测；④ 靶点活性（ChEMBL/BindingDB）。这些需湿实验或人工下载专库，非免费 API 可自动完成。

## 必须同步给对接方的方法学边界

- **源库 SDF 无任何实验值**：理化/相似度是 RDKit 确定性计算；ADMET 毒性是 admet-ai 的 ML 预测（非实测）。
- **"降脂活性"是相似度代理（类比假设）**，不是对候选自身靶点的直接活性预测。
- 全流程仅为**优先级排序**，最终需湿实验确认。
- hERG 自动验证覆盖极低且偏向 inactive，**灵敏度未经验证**；建议到 admetSAR 3.0 / ADMETNet
  拉取含 active 样本的实验 hERG/AMES/DILI 集做完整一致率评估（见 `VALIDATION_PLAN.md` 第 5 节）。

## 快速复跑

```bash
# 重新筛选
python masld_screener/masld_screener.py "T001 TargetMol现货产品22966.sdf" -o masld_screener/T001_ranked.csv
# 重新做 hERG 验证（默认快速模式）
python masld_screener/validate_predictions.py -i masld_screener/T001_ranked.csv -o masld_screener/validation_report.csv
# 额外做 AMES/DILI 验证（较慢，可能超时）
python masld_screener/validate_predictions.py -i masld_screener/T001_ranked.csv -o masld_screener/validation_report.csv --deep
# 溶解度/生物利用度实验比对（需 exp_datasets/ 下的 AqSolDB 与 Bioavailability_Ma 原始文件）
python masld_screener/validate_solubility_bioavail.py -i masld_screener/T001_ranked.csv -o masld_screener/solubility_bioavail_report.csv
# 逐一 PubChem 身份核验（写回 known_flag/known_as/pubchem_cid/pubchem_name）
python masld_screener/confirm_identity.py -i masld_screener/T001_ranked.csv -o masld_screener/T001_ranked.csv
# physchem 双源对账（RDKit vs PubChem），输出 physchem_crosscheck.csv
python masld_screener/crosscheck_physchem.py
# 对账结果分级标注（写回 severity/interpretation 列）
python masld_screener/finalize_crosscheck.py
```
