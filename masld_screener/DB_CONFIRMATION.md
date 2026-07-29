# MASLD Screener — 数据库交叉核验报告 (DB_CONFIRMATION)

- 生成日期：2026-07-28
- 核验来源：PubChem PUG-REST（结构 / 分子式 / 分子量）、公开文献与临床试验库（靶点 / 临床阶段）
- 目的：应队友要求，对“数据库里能查到的”事实（参考分子结构/靶点/阶段、候选分子身份）做外部交叉确认。

---

## 1. 参考分子（9 个已知药）结构核验

| # | 分子 | 与我方 SMILES 比对 | PubChem 分子式 | 分子量 | 结论 |
|---|------|-------------------|---------------|--------|------|
| 1 | Resmetirom | 公式/分子量一致 | C17H12Cl2N6O4 | 435.2 | ✅ 一致 |
| 2 | VK2809 | 名称解析→公式/分子量一致 | C28H32ClO5P | 515.0 | ✅ 一致 |
| 3 | OCA | 与 PubChem 异构 SMILES 完全相同 | C26H44O4 | 420.6 | ✅ 精确匹配 |
| 4 | Cilofexor | 与 PubChem 异构 SMILES 完全相同 | C28H22Cl3N3O5 | 586.8 | ✅ 精确匹配 |
| 5 | Lanifibranor | 公式/分子量一致 | C19H15ClN2O4S2 | 434.9 | ✅ 一致 |
| 6 | Firsocostat | 与 PubChem 异构 SMILES 完全相同 | C28H31N3O8S | 569.6 | ✅ 精确匹配 |
| 7 | Denifanstat | 与 PubChem 异构 SMILES 完全相同 | C27H29N5O | 439.6 | ✅ 精确匹配 |
| 8 | BI-1467335 | 游离碱（去 HCl 盐） | C15H21FN2O2（游离碱） | ~300 | ✅ 一致（设计性去盐） |
| 9 | Aramchol | 与 PubChem 异构 SMILES 完全相同 | C44H79NO5 | 702.1 | ✅ 精确匹配 |

**结论**：9/9 参考分子结构正确。OCA / Cilofexor / Firsocostat / Denifanstat / Aramchol 与我方 SMILES 逐字符一致；其余 4 个公式与分子量一致。BI-1467335 按设计去掉盐酸盐反离子（避免指纹引入氯离子片段），PubChem 条目含 `.Cl` 盐型，连接性相同。

---

## 2. 参考分子 靶点 / 临床阶段 核验

| 分子 | 靶点 | 我方阶段 | 核验结论 |
|------|------|---------|---------|
| Resmetirom | THR-β | approved | ✅ 2024 获批（Rezdiffra） |
| VK2809 | THR-β | Phase 2 | ✅ 肝选择性 THR-β 激动剂，NASH Phase 2/3 |
| OCA | FXR | approved (PBC) | ✅ 获批 PBC；MASH Phase 3（REGENERATE） |
| Cilofexor | FXR | **Phase 3** | ⚠️ 修正：NASH 为 **Phase 2**（已发表 Phase 2 RCT NCT02854605）；另有 PSC Phase 3。Gilead 的 cilofexor+firsocostat MASH 组合已于 2025 年终止开发 |
| Lanifibranor | PPAR α/γ/δ | Phase 3 | ✅ NATiV3 Phase 3 进行中 |
| Firsocostat | ACC | Phase 2 | ✅ 另：MASH 组合 2025 终止 |
| Denifanstat | FASN | Phase 3 | ✅ TVB-2640，FASCINATE-2 Phase 2b/3 |
| BI-1467335 | AOC3/VAP-1 | Phase 2 | ✅ Phase IIa（NCT03166735）达靶 engagement；后因 MAO-B DDI 风险终止 |
| Aramchol | SCD1 | Phase 2 | ✅ Phase 2b/3 |

**结论**：靶点全部正确。阶段需修正 1 处（Cilofexor → NASH Phase 2）。另：Cilofexor（MASH 组合）、BI-1467335、Firsocostat（MASH 组合）已终止开发——作为相似度锚点仍有效，但作为“在研对标”需注意。

---

## 3. Top 候选身份核验（CAS → PubChem）

| Rank | ID | CAS | PubChem 身份 | 我方 MW | SDF MW | 命中参考 | 说明 |
|------|----|-----|------------|--------|--------|---------|------|
| 1 | T7395 | 927961-18-0 | **Lanifibranor** | 434.9 | 434.9 | Lanifibranor (TC 1.0) | 即参考药本体 |
| 2 | T7331 | 2955-27-3 | Ursocholic acid | 408.6 | 408.6 | Aramchol (0.69) | 胆汁酸 |
| 3 | T19189 | 2464-18-8 | Allocholic acid | 408.6 | 408.6 | Aramchol (0.69) | 胆汁酸 |
| 4 | T0700 | 128-13-2 | Ursodeoxycholic acid (UDCA) | 392.6 | 392.6 | OCA (0.67) | 内源 FXR 配体 |
| 5 | T0847 | 474-25-9 | Chenodeoxycholic acid (CDCA) | 392.6 | 392.6 | OCA (0.67) | 内源 FXR 配体 |
| 6 | T5072 | 1192657-83-2 | Glycocholic acid | 465.6* | 483.6† | Aramchol (0.62) | 甘氨酸结合胆汁酸；†含水 |
| 7 | T13522 | 4651-67-6 | 7-Ketolithocholic acid | 390.6 | 390.6 | OCA (0.61) | 胆汁酸衍生物 |
| 8 | TMIH-0148 | 116380-66-6 | Cholic acid-2,2,4,4-D4 | 412.6 | 412.6 | Aramchol (0.60) | 氘代胆酸 |
| 9 | T1789 | 459789-99-2 | **Obeticholic acid (OCA)** | 420.6 | 420.6 | OCA (TC 1.0) | 即参考药本体 |
| 10 | T3595 | 920509-32-6 | **Resmetirom** | 435.2 | 435.2 | Resmetirom (TC 1.0) | 即参考药本体 |
| 11 | T7184 | 1434635-54-7 | **Firsocostat (GS-0976)** | 569.6 | 569.6 | Firsocostat (TC 1.0) | 即参考药本体 |
| 12 | TQ0243 | 1434639-57-2 | Firsocostat 类似物 | 568.6 | 568.6 | Firsocostat (0.89) | 噻吩并嘧啶酰胺类似物 |
| 13 | T5234 | 64480-66-6 | Glycoursodeoxycholic acid | 449.6 | 449.6 | OCA (0.55) | 甘氨酸结合 UDCA |
| 14 | TWA2417 | 145-42-6 | Sodium Taurocholate | 514.7* | 537.7† | OCA (0.54) | 牛磺酸结合胆酸；†钠盐 |

> * 我方 SMILES 计算的游离酸 MW；† 源库标注为盐/水合物（Na 盐或 hydrate），差值来自抗衡离子/结晶水，非结构错误。

**结论**：Top 候选身份经 CAS → PubChem 全数确认，无身份存疑项。

---

## 4. 关键发现与给队友的提示

1. **方法学验证成功**：4 个 Top 候选（T7395 / T1789 / T3595 / T7184）以 TC=1.0 精确命中对应参考药，证明相似性检索逻辑正确、参考集可信。
2. **Top 命中富集于“已知药 + 内源胆汁酸”**：前 14 名里 4 个是临床药本体，10 个是胆汁酸及其结合物（UDCA / CDCA / 胆酸 / 鹅去氧胆酸 / Glyco- / Tauro- 结合物）。因 OCA(FXR) 与 Aramchol(SCD1/胆汁酸偶联) 的相似锚点本身富集胆酸骨架，**结构新颖性有限，IP 风险需评估**。
3. **真正新颖的 chemotype 较少**：TQ0243（Firsocostat 类似物）是少数接近但非本体的结构；若队友要“新骨架”，建议把参考药本体及其直接类似物/胆汁酸在后续对接中视情况剔除或降权。
4. **盐/水合物处理一致**：候选 MW 与我方计算值、PubChem 自由酸值一致；源库标注的盐/水合物差异已识别，不影响结构身份。
5. **参考药开发状态更新**：Cilofexor（NASH Phase 2，组合已终止）、BI-1467335（已终止）、Firsocostat（组合已终止）——作为相似度锚点仍有效，但作为“在研对标”需注意。

---

## 5. 最终结论

- 参考集 **9/9 结构正确**，靶点正确，1 处阶段标注已修正（Cilofexor）。
- 候选身份经 CAS → PubChem **全数确认**，无身份存疑项。
- 本核验确认的是**结构 / 身份 / 标注正确性**，**仍非“预测值 vs 实验值”的置信度标定**（源库无实验活性数据）。若需真正置信度，请用 Top 候选（已知药/胆汁酸）的公开 Assay 数据做对照（参见上一轮说明）。
