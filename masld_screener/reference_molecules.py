# -*- coding: utf-8 -*-
"""
reference_molecules.py — MASLD 临床阶段已知活性小分子参考集。

数据来源：PubChem PUG-REST（IsomericSMILES），检索日期 2026-07-28。
校验方式：每个 SMILES 经 RDKit 解析后计算的分子量(MW)与分子式(Formula)
          与 PubChem 报告值逐一比对，全部一致后录入。

说明：
  - aramchol 的脂酰链为花生酰基（arachidyl, C20），已据此修正链长。
  - BI-1467335 以游离碱形式存储（去除盐酸盐抗衡离子），以保证指纹不含氯离子片段。
  - 数据库核验（2026-07-28，见 DB_CONFIRMATION.md）：
      9/9 参考分子 SMILES 经 PubChem 结构比对一致；靶点全部正确。
      Cilofexor 的 NASH 阶段为 Phase 2（另有 PSC Phase 3）。
      部分参考药（Cilofexor MASH 组合、BI-1467335、Firsocostat MASH 组合）
      已终止开发，此处仅作相似度检索锚点，不代表在研状态。

每个分子以字典存储：{"smiles": "...", "target": "...", "stage": "..."}
"""

REFERENCE_MOLECULES = [
    {
        "name": "Resmetirom",
        "smiles": "CC(C)C1=CC(=NNC1=O)OC2=C(C=C(C=C2Cl)N3C(=O)NC(=O)C(=N3)C#N)Cl",
        "target": "THR-beta",
        "stage": "approved",
    },
    {
        "name": "VK2809",
        "smiles": "CC1=CC(=CC(=C1CC2=CC(=C(C=C2)O)C(C)C)C)OC[P@@]3(=O)OCC[C@H](O3)C4=CC(=CC=C4)Cl",
        "target": "THR-beta",
        "stage": "Phase 2",
    },
    {
        "name": "Obeticholic acid (OCA)",
        "smiles": "CC[C@@H]1[C@@H]2C[C@@H](CC[C@@]2([C@H]3CC[C@]4([C@H]([C@@H]3[C@@H]1O)CC[C@@H]4[C@H](C)CCC(=O)O)C)C)O",
        "target": "FXR",
        "stage": "approved (PBC)",
    },
    {
        "name": "Cilofexor (GS-9674)",
        "smiles": "C1CC1C2=C(C(=NO2)C3=C(C=CC=C3Cl)Cl)COC4=CC(=C(C=C4)C5(CN(C5)C6=NC=CC(=C6)C(=O)O)O)Cl",
        "target": "FXR",
        "stage": "Phase 2 (NASH); PSC Phase 3",
    },
    {
        "name": "Lanifibranor (IVA-337)",
        "smiles": "C1=CC2=C(C=C1S(=O)(=O)N3C4=C(C=C(C=C4)Cl)C=C3CCCC(=O)O)SC=N2",
        "target": "PPAR-alpha/gamma/delta",
        "stage": "Phase 3",
    },
    {
        "name": "Firsocostat (GS-0976)",
        "smiles": "CC1=C(SC2=C1C(=O)N(C(=O)N2C[C@@H](C3=CC=CC=C3OC)OC4CCOCC4)C(C)(C)C(=O)O)C5=NC=CO5",
        "target": "ACC",
        "stage": "Phase 2",
    },
    {
        "name": "Denifanstat (TVB-2640)",
        "smiles": "CC1=CC(=C(C=C1C(=O)N2CCC(CC2)C3=CC=C(C=C3)C#N)C4=NNC(=N4)C)C5CCC5",
        "target": "FASN",
        "stage": "Phase 3",
    },
    {
        "name": "BI-1467335",
        "smiles": "CC(C)(C)NC(=O)C1=CC=C(C=C1)OC/C(=C/F)/CN",
        "target": "AOC3/VAP-1",
        "stage": "Phase 2",
    },
    {
        "name": "Aramchol",
        "smiles": "CCCCCCCCCCCCCCCCCCCC(=O)N[C@H]1CC[C@]2([C@@H](C1)C[C@H]([C@@H]3[C@@H]2C[C@H]([C@]4([C@H]3CC[C@@H]4[C@H](C)CCC(=O)O)C)O)O)C",
        "target": "SCD1",
        "stage": "Phase 2",
    },
]


def get_reference_molecules():
    """返回参考分子列表的副本，避免外部修改原始数据。"""
    return [dict(m) for m in REFERENCE_MOLECULES]


if __name__ == "__main__":
    for m in REFERENCE_MOLECULES:
        print(f"{m['name']:28s} {m['target']:22s} {m['stage']}")
