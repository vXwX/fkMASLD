from rdkit import Chem
from rdkit.Chem import Descriptors
import csv

sdf_path = '/home/gaoxiangxu/Competition/AI4S_202607/molecule/T001 Jan-2026/T001 TargetMol现货产品22966.sdf'
csv_ok = '/home/gaoxiangxu/Competition/AI4S_202607/molecule/T001_TargetMol_SMILES.csv'
csv_fail = '/home/gaoxiangxu/Competition/AI4S_202607/molecule/T001_TargetMol_SMILES_failed.csv'
csv_dedup = '/home/gaoxiangxu/Competition/AI4S_202607/molecule/T001_TargetMol_SMILES_dedup_report.csv'

suppl = Chem.SDMolSupplier(sdf_path)
smiles_to_ids = {}

with open(csv_ok, 'w', newline='') as fok, open(csv_fail, 'w', newline='') as ffail, open(csv_dedup, 'w', newline='') as fdup:
    wok = csv.writer(fok)
    wfail = csv.writer(ffail)
    wdup = csv.writer(fdup)
    wok.writerow(['ID', 'SMILES', 'MolWt'])
    wfail.writerow(['Index', 'SDF_Line'])
    wdup.writerow(['SMILES', 'Duplicate_IDs', 'Count'])
    for idx, mol in enumerate(suppl, 1):
        if mol is None:
            wfail.writerow([idx, f'Entry #{idx}'])
            continue
        smiles = Chem.MolToSmiles(mol)
        mol_id = mol.GetProp('ID') if mol.HasProp('ID') else ''
        smiles_to_ids.setdefault(smiles, []).append(mol_id)

    smiles_seen = set()
    for smiles, ids in smiles_to_ids.items():
        if len(ids) > 1:
            wdup.writerow([smiles, ';'.join(ids), len(ids)])
        if smiles in smiles_seen:
            continue
        smiles_seen.add(smiles)
        molwt = Chem.Descriptors.MolWt(Chem.MolFromSmiles(smiles))
        wok.writerow([ids[0], smiles, f'{molwt:.2f}'])