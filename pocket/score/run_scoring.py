import csv
import json
import os
import pickle
import subprocess
import numpy as np
from biopandas.pdb import PandasPdb
import lmdb

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
APO_DIR = os.path.join(PROJECT_DIR, 'protein', 'pdb_apo')
POCKET_DIR = os.path.join(PROJECT_DIR, 'pocket', 'pocket_output')
SCORE_DIR = os.path.join(SCRIPT_DIR, 'results')
DRUGCLIP_DIR = os.path.expanduser('~/Codeshub/DrugCLIP')
WEIGHT_PATH = os.path.join(DRUGCLIP_DIR, 'checkpoint_best.pt')
MOL_LMDB = os.path.join(SCRIPT_DIR, 'mols.lmdb')
ID_MAP = os.path.join(SCRIPT_DIR, 'mol_id_map.csv')
MOL_CSV = os.path.join(PROJECT_DIR, 'molecule', 'T001_TargetMol_SMILES.csv')

POCKET_RADIUS = 8.0
POCKET_ATOM_LIST = ["C", "N", "O", "S", "H", "F", "Cl", "Br", "I", "P", "B", "Si", "Se", "c", "n", "o", "s"]
MOL_ATOM_LIST = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "B", "Si", "Se", "c", "n", "o", "s"]

GENES = ['THRB', 'NR1H4', 'PPARA', 'PPARG', 'PPARD', 'ACACA',
         'FASN', 'AOC3', 'SCD1', 'SCD', 'NR1H3', 'MLXIPL',
         'XPR1', 'GPAM', 'RXRA', 'THRSP', 'RXRB', 'RXRG']
SKIP = ['DGAT2', 'ELOVL6']


def get_pocket_center(pocket_pdb):
    coords = []
    with open(pocket_pdb) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except:
                    pass
    if not coords:
        return None
    return np.mean(coords, axis=0)


def get_element_from_atom_name(atom_name):
    elem = atom_name.strip()
    if elem and elem[0].isalpha():
        e = elem[0].upper()
        if len(elem) > 1 and elem[1].islower():
            return elem[:2].capitalize()
        return e
    return 'C'


def extract_pocket_atoms(apo_pdb, center):
    ppdb = PandasPdb().read_pdb(apo_pdb)
    atom_df = ppdb.df['ATOM']
    coords = atom_df[['x_coord', 'y_coord', 'z_coord']].values
    dists = np.linalg.norm(coords - center, axis=1)
    mask = dists < POCKET_RADIUS
    if mask.sum() == 0:
        mask = dists < 12.0
    return {
        'atom_types': [get_element_from_atom_name(n) for n in atom_df.loc[mask, 'atom_name']],
        'coords': coords[mask],
    }


def write_lmdb(entries, lmdb_path):
    env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False,
                    readahead=False, meminit=False, map_size=1099511627776)
    with env.begin(write=True) as txn:
        for i, e in enumerate(entries):
            txn.put(str(i).encode('ascii'), pickle.dumps(e))
    env.close()


def write_dict(atom_list, dict_path):
    with open(dict_path, 'w') as f:
        for atom in atom_list:
            f.write(f'{atom} 1\n')


def run_drugclip(mol_lmdb, pocket_lmdb, emb_dir, out_dir):
    cmd = [
        'python', os.path.join(DRUGCLIP_DIR, 'unimol', 'retrieval.py'),
        '--user-dir', os.path.join(DRUGCLIP_DIR, 'unimol'),
        os.path.join(DRUGCLIP_DIR, 'data'),
        '--valid-subset', 'test',
        '--results-path', out_dir,
        '--num-workers', '4',
        '--batch-size', '48',
        '--task', 'drugclip',
        '--loss', 'in_batch_softmax',
        '--arch', 'drugclip',
        '--max-pocket-atoms', '256',
        '--fp16',
        '--seed', '1',
        '--path', WEIGHT_PATH,
        '--mol-path', mol_lmdb,
        '--pocket-path', pocket_lmdb,
        '--emb-dir', emb_dir,
    ]
    env = os.environ.copy()
    env['CUDA_VISIBLE_DEVICES'] = '0'
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=86400, env=env)
    if result.returncode != 0:
        print(f'  DrugCLIP ERROR: {result.stderr[-300:]}')
        return False
    return True


def process_gene(gene, smiles_map):
    print(f'\n{"="*60}')
    print(f'{gene}...')
    out_dir = os.path.join(SCORE_DIR, gene)
    os.makedirs(out_dir, exist_ok=True)

    pocket_dir = os.path.join(POCKET_DIR, gene, 'pocket_detected')
    apo_path = os.path.join(APO_DIR, f'{gene}_apo.pdb')
    if not os.path.exists(pocket_dir) or not os.path.exists(apo_path):
        print(f'  Skip: missing data')
        return

    pocket_pdbs = sorted([f for f in os.listdir(pocket_dir) if f.endswith('.pdb')])

    # Build pocket LMDB
    pocket_entries = []
    for pp in pocket_pdbs:
        pid = pp.replace('.pdb', '')
        center = get_pocket_center(os.path.join(pocket_dir, pp))
        if center is None:
            continue
        pocket_atoms = extract_pocket_atoms(apo_path, center)
        if len(pocket_atoms['atom_types']) < 10:
            continue
        pocket_entries.append({
            'pocket_atoms': pocket_atoms['atom_types'],
            'pocket_coordinates': [pocket_atoms['coords'][i] for i in range(len(pocket_atoms['coords']))],
            'pocket': pid,
            'pocket_index': 0,
        })

    if not pocket_entries:
        print(f'  No valid pockets')
        return
    print(f'  Pockets: {len(pocket_entries)}')

    pocket_lmdb = os.path.join(out_dir, 'pockets.lmdb')
    emb_dir = os.path.join(out_dir, 'emb')
    write_lmdb(pocket_entries, pocket_lmdb)
    write_dict(POCKET_ATOM_LIST, os.path.join(out_dir, 'dict_pkt.txt'))
    write_dict(MOL_ATOM_LIST, os.path.join(out_dir, 'dict_mol.txt'))

    # Run DrugCLIP
    print(f'  Running DrugCLIP...')
    if not run_drugclip(MOL_LMDB, pocket_lmdb, emb_dir, out_dir):
        return

    # Parse results
    ranked_file = os.path.join(emb_dir, 'ranked_compounds.txt')
    if os.path.exists(ranked_file):
        scores = []
        with open(ranked_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    scores.append((parts[0], float(parts[1])))
        scores.sort(key=lambda x: -x[1])

        smi_to_id = {v: k for k, v in smiles_map.items()}

        with open(os.path.join(out_dir, 'scores.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['Rank', 'ID', 'SMILES', 'Score'])
            for rank, (smi, score) in enumerate(scores, 1):
                mol_id = smi_to_id.get(smi, '')
                w.writerow([rank, mol_id, smi, f'{score:.4f}'])
        print(f'  Scores: {len(scores)} molecules ranked')
    else:
        print(f'  No ranked_compounds.txt')


def main():
    os.makedirs(SCORE_DIR, exist_ok=True)
    if not os.path.exists(MOL_LMDB):
        print(f'Error: {MOL_LMDB} not found. Run prepare_mols.py first')
        return

    # Load ID map
    id_map = {}
    if os.path.exists(ID_MAP):
        with open(ID_MAP) as f:
            for row in csv.DictReader(f):
                id_map[int(row['Index'])] = row['ID']
    print(f'Loaded {len(id_map)} ID mappings')

    # Load SMILES map
    smiles_map = {}
    if os.path.exists(MOL_CSV):
        with open(MOL_CSV) as f:
            for row in csv.DictReader(f):
                smiles_map[row['ID']] = row['SMILES']
    print(f'Loaded {len(smiles_map)} SMILES')

    for gene in GENES:
        if gene in SKIP:
            continue
        process_gene(gene, smiles_map)
    print(f'\nAll done.')


if __name__ == '__main__':
    main()