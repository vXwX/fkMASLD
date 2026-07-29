import csv
import os
import pickle
import signal
import sys
from multiprocessing import Pool, cpu_count
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
MOL_CSV = os.path.join(PROJECT_DIR, 'molecule', 'T001_TargetMol_SMILES.csv')
OUTPUT = os.path.join(SCRIPT_DIR, 'mols.lmdb')
ID_MAP = os.path.join(SCRIPT_DIR, 'mol_id_map.csv')

NUM_CONFS = 5


class TimeoutError(Exception):
    pass


def timeout_handler(signum, frame):
    raise TimeoutError


def gen_conformation(smiles):
    signal.alarm(120)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        signal.alarm(0)
        return None
    try:
        mol = Chem.AddHs(mol)
        params = AllChem.EmbedMultipleConfs(
            mol, numConfs=NUM_CONFS, numThreads=1,
            pruneRmsThresh=1, maxAttempts=200, useRandomCoords=True)
        if not params:
            signal.alarm(0)
            return None
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1)
        mol = Chem.RemoveHs(mol)
    except:
        signal.alarm(0)
        return None
    signal.alarm(0)
    if mol.GetNumConformers() == 0:
        return None
    return mol


def mol_to_entry(mol, mol_id, smiles):
    coords = [np.array(mol.GetConformer(i).GetPositions())
              for i in range(mol.GetNumConformers())]
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    return {'atoms': atoms, 'coordinates': coords, 'smi': smiles, 'mol_id': mol_id}


def worker(args):
    mol_id, smiles = args
    signal.signal(signal.SIGALRM, timeout_handler)
    try:
        mol = gen_conformation(smiles)
    except TimeoutError:
        return None
    if mol is None:
        return None
    return mol_to_entry(mol, mol_id, smiles)


def main():
    molecules = []
    with open(MOL_CSV) as f:
        for row in csv.DictReader(f):
            molecules.append((row['ID'], row['SMILES']))
    total = len(molecules)
    print(f'Molecules: {total}', flush=True)

    print('Generating conformations...', flush=True)
    entries = []
    id_map = []
    with Pool(processes=cpu_count()) as pool:
        for result in tqdm(pool.imap_unordered(worker, molecules, chunksize=200), total=total, desc='Conformations'):
            if result:
                entries.append(result)
                id_map.append(result['mol_id'])
    print(f'Valid: {len(entries)}/{total}')

    import lmdb
    env = lmdb.open(OUTPUT, subdir=False, readonly=False, lock=False,
                    readahead=False, meminit=False, map_size=1099511627776)
    with env.begin(write=True) as txn:
        for i, e in enumerate(entries):
            txn.put(str(i).encode('ascii'), pickle.dumps(e))
    env.close()
    print(f'Saved: {OUTPUT}')

    with open(ID_MAP, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Index', 'ID'])
        for idx, mid in enumerate(id_map):
            w.writerow([idx, mid])
    print(f'Saved: {ID_MAP}')


if __name__ == '__main__':
    main()