import csv
import json
import time
import os
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
SDF_CSV = os.path.join(PROJECT_DIR, 'molecule', 'T001_TargetMol_SMILES.csv')
OUTPUT = os.path.join(SCRIPT_DIR, 'properties.csv')
PROGRESS_FILE = os.path.join(SCRIPT_DIR, 'progress_properties.json')

REQUEST_DELAY = 0.25
MAX_RETRIES = 5
RETRY_DELAY = 15

COLUMNS = ['ID', 'SMILES', 'CID', 'MolecularFormula', 'MolecularWeight',
           'XLogP', 'TPSA', 'Complexity', 'HBondDonorCount',
           'HBondAcceptorCount', 'RotatableBondCount']


def query_pubchem(url, retries=MAX_RETRIES):
    for i in range(retries):
        try:
            req = urllib.request.Request(url)
            resp = urllib.request.urlopen(req, timeout=30)
            return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code == 503:
                print(f'    [Retry {i+1}/{retries}] Server busy (503), waiting {RETRY_DELAY}s...')
                time.sleep(RETRY_DELAY)
            else:
                print(f'    [Error] HTTP {e.code}')
                return None
        except Exception as e:
            print(f'    [Retry {i+1}/{retries}] Error: {e}')
            time.sleep(RETRY_DELAY)
    return None


def smiles_to_cid(smiles):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON'
    data = query_pubchem(url)
    if data and 'IdentifierList' in data:
        cids = data['IdentifierList']['CID']
        return cids[0] if cids else None
    return None


def get_properties(cid):
    props = ','.join(['MolecularFormula', 'MolecularWeight', 'XLogP',
                      'TPSA', 'Complexity', 'HBondDonorCount',
                      'HBondAcceptorCount', 'RotatableBondCount'])
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{props}/JSON'
    data = query_pubchem(url)
    if data and 'PropertyTable' in data:
        props_list = data['PropertyTable']['Properties']
        return props_list[0] if props_list else None
    return None


def load_progress():
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE) as f:
            return json.load(f)
    return {'processed': [], 'failed': []}


def save_progress(progress):
    with open(PROGRESS_FILE, 'w') as f:
        json.dump(progress, f)


def main():
    molecules = []
    with open(SDF_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            molecules.append(row)
    print(f'共 {len(molecules)} 个分子待处理')

    progress = load_progress()
    processed_set = set(progress['processed'])
    print(f'已处理 {len(processed_set)} 个，跳过已处理的')

    outfile = open(OUTPUT, 'a', newline='')
    writer = csv.DictWriter(outfile, fieldnames=COLUMNS)
    if os.path.getsize(OUTPUT) == 0:
        writer.writeheader()

    success_count = 0
    fail_count = 0

    try:
        for idx, mol in enumerate(molecules):
            mol_id = mol['ID']
            smiles = mol['SMILES']

            if mol_id in processed_set:
                continue

            print(f'[{idx+1}/{len(molecules)}] {mol_id}: {smiles[:50]}...')

            time.sleep(REQUEST_DELAY)
            cid = smiles_to_cid(smiles)
            if cid is None:
                print(f'  -> CID not found')
                progress['failed'].append({'ID': mol_id, 'SMILES': smiles, 'reason': 'CID not found'})
                fail_count += 1
                save_progress(progress)
                continue

            time.sleep(REQUEST_DELAY)
            props = get_properties(cid)
            if props is None:
                print(f'  -> No properties for CID {cid}')
                progress['failed'].append({'ID': mol_id, 'SMILES': smiles, 'CID': cid, 'reason': 'No properties'})
                fail_count += 1
                save_progress(progress)
                continue

            out_row = {'ID': mol_id, 'SMILES': smiles, 'CID': cid}
            for k in ['MolecularFormula', 'MolecularWeight', 'XLogP',
                       'TPSA', 'Complexity', 'HBondDonorCount',
                       'HBondAcceptorCount', 'RotatableBondCount']:
                out_row[k] = props.get(k, '')
            writer.writerow(out_row)

            progress['processed'].append(mol_id)
            processed_set.add(mol_id)
            success_count += 1

            if success_count % 100 == 0:
                save_progress(progress)
                outfile.flush()
                print(f'\n--- 进度保存: {success_count} 成功, {fail_count} 失败 ---\n')

    except KeyboardInterrupt:
        print('\n手动中断，保存进度...')
    finally:
        save_progress(progress)
        outfile.close()

    print(f'\n完成! 成功: {success_count}, 失败: {fail_count}')
    print(f'理化性质数据: {OUTPUT}')


if __name__ == '__main__':
    main()