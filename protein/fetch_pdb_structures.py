import csv
import os
import time
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, 'pdb_files')
APO_DIR = os.path.join(SCRIPT_DIR, 'pdb_apo')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(APO_DIR, exist_ok=True)

GENES = ['THRB', 'NR1H4', 'PPARA', 'PPARG', 'PPARD', 'ACACA',
         'FASN', 'AOC3', 'SCD1', 'SCD', 'NR1H3', 'MLXIPL',
         'DGAT2', 'XPR1', 'ELOVL6', 'GPAM', 'RXRA', 'THRSP',
         'RXRB', 'RXRG']

UNIPROT_MAP = {
    'THRB': 'P10828', 'NR1H4': 'Q96RI1', 'PPARA': 'Q07869',
    'PPARG': 'P37231', 'PPARD': 'Q03181', 'ACACA': 'Q13085',
    'FASN': 'P49327', 'AOC3': 'Q16853', 'SCD1': 'O00767',
    'SCD': 'O00767', 'NR1H3': 'Q13133', 'MLXIPL': 'Q9NP71',
    'DGAT2': 'Q96PD7', 'XPR1': 'Q9UBH6', 'ELOVL6': 'Q9H5J4',
    'GPAM': 'Q9HCL2', 'RXRA': 'P19793', 'THRSP': 'Q9UKT5',
    'RXRB': 'P28702', 'RXRG': 'P48443',
}

REQUEST_DELAY = 0.5


def search_pdb_by_uniprot(uniprot_id):
    url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{uniprot_id}'
    try:
        req = urllib.request.Request(url)
        resp = urllib.request.urlopen(req, timeout=30)
        data = json.loads(resp.read())
        return data
    except urllib.error.HTTPError as e:
        if e.code != 404:
            print(f'  HTTP {e.code}')
        return None
    except Exception as e:
        print(f'  Error: {e}')
        return None


import json


def strip_to_apo(pdb_content):
    lines = []
    for line in pdb_content.splitlines():
        if line.startswith('ATOM'):
            lines.append(line)
        elif line.startswith('TER'):
            lines.append(line)
        elif line.startswith('END'):
            lines.append(line)
    return '\n'.join(lines)


def search_pdb_by_uniprot_alt(uniprot_id):
    """RCSB Search API: find PDB entries by UniProt ID"""
    url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 100
            },
            "sort": [
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
            ],
            "scoring_strategy": "combined"
        }
    }
    try:
        req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                     headers={'Content-Type': 'application/json'})
        resp = urllib.request.urlopen(req, timeout=30)
        return json.loads(resp.read())
    except Exception as e:
        print(f'  Search error: {e}')
        return None


def get_entry_info(pdb_id):
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
    try:
        req = urllib.request.Request(url)
        resp = urllib.request.urlopen(req, timeout=30)
        return json.loads(resp.read())
    except Exception:
        return None


def get_best_pdb(uniprot_id):
    data = search_pdb_by_uniprot_alt(uniprot_id)
    if not data:
        return None

    result_list = data.get('result_set', [])
    if not result_list:
        return None

    candidates = []
    for item in result_list:
        pdb_id = item.get('identifier', '')
        if not pdb_id:
            continue
        time.sleep(0.2)
        info = get_entry_info(pdb_id)
        if not info:
            continue

        resolution = info.get('rcsb_entry_info', {}).get('resolution_combined', [99.0])
        resolution = resolution[0] if resolution else 99.0
        title = info.get('struct', {}).get('title', '')
        keywords = info.get('struct_keywords', {}).get('pdbx_keywords', '')
        deposit_date = info.get('rcsb_accession_info', {}).get('deposit_date', '')
        structure_method = info.get('exptl', [{}])[0].get('method', '') if info.get('exptl') else ''

        has_mutation = any(kw in (title + ' ' + keywords).lower()
                           for kw in ['mutant', 'mutation', 'variant'])

        candidates.append({
            'pdb_id': pdb_id.upper(),
            'resolution': resolution,
            'has_mutation': has_mutation,
            'title': title,
            'method': structure_method,
            'deposit_date': deposit_date,
        })

    if not candidates:
        return None

    candidates.sort(key=lambda x: (
        x['has_mutation'],
        x['resolution'] if isinstance(x['resolution'], (int, float)) and x['resolution'] < 50 else 99.0,
        x['deposit_date'] or '9999',
    ))
    return candidates[0]


def download_pdb(pdb_id, out_path):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    try:
        req = urllib.request.Request(url)
        resp = urllib.request.urlopen(req, timeout=60)
        content = resp.read()
        with open(out_path, 'wb') as f:
            f.write(content)
        return content
    except Exception as e:
        print(f'  Download failed: {e}')
        return None


def main():
    results_log = []

    for i, gene in enumerate(GENES, 1):
        print(f'\n[{i}/{len(GENES)}] {gene} (UniProt: {UNIPROT_MAP[gene]})...')

        time.sleep(REQUEST_DELAY)
        best = get_best_pdb(UNIPROT_MAP[gene])

        if not best:
            print(f'  No structures found')
            results_log.append({'gene': gene, 'uniprot': UNIPROT_MAP[gene],
                                'pdb_id': '', 'resolution': '', 'method': '',
                                'title': '', 'status': 'no_structure'})
            continue

        pdb_id = best['pdb_id']
        resolution = best['resolution']
        title = best['title'][:80]
        print(f'  Best: {pdb_id}  resolution={resolution:.2f}A  {best["method"]}')
        print(f'  Title: {title}')

        out_path = os.path.join(OUT_DIR, f'{gene}_{pdb_id}.pdb')
        content = download_pdb(pdb_id, out_path)
        if content:
            apo_path = os.path.join(APO_DIR, f'{gene}_apo.pdb')
            apo_content = strip_to_apo(content.decode())
            with open(apo_path, 'w') as f:
                f.write(apo_content)
            atom_count = apo_content.count('\nATOM') + 1
            print(f'  Saved apo: {apo_path} ({atom_count} atoms)')
        status = 'downloaded' if content else 'download_failed'

        results_log.append({
            'gene': gene,
            'uniprot': UNIPROT_MAP[gene],
            'pdb_id': pdb_id,
            'resolution': f'{resolution:.2f}' if isinstance(resolution, (int, float)) and resolution < 50 else 'N/A',
            'method': best['method'],
            'title': title,
            'status': status,
        })

    log_path = os.path.join(SCRIPT_DIR, 'pdb_download_log.csv')
    with open(log_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gene', 'uniprot', 'pdb_id', 'resolution',
                                           'method', 'title', 'status'])
        w.writeheader()
        w.writerows(results_log)

    print(f'\n{"="*60}')
    print(f'Log: {log_path}')
    print(f'PDB files: {OUT_DIR}/')
    print()
    for r in results_log:
        icon = 'OK' if r['status'] == 'downloaded' else '--'
        print(f'  [{icon}] {r["gene"]:8s} {r["pdb_id"]:6s} {r["resolution"]:6s}  {r["method"][:20]}')


if __name__ == '__main__':
    main()