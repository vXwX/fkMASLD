"""
基于 STRING API 查询 MASLD 靶点基因的蛋白-蛋白互作网络，扩展相关蛋白。
用法: conda run -n tmp python protein/relpro/find_related_proteins.py
"""
import csv
import json
import os
import time
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = SCRIPT_DIR

SEED_GENES = ['THRB', 'NR1H4', 'PPARA', 'PPARG', 'PPARD',
              'ACACA', 'FASN', 'AOC3', 'SCD1']

SPECIES = 9606
REQUEST_DELAY = 1.0


def query_string(method, params):
    url = f'https://string-db.org/api/json/{method}'
    data = urllib.parse.urlencode(params).encode()
    for retry in range(5):
        try:
            req = urllib.request.Request(url, data=data)
            resp = urllib.request.urlopen(req, timeout=60)
            return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            print(f'  [Retry {retry+1}/5] HTTP {e.code}, waiting...')
            time.sleep(5)
        except Exception as e:
            print(f'  [Retry {retry+1}/5] Error: {e}')
            time.sleep(5)
    return None


def main():
    import urllib.parse

    print(f'Seed genes: {SEED_GENES}')
    params = {
        'identifiers': '%0d'.join(SEED_GENES),
        'species': SPECIES,
    }

    # Step 1: get protein IDs
    print('\n=== Getting STRING IDs ===')
    result = query_string('get_string_ids', {**params, 'limit': 1})
    if not result:
        print('Failed to get STRING IDs')
        return
    string_ids = [item['stringId'] for item in result]
    print(f'Got {len(string_ids)} STRING IDs')

    # Step 2: get interaction network
    print('\n=== Getting interaction network ===')
    net_params = {
        'identifiers': '%0d'.join(string_ids),
        'species': SPECIES,
        'required_score': 400,
        'add_nodes': 10,
    }
    network = query_string('network', net_params)
    if not network:
        print('Failed to get network')
        return

    print(f'Got {len(network)} interactions')
    with open(os.path.join(OUT_DIR, 'string_network.json'), 'w') as f:
        json.dump(network, f, indent=2)

    # Step 3: get all node info
    all_pref_ids = set()
    for edge in network:
        all_pref_ids.add(edge['preferredName_A'])
        all_pref_ids.add(edge['preferredName_B'])

    print(f'\n=== Extended proteins ({len(all_pref_ids)} total) ===')
    params_info = {
        'identifiers': '%0d'.join(all_pref_ids),
        'species': SPECIES,
    }
    info = query_string('get_string_ids', params_info)
    if info:
        with open(os.path.join(OUT_DIR, 'string_protein_info.json'), 'w') as f:
            json.dump(info, f, indent=2)

    # Step 4: build edge list for extended network
    edges = []
    nodes = {}
    for pref_id in all_pref_ids:
        nodes[pref_id] = {'is_seed': pref_id in SEED_GENES, 'interactors': 0}
    for edge in network:
        a, b = edge['preferredName_A'], edge['preferredName_B']
        score = edge['score']
        edges.append({'protein_A': a, 'protein_B': b, 'score': score})
        if a in nodes:
            nodes[a]['interactors'] += 1
        if b in nodes:
            nodes[b]['interactors'] += 1

    # Sort non-seed proteins by interaction count
    extensions = sorted(
        [(g, info) for g, info in nodes.items() if not info['is_seed']],
        key=lambda x: -x[1]['interactors']
    )

    # Save edge list
    with open(os.path.join(OUT_DIR, 'string_edges.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Protein_A', 'Protein_B', 'Score'])
        for e in sorted(edges, key=lambda x: -x['score']):
            w.writerow([e['protein_A'], e['protein_B'], e['score']])

    # Save extended proteins
    with open(os.path.join(OUT_DIR, 'extended_proteins.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Gene', 'Is_Seed', 'Interactor_Count'])
        for g, info in sorted(nodes.items(), key=lambda x: (-x[1]['is_seed'], -x[1]['interactors'])):
            w.writerow([g, 'Yes' if info['is_seed'] else 'No', info['interactors']])

    print(f'\nExtended proteins (non-seed, by interaction count):')
    for g, info in extensions[:20]:
        print(f'  {g:20s}  interactors: {info["interactors"]}')

    print(f'\nOutput files:')
    print(f'  {OUT_DIR}/string_network.json')
    print(f'  {OUT_DIR}/string_protein_info.json')
    print(f'  {OUT_DIR}/string_edges.csv')
    print(f'  {OUT_DIR}/extended_proteins.csv')


if __name__ == '__main__':
    main()