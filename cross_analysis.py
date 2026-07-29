import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = SCRIPT_DIR
RANK_CSV = os.path.join(PROJECT_DIR, 'rank', 'ranked_molecules.csv')
SCORE_CSV = os.path.join(PROJECT_DIR, 'pocket', 'score', 'results', 'all_scores_summary.csv')
OUTPUT_ALL = os.path.join(PROJECT_DIR, 'top_candidates.csv')
OUTPUT_UNIQ = os.path.join(PROJECT_DIR, 'top_candidates_unique.csv')


def main():
    ranked = {}
    with open(RANK_CSV) as f:
        for row in csv.DictReader(f):
            ranked[row['ID']] = row

    best = {}
    with open(SCORE_CSV) as f:
        for row in csv.DictReader(f):
            mid = row['ID']
            score = float(row['Score'])
            r = ranked.get(mid)
            if not r:
                continue
            act = r['Activity_Level']
            tox = r['Toxicity_Level']
            if act not in ('高活性', '活性') or tox not in ('低毒性', '毒性'):
                continue
            if mid not in best or score > best[mid]['DrugCLIP_Score']:
                best[mid] = {
                    'ID': mid, 'SMILES': row['SMILES'],
                    'Best_Gene': row['Gene'], 'DrugCLIP_Score': score,
                    'Activity_Level': act, 'Toxicity_Level': tox,
                    'Active_Targets': r.get('Active_Targets', ''),
                    'Active_Accessions': r.get('Active_Accessions', ''),
                }

    SORT_ORDER = {
        ('高活性', '低毒性'): 0,
        ('高活性', '毒性'): 1,
        ('活性', '低毒性'): 2,
        ('活性', '毒性'): 3,
    }

    def sort_key(x):
        return (SORT_ORDER.get((x['Activity_Level'], x['Toxicity_Level']), 9), -x['DrugCLIP_Score'])

    uniq = sorted(best.values(), key=sort_key)
    cols = ['ID', 'SMILES', 'Best_Gene', 'DrugCLIP_Score',
            'Activity_Level', 'Toxicity_Level',
            'Active_Targets', 'Active_Accessions']

    with open(OUTPUT_UNIQ, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(uniq)

    with open(OUTPUT_ALL, 'w', newline='') as f:
        cols_all = cols[:3] + ['Gene'] + cols[3:]
        w = csv.DictWriter(f, fieldnames=cols_all)
        w.writeheader()
        all_rows = []
        with open(SCORE_CSV) as sf:
            for row in csv.DictReader(sf):
                mid = row['ID']
                r = ranked.get(mid)
                if not r:
                    continue
                act = r['Activity_Level']
                tox = r['Toxicity_Level']
                if act in ('高活性', '活性') and tox in ('低毒性', '毒性'):
                    all_rows.append({
                        'ID': mid, 'SMILES': row['SMILES'], 'Gene': row['Gene'],
                        'Best_Gene': best.get(mid, {}).get('Best_Gene', ''),
                        'DrugCLIP_Score': row['Score'],
                        'Activity_Level': act, 'Toxicity_Level': tox,
                        'Active_Targets': r.get('Active_Targets', ''),
                        'Active_Accessions': r.get('Active_Accessions', ''),
                    })
        all_rows.sort(key=lambda x: (
            SORT_ORDER.get((x['Activity_Level'], x['Toxicity_Level']), 9),
            -float(x['DrugCLIP_Score']),
        ))
        w.writerows(all_rows)

    print(f'All pairs: {sum(1 for _ in open(OUTPUT_ALL)) - 1}')
    print(f'Unique: {len(uniq)}')
    print(f'Saved: {OUTPUT_ALL}')
    print(f'Saved: {OUTPUT_UNIQ}')
    print()
    print(f'Top 20 unique:')
    print(f'{"Rank":<6} {"ID":<12} {"Best_Gene":<10} {"Score":<8} {"Act":<8} {"Tox":<8} {"SMILES":<40}')
    print('-' * 97)
    for i, c in enumerate(uniq[:20], 1):
        print(f'{i:<6} {c["ID"]:<12} {c["Best_Gene"]:<10} {c["DrugCLIP_Score"]:<8.4f} {c["Activity_Level"]:<8} {c["Toxicity_Level"]:<8} {c["SMILES"][:38]}')

    dist = {}
    for c in uniq:
        key = f"{c['Activity_Level']} x {c['Toxicity_Level']}"
        dist[key] = dist.get(key, 0) + 1
    print(f'\nDistribution:')
    for k, v in sorted(dist.items()):
        print(f'  {k}: {v}')


if __name__ == '__main__':
    main()