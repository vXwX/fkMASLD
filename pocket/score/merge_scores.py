import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCORES_DIR = os.path.join(SCRIPT_DIR, 'results')
OUTPUT = os.path.join(SCRIPT_DIR, 'scores', 'all_scores_summary.csv')


def main():
    all_rows = []
    genes = sorted([d for d in os.listdir(SCORES_DIR) if os.path.isdir(os.path.join(SCORES_DIR, d))])

    for gene in genes:
        csv_path = os.path.join(SCORES_DIR, gene, 'scores.csv')
        if not os.path.exists(csv_path):
            continue
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                all_rows.append({
                    'Gene': gene,
                    'Rank': row['Rank'],
                    'ID': row['ID'],
                    'SMILES': row['SMILES'],
                    'Score': row['Score'],
                })

    all_rows.sort(key=lambda x: (-float(x['Score']), x['Gene'], int(x['Rank'])))

    with open(OUTPUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['Gene', 'Rank', 'ID', 'SMILES', 'Score'])
        w.writeheader()
        w.writerows(all_rows)

    print(f'Total: {len(all_rows)} entries, {len(genes)} genes')
    print(f'Saved: {OUTPUT}')
    print(f'\nTop 10 across all proteins:')
    print(f'{"Rank":<6} {"Gene":<8} {"ID":<10} {"Score":<8} {"SMILES":<40}')
    print('-' * 75)
    for i, r in enumerate(all_rows[:10], 1):
        print(f'{i:<6} {r["Gene"]:<8} {r["ID"]:<10} {r["Score"]:<8} {r["SMILES"][:38]}')


if __name__ == '__main__':
    main()