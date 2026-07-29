import csv
import os
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(SCRIPT_DIR, 'ranked_molecules.csv')
OUTPUT = os.path.join(SCRIPT_DIR, 'ranked_target_summary.csv')


def main():
    accession_counter = Counter()
    gene_counter = Counter()

    with open(INPUT) as f:
        for row in csv.DictReader(f):
            targets = row.get('Active_Targets', '').strip()
            accs = row.get('Active_Accessions', '').strip()
            if targets:
                for t in targets.split(';'):
                    t = t.strip()
                    if t:
                        gene_counter[t] += 1
            if accs:
                for a in accs.split(';'):
                    a = a.strip()
                    if a:
                        accession_counter[a] += 1

    print(f'Unique gene IDs: {len(gene_counter)}')
    print(f'Unique RefSeq accessions: {len(accession_counter)}')

    with open(OUTPUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Rank', 'Type', 'ID', 'Molecule_Count'])
        for rank, (acc, count) in enumerate(accession_counter.most_common(), 1):
            w.writerow([rank, 'RefSeq', acc, count])
        for rank, (gene, count) in enumerate(gene_counter.most_common(), 1):
            w.writerow([rank, 'GeneID', gene, count])

    print(f'\nTop 20 RefSeq accessions:')
    print(f'{"Rank":<6} {"Accession":<16} {"Count":<8}')
    print('-' * 32)
    for rank, (acc, count) in enumerate(accession_counter.most_common(20), 1):
        print(f'{rank:<6} {acc:<16} {count:<8}')

    print(f'\nTop 20 Gene IDs:')
    print(f'{"Rank":<6} {"GeneID":<16} {"Count":<8}')
    print('-' * 32)
    for rank, (gene, count) in enumerate(gene_counter.most_common(20), 1):
        print(f'{rank:<6} {gene:<16} {count:<8}')

    print(f'\nSaved: {OUTPUT}')


if __name__ == '__main__':
    main()