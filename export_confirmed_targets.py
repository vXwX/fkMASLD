import csv

our_symbols = {
    'THRB': '7068', 'NR1H4': '9971', 'PPARA': '5465',
    'PPARG': '5468', 'PPARD': '5467', 'ACACA': '31',
    'FASN': '2194', 'AOC3': '26', 'SCD1': '6319',
    'SCD': '6319', 'NR1H3': '10062', 'MLXIPL': '51085',
    'XPR1': '84333', 'GPAM': '57678', 'RXRA': '6256',
    'THRSP': '57103', 'RXRB': '6257', 'RXRG': '6258',
}

gene_id_to_symbol = {v: k for k, v in our_symbols.items()}
our_ids = set(our_symbols.values())

matches = []
with open('/home/gaoxiangxu/Competition/AI4S_202607/top_candidates_unique.csv') as f:
    for row in csv.DictReader(f):
        targets = row.get('Active_Targets', '').strip()
        if not targets:
            continue
        tids = set(t.strip() for t in targets.split(';') if t.strip())
        overlap = tids & our_ids
        if not overlap:
            continue
        known_genes = ';'.join(gene_id_to_symbol.get(eid, eid) for eid in sorted(overlap))
        matches.append({
            'ID': row['ID'],
            'SMILES': row['SMILES'],
            'Best_Gene': row['Best_Gene'],
            'DrugCLIP_Score': row['DrugCLIP_Score'],
            'Activity_Level': row['Activity_Level'],
            'Toxicity_Level': row['Toxicity_Level'],
            'Known_Target_Genes': known_genes,
            'Active_Targets': targets,
        })

matches.sort(key=lambda x: -float(x['DrugCLIP_Score']))

with open('/home/gaoxiangxu/Competition/AI4S_202607/top_candidates_confirmed_targets.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=[
        'ID', 'SMILES', 'Best_Gene', 'DrugCLIP_Score',
        'Activity_Level', 'Toxicity_Level',
        'Known_Target_Genes', 'Active_Targets'])
    w.writeheader()
    w.writerows(matches)

print(f'Saved: {len(matches)} molecules')