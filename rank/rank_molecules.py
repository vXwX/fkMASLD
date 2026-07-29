import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
ACTIVITY_CSV = os.path.join(PROJECT_DIR, 'data', 'activity.csv')
TOXICITY_CSV = os.path.join(PROJECT_DIR, 'data', 'toxicity.csv')
OUTPUT = os.path.join(SCRIPT_DIR, 'ranked_molecules.csv')

HIGH_TOX_TARGETS = {'hERG', 'QT prolongation', 'cardiotox', 'hepatotox',
                    'nephrotox', 'neurotox', 'genotox', 'death', 'lethal',
                    'carcinogen', 'mutagen', 'teratogen'}


def count_active_by_target(csv_path, outcome_col, target_col, accession_col):
    result = {}
    try:
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                mol_id = row.get('ID', '')
                outcome = row.get(outcome_col, '')
                target = row.get(target_col, '')
                accession = row.get(accession_col, '')
                if mol_id not in result:
                    result[mol_id] = {'active_count': 0, 'total_count': 0,
                                      'active_targets': set(), 'active_accessions': set()}
                result[mol_id]['total_count'] += 1
                if outcome.lower() == 'active':
                    result[mol_id]['active_count'] += 1
                    if target:
                        result[mol_id]['active_targets'].add(target)
                    if accession:
                        result[mol_id]['active_accessions'].add(accession)
    except FileNotFoundError:
        pass
    return result


def classify_activity(mol_id, activity_data):
    data = activity_data.get(mol_id)
    if not data or data['total_count'] == 0:
        return '低活性'
    active_ratio = data['active_count'] / data['total_count']
    if active_ratio > 0.3 and len(data['active_targets']) >= 2:
        return '高活性'
    if active_ratio > 0.05 or len(data['active_targets']) >= 1:
        return '活性'
    return '低活性'


def classify_toxicity(mol_id, toxicity_data):
    data = toxicity_data.get(mol_id)
    if not data or data['total_count'] == 0:
        return '低毒性'
    has_high_risk = any(
        kw in name.lower()
        for name in data.get('assay_names', set())
        for kw in HIGH_TOX_TARGETS
    )
    if has_high_risk:
        return '高毒性'
    if data['active_count'] > 0:
        return '毒性'
    return '低毒性'


def collect_toxicity_details(csv_path):
    result = {}
    try:
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                mol_id = row.get('ID', '')
                outcome = row.get('Activity_Outcome', '')
                assay_name = row.get('Assay_Name', '')
                if mol_id not in result:
                    result[mol_id] = {'active_count': 0, 'total_count': 0,
                                      'assay_names': set(), 'targets': set()}
                result[mol_id]['total_count'] += 1
                result[mol_id]['assay_names'].add(assay_name)
                if outcome.lower() == 'active':
                    result[mol_id]['active_count'] += 1
                    target = row.get('Target_GeneID', '') or row.get('Target_Accession', '')
                    if target:
                        result[mol_id]['targets'].add(target)
    except FileNotFoundError:
        pass
    return result


def main():
    activity_data = count_active_by_target(ACTIVITY_CSV, 'Activity_Outcome', 'Target_GeneID', 'Target_Accession')
    toxicity_data = collect_toxicity_details(TOXICITY_CSV)

    all_ids = set(activity_data.keys()) | set(toxicity_data.keys())
    if not all_ids:
        print('No data found. Make sure activity.csv and toxicity.csv exist.')
        return

    ranked = []
    for mol_id in all_ids:
        act_class = classify_activity(mol_id, activity_data)
        tox_class = classify_toxicity(mol_id, toxicity_data)
        act_info = activity_data.get(mol_id, {})
        tox_info = toxicity_data.get(mol_id, {})

        rank_score = {
            ('高活性', '低毒性'): 1,
            ('高活性', '毒性'): 2,
            ('高活性', '高毒性'): 3,
            ('活性', '低毒性'): 4,
            ('活性', '毒性'): 5,
            ('活性', '高毒性'): 6,
            ('低活性', '低毒性'): 7,
            ('低活性', '毒性'): 8,
            ('低活性', '高毒性'): 9,
        }.get((act_class, tox_class), 9)

        ranked.append({
            'ID': mol_id,
            'Activity_Level': act_class,
            'Toxicity_Level': tox_class,
            'Rank_Score': rank_score,
            'Active_Count': act_info.get('active_count', 0),
            'Total_Assay_Count': act_info.get('total_count', 0) + tox_info.get('total_count', 0),
            'Active_Targets': ';'.join(act_info.get('active_targets', set())),
            'Active_Accessions': ';'.join(act_info.get('active_accessions', set())),
            'Toxicity_Targets': ';'.join(tox_info.get('targets', set())),
        })

    ranked.sort(key=lambda x: (x['Rank_Score'], -x['Active_Count']))

    with open(OUTPUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'ID', 'Activity_Level', 'Toxicity_Level', 'Rank_Score',
            'Active_Count', 'Total_Assay_Count',
            'Active_Targets', 'Active_Accessions', 'Toxicity_Targets'
        ])
        writer.writeheader()
        writer.writerows(ranked)

    print(f'Done: {OUTPUT}')
    print(f'Total molecules: {len(ranked)}')
    print('\nDistribution:')
    dist = {}
    for r in ranked:
        key = f"{r['Activity_Level']} x {r['Toxicity_Level']}"
        dist[key] = dist.get(key, 0) + 1
    for k, v in sorted(dist.items(), key=lambda x: list(dist.keys()).index(x[0])):
        print(f'  {k}: {v}')


if __name__ == '__main__':
    main()