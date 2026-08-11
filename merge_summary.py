import os

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_CSV = os.path.join(SCRIPT_DIR, 'T001_physchem_pass_with_sim.csv')
DRUGCLIP_CSV = os.path.join(SCRIPT_DIR, 'pocket', 'score', 'results', 'all_scores_summary.csv')
GIN_CSV = os.path.join(SCRIPT_DIR, 'gin', 'results', 'all_scores.csv')
OUT_CSV = os.path.join(SCRIPT_DIR, 'Top6000.csv')


def main():
    sim = pd.read_csv(SIM_CSV)

    dc = pd.read_csv(DRUGCLIP_CSV)
    dc_best = dc.loc[dc.groupby('ID')['Score'].idxmax(), ['ID', 'Gene', 'Score']]
    dc_best = dc_best.rename(columns={'Gene': 'DrugCLIP_Target', 'Score': 'DrugCLIP_Score'})

    gin = pd.read_csv(GIN_CSV)
    gin['active_score'] = gin['P_act'] * gin['STRENGTH']
    gin['toxicity_score'] = gin['P_tox_1'] + gin['P_tox_2']
    gin['recommend_score'] = gin['active_score'] * (1.0 - gin['toxicity_score'])
    gin = gin[['ID', 'active_score', 'toxicity_score', 'recommend_score']]

    out = sim[['name', 'smiles', 'max_tc']].merge(dc_best, how='left', left_on='name', right_on='ID')
    out = out.merge(gin, how='left', left_on='name', right_on='ID')
    out = out.drop(columns=['ID_x', 'ID_y']).rename(columns={'name': 'name'})

    out = out.rename(columns={
        'name': 'name',
        'smiles': 'smiles',
        'max_tc': 'max_tc',
        'DrugCLIP_Score': 'DrugCLIP_Score',
        'DrugCLIP_Target': 'DrugCLIP_Target',
        'active_score': 'active_score',
        'toxicity_score': 'toxicity_score',
        'recommend_score': 'recommend_score',
    })

    out.to_csv(OUT_CSV, index=False, float_format='%.6f')
    print(f'Saved: {OUT_CSV} ({len(out)} rows)')
    print(out.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
