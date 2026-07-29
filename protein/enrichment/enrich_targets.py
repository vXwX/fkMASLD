"""
MASLD 参考药物靶点富集分析
用法: conda run -n tmp python protein/enrichment/enrich_targets.py
"""
import csv
import os

import gseapy as gp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = SCRIPT_DIR

TARGET_GENES = {
    'THR-beta': 'THRB',
    'FXR': 'NR1H4',
    'PPAR-alpha/gamma/delta': ['PPARA', 'PPARG', 'PPARD'],
    'ACC': 'ACACA',
    'FASN': 'FASN',
    'AOC3/VAP-1': 'AOC3',
    'SCD1': 'SCD1',
}


def main():
    genes = set()
    target_info = []
    for target, gene in TARGET_GENES.items():
        if isinstance(gene, list):
            genes.update(gene)
            for g in gene:
                target_info.append((target, g))
        else:
            genes.add(gene)
            target_info.append((target, gene))

    genes = list(genes)
    print(f'Input genes ({len(genes)}): {genes}')
    print()

    # Enrichr KEGG
    print('=== KEGG Pathway ===')
    try:
        enr_kegg = gp.enrichr(
            gene_list=genes,
            gene_sets='KEGG_2021_Human',
            organism='human',
            outdir=None,
            no_plot=True,
        )
        results = enr_kegg.results
        results.to_csv(os.path.join(OUT_DIR, 'enrichment_kegg.csv'), index=False)
        if not results.empty:
            for _, r in results.head(10).iterrows():
                print(f"  {r['Term']:50s} p={r['P-value']:.2e}  overlap={r['Overlap']}")
        else:
            print('  No significant enrichment.')
    except Exception as e:
        print(f'  Error: {e}')
    print()

    # Enrichr GO Biological Process
    print('=== GO Biological Process ===')
    try:
        enr_go = gp.enrichr(
            gene_list=genes,
            gene_sets='GO_Biological_Process_2023',
            organism='human',
            outdir=None,
            no_plot=True,
        )
        results = enr_go.results
        results.to_csv(os.path.join(OUT_DIR, 'enrichment_go_bp.csv'), index=False)
        if not results.empty:
            for _, r in results.head(10).iterrows():
                print(f"  {r['Term']:50s} p={r['P-value']:.2e}  overlap={r['Overlap']}")
        else:
            print('  No significant enrichment.')
    except Exception as e:
        print(f'  Error: {e}')
    print()

    # Enrichr Reactome
    print('=== Reactome ===')
    try:
        enr_react = gp.enrichr(
            gene_list=genes,
            gene_sets='Reactome_2022',
            organism='human',
            outdir=None,
            no_plot=True,
        )
        results = enr_react.results
        results.to_csv(os.path.join(OUT_DIR, 'enrichment_reactome.csv'), index=False)
        if not results.empty:
            for _, r in results.head(10).iterrows():
                print(f"  {r['Term']:50s} p={r['P-value']:.2e}  overlap={r['Overlap']}")
        else:
            print('  No significant enrichment.')
    except Exception as e:
        print(f'  Error: {e}')
    print()

    # 保存靶点-基因映射
    with open(os.path.join(OUT_DIR, 'target_gene_map.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Target', 'Gene'])
        w.writerows(target_info)

    print('Done. Files saved to:', OUT_DIR)


if __name__ == '__main__':
    main()