import csv
import os
import mygene

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(SCRIPT_DIR, 'ranked_target_summary.csv')
OUTPUT = os.path.join(SCRIPT_DIR, 'ranked_target_proteins.csv')

mg = mygene.MyGeneInfo()


def lookup(query):
    try:
        return mg.query(query, fields='symbol,name,uniprot', species='human', size=1)
    except:
        return None


def main():
    accessions = []
    gene_ids = []
    with open(INPUT) as f:
        reader = csv.DictReader(f)
        for row in reader:
            typ = row['Type']
            id_val = row['ID']
            count = int(row['Molecule_Count'])
            if typ == 'RefSeq':
                accessions.append((id_val, count))
            else:
                gene_ids.append((id_val, count))

    print(f'Mapping {len(gene_ids)} Gene IDs, {len(accessions)} Accessions...')

    gene_map = {}
    batch = [gid for gid, _ in gene_ids]
    if batch:
        results = mg.querymany(batch, scopes='entrezgene', fields='symbol,name,uniprot',
                               species='human', returnall=True)
        for item in results.get('out', []):
            q = item.get('query', '')
            symbol = item.get('symbol', '')
            name = item.get('name', '')
            uniprot = item.get('uniprot', {})
            if isinstance(uniprot, dict):
                uniprot = uniprot.get('Swiss-Prot', '')
            gene_map[q] = (symbol, name, uniprot)

    acc_map = {}
    batch = [acc for acc, _ in accessions]
    if batch:
        results = mg.querymany(batch, scopes='refseq,uniprot,accession',
                               fields='symbol,name,uniprot,entrezgene',
                               species='human', returnall=True)
        for item in results.get('out', []):
            q = item.get('query', '')
            symbol = item.get('symbol', '')
            name = item.get('name', '')
            uniprot = item.get('uniprot', {})
            if isinstance(uniprot, dict):
                uniprot = uniprot.get('Swiss-Prot', '')
            entrez = item.get('entrezgene', '')
            acc_map[q] = (symbol, name, uniprot, entrez)

    with open(OUTPUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Rank', 'Original_ID', 'Type', 'Gene_Symbol', 'Protein_Name',
                     'UniProt_ID', 'Entrez_GeneID', 'Molecule_Count'])

        rank = 0
        for acc, count in accessions:
            rank += 1
            info = acc_map.get(acc, ('', '', '', ''))
            w.writerow([rank, acc, 'RefSeq', info[0], info[1], info[2], info[3], count])

        for gid, count in gene_ids:
            rank += 1
            info = gene_map.get(gid, ('', '', ''))
            w.writerow([rank, gid, 'GeneID', info[0], info[1], info[2], gid, count])

    print(f'\nSaved: {OUTPUT}')


if __name__ == '__main__':
    main()