"""
PubChem 活性和毒性数据获取脚本
用法: conda run -n tmp python data/fetch_pubchem.py
"""
import csv
import json
import time
import os
import urllib.request
import urllib.error
from datetime import datetime

# ========== 配置 ==========
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
SDF_CSV = os.path.join(PROJECT_DIR, 'molecule', 'T001_TargetMol_SMILES.csv')
OUT_DIR = SCRIPT_DIR
ACTIVITY_CSV = os.path.join(OUT_DIR, 'activity.csv')
TOXICITY_CSV = os.path.join(OUT_DIR, 'toxicity.csv')
PROGRESS_FILE = os.path.join(OUT_DIR, 'progress.json')

REQUEST_DELAY = 0.25  # 每次请求间隔（秒），5次/秒
MAX_RETRIES = 5
RETRY_DELAY = 15  # 重试等待时间

# 毒性相关关键词（出现在 Assay Name 中）
TOXICITY_KEYWORDS = [
    'toxic', 'tox', 'death', 'lethal', 'ld50', 'lc50',
    'cell death', 'viability', 'necrosis', 'apoptosis',
    'carcinogen', 'mutagen', 'ames', 'hepatotox', 'nephrotox',
    'cardiotox', 'neurotox', 'cytotox', 'cytotoxic',
    'mitochondrial toxicity', 'genotox', 'teratogen',
    'acute toxicity', 'chronic toxicity', 'subchronic',
    'organ toxicity', 'liver toxicity', 'kidney toxicity',
    'heart toxicity', 'brain toxicity', 'lung toxicity',
    'sperm', 'embryo', 'developmental', 'reproductive',
    'hERG', 'channel block', 'QT prolongation',
]

# BioAssay 列名
ASSAY_COLUMNS = [
    'ID', 'SMILES', 'CID', 'AID', 'Activity_Outcome',
    'Target_Accession', 'Target_GeneID', 'Activity_Value_uM',
    'Activity_Name', 'Assay_Name', 'Assay_Type', 'PubMed_ID'
]


def query_pubchem(url, retries=MAX_RETRIES):
    """带重试的 PubChem REST API 查询"""
    for i in range(retries):
        try:
            req = urllib.request.Request(url)
            resp = urllib.request.urlopen(req, timeout=30)
            return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code == 503:
                print(f'    [Retry {i+1}/{retries}] Server busy (503), waiting {RETRY_DELAY}s...')
                time.sleep(RETRY_DELAY)
            else:
                print(f'    [Error] HTTP {e.code}: {e.reason}')
                return None
        except Exception as e:
            print(f'    [Retry {i+1}/{retries}] Error: {e}')
            time.sleep(RETRY_DELAY)
    return None


def smiles_to_cid(smiles):
    """SMILES → CID"""
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON'
    data = query_pubchem(url)
    if data and 'IdentifierList' in data:
        cids = data['IdentifierList']['CID']
        return cids[0] if cids else None
    return None


def get_bioassay_summary(cid):
    """CID → BioAssay Summary"""
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON'
    data = query_pubchem(url)
    if data and 'Table' in data:
        return data['Table']
    return None


def is_toxicity_assay(assay_name):
    """判断是否为毒性相关实验"""
    name_lower = assay_name.lower()
    return any(kw in name_lower for kw in TOXICITY_KEYWORDS)


def parse_assay_row(cells, columns):
    """解析一行 Assay 数据"""
    row = {}
    for i, col in enumerate(columns):
        row[col] = cells[i] if i < len(cells) else ''
    return row


def load_progress():
    """加载进度"""
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE) as f:
            return json.load(f)
    return {'processed': [], 'failed': []}


def save_progress(progress):
    """保存进度"""
    with open(PROGRESS_FILE, 'w') as f:
        json.dump(progress, f)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 读取 SMILES 列表
    molecules = []
    with open(SDF_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            molecules.append(row)
    print(f'共 {len(molecules)} 个分子待处理')

    # 加载进度
    progress = load_progress()
    processed_set = set(progress['processed'])
    print(f'已处理 {len(processed_set)} 个，跳过已处理的')

    # 打开输出文件（追加模式）
    activity_file = open(ACTIVITY_CSV, 'a', newline='')
    toxicity_file = open(TOXICITY_CSV, 'a', newline='')
    act_writer = csv.DictWriter(activity_file, fieldnames=ASSAY_COLUMNS)
    tox_writer = csv.DictWriter(toxicity_file, fieldnames=ASSAY_COLUMNS)

    # 如果文件为空则写入表头
    if os.path.getsize(ACTIVITY_CSV) == 0:
        act_writer.writeheader()
    if os.path.getsize(TOXICITY_CSV) == 0:
        tox_writer.writeheader()

    success_count = 0
    fail_count = 0

    try:
        for idx, mol in enumerate(molecules):
            mol_id = mol['ID']
            smiles = mol['SMILES']

            # 跳过已处理的
            if mol_id in processed_set:
                continue

            print(f'[{idx+1}/{len(molecules)}] {mol_id}: {smiles[:50]}...')

            # Step 1: SMILES → CID
            time.sleep(REQUEST_DELAY)
            cid = smiles_to_cid(smiles)
            if cid is None:
                print(f'  → CID not found, skipping')
                progress['failed'].append({'ID': mol_id, 'SMILES': smiles, 'reason': 'CID not found'})
                fail_count += 1
                save_progress(progress)
                continue

            print(f'  → CID: {cid}')

            # Step 2: CID → BioAssay Summary
            time.sleep(REQUEST_DELAY)
            table = get_bioassay_summary(cid)
            if table is None:
                print(f'  → No assay data, skipping')
                progress['failed'].append({'ID': mol_id, 'SMILES': smiles, 'CID': cid, 'reason': 'No assay data'})
                fail_count += 1
                save_progress(progress)
                continue

            rows = table.get('Row', [])
            cols = table.get('Columns', {}).get('Column', [])
            print(f'  → {len(rows)} assay records')

            # Step 3: 解析并分类
            activity_count = 0
            toxicity_count = 0

            for row in rows:
                cells = row.get('Cell', [])
                parsed = parse_assay_row(cells, cols)

                # 构建输出行
                out_row = {
                    'ID': mol_id,
                    'SMILES': smiles,
                    'CID': cid,
                    'AID': parsed.get('AID', ''),
                    'Activity_Outcome': parsed.get('Activity Outcome', ''),
                    'Target_Accession': parsed.get('Target Accession', ''),
                    'Target_GeneID': parsed.get('Target GeneID', ''),
                    'Activity_Value_uM': parsed.get('Activity Value [uM]', ''),
                    'Activity_Name': parsed.get('Activity Name', ''),
                    'Assay_Name': parsed.get('Assay Name', ''),
                    'Assay_Type': parsed.get('Assay Type', ''),
                    'PubMed_ID': parsed.get('PubMed ID', ''),
                }

                # 写入活性数据
                act_writer.writerow(out_row)
                activity_count += 1

                # 判断是否为毒性实验
                if is_toxicity_assay(parsed.get('Assay Name', '')):
                    tox_writer.writerow(out_row)
                    toxicity_count += 1

            print(f'  → Activity: {activity_count}, Toxicity: {toxicity_count}')

            # 更新进度
            progress['processed'].append(mol_id)
            processed_set.add(mol_id)
            success_count += 1

            # 每 100 个分子保存一次进度
            if success_count % 100 == 0:
                save_progress(progress)
                activity_file.flush()
                toxicity_file.flush()
                print(f'\n--- 进度保存: {success_count} 成功, {fail_count} 失败 ---\n')

    except KeyboardInterrupt:
        print('\n\n手动中断，保存进度...')
    finally:
        save_progress(progress)
        activity_file.close()
        toxicity_file.close()

    print(f'\n完成! 成功: {success_count}, 失败: {fail_count}')
    print(f'活动数据: {ACTIVITY_CSV}')
    print(f'毒性数据: {TOXICITY_CSV}')


if __name__ == '__main__':
    main()
