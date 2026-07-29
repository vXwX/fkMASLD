import os
import subprocess
import sys

D3POCKETS_DIR = os.path.expanduser('~/D3Pockets/D3Pockets_2.7.2')
APO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'protein', 'pdb_apo')
OUT_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pocket_output')
RADIUS_FILE = os.path.join(D3POCKETS_DIR, 'data', 'atomic_radius.json')

CONDA_ENV = 'D3Pockets2J'

# Skip these (no structure or membrane protein)
SKIP = ['DGAT2', 'ELOVL6']


def run_d3pockets(pdb_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cmd = [
        'conda', 'run', '-n', CONDA_ENV, 'python', os.path.join(D3POCKETS_DIR, 'main.py'),
        '-t', 'Pockets',
        '-p', pdb_path,
        '-o', out_dir,
        '--radiusfile', RADIUS_FILE,
        '--bradius', '4.0',
        '--sradius', '2.0',
        '--stepsize', '1.0',
        '--batchsize', '1000',
    ]
    print(f'Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd=D3POCKETS_DIR)
    if result.returncode != 0:
        print(f'  ERROR: {result.stderr[-500:]}')
    else:
        print(f'  OK')
    return result.returncode == 0


def main():
    apo_files = sorted([f for f in os.listdir(APO_DIR) if f.endswith('_apo.pdb')])
    print(f'Found {len(apo_files)} apo protein files')

    success = []
    failed = []

    for fname in apo_files:
        gene = fname.split('_apo')[0]

        pdb_path = os.path.join(APO_DIR, fname)
        out_dir = os.path.join(OUT_BASE, gene)
        if os.path.exists(os.path.join(out_dir, 'pocket_detected')):
            existing = [f for f in os.listdir(os.path.join(out_dir, 'pocket_detected')) if f.endswith('.pdb')]
            if existing:
                print(f'\n[{gene}] Already processed ({len(existing)} pockets), skip')
                success.append(gene)
                continue

        print(f'\n[{gene}] Processing {fname}...')
        ok = run_d3pockets(pdb_path, out_dir)
        if ok:
            success.append(gene)
        else:
            failed.append(gene)

    print(f'\n{"="*50}')
    print(f'Done. Success: {len(success)}, Failed: {len(failed)}')
    if failed:
        print(f'Failed: {failed}')
    for g in success:
        pkt_dir = os.path.join(OUT_BASE, g, 'pockets')
        if os.path.exists(pkt_dir):
            count = len([f for f in os.listdir(pkt_dir) if f.endswith('.pdb')])
            print(f'  {g}: {count} pockets')


if __name__ == '__main__':
    main()