"""Prepare portable terrain/configuration inputs for the published article case."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    cfg=json.loads((root/'configs/exponential_20260916.json').read_text(encoding='utf8'))
    g=cfg['preparation']['grid']
    gok=np.load(root/'data/article/gok.npy',allow_pickle=False)
    mask=np.load(root/'data/article/below_gok.npy',allow_pickle=False)
    nxy=g['nx']*g['ny']
    if gok.size!=nxy or not np.isfinite(gok).all():raise ValueError('Invalid sampled terrain')
    z=g['zmn']+np.arange(g['nz'])*g['zsiz']
    expected=((z[:,None]+g['zsiz']/2)<=gok.reshape(1,nxy)).ravel()
    if not np.array_equal(mask,expected):raise ValueError('Terrain mask/grid mismatch')
    out=args.work_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(root/'configs/exponential_20260916.json',out/'exponential.json')
    for name in ['gok.npy','below_gok.npy','terrain.json','fit_record.json']:
        shutil.copy2(root/'data/article'/name,out/name)
    print(f'Prepared {out}: {int(mask.sum())} below-ground cells; no simulations executed.')

if __name__=='__main__':main()
