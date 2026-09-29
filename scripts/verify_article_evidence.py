"""Check packaged article evidence without GSLIB or the private full run."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

def read_json(path):return json.loads(path.read_text(encoding='utf-8-sig'))

def main():
    root=Path(__file__).resolve().parents[1]
    cfg=read_json(root/'configs/exponential_20260916.json')
    assert cfg==read_json(root/'pics/data/exponential.json'),'Article config mismatch'
    g=cfg['preparation']['grid'];nxy=g['nx']*g['ny']
    assert (g['nx'],g['ny'],g['nz'])==(100,100,200)
    terrain=np.load(root/'data/article/gok.npy',allow_pickle=False)
    mask=np.load(root/'data/article/below_gok.npy',allow_pickle=False)
    assert terrain.size==nxy and np.isfinite(terrain).all()
    expected=((g['zmn']+np.arange(g['nz'])[:,None]*g['zsiz']+g['zsiz']/2)<=terrain.reshape(1,nxy)).ravel()
    np.testing.assert_array_equal(mask,expected)
    assert mask.sum()==857561
    with (root/'pics/data/fold_accuracy.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    total=correct=0
    for row in rows:
        fold=int(row['fold'])
        matrix=np.loadtxt(root/f'pics/data/fold_{fold}_confusion.csv',delimiter=',',skiprows=1,usecols=(1,2,3),dtype=int)
        n=int(matrix.sum());hits=int(np.trace(matrix))
        assert n==int(row['N'])
        assert abs(hits/n-float(row['accuracy']))<1e-12
        assert int(row['shared_train_test_voxels'])==0
        for i,name in enumerate(['clay','sand','gravel']):
            tp=matrix[i,i];fn=matrix[i,:].sum()-tp;fp=matrix[:,i].sum()-tp;tn=n-tp-fn-fp
            assert abs((tp+tn)/n-float(row[name+'_indicator_acc']))<1e-12
            assert abs(tp/(tp+fn)-float(row[name+'_recall']))<1e-12
        total+=n;correct+=hits
    assert total==29961 and abs(correct/total-0.8198324488501719)<1e-12
    with (root/'pics/data/borehole_folds.csv').open(encoding='utf-8-sig',newline='') as f:
        assignments=list(csv.DictReader(f))
    by_id={};by_group={}
    for row in assignments:
        for mapping,key in [(by_id,row['borehole_id']),(by_group,row['group'])]:
            mapping.setdefault(key,set()).add(row['fold'])
    assert all(len(v)==1 for v in by_id.values()) and all(len(v)==1 for v in by_group.values())
    with (root/'pics/data/mean_entropy_by_realizations.csv').open(encoding='utf-8-sig',newline='') as f:stability=list(csv.DictReader(f))
    audit=read_json(root/'docs/article/source_audit.json')
    assert [int(r['number_of_realizations']) for r in stability]==[20,40,60,80,100]
    assert all(int(r['voxel_count'])==857561 for r in stability)
    np.testing.assert_allclose([float(r['mean_entropy']) for r in stability],audit['entropy']['means'],atol=1e-14,rtol=0)
    for file,digest in read_json(root/'docs/core_hashes.json').items():
        assert hashlib.sha256((root/file).read_bytes()).hexdigest()==digest,file
    vtk=root/'vtk/model_gok_exponential.vtk'
    assert hashlib.sha256(vtk.read_bytes()).hexdigest()==read_json(root/'pics/provenance.json')['vtk_sha256']
    print(f'PASS: grid/mask, 5 confusion matrices and metrics, fold grouping, stability tables, core hashes and VTK; accuracy={correct/total:.12f}; N={total}.')
    print('Scope: packaged evidence consistency; not a fresh simulation or independent statistical validation.')

if __name__=='__main__':main()
