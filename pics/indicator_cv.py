"""Five-fold whole-borehole validation of the exponential hierarchical model."""
import argparse, json, os, pickle, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ap=argparse.ArgumentParser();ap.add_argument('--fold',type=int)
for name in ['work-dir','repo','gslib']:ap.add_argument('--'+name,type=Path,required=True)
a=ap.parse_args();P=a.work_dir.resolve();REPO=a.repo.resolve();RUN=P/'run';EXE=a.gslib.resolve()
sys.path.insert(0,str(REPO/'scripts'))
from workflow_support import U,load_config,simulate_gaussian,simulate_indicators
cfg=load_config(P/'exponential.json');g=cfg['preparation']['grid'];nx,ny,nz=g['nx'],g['ny'],g['nz']
os.chdir(RUN);out=Path('cv');out.mkdir(exist_ok=True)

if a.fold is None:
    s=pd.read_csv('01/survey_legacy_selected.csv');ids=s.ObjektIDObjektID.to_numpy().astype(str)
    xy=s[['xcoord','ycoord']].to_numpy();e=s.Ansatzhoeh.to_numpy()
    shifted=xy+(e-(np.ceil(e/.5)*.5-.5))[:,None]
    parent=np.arange(len(s))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def merge(i,j):parent[root(i)]=root(j)
    seen={}
    # Co-group identical IDs and every original/composited XY grid column.
    for i,ID in enumerate(ids):
        keys=[('id',ID)]
        for v in [xy[i],shifted[i]]:
            keys.append(('cell',int(np.ceil((v[0]-g['xmn'])/g['xsiz']+.5)-1),int(np.ceil((v[1]-g['ymn'])/g['ysiz']+.5)-1)))
        for key in keys:
            if key in seen:merge(i,seen[key])
            else:seen[key]=i
    groups=np.array([root(i) for i in range(len(s))]);unique=np.unique(groups)
    np.random.RandomState(12345).shuffle(unique)
    foldmap={key:i%5 for i,key in enumerate(unique)};folds=np.array([foldmap[k] for k in groups])
    assignments=pd.DataFrame({'borehole_id':ids,'x':xy[:,0],'y':xy[:,1],'group':groups,'fold':folds})
    assignments.to_csv(out/'borehole_folds.csv',index=False)
    assert assignments.groupby('borehole_id').fold.nunique().max()==1
    tree=cKDTree(shifted)
    cat=np.load('01/Cl_Sa_points.npy');full=np.load('01/tertiary_points.npy')
    full=full[np.isin(full[:,3],[1,2,3])]
    dc,ic=tree.query(cat[:,:2]);dt,it=tree.query(full[:,:2]);assert max(dc.max(),dt.max())<1e-5
    qb=np.load('01/qbasis.npy',allow_pickle=True);v=qb[:,1:].astype(float);valid=np.isfinite(v).all(axis=1)
    qb=qb[valid];v=v[valid];idfold=dict(zip(ids,folds));qfold=np.array([idfold[ID] for ID in qb[:,0]])
    np.savez(out/'inputs.npz',cat=cat,catfold=folds[ic],test=full,testfold=folds[it],testid=ids[it],qb=v,qfold=qfold)
    def run_fold(i):
        with (out/f'fold_{i}.log').open('w') as log:
            subprocess.run([sys.executable,'-B','-u',str(Path(__file__).resolve()),'--fold',str(i),'--work-dir',str(P),'--repo',str(REPO),'--gslib',str(EXE)],stdout=log,stderr=subprocess.STDOUT,check=True)
        print('CV fold complete:',i+1,flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run_fold,range(5)))
    metrics=pd.concat([pd.read_csv(out/f'fold_{i}/metrics.csv') for i in range(5)],ignore_index=True)
    metrics.to_csv(out/'fold_accuracy.csv',index=False)
    (out/'scope.json').write_text(json.dumps(dict(k=5,realizations_per_fold=100,grouping='Entire borehole IDs; neighbouring holes sharing original or composite XY grid columns grouped together',
        seed=12345,preprocessing='Training-only original save_train (Quantile transform); training-only categorical conditioning',
        model='Fixed exponential variograms; conditional CV, not nested variogram selection',
        target='Three-class hierarchical lithology below GOK, plus binary SISIM accuracy on clay/sand observations',
        indicator_accuracy='One-vs-rest (TP + TN) / N; class recall also reported',weighting='Per observed composite; borehole-macro accuracy separately recorded'),indent=2))
    print('Five-fold indicator CV complete',flush=True)
    sys.exit()

i=a.fold;folder=out/f'fold_{i}'
if folder.exists():
    for _ in range(1200):
        if (folder/'metrics.csv').exists():
            print('Completed fold reused:',i+1);sys.exit()
        time.sleep(1)
    raise RuntimeError('Existing fold did not complete; inspect its log before resuming')
folder.mkdir(exist_ok=False)
data=np.load(out/'inputs.npz');test=data['test'][data['testfold']==i];testids=data['testid'][data['testfold']==i]
traincat=data['cat'][data['catfold']!=i];trainqb=data['qb'][data['qfold']!=i]
loc,inside=U.checkgrid(*test[:,:3].T,g,XY_option=False)
mask=np.load(P/'below_gok.npy');inside &= mask[np.clip(loc,0,len(mask)-1)]
test=test[inside];testids=testids[inside];loc=loc[inside].astype(int)
trainloc,trainin=U.checkgrid(*traincat[:,:3].T,g,XY_option=False)
assert not set(trainloc[trainin]) & set(loc)
fitted=U.save_train([trainqb],str(folder/'boundary_'),1,g,moving_window=True,radius_xy=[300,300],tolerance=20,trans_method='Johnson')[0]
with (folder/'transform.pkl').open('wb') as f:pickle.dump(fitted,f)
params=simulate_gaussian(cfg,cfg['variograms']['continuous'],folder/'boundary_0.out',folder/'sgsim',EXE)
method,transform,trend=fitted
base=(U.backtransform(U.readsgsim(str(folder/'sgsim.out')),method,transform)+np.tile(trend,100)).reshape(100,nx*ny)
U.PANDASDF2GSLIBGeoEAS(pd.DataFrame(traincat,columns=['X','Y','Z','ST']),str(folder/'tertiary.gslib'))
sim=simulate_indicators(cfg,cfg['variograms']['categorical'],folder/'tertiary.gslib',folder/'sisim',EXE)
binary=np.asarray(sim[:,loc]);zz=g['zmn']+(loc//(nx*ny))*g['zsiz']
combined=np.where(zz[None,:]>=base[:,loc%(nx*ny)],3,binary+1)
pred,_,pc,ps,pg,_,_=U.evaluate_nmodel(combined,100)
truth=test[:,3].astype(int);prob=np.stack([pc,ps,pg],axis=1)
result=pd.DataFrame({'borehole_id':testids,'x':test[:,0],'y':test[:,1],'z':test[:,2],'voxel':loc,
    'observed':truth,'predicted':pred,'p_clay':pc,'p_sand':ps,'p_gravel':pg})
result.to_csv(folder/'predictions.csv',index=False)
matrix=np.zeros((3,3),int);np.add.at(matrix,(truth-1,pred.astype(int)-1),1)
pd.DataFrame(matrix,index=['clay','sand','gravel'],columns=['clay','sand','gravel']).to_csv(folder/'confusion.csv')
m=dict(fold=i+1,N=len(truth),test_boreholes=len(np.unique(testids)),accuracy=float((pred==truth).mean()),
       borehole_macro_accuracy=float(result.assign(correct=pred==truth).groupby('borehole_id').correct.mean().mean()),shared_train_test_voxels=0)
for code,name in enumerate(['clay','sand','gravel'],1):
    obs=truth==code;est=pred==code
    m[name+'_indicator_acc']=float((obs==est).mean())
    m[name+'_recall']=float(est[obs].mean()) if obs.any() else None
    m[name+'_precision']=float(obs[est].mean()) if est.any() else None
    m[name+'_N']=int(obs.sum())
two=truth<3;binpred=(binary.mean(axis=0)>.5).astype(int)+1
m['sisim_binary_acc']=float((binpred[two]==truth[two]).mean())
pd.DataFrame([m]).to_csv(folder/'metrics.csv',index=False)
print(json.dumps(m),flush=True)
