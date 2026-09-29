from pathlib import Path
import argparse, copy, hashlib, json, shutil, sys
import numpy as np
import rasterio

parser=argparse.ArgumentParser()
for name in ['work-dir','repo','source-run','terrain']:parser.add_argument('--'+name,type=Path,required=True)
args=parser.parse_args()
ROOT=args.work_dir.resolve();ROOT.mkdir(parents=True,exist_ok=True)
REPO=args.repo.resolve();OLD=args.source_run.resolve()
RUN=ROOT/'run';RUN.mkdir(exist_ok=True)
sys.path.insert(0,str(REPO/'scripts'))
import UTILITIES_Wang as U
cfg=json.loads((REPO/'configs/example_1km.json').read_text(encoding='utf-8'))
fits=json.loads((OLD/'02/fitting_candidates.json').read_text())
def fitted(case,i):
    v=next(f['legacy_fit'] for f in fits if f['case']==case and f['direction']==i)
    assert v[-1]=='Exponential'
    return copy.deepcopy(v)
def model(summary,vertical):
    major,minor=summary['ranges'][:2]
    return dict(nst=1,c0=summary['nugget'],nest=[[2,summary['sill']-summary['nugget'],0,0,0,major,minor,vertical]])
c=U.summarize_variograms_2D(fitted('continuous',1),fitted('continuous',0),False)
v=fitted('categorical_vertical',0)
variance=float(np.var(np.load(OLD/'01/Cl_Sa_points.npy')[:,3]))
v[0]*=variance;v[2]*=variance
t=U.summarize_variograms(fitted('categorical_horizontal',1),fitted('categorical_horizontal',0),v,False)
assert fitted('continuous',1)[1]>fitted('continuous',0)[1]
assert fitted('categorical_horizontal',1)[1]>fitted('categorical_horizontal',0)[1]
cfg['profile']='real_munich_exponential_20260916'
cfg['variograms']={'continuous':model(c,6),'categorical':model(t,t['ranges'][2])}
cfg['validation']['categorical_folds']=5
config=ROOT/'exponential.json';config.write_text(json.dumps(cfg,indent=2))
digest=hashlib.sha256(config.read_bytes()).hexdigest()
for stage in ['01','02']:
    shutil.copytree(OLD/stage,RUN/stage)
    path=RUN/stage/'run_manifest.json';m=json.loads(path.read_text(encoding='utf-8'))
    m.update(reused_from=str(OLD/stage),original_config_digest=m['config_digest'],config_digest=digest,
             reuse_scope='Unchanged preparation and experimental variograms; theoretical simulation models replaced explicitly by exponential summaries')
    path.write_text(json.dumps(m,indent=2))
(RUN/'02/adopted_variograms.json').write_text(json.dumps({'models':cfg['variograms'],'source':'Original exponential fits and original summarize_variograms functions; vertical indicator sill converted from standardized units'},indent=2))
(ROOT/'fit_record.json').write_text(json.dumps(dict(individual_fits=fits,indicator_variance=variance,
    continuous_summary=c,indicator_summary=t,models=cfg['variograms']),indent=2))
g=cfg['preparation']['grid'];nx,ny,nz=g['nx'],g['ny'],g['nz']
x=g['xmn']+np.arange(nx)*g['xsiz'];y=g['ymn']+np.arange(ny)*g['ysiz'];xx,yy=np.meshgrid(x,y)
xy=np.c_[xx.ravel(),yy.ravel()];elev=np.full(nx*ny,np.nan);sources=[]
folder=args.terrain.resolve()
for f in folder.glob('*.tif'):
    with rasterio.open(f) as d:
        b=d.bounds;use=(xy[:,0]>=b.left)&(xy[:,0]<b.right)&(xy[:,1]>=b.bottom)&(xy[:,1]<b.top)
        if not use.any():continue
        assert d.crs.to_epsg()==25832
        samples=np.ma.concatenate([s for s in d.sample(xy[use],masked=True)]).astype(float)
        elev[use]=samples.filled(np.nan)
        sources.append({'name':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
assert np.isfinite(elev).all()
np.save(ROOT/'gok.npy',elev)
# Hide any cell whose top face exceeds sampled GOK (stricter than centre-only).
z=g['zmn']+np.arange(nz)*g['zsiz']
mask=(z[:,None]+g['zsiz']/2)<=elev[None,:]
np.save(ROOT/'below_gok.npy',mask.ravel())
(ROOT/'terrain.json').write_text(json.dumps(dict(crs='EPSG:25832',rotation=0,sources=sources,
    sample='Nearest native 1 m raster cell at each 10 m grid centre',
    inclusion='Voxel top elevation <= sampled GOK',gok_min=float(elev.min()),gok_max=float(elev.max()),
    retained_cells=int(mask.sum()),excluded_cells=int((~mask).sum())),indent=2))
print('Prepared exponential profile and ground mask:',int(mask.sum()),'subsurface cells')
