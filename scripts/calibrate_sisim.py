"""Independent, resumable SISIM variogram calibration; original run is read-only.

Uses the existing GSLIB parameter writer and exact GAM-equivalent estimator from
validation_model.py. Never changes SGSIM, search settings, data or original files.
Commands: screen --models candidates.json --nsim 2; cv --models selected.json --nsim 100.
Each output directory is immutable: changed model/settings require a new directory.
This adapter targets the saved 100 x 100 x 200 grid with 10 x 10 x 0.5 m cells.
The extra Neogene-domain analysis is reported separately from the full-grid GAM.
Use report followed by audit after 100 realizations and all five CV folds finish.
"""
from __future__ import annotations
import argparse, copy, hashlib, json, pickle, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from validation_model import (axial_variogram, read_sg, theoretical, reliability, PROB_COLS,
                              cluster_bootstrap_weights, weighted_ratio_ci)

def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding='utf-8')

def experiments(run):
    variance = float(np.var(np.load(run/'01/Cl_Sa_points.npy')[:,3]))
    result = {}
    for direction, filename, step in [('EW','categorical_horizontal_direction_0.csv',10),
                                      ('NS','categorical_horizontal_direction_1.csv',10),
                                      ('Z','categorical_vertical_direction_0.csv',.5)]:
        df = pd.read_csv(run/'02'/filename)
        df = df[(df.npairs>0)&(df.avsepdist>0)].sort_values('avsepdist')
        h = np.arange(1, 16 if direction!='Z' else 31)*step
        gamma = df.svargval.to_numpy()*(variance if direction=='Z' else 1)
        # Compare the fixed experimental curve interpolated to actual grid lags;
        # do not silently extrapolate or claim GAM/GAMV pair selection identical.
        valid=(h>=df.avsepdist.min())&(h<=df.avsepdist.max())
        result[direction]=(h[valid],np.interp(h[valid],df.avsepdist,gamma),df,variance)
    return result

def simulate(args, model, folder, conditioning, nsim, seed):
    if any(row[0] not in [1,2] for row in model['nest']):
        raise ValueError('This SISIM runner permits spherical/exponential indicator models only')
    folder.mkdir(parents=True,exist_ok=True)
    template=json.loads((args.run/'03/sisim_parameters_template.json').read_text())
    p=copy.deepcopy(template);p['nsim']=nsim;p['seed']=seed
    p['variogram']=[[model['nst'],model['c0'],model['nest']] for _ in range(2)]
    signature={'parameters':p,'conditioning_sha256':hashlib.sha256(conditioning.read_bytes()).hexdigest()}
    if (folder/'signature.json').exists():
        if json.loads((folder/'signature.json').read_text())!=signature:
            raise ValueError(f'Changed settings in existing output: {folder}')
    else: save(folder/'signature.json',signature)
    if (folder/'realizations.npy').exists():
        return np.load(folder/'realizations.npy',mmap_mode='r')
    shutil.copy2(conditioning,folder/'data.dat')
    from GSLIB_parfl_011 import GSLIB_SISIM
    GSLIB_SISIM(str(args.exe),str(folder/'sisim.par'),'data.dat','sisim.dbg','sisim.out',p)
    start=time.time()
    with (folder/'console.log').open('w') as log:
        subprocess.run([str(args.exe),'sisim.par'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,check=True)
    g=p['grid_params'];cells=g['nx']*g['ny']*g['nz']
    values=read_sg(folder/'sisim.out',cells)
    if values.shape!=(nsim,cells) or not np.isin(values,[0,1]).all():
        raise ValueError(f'Invalid SISIM output {folder}: {values.shape}')
    np.save(folder/'realizations.npy',values.astype(np.uint8))
    save(folder/'timing.json',{'seconds':time.time()-start,'nsim':nsim})
    # Parsed binary is lossless; discard only the bulky ASCII generated here.
    (folder/'sisim.out').unlink()
    print(f'Completed {folder.name}: {nsim} realizations in {time.time()-start:.1f}s',flush=True)
    return np.load(folder/'realizations.npy',mmap_mode='r')

def masked_variogram(field,mask,axis,nlag):
    """Exact GAM pair estimator, with invalid cells excluded at BOTH endpoints."""
    result=[]
    for lag in range(1,nlag+1):
        left=[slice(None)]*3;right=left.copy()
        left[axis]=slice(None,-lag);right[axis]=slice(lag,None)
        left,right=tuple(left),tuple(right)
        valid=mask[left]&mask[right]
        result.append(.5*np.count_nonzero((field[left]!=field[right])&valid)/np.count_nonzero(valid))
    return result

def domain_metrics(args,name,a,folder):
    """Additional scope diagnostic, NEVER substitutes for complete-grid metrics.

    Raw 0/1 SISIM values, restricted below the unchanged paired SGSIM boundary and
    terrain. These are the SISIM values actually used in the final Neogene model.
    Original experimental targets are kept fixed; their sampling support differs.
    """
    g=json.loads((args.run/'03/sisim_parameters_template.json').read_text())['grid_params']
    bases=np.load(args.run/'03/qbasis_realizations.npy',mmap_mode='r')
    terrain=np.load(args.run.parent/'below_gok.npy').reshape(g['nz'],g['ny'],g['nx'])
    z=g['zmn']+g['zsiz']*np.arange(g['nz'])
    curves={d:[] for d in ['EW','NS','Z']}
    for i,v in enumerate(a):
        field=v.reshape(g['nz'],g['ny'],g['nx'])
        mask=(z[:,None,None]<bases[i].reshape(g['ny'],g['nx'])[None,:,:])&terrain
        for direction,axis in [('EW',2),('NS',1),('Z',0)]:
            curves[direction].append(masked_variogram(field,mask,axis,30))
    metrics={'name':name,'nsim':len(a),'scope':'Raw SISIM restricted to actual Neogene region'};rows=[]
    for direction in curves:
        values=np.array(curves[direction]);np.save(folder/f'neogene_{direction}_curves.npy',values)
        step=.5 if direction=='Z' else 10;h=np.arange(1,31)*step
        x,y,_,_=experiments(args.run)[direction];idx=np.round(x/step).astype(int)-1
        mean=values.mean(axis=0);lo,hi=np.quantile(values,[.05,.95],axis=0)
        metrics[direction+'_MAE']=float(np.mean(abs(mean[idx]-y)))
        metrics[direction+'_bias']=float(np.mean(mean[idx]-y))
        metrics[direction+'_first_error']=float(mean[0]-y[0])
        metrics[direction+'_envelope_fraction']=float(((y>=lo[idx])&(y<=hi[idx])).mean())
        for j,lag in enumerate(h):
            rows.append(dict(direction=direction,distance=lag,mean=mean[j],p05=lo[j],p95=hi[j],experimental=np.interp(lag,x,y) if x.min()<=lag<=x.max() else np.nan))
    metrics['objective']=float(np.mean([metrics[d+'_MAE'] for d in ['EW','NS','Z']]))
    save(folder/'neogene_metrics.json',metrics);pd.DataFrame(rows).to_csv(folder/'neogene_curves.csv',index=False)
    print('Neogene scope: '+json.dumps(metrics),flush=True)
    return metrics

def spatial(args, name, model):
    folder=args.output/name
    if args.nsim<=5:
        a=simulate(args,model,folder,args.run/'01/tertiary.gslib',args.nsim,69069)
    else:
        parts=[simulate(args,model,folder/f'part_{part:02}',args.run/'01/tertiary.gslib',
                        min(5,args.nsim-start),69069+11*part)
               for part,start in enumerate(range(0,args.nsim,5))]
        a=np.concatenate(parts)
        np.save(folder/'realizations.npy',a)
    return spatial_metrics(args, name, a, folder)

def spatial_metrics(args,name,a,folder):
    g=json.loads((args.run/'03/sisim_parameters_template.json').read_text())['grid_params']
    exp=experiments(args.run);rows=[];metrics={'name':name,'nsim':len(a)}
    for direction,axis,step in [('EW',2,10),('NS',1,10),('Z',0,.5)]:
        curves=np.array([axial_variogram(v.reshape(g['nz'],g['ny'],g['nx']),axis,30)[1] for v in a])
        h=np.arange(1,curves.shape[1]+1)*step
        target_h,target,_,variance=exp[direction]
        idx=np.round(target_h/step).astype(int)-1
        mean=curves.mean(axis=0);lo,hi=np.quantile(curves,[.05,.95],axis=0)
        delta=mean[idx]-target
        metrics[direction+'_MAE']=float(np.abs(delta).mean())
        metrics[direction+'_bias']=float(delta.mean())
        metrics[direction+'_RMSE']=float(np.sqrt(np.mean(delta**2)))
        metrics[direction+'_envelope_fraction']=float(((target>=lo[idx])&(target<=hi[idx])).mean())
        for j,x in enumerate(h):
            rows.append(dict(direction=direction,distance=x,mean=mean[j],p05=lo[j],p95=hi[j],
                             experimental=np.interp(x,target_h,target) if target_h.min()<=x<=target_h.max() else np.nan))
        np.save(folder/f'{direction}_curves.npy',curves)
    metrics['objective']=float(np.mean([metrics[d+'_MAE'] for d in ['EW','NS','Z']]))
    metrics['mean_sand_proportion']=float(a.mean())
    save(folder/'metrics.json',metrics);pd.DataFrame(rows).to_csv(folder/'curves.csv',index=False)
    domain_metrics(args,name,a,folder)
    print(json.dumps(metrics),flush=True)
    return metrics

def cv_fold(args,name,model,fold):
    folder=args.output/name/f'fold_{fold}';source=args.run/f'cv/fold_{fold}'
    truth=pd.read_csv(source/'predictions.csv',dtype={'borehole_id':str})
    # Fixed saved fold SGSIM (training-only transform/trend), reused identically.
    with (source/'transform.pkl').open('rb') as f: method,transform,trend=pickle.load(f)
    if method!='Quantile':raise ValueError(method)
    g=json.loads((args.run/'03/sisim_parameters_template.json').read_text())['grid_params']
    xy=truth.voxel.to_numpy()%10000
    gaussian=read_sg(source/'sgsim.out',10000)
    bases=transform.inverse_transform(gaussian[:,xy].reshape(-1,1)).reshape(len(gaussian),-1)+np.asarray(trend)[xy]
    parts=[]
    for part,start in enumerate(range(0,args.nsim,5)):
        count=min(5,args.nsim-start)
        a=simulate(args,model,folder/f'part_{part:02}',source/'tertiary.gslib',count,69069+11*part)
        binary=np.asarray(a[:,truth.voxel])
        z=g['zmn']+(truth.voxel.to_numpy()//10000)*g['zsiz']
        parts.append(np.where(z[None,:]>=bases[start:start+count],3,binary+1))
    combined=np.concatenate(parts)
    probs=np.stack([(combined==c).mean(axis=0) for c in [1,2,3]],axis=1)
    # Preserve the original evaluate_nmodel tie-breaking exactly.
    pc,ps,pg=probs.T
    clay=(pc>=ps)&(pc>pg);sand=(ps>pc)&(ps>=pg)
    pred=np.where(clay,1,np.where(sand,2,3))
    result=truth[['borehole_id','x','y','z','voxel','observed']].copy()
    result['predicted']=pred
    for j,key in enumerate(['p_clay','p_sand','p_gravel']):result[key]=probs[:,j]
    result.to_csv(folder/'predictions.csv',index=False)
    print(f'CV completed {name} fold {fold}',flush=True)

def report(args,models):
    import matplotlib.pyplot as plt
    if len(models)!=1:raise ValueError('Report requires one frozen candidate')
    name,model=next(iter(models.items()))
    spatial_dir=args.spatial_dir/name
    base=args.output/'baseline_spatial';base.mkdir(exist_ok=True)
    if (base/'metrics.json').exists() and (base/'neogene_metrics.json').exists():
        baseline=json.loads((base/'metrics.json').read_text())
    else:
        baseline=spatial_metrics(args,'baseline',np.load(args.run/'03/sisim/realizations.npy',mmap_mode='r'),base)
    candidate=json.loads((spatial_dir/'metrics.json').read_text())
    if candidate['nsim']!=100 or baseline['nsim']!=100:
        raise ValueError('The final report requires 100 realizations for both spatial comparisons')
    for fold in range(5):
        count=sum(len(np.load(args.cv_dir/name/f'fold_{fold}/part_{j:02}/realizations.npy',mmap_mode='r')) for j in range(20))
        if count!=100:raise ValueError(f'CV fold {fold} does not contain the required 100 realizations')
    exp=experiments(args.run)
    for scope,prefix in [('Complete grid',''),('Actual Neogene region','neogene_')]:
        fig,axes=plt.subplots(1,3,figsize=(14,4.4))
        for direction,ax in zip(['EW','NS','Z'],axes):
            x,y,df,var=exp[direction]
            limit=150 if direction!='Z' else 15
            ax.scatter(df.avsepdist,df.svargval*(var if direction=='Z' else 1),color='black',s=16,label='Experimental GAMV',zorder=5)
            for path,label,color in [(base,'Original SISIM','#9b798b'),(spatial_dir,'Calibrated SISIM','#087d8b')]:
                curves=np.load(path/f'{prefix}{direction}_curves.npy')
                h=np.arange(1,curves.shape[1]+1)*(10 if direction!='Z' else .5)
                lo,hi=np.quantile(curves,[.05,.95],axis=0)
                ax.fill_between(h,lo,hi,color=color,alpha=.20)
                ax.plot(h,curves.mean(axis=0),color=color,lw=2,label=label)
            h=np.linspace(.001,limit,400)
            ax.plot(h,theoretical(h,model,{'EW':'E-W','NS':'N-S','Z':'Vertical'}[direction]),'--',color='#d48a1b',lw=1.5,label='Calibrated input')
            ax.set(xlim=(0,limit),ylim=(0,.31),xlabel='Separation [m]',ylabel='Semivariance',title=direction)
            ax.grid(alpha=.18)
            if direction=='EW':ax.legend(fontsize=8)
        fig.suptitle(f'Raw binary SISIM | {scope} | mean and 90% realization envelope')
        fig.tight_layout();fig.savefig(args.output/f'{prefix}spatial_reproduction.png',dpi=180)
        fig.savefig(args.output/f'{prefix}spatial_reproduction.pdf');plt.close(fig)
    assignment=pd.read_csv(args.run/'cv/borehole_folds.csv',dtype={'borehole_id':str}).drop_duplicates('borehole_id').set_index('borehole_id')
    datasets={}
    for label,root in [('Original',args.run/'cv'),('Calibrated',args.cv_dir/name)]:
        parts=[]
        for fold in range(5):
            df=pd.read_csv(root/f'fold_{fold}/predictions.csv',dtype={'borehole_id':str})
            df['fold']=fold;df['group']=assignment.reindex(df.borehole_id)['group'].to_numpy()
            parts.append(df)
        datasets[label]=pd.concat(parts,ignore_index=True)
    pd.testing.assert_frame_equal(datasets['Original'][['borehole_id','voxel','observed','fold']],datasets['Calibrated'][['borehole_id','voxel','observed','fold']])
    scores={};fig,(ax,hist)=plt.subplots(2,1,figsize=(6.5,7),gridspec_kw={'height_ratios':[3,1]})
    ax.plot([0,1],[0,1],'--',color='black',lw=1,label='Ideal')
    for (label,df),color in zip(datasets.items(),['#9b798b','#087d8b']):
        prob=df[PROB_COLS].to_numpy();truth=df.observed.to_numpy();correct=df.predicted.to_numpy()==truth
        confidence=prob.max(axis=1);groups=df.group.to_numpy();folds=df.fold.to_numpy()
        target=truth[:,None]==np.array([1,2,3])
        score={'N':len(df),'accuracy':float(correct.mean()),'mean_confidence':float(confidence.mean()),'multiclass_brier':float(np.mean(np.sum((prob-target)**2,axis=1)))}
        for bins in [5,10,20,50]:
            tab,ece=reliability(confidence,correct,groups,folds,bins,1000,12345)
            tab.to_csv(args.output/f'calibration_{label.lower()}_{bins}_bins.csv',index=False)
            score[f'ECE_{bins}_bins']=ece
            score[f'unweighted_bin_MAE_{bins}_bins']=float(np.abs(tab.mean_probability-tab.observed_frequency).mean())
            if bins==10:
                ax.plot(tab.mean_probability,tab.observed_frequency,'o-',color=color,lw=1.7,ms=4,label=f'{label}: ECE={ece:.3f}')
                ax.vlines(tab.mean_probability,tab.ci_low,tab.ci_high,color=color,lw=1,alpha=.7)
                hist.step((tab.lower+tab.upper)/2,tab['count'],where='mid',color=color,label=label)
        scores[label]=score
        df.to_csv(args.output/f'predictions_{label.lower()}.csv',index=False)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='Mean predicted confidence',ylabel='Fraction correctly classified',title='Whole-model out-of-fold calibration')
    ax.grid(alpha=.2);ax.legend(fontsize=9)
    hist.set(xlim=(0,1),xlabel='Confidence bin',ylabel='Composites');hist.legend(fontsize=8)
    fig.suptitle('Five whole-borehole folds | 95% group-bootstrap intervals',fontsize=11)
    fig.tight_layout();fig.savefig(args.output/'overall_calibration.png',dpi=180)
    fig.savefig(args.output/'overall_calibration.pdf');plt.close(fig)
    save(args.output/'comparison.json',{'spatial_original':baseline,'spatial_calibrated':candidate,'calibration':scores,
        'neogene_original':json.loads((base/'neogene_metrics.json').read_text()),'neogene_calibrated':json.loads((spatial_dir/'neogene_metrics.json').read_text()),
        'calibration_definition':'Top-label confidence versus observed correct fraction; one overall three-class curve. Not the single sand probability used by Wiegel.',
        'scope':'Spatially selected full-data variogram; fixed-variogram whole-borehole CV, not nested selection validation. SGSIM unchanged.',
        'spatial_target':'Fixed experimental curve at GAM axial distances, no independent curve rescaling; vertical GAMV de-standardized by conditioning indicator variance.'})

def ensemble_audit(args,models):
    """Paired Monte Carlo stability and binary diagnostic, after selection only."""
    name=next(iter(models));candidate=args.cv_dir/name
    mapping=pd.read_csv(args.run/'cv/borehole_folds.csv',dtype={'borehole_id':str}).drop_duplicates('borehole_id').set_index('borehole_id')
    records=[];paired={};binary_records=[]
    def ece(p,y,bins=10):
        idx=np.searchsorted(np.linspace(0,1,bins+1)[1:-1],p)
        return float(sum(abs(np.sum(y[idx==j]-p[idx==j])) for j in range(bins))/len(p))
    for count in args.audit_counts:
        frames={'Original':[],'Calibrated':[]}
        for fold in range(5):
            source=args.run/f'cv/fold_{fold}';truth=pd.read_csv(source/'predictions.csv',dtype={'borehole_id':str})
            loc=truth.voxel.to_numpy();xy=loc%10000
            with (source/'transform.pkl').open('rb') as f:method,transform,trend=pickle.load(f)
            if method!='Quantile':raise ValueError(method)
            gaussian=read_sg(source/'sgsim.out',10000)[:count,xy]
            bases=transform.inverse_transform(gaussian.reshape(-1,1)).reshape(count,-1)+np.asarray(trend)[xy]
            g=json.loads((args.run/'03/sisim_parameters_template.json').read_text())['grid_params']
            above=(g['zmn']+(loc//10000)*g['zsiz'])[None,:]>=bases
            for label in frames:
                if label=='Original':binary=np.load(source/'sisim/realizations.npy',mmap_mode='r')[:count,loc]
                else:binary=np.concatenate([np.load(candidate/f'fold_{fold}/part_{j:02}/realizations.npy',mmap_mode='r')[:,loc] for j in range(count//5)])
                combined=np.where(above,3,binary+1)
                prob=np.stack([(combined==k).mean(axis=0) for k in [1,2,3]],axis=1)
                pc,ps,pg=prob.T;pred=np.where((pc>=ps)&(pc>pg),1,np.where((ps>pc)&(ps>=pg),2,3))
                observed=truth.observed.to_numpy();target=observed[:,None]==np.array([1,2,3])
                if count==100:
                    reference=pd.read_csv((source if label=='Original' else candidate/f'fold_{fold}')/'predictions.csv')
                    np.testing.assert_array_equal(prob,reference[PROB_COLS].to_numpy())
                    np.testing.assert_array_equal(pred,reference.predicted)
                frames[label].append(pd.DataFrame(dict(fold=fold,group=mapping.reindex(truth.borehole_id)['group'].to_numpy(),
                    confidence=prob.max(axis=1),correct=pred==observed,brier=np.sum((prob-target)**2,axis=1),
                    binary_p_sand=binary.mean(axis=0),observed=observed,zero_true_prob=prob[np.arange(len(prob)),observed-1]==0)))
        for label,parts in frames.items():
            d=pd.concat(parts,ignore_index=True);binary=d[d.observed<3]
            records.append({'model':label,'realizations':count,'accuracy':float(d.correct.mean()),
                'ECE_10':ece(d.confidence.to_numpy(),d.correct.to_numpy()),'Brier':float(d.brier.mean()),
                'confidence':float(d.confidence.mean()),'zero_true_probability':int(d.zero_true_prob.sum()),
                'binary_sisim_ECE_10':ece(binary.binary_p_sand.to_numpy(),(binary.observed==2).to_numpy()),
                'binary_sisim_Brier':float(np.mean((binary.binary_p_sand-(binary.observed==2))**2))})
            if count==max(args.audit_counts):paired[label]=d
            if count==100:
                q=binary.binary_p_sand.to_numpy();y=(binary.observed==2).to_numpy()
                idx=np.searchsorted(np.linspace(0,1,51)[1:-1],q)
                for b in range(50):
                    use=idx==b
                    if use.any():binary_records.append(dict(model=label,bin=b,count=int(use.sum()),predicted=float(q[use].mean()),observed=float(y[use].mean())))
    pd.DataFrame(records).to_csv(args.output/'ensemble_stability.csv',index=False)
    pd.DataFrame(binary_records).to_csv(args.output/'binary_sisim_paper_style_calibration.csv',index=False)
    if not paired:
        print(pd.DataFrame(records).to_string(index=False));return
    old,new=paired['Original'],paired['Calibrated']
    idx,weights=cluster_bootstrap_weights(old.group.to_numpy(),old.fold.to_numpy(),2000,12345)
    paired_scores={}
    for key in ['correct','brier','confidence']:
        delta=new[key].astype(float).to_numpy()-old[key].astype(float).to_numpy()
        ci=weighted_ratio_ci(delta,np.ones(len(delta)),idx,weights)
        paired_scores[key]={'difference_new_minus_old':float(delta.mean()),'group_bootstrap_95_interval':ci.tolist()}
    # ECE is nonlinear: recompute each model's absolute bin gaps inside every
    # SAME grouped bootstrap sample, then compare the two scores.
    denominator=weights@np.bincount(idx,minlength=weights.shape[1])
    boot_ece={}
    for label,d in [('Original',old),('Calibrated',new)]:
        bins=np.searchsorted(np.linspace(0,1,11)[1:-1],d.confidence.to_numpy())
        gaps=np.zeros((weights.shape[1],10))
        np.add.at(gaps,(idx,bins),d.correct.to_numpy()-d.confidence.to_numpy())
        boot_ece[label]=np.abs(weights@gaps).sum(axis=1)/denominator
    delta_ece=boot_ece['Calibrated']-boot_ece['Original']
    paired_scores['ECE_10']={'difference_new_minus_old':ece(new.confidence.to_numpy(),new.correct.to_numpy())-ece(old.confidence.to_numpy(),old.correct.to_numpy()),
        'group_bootstrap_95_interval':np.quantile(delta_ece,[.025,.975]).tolist()}
    save(args.output/'paired_differences.json',paired_scores)
    print(pd.DataFrame(records).to_string(index=False));print(json.dumps(paired_scores,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['screen','cv','baseline','report','audit'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--exe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--models',type=Path);p.add_argument('--nsim',type=int,default=2)
    p.add_argument('--workers',type=int,default=2);p.add_argument('--folds',type=int,nargs='+',default=[0,1,2,3,4])
    p.add_argument('--spatial-dir',type=Path);p.add_argument('--cv-dir',type=Path)
    p.add_argument('--audit-counts',type=int,nargs='+',default=[20,50,100])
    args=p.parse_args();args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(args.repo/'scripts'))
    grid=json.loads((args.run/'03/sisim_parameters_template.json').read_text())['grid_params']
    if [grid[k] for k in ['nx','ny','nz','xsiz','ysiz','zsiz']] != [100,100,200,10,10,.5]:
        raise ValueError('This independent adapter is for the saved 100 x 100 x 200 project grid only')
    if not 1<=args.nsim<=100:raise ValueError('nsim must be in 1..100 to pair with saved SGSIM')
    if any(n<5 or n>100 or n%5 for n in args.audit_counts):
        raise ValueError('Audit counts must be multiples of five in 5..100')
    if args.stage=='baseline':
        spatial_metrics(args,'baseline',np.load(args.run/'03/sisim/realizations.npy',mmap_mode='r')[:args.nsim],args.output)
        return
    models=json.loads(args.models.read_text(encoding='utf-8-sig'))
    if args.stage=='report':
        report(args,models);return
    if args.stage=='audit':
        ensemble_audit(args,models);return
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        if args.stage=='screen':
            results=list(pool.map(lambda item:spatial(args,*item),models.items()))
            pd.DataFrame(results).sort_values('objective').to_csv(args.output/'ranking.csv',index=False)
        else:
            jobs=[(name,m,fold) for name,m in models.items() for fold in args.folds]
            list(pool.map(lambda item:cv_fold(args,*item),jobs))

if __name__=='__main__':main()
