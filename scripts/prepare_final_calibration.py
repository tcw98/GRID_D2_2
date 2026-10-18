"""Audit response_step_50 CV and extract its binary calibration plot tables.

--init freezes existing inputs and copies unchanged SGSIM validation data only.
Without --init, refuses to finalize until all five folds have 100 realizations.
"""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[1]
DATA=R/'pics/final_20261018/data'
RUN=CV=HERE=SGSAVED=EXE=None
MODEL=R/'docs/final_20261018/selected_variogram.json'
ADAPTER=R/'scripts/calibrate_sisim.py'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def save(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False),encoding='utf-8')

def inputs():
    p=[RUN/'cv/borehole_folds.csv',RUN/'cv/inputs.npz',RUN/'03/sisim_parameters_template.json',
       MODEL,ADAPTER,
       EXE]
    for i in range(5):p.extend(RUN/f'cv/fold_{i}/{name}' for name in ['tertiary.gslib','sgsim.out','transform.pkl','predictions.csv'])
    for name in ['continuous_heldout_predictions.csv','continuous_quantile_calibration.csv','continuous_interval_coverage.csv']:
        p.append(SGSAVED/name)
    return {str(p):sha(p) for p in p}

def main():
    global RUN,CV,HERE,SGSAVED,EXE
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--init',action='store_true')
    for name in ['run','cv-root','sgsim-validation-dir','audit-dir','exe']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    RUN=a.run.resolve();CV=a.cv_root.resolve();SGSAVED=a.sgsim_validation_dir.resolve()
    HERE=a.audit_dir.resolve();HERE.mkdir(parents=True,exist_ok=True);EXE=a.exe.resolve()
    if a.init:
        protocol=dict(candidate='response_step_50',folds=5,realizations_per_fold=100,
            input_hashes=inputs(),model=json.loads((MODEL).read_text(encoding='utf-8'))['response_step_50'],
            grouping='Unchanged whole-borehole/group folds; source training-only data and SGSIM reused.',
            binary_target='Raw SISIM sand probability at withheld observed Neogene samples (observed classes 1 or 2), not final hierarchical sand probability.',
            bins='Fixed 50 bins of width .02, intervals (lower,upper], p=0 included in first bin. Empty bins remain NaN. No smoothing.',
            sgsim='Unchanged 20261015 continuous out-of-fold predictions. ECDF of saved mid-PIT, same 100-member ensembles; no new SGSIM.',
            plot='One final SISIM 2% staircase plus bisector; one SGSIM cumulative-rank curve plus bisector. No old-model overlay.',
            reference='Wiegel, Peña-Olarte & Cudmani 2024, doi:10.1002/gete.202400011, Figs. 7a and 9. Our CV remains five-fold for both.',
            selection='Variogram fixed before this CV; not adjusted using its calibration curve. Fixed full-data variogram CV is conditional, not nested model-selection validation.')
        dest=HERE/'cv_protocol.json'
        if dest.exists():assert json.loads(dest.read_text(encoding='utf-8'))==protocol
        save(dest,protocol)
        for name in ['continuous_heldout_predictions.csv','continuous_quantile_calibration.csv','continuous_interval_coverage.csv']:
            shutil.copy2(SGSAVED/name,DATA/name)
        print('Input hashes, fixed CV protocol and unchanged SGSIM data ready');return
    protocol=json.loads((HERE/'cv_protocol.json').read_text(encoding='utf-8'))
    assert inputs()==protocol['input_hashes'],'Protected source changed'
    model=protocol['model'];expected_variogram=[[4,0,model['nest']]]*2
    template=json.loads((RUN/'03/sisim_parameters_template.json').read_text(encoding='utf-8'))
    assignment=pd.read_csv(RUN/'cv/borehole_folds.csv',dtype={'borehole_id':str}).drop_duplicates('borehole_id').set_index('borehole_id')
    frames=[];foldrows=[];batchrows=[];full=[]
    for i in range(5):
        folder=CV/f'fold_{i}';source=RUN/f'cv/fold_{i}'
        if not (folder/'predictions.csv').exists():raise ValueError(f'Fold {i+1} incomplete')
        truth=pd.read_csv(source/'predictions.csv',dtype={'borehole_id':str})
        prediction=pd.read_csv(folder/'predictions.csv',dtype={'borehole_id':str})
        pd.testing.assert_frame_equal(truth[['borehole_id','voxel','observed']],prediction[['borehole_id','voxel','observed']])
        selected=truth.observed.isin([1,2]);binary=truth.loc[selected,['borehole_id','x','y','z','voxel','observed']].copy()
        binary['fold']=i+1;binary['group']=assignment.reindex(binary.borehole_id)['group'].to_numpy()
        assert not binary.group.isna().any()
        voxels=binary.voxel.to_numpy();success=np.zeros(len(binary),dtype=np.int64);count=0
        for j in range(20):
            part=folder/f'part_{j:02}';raw=part/'realizations.npy'
            if not (part/'timing.json').exists():raise ValueError(f'Fold {i+1} batch {j} incomplete')
            timing=json.loads((part/'timing.json').read_text(encoding='utf-8'));assert timing['nsim']==5
            signature=json.loads((part/'signature.json').read_text(encoding='utf-8'));p=signature['parameters']
            assert p['variogram']==expected_variogram and p['nsim']==5 and p['seed']==69069+11*j
            baseline=p.copy();baseline['variogram']=template['variogram'];baseline['nsim']=template['nsim'];baseline['seed']=template['seed']
            assert baseline==template,'Other native parameters changed'
            assert signature['conditioning_sha256']==sha(source/'tertiary.gslib')==sha(part/'data.dat')
            assert 'Stop - Program terminated.' in (part/'console.log').read_text(errors='replace')
            field=np.load(raw,mmap_mode='r');assert field.shape==(5,2000000) and field.dtype==np.uint8
            assert field.min()==0 and field.max()==1
            success+=field[:,voxels].sum(axis=0,dtype=np.int64);count+=5
            batchrows.append(dict(fold=i+1,part=j,realizations=5,seconds=timing['seconds'],binary_sha256=sha(raw),signature_sha256=sha(part/'signature.json')))
        assert count==100
        binary['sand_count']=success;binary['realizations']=count;binary['p_sand']=success/count
        binary['is_sand']=(binary.observed==2).astype(int)
        binary['predicted']=np.where(binary.p_sand>.5,2,1)  # exact original clay tie rule
        frames.append(binary)
        probabilities=prediction[['p_clay','p_sand','p_gravel']].to_numpy()
        np.testing.assert_allclose(probabilities.sum(axis=1),1,rtol=0,atol=np.finfo(float).eps)
        foldrows.append(dict(fold=i+1,realizations=count,test_N=len(truth),binary_N=len(binary),
            binary_accuracy=float((binary.predicted==binary.observed).mean()),binary_brier=float(np.mean((binary.p_sand-binary.is_sand)**2)),
            overall_accuracy=float((prediction.predicted==prediction.observed).mean())))
        prediction['fold']=i+1;full.append(prediction)
    data=pd.concat(frames,ignore_index=True);assert len(data)==22313
    assert data.observed.value_counts().to_dict()=={2:12876,1:9437}
    assert data.groupby('borehole_id').fold.nunique().max()==1
    assert data.groupby('group').fold.nunique().max()==1
    data.to_csv(DATA/'binary_sisim_cv_predictions.csv',index=False)
    table=[];edges=np.linspace(0,1,51);bins=np.searchsorted(edges[1:-1],data.p_sand.to_numpy())
    for b in range(50):
        d=data[bins==b]
        table.append(dict(bin=b,lower=edges[b],upper=edges[b+1],count=len(d),sand_count=int(d.is_sand.sum()),
            mean_probability=float(d.p_sand.mean()) if len(d) else np.nan,
            observed_frequency=float(d.is_sand.mean()) if len(d) else np.nan))
    table=pd.DataFrame(table);assert table['count'].sum()==len(data)
    table.to_csv(DATA/'binary_sisim_calibration_2pct.csv',index=False)
    pd.DataFrame(foldrows).to_csv(DATA/'final_cv_fold_metrics.csv',index=False)
    pd.DataFrame(batchrows).to_csv(HERE/'native_batches_verified.csv',index=False)
    pd.concat(full,ignore_index=True).to_csv(HERE/'final_hierarchical_cv_predictions.csv',index=False)
    valid=table['count']>0;gap=abs(table.mean_probability-table.observed_frequency)
    continuous=pd.read_csv(DATA/'continuous_heldout_predictions.csv')
    assert len(continuous)==664 and continuous.ensemble_size.eq(100).all()
    for name in ['continuous_heldout_predictions.csv','continuous_quantile_calibration.csv','continuous_interval_coverage.csv']:
        assert sha(DATA/name)==sha(SGSAVED/name)
    qt=pd.read_csv(DATA/'continuous_quantile_calibration.csv')
    summary=dict(status='complete',candidate='response_step_50',folds=5,realizations_per_fold=[r['realizations'] for r in foldrows],
        actual_new_sisim_realizations=500,actual_native_batches=100,total_heldout_N=sum(r['test_N'] for r in foldrows),
        binary_N=len(data),binary_boreholes=int(data.borehole_id.nunique()),binary_groups=int(data.group.nunique()),
        binary_accuracy=float((data.predicted==data.observed).mean()),binary_brier=float(np.mean((data.p_sand-data.is_sand)**2)),
        binary_ECE_2pct=float(np.sum(table['count']*gap.fillna(0))/len(data)),binary_unweighted_bin_MAE=float(gap[valid].mean()),
        binary_empty_bins=int((~valid).sum()),binary_min_nonempty_bin_count=int(table.loc[valid,'count'].min()),
        continuous_N=len(continuous),continuous_folds=int(continuous.fold.nunique()),sgsim_rerun=False,
        continuous_saved_quantile_MAE=float(abs(qt.nominal-qt.observed).mean()),
        selected_variogram_sha256=sha(MODEL),
        cv_protocol_sha256=sha(HERE/'cv_protocol.json'),protected_input_hashes_unchanged=True,
        interpretation='Fixed full-data variogram, whole-borehole/group five-fold CV; conditional validation, not nested variogram selection.',
        binary_scope=protocol['binary_target'],plot_bins=protocol['bins'],reference=protocol['reference'])
    save(DATA/'final_cv_summary.json',summary);save(HERE/'summary.json',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
