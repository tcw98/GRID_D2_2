"""Hard conditioning, variogram reproduction, original k-fold metrics and checks."""
import matplotlib.pyplot as plt
from workflow_support import *


def conditioning(points,values,simulations,grid,output,tolerance,continuous=False):
    z=np.full(len(points),grid["zmn"]) if continuous else points[:,2]
    loc,inside=U.checkgrid(points[:,0],points[:,1],z,grid,XY_option=continuous)
    frame=pd.DataFrame(dict(source_point=np.arange(len(points)),x=points[:,0],y=points[:,1],z=points[:,2],
                            observed=values,in_grid=inside,voxel=loc))
    frame["same_voxel_conflict"]=False
    frame["fraction_reproduced"]=np.nan
    frame["mean_absolute_difference"]=np.nan
    for voxel,group in frame.loc[inside].groupby("voxel"):
        # No conflict resolution here: report the original data collision.
        conflict=np.ptp(group.observed)>tolerance
        frame.loc[group.index,"same_voxel_conflict"]=conflict
        for idx,row in group.iterrows():
            diff=np.abs(simulations[:,int(voxel)]-row.observed)
            frame.loc[idx,"fraction_reproduced"]=np.mean(diff<=tolerance)
            frame.loc[idx,"mean_absolute_difference"]=np.mean(diff)
    frame.to_csv(output,index=False)
    unambiguous=frame.loc[inside & ~frame.same_voxel_conflict]
    return dict(in_grid=int(inside.sum()),outside_grid=int((~inside).sum()),conflicting_points=int(frame.same_voxel_conflict.sum()),
                unambiguous_all_reproduced=bool(len(unambiguous)>0 and unambiguous.fraction_reproduced.eq(1).all()),
                tolerance=tolerance)


def reproduction(cfg,gslib,out):
    folder=out/"variogram_reproduction";folder.mkdir()
    for name,cat in [("continuous",False),("categorical",True)]:
        grid=copy.deepcopy(cfg["preparation"]["grid"])
        if not cat: grid["nz"]=1
        experimental=[pd.read_csv(f"02/{'categorical_horizontal' if cat else 'continuous'}_direction_{i}.csv") for i in range(2)]
        fig,axes=plt.subplots(1,2,figsize=(11,4))
        for i,t in enumerate(experimental): axes[i].plot(t.avsepdist,t.svargval,"o",label="input experimental")
        for realization in cfg["validation"]["reproduction_realizations"]:
            data=f"03/sisim/sisim_{(realization-1)//5}.out" if cat else "03/sgsim.out"
            number=(realization-1)%5+1 if cat else realization
            # GAM directions match GAMV east then north; original regular-grid writer.
            params=dict(nvar=1,ivar1=1,tmin=-1.e21,tmax=1.e21,nsim=number,grid_params=grid,
                ndir=2,nlag=30 if cat else 50,xdir=[[1,0,0],[0,1,0]],standardize=0,nvarg=1,ivtail=1,ivhead=1,ivtype=1)
            stem=folder/f"{name}_{realization}"
            G.GSLIB_GAM(str(executable(gslib,"gam")),str(stem)+".par",data,str(stem)+".out",params)
            run_gslib(gslib,"gam",str(stem)+".par")
            tables=U.read_GAMV_GAM(str(stem)+".out",params)
            for i,t in enumerate(tables):
                t.to_csv(str(stem)+f"_direction_{i}.csv",index=False)
                good=t.npairs>0;axes[i].plot(t.avsepdist[good],t.svargval[good],label=f"realization {realization}")
        for i,ax in enumerate(axes): ax.set(title=["E-W","N-S"][i],xlabel="Distance [m]",ylabel="Semivariance");ax.legend(fontsize=7)
        fig.tight_layout();fig.savefig(folder/f"{name}.png",dpi=150);plt.close(fig)
        if cat:
            experimental=pd.read_csv("02/categorical_vertical_direction_0.csv")
            fig,ax=plt.subplots();ax.plot(experimental.avsepdist,experimental.svargval,"o",label="input vertical (standardized)")
            for realization in cfg["validation"]["reproduction_realizations"]:
                params.update(nsim=(realization-1)%5+1,ndir=1,xdir=[[0,0,1]],standardize=1)
                data=f"03/sisim/sisim_{(realization-1)//5}.out"
                stem=folder/f"categorical_vertical_{realization}"
                G.GSLIB_GAM(str(executable(gslib,"gam")),str(stem)+".par",data,str(stem)+".out",params)
                run_gslib(gslib,"gam",str(stem)+".par")
                table=U.read_GAMV_GAM(str(stem)+".out",params)[0]
                table.to_csv(str(stem)+".csv",index=False)
                keep=table.npairs>0;ax.plot(table.avsepdist[keep],table.svargval[keep],label=f"realization {realization}")
            ax.set(xlabel="Vertical distance [m]",ylabel="Standardized semivariance");ax.legend(fontsize=7)
            fig.tight_layout();fig.savefig(folder/"categorical_vertical.png",dpi=150);plt.close(fig)


def cross_validation(cfg,gslib,out,models):
    folder=out/"cross_validation";folder.mkdir()
    grid=cfg["preparation"]["grid"]
    k=cfg["validation"]["continuous_folds"]
    np.random.seed(cfg["validation"]["split_seed"])
    coords=np.load("01/coords.npy")
    train,_,test,_=U.split_array(k,coords.copy(),coords_t=None)
    # Preserve save_train's ACTUAL behaviour: Quantile, +/-10 residual filter,
    # .out training files, even though the legacy caller requests Johnson.
    prefix=str(folder/"qbasis_")
    transform_parameters=U.save_train(train,prefix,k,grid,moving_window=True,
        radius_xy=cfg["preparation"]["radius_xy"],tolerance=cfg["preparation"]["tolerance"],trans_method="Johnson")
    for i in range(k):
        np.save(folder/f"continuous_test_{i}.npy",test[i])
        np.save(folder/f"continuous_train_{i}.npy",train[i])
        params=simulate_gaussian(cfg,models["continuous"],prefix+str(i)+".out",folder/f"sgsim_{i}",gslib)
    errors,quantiles,widths=U.k_fold_cross_validation(k,test,str(folder/"sgsim_"),params,transform_parameters)
    flat=np.concatenate(errors)
    if not len(flat) or not np.isfinite(flat).all(): raise ValueError("No valid continuous held-out errors")
    np.save(folder/"continuous_quantile_ranks.npy",quantiles)
    np.save(folder/"continuous_interval_widths.npy",widths)
    rows=[dict(fold=i,evaluated_points=len(e),legacy_ensemble_MAE=float(np.mean(e)) if len(e) else None) for i,e in enumerate(errors)]
    pd.DataFrame(rows).to_csv(folder/"continuous_folds.csv",index=False)
    # The existing binary CV reader assumes 100 realizations, no per-row padding.
    if grid["nx"]*grid["ny"]*grid["nz"]%8: raise ValueError("Legacy binary CV requires voxel count divisible by 8")
    kcat=cfg["validation"]["categorical_folds"]
    np.random.seed(cfg["validation"]["split_seed"])
    points=np.load("01/Cl_Sa_points.npy")
    traincat,testcat=U.split_array_3D(kcat,points)
    U.save_train_3D(traincat,kcat,str(folder/"tertiary_"))
    overlaps=[]
    for i in range(kcat):
        np.save(folder/f"categorical_test_{i}.npy",testcat[i]);np.save(folder/f"categorical_train_{i}.npy",traincat[i])
        a,ma=U.checkgrid(*traincat[i][:,:3].T,grid,XY_option=False)
        b,mb=U.checkgrid(*testcat[i][:,:3].T,grid,XY_option=False)
        overlaps.append(dict(fold=i,shared_voxels=len(set(a[ma])&set(b[mb])),test_in_grid=int(mb.sum()),
                             shared_xy=len(set(map(tuple,traincat[i][:,:2]))&set(map(tuple,testcat[i][:,:2])))))
        sims=simulate_indicators(cfg,models["categorical"],folder/f"tertiary_{i}.gslib",folder/f"sis_fold_{i}",gslib)
        U.write_boolean_array_to_binary(sims,str(folder/f"fold_{i}binary.binary"))
    pd.DataFrame(overlaps).to_csv(folder/"categorical_fold_overlap.csv",index=False)
    params=sisim_params(cfg,models["categorical"])
    correct,wrong,complete=U.k_fold_cross_validation_3D(kcat,testcat,str(folder/"fold_"),params)
    if not len(correct)+len(wrong): raise ValueError("No categorical held-out samples in grid")
    metrics=U.calculate_classification_metrics_k_fold(kcat,testcat,str(folder/"fold_"),grid)
    np.savez(folder/"categorical_calibration.npz",correct=correct,wrong=wrong,confidence=complete)
    summary=dict(continuous_legacy_ensemble_MAE=float(np.mean(flat)),continuous_evaluated=len(flat),
        categorical=metrics,split_seed=cfg["validation"]["split_seed"],
        preprocessing="Fold-specific legacy save_train; actual Quantile transform. Variograms fixed from adopted model (not nested refitting).",
        scope="Validation of boundary and indicator submodels; not joint hierarchical cross-validation",
        leakage_diagnostics=overlaps)
    dump(folder/"metrics.json",summary)
    return summary


def main():
    p=parser(__doc__,True);p.add_argument("--cv",action="store_true",help="Execute original k-fold simulations (100 realizations per fold)")
    args=p.parse_args()
    with stage(args,"05",needs=("01","02","03","04")) as (cfg,out,manifest):
        grid=cfg["preparation"]["grid"];n=grid["nx"]*grid["ny"]*grid["nz"]
        model=np.load("03/combined_realizations.npy",mmap_mode="r")
        arrays={k:np.load(f"04/{k}.npy",mmap_mode="r") for k in ["p_clay","p_sand","p_gravel","s_type","entropy","var"]}
        checks=dict(realizations_shape=model.shape==(100,n),valid_classes=True,probability_sum=True,
                    probability_range=True,finite_uq=True,legacy_entropy_bounds=True,variability_definition=True)
        for start in range(0,n,20000):
            sl=slice(start,start+20000);probs=np.stack([arrays[k][sl] for k in ["p_clay","p_sand","p_gravel"]])
            checks["valid_classes"] &= bool(np.isin(model[:,sl],[1,2,3]).all())
            checks["probability_sum"] &= bool(np.allclose(probs.sum(axis=0),1,rtol=0,atol=1.e-14))
            checks["probability_range"] &= bool(((probs>=0)&(probs<=1)).all())
            checks["finite_uq"] &= all(bool(np.isfinite(a[sl]).all()) for a in arrays.values())
            h=arrays["entropy"][sl]
            checks["legacy_entropy_bounds"] &= bool(((h>=0)&(h<=np.log(3)+.003)).all())
            checks["variability_definition"] &= bool(np.array_equal(arrays["var"][sl],1-probs.max(axis=0)))
        c=np.load("04/continuous_statistics.npz")
        checks["quantile_order"]=bool(((c["P10"]<=c["P50"])&(c["P50"]<=c["P90"])).all())
        dump(out/"output_consistency.json",checks)
        if not all(checks.values()): raise AssertionError("Output consistency failed")
        points=np.load("01/coords_t.npy");keep=(points[:,2]>=-3)&(points[:,2]<=3)
        grid2=grid.copy();grid2["nz"]=1
        hard={"continuous_gaussian":conditioning(points[keep],points[keep,2],np.load("03/sgsim_gaussian.npy",mmap_mode="r"),grid2,out/"hard_continuous.csv",1.e-3,True)}
        cs=np.load("01/Cl_Sa_points.npy")
        hard["categorical"]=conditioning(cs,cs[:,3],np.load("03/sisim/realizations.npy",mmap_mode="r"),grid,out/"hard_categorical.csv",0)
        composite=np.load("01/tertiary_points.npy");composite=composite[np.isin(composite[:,3],[1,2,3])]
        hard["final_hierarchy"]=conditioning(composite,composite[:,3],model,grid,out/"hard_combined.csv",0)
        dump(out/"hard_data_summary.json",hard)
        reproduction(cfg,args.gslib_dir,out)
        manifest["cross_validation"]="executed" if args.cv else "not_requested"
        if args.cv: manifest["cv_metrics"]=cross_validation(cfg,args.gslib_dir,out,cfg["variograms"])
        manifest["hard_conditioning"]=hard
        manifest["assessment"]="Diagnostics are reported, not repaired; execution complete does not imply every geological check passed"


if __name__=="__main__": main()
