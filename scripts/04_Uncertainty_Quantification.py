"""Original class probabilities, mode, entropy and variability; continuous summaries.

Entropy is deliberately UNCHANGED: natural log, zero probabilities replaced by
0.0001, NO normalization. Keep original category-tie behaviour and fs_type.
"""
import matplotlib.pyplot as plt
import pyvista as pv
from workflow_support import *


NAMES=("s_type","fs_type","p_clay","p_sand","p_gravel","entropy","var")


def legacy_uq(model,folder,chunk):
    arrays={name:np.lib.format.open_memmap(folder/f"{name}.npy",mode="w+",dtype=float,shape=(model.shape[1],)) for name in NAMES}
    for start in range(0,model.shape[1],chunk):
        end=min(start+chunk,model.shape[1])
        result=U.evaluate_nmodel(model[:,start:end],len(model))
        for name,values in zip(NAMES,result): arrays[name][start:end]=values
    for a in arrays.values(): a.flush()
    return arrays


def main():
    args=parser(__doc__).parse_args()
    with stage(args,"04",needs=("01","03")) as (cfg,out,manifest):
        model=np.load("03/combined_realizations.npy",mmap_mode="r")
        if len(model)!=100: raise ValueError("Expected the complete 100-realization ensemble")
        chunk=cfg["uq"]["chunk_cells"]
        arrays=legacy_uq(model,out,chunk)
        q=np.load("03/qbasis_realizations.npy",mmap_mode="r")
        continuous=dict(mean=np.mean(q,axis=0),variance=np.var(q,axis=0))
        continuous.update(zip(["P10","P50","P90"],np.percentile(q,[10,50,90],axis=0)))
        np.savez(out/"continuous_statistics.npz",**continuous)
        rows=[]
        counts=cfg["uq"]["stability_counts"]
        for low,high in zip(counts[:-1],counts[1:]):
            sums=np.zeros(3);n=0
            for start in range(0,model.shape[1],chunk):
                a=U.evaluate_nmodel(model[:low,start:start+chunk],low)
                b=U.evaluate_nmodel(model[:high,start:start+chunk],high)
                size=len(a[0]);n+=size
                sums+=np.array([np.abs(a[5]-b[5]).sum(),sum(np.abs(a[i]-b[i]).sum() for i in [2,3,4])/3,np.abs(a[6]-b[6]).sum()])
            rows.append(dict(N=low,N_next=high,mean_abs_entropy_change=sums[0]/n,mean_abs_probability_change=sums[1]/n,mean_abs_variability_change=sums[2]/n,
                             mean_abs_P50_change=float(np.mean(np.abs(np.median(q[:low],axis=0)-np.median(q[:high],axis=0))))))
        frame=pd.DataFrame(rows);frame.to_csv(out/"realization_stability.csv",index=False)
        fig,ax=plt.subplots()
        for name in ["mean_abs_entropy_change","mean_abs_probability_change","mean_abs_variability_change"]:
            ax.plot(frame.N_next,frame[name],"o-",label=name)
        ax.set(xlabel="Larger ensemble size",ylabel="Mean absolute change");ax.legend(fontsize=7);fig.tight_layout();fig.savefig(out/"stability.png",dpi=150);plt.close(fig)
        grid=cfg["preparation"]["grid"]
        # Geometry-only API compatibility with PyVista 0.43 (UniformGrid -> ImageData).
        mesh=pv.ImageData(dimensions=(grid["nx"]+1,grid["ny"]+1,grid["nz"]+1),
            spacing=(grid["xsiz"],grid["ysiz"],grid["zsiz"]),
            origin=(grid["xmn"]-grid["xsiz"]/2,grid["ymn"]-grid["ysiz"]/2,grid["zmn"]-grid["zsiz"]/2))
        for name,values in arrays.items(): mesh.cell_data[name]=np.asarray(values)
        mesh.save(out/"model_uq_rotated_frame.vti")
        manifest.update(entropy="legacy natural-log entropy, zero replacement 0.0001, not normalized",
                        terrain_mask="not supplied; volume includes above-ground voxels",continuous_quantiles="NumPy linear percentiles, absolute base elevation",
                        stability="fixed prefixes of the SAME ensemble; diagnostic, not proof of geological accuracy")


if __name__=="__main__": main()
