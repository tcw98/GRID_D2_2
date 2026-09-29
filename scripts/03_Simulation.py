"""SGSIM + SISIM, matched realizations and original Munich hierarchical combination."""
import pickle
from workflow_support import *


def main():
    args=parser(__doc__,True).parse_args()
    with stage(args,"03",needs=("01","02")) as (cfg,out,manifest):
        models=json.loads(Path("02/adopted_variograms.json").read_text())["models"]
        params=simulate_gaussian(cfg,models["continuous"],"01/qbasis_no_check_300.gslib",out/"sgsim",args.gslib_dir)
        with Path("01/transformation_parameters.pkl").open("rb") as f: method,trans=pickle.load(f)
        trend=np.load("01/grid_data.npy")
        gaussian,detrended,elevations=U.post_processing_2D(str(out/"sgsim.out"),trend,params,method,trans)
        grid=cfg["preparation"]["grid"]
        nxy=grid["nx"]*grid["ny"];ncell=nxy*grid["nz"]
        if len(elevations)!=100*nxy or not np.isfinite(elevations).all():
            raise ValueError("Invalid SGSIM shape or nonfinite back-transformed elevation")
        for name,values in [("sgsim_gaussian",gaussian),("sgsim_detrended",detrended),("qbasis_realizations",elevations)]:
            np.save(out/f"{name}.npy",values.reshape(100,nxy))
        sisim=simulate_indicators(cfg,models["categorical"],"01/tertiary.gslib",out/"sisim",args.gslib_dir)
        combined=np.lib.format.open_memmap(out/"combined_realizations.npy",mode="w+",dtype=np.uint8,shape=(100,ncell))
        bases=elevations.reshape(100,nxy)
        for i in range(100):
            # Pure per-realization functions: same pairing and x-fastest ordering.
            qb3d=U.to_3D(bases[i],grid,1)
            result=U.combine_sim_ohnelm(sisim[i:i+1],qb3d,1)[0]
            if not np.isin(result,[1,2,3]).all(): raise ValueError("Unexpected combined class")
            combined[i]=result
        combined.flush()
        dump(out/"sgsim_parameters.json",params)
        dump(out/"sisim_parameters_template.json",sisim_params(cfg,models["categorical"]))
        manifest.update(realizations=100,shape=[100,ncell],organization="20 parts x 5; numeric part order; x fastest, then y, then z",
                        pairing="SGSIM i with SISIM i; no Cartesian product",classes={"1":"clay group","2":"sand","3":"Quaternary/gravel"})
        manifest["executables"]={name:digest(executable(args.gslib_dir,name)) for name in ["sgsim","sisim"]}


if __name__=="__main__": main()
