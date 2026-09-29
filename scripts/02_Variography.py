"""Legacy VARMAP/GAMV, directional fitting and explicit parameter sensitivity."""
import copy
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from workflow_support import *


def main():
    p=parser(__doc__,True)
    p.add_argument("--sweep",action="store_true",help="Run the original 25 lag/tolerance combinations for continuous data")
    args=p.parse_args()
    with stage(args,"02",needs=("01",)) as (cfg,out,manifest):
        cases=[("continuous",False,False,"01/qbasis_no_check_300.gslib"),
               ("categorical_horizontal",True,False,"01/tertiary.gslib"),
               ("categorical_vertical",True,True,"01/tertiary.gslib")]
        fitting=[]
        for name,cat,vertical,data in cases:
            params=gamv_params(cat,vertical)
            tables=read_variogram(args.gslib_dir,data,out/name,params)
            dump(out/f"{name}_settings.json",params)
            fig,ax=plt.subplots(figsize=(8,5))
            for i,table in enumerate(tables):
                ax.plot(table.avsepdist,table.svargval,"o",label=f"direction {i}")
                keep=(table.npairs>0)&np.isfinite(table.avsepdist)&np.isfinite(table.svargval)
                try:
                    if keep.sum()<4: raise ValueError("Fewer than four populated lag bins")
                    fit=U.fit_2_curve(ax,table.avsepdist[keep].to_numpy(),table.svargval[keep].to_numpy(),["blue","orange"][i%2],label=f"candidate {i}")
                    fitting.append(dict(case=name,direction=i,legacy_fit=fit,adopted=False))
                except (ValueError,RuntimeError,TypeError) as exc:
                    fitting.append(dict(case=name,direction=i,status="fit_unavailable",reason=str(exc),adopted=False))
            if not vertical:
                model=cfg["variograms"]["categorical" if cat else "continuous"]
                for label,fn in [("major",U.MAKENESTEDVARGM_Major),("minor",U.MAKENESTEDVARGM_Minor)]:
                    curve=fn(model["nst"],model["c0"],model["nest"],500,1)
                    ax.plot(curve.x,curve.y,label=f"adopted legacy {label}")
            ax.set(xlabel="Lag distance [m]",ylabel="Semivariance",title=name)
            ax.legend();fig.tight_layout();fig.savefig(out/f"{name}.png",dpi=150);plt.close(fig)
        dump(out/"fitting_candidates.json",fitting)
        for name,cat,data in [("continuous",False,cases[0][3]),("categorical",True,cases[1][3])]:
            params=varmap_params(cat)
            stem=out/f"varmap_{name}"
            G.GSLIB_VARMAP(str(executable(args.gslib_dir,"varmap")),str(stem)+".par",data,str(stem)+".out",params)
            run_gslib(args.gslib_dir,"varmap",str(stem)+".par")
            # read_VARMAP overwrites overflow markers; retain original bytes.
            shutil.copyfile(str(stem)+".out",str(stem)+"_raw.out")
            raw=Path(str(stem)+"_raw.out").read_text()
            table=U.read_VARMAP(str(stem)+".out")
            pd.DataFrame({"variogram":table}).to_csv(str(stem)+".csv",index=False)
            dump(str(stem)+"_settings.json",dict(parameters=params,legacy_overflow_replacements=raw.count("**********")))
            values=np.asarray(table).reshape(2*params["nzlag"]+1,2*params["nylag"]+1,2*params["nxlag"]+1)
            fig,ax=plt.subplots();img=ax.imshow(values[params["nzlag"]],origin="lower");fig.colorbar(img,ax=ax,label="Semivariance")
            ax.set(title=f"{name}: central XY variogram map",xlabel="X lag index",ylabel="Y lag index")
            fig.tight_layout();fig.savefig(str(stem)+".png",dpi=150);plt.close(fig)
            if cat:
                for plane,section in [("XZ",values[:,params["nylag"],:]),("YZ",values[:,:,params["nxlag"]])]:
                    fig,ax=plt.subplots();img=ax.imshow(section,origin="lower",aspect="auto")
                    fig.colorbar(img,ax=ax,label="Semivariance");ax.set(title=f"{name}: central {plane} variogram map",xlabel=plane[0]+" lag index",ylabel="Z lag index")
                    fig.tight_layout();fig.savefig(str(stem)+f"_{plane}.png",dpi=150);plt.close(fig)
        sweep=[]
        if args.sweep:
            (out/"sensitivity").mkdir()
            for xlag in range(5,30,5):
                for lagtol in range(3,18,3):
                    params=gamv_params()
                    params["xlag"]=xlag;params["lagtol"]=lagtol
                    for direction in params["xdir"]: direction[2]=xlag*3
                    tables=read_variogram(args.gslib_dir,cases[0][3],out/"sensitivity"/f"lag_{xlag}_tol_{lagtol}",params)
                    for i,table in enumerate(tables):
                        sweep.append(dict(xlag=xlag,lagtol=lagtol,bandwidth=xlag*3,direction=i,pairs=int(table.npairs.sum())))
            pd.DataFrame(sweep).to_csv(out/"sensitivity_summary.csv",index=False)
        dump(out/"adopted_variograms.json",dict(models=cfg["variograms"],source="explicit legacy branch, not automatically replaced by fit candidates"))
        manifest.update(sensitivity_executed=bool(args.sweep),fitting_note="Legacy fit_2_curve model-selection/plot behaviour retained; candidates do not change simulation parameters")
        manifest["executables"]={name:digest(executable(args.gslib_dir,name)) for name in ["gamv","varmap"]}


if __name__=="__main__": main()
