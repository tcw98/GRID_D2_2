"""Compare saved outputs with direct original functions (no new simulations).

Intended for small regression runs: the original all-at-once UQ and combination
would require much more memory than the streaming production workflow.
"""
import pickle
from workflow_support import *


def main():
    p=parser(__doc__);args=p.parse_args()
    run=args.run_dir.resolve();cfg=load_config(args.config.resolve())
    grid=cfg["preparation"]["grid"]
    if grid["nx"]*grid["ny"]*grid["nz"]>100000:
        raise ValueError("Run this all-at-once regression on a smoke grid (<=100000 voxels)")
    checks={}
    def compare(name,actual,expected):
        a,b=np.asarray(actual),np.asarray(expected)
        checks[name]=dict(passed=bool(np.array_equal(a,b,equal_nan=True)),shape=list(a.shape),
                          max_absolute_difference=float(np.max(np.abs(a.astype(float)-b.astype(float)))))
    with (run/"regression.log").open("w",encoding="utf-8") as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        raw=U.readsgsim(str(run/"03/sgsim.out"))
        compare("sgsim_reader",np.load(run/"03/sgsim_gaussian.npy"),raw.reshape(100,-1))
        params=json.loads((run/"03/sgsim_parameters.json").read_text())
        with (run/"01/transformation_parameters.pkl").open("rb") as f: method,trans=pickle.load(f)
        _,_,bases=U.post_processing_2D(str(run/"03/sgsim.out"),np.load(run/"01/grid_data.npy"),params,method,trans)
        compare("sgsim_backtransform_and_trend",np.load(run/"03/qbasis_realizations.npy"),bases.reshape(100,-1))
        indicators=U.read_sisim(str(run/"03/sisim/sisim_"))
        compare("sisim_20_by_5_reader",np.load(run/"03/sisim/realizations.npy"),indicators)
        reference=U.combine_sim_ohnelm(indicators,U.to_3D(bases,grid,100),100)
        compare("hierarchical_combination",np.load(run/"03/combined_realizations.npy"),reference)
        expected=U.evaluate_nmodel(reference,100)
        for name,value in zip(["s_type","fs_type","p_clay","p_sand","p_gravel","entropy","var"],expected):
            compare("uq_"+name,np.load(run/f"04/{name}.npy"),value)
    result=dict(reference="Direct original functions reading the same actual GSLIB output files",checks=checks,passed=all(r["passed"] for r in checks.values()))
    dump(run/"numerical_regression.json",result)
    if not result["passed"]: raise AssertionError("Original-function regression mismatch")
    print(f"All {len(checks)} direct original-function regression checks passed")


if __name__=="__main__": main()
