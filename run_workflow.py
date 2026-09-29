"""Run the five standalone scripts in order; no package installation required."""
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",type=Path,default=root/"configs/example_1km.json")
    p.add_argument("--wd",type=Path,default=root/"data/real_munich")
    p.add_argument("--run-dir",type=Path,required=True)
    p.add_argument("--gslib-dir",type=Path,default=os.environ.get("GSLIB_DIR"))
    p.add_argument("--from-stage",type=int,choices=range(1,6),default=1)
    p.add_argument("--to-stage",type=int,choices=range(1,6),default=5)
    p.add_argument("--skip-cv",action="store_true")
    p.add_argument("--sweep",action="store_true")
    p.add_argument("--verify-preparation",action="store_true")
    a=p.parse_args()
    if a.from_stage>a.to_stage: p.error("from-stage must not exceed to-stage")
    scripts=["01_Data_Quality_and_Preparation.py","02_Variography.py","03_Simulation.py","04_Uncertainty_Quantification.py","05_Validation.py"]
    for i in range(a.from_stage,a.to_stage+1):
        cmd=[sys.executable,"-B",str(root/"scripts"/scripts[i-1]),"--config",str(a.config.resolve())]
        if i==1:
            cmd += ["--wd",str(a.wd.resolve()),"--output",str(a.run_dir.resolve()/"01")]
            if a.verify_preparation: cmd += ["--verify"]
        else: cmd += ["--run-dir",str(a.run_dir.resolve())]
        if i in [2,3,5] and a.gslib_dir: cmd += ["--gslib-dir",str(a.gslib_dir.resolve())]
        if i==2 and a.sweep: cmd += ["--sweep"]
        if i==5 and not a.skip_cv: cmd += ["--cv"]
        subprocess.run(cmd,check=True)


if __name__=="__main__": main()
