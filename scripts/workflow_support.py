"""Small I/O and ordered-parameter helpers; numerical kernels stay in UTILITIES."""
from __future__ import annotations
import argparse
import contextlib
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import UTILITIES_Wang as U
import GSLIB_parfl_011 as G

ROOT = Path(__file__).resolve().parents[1]
GRID_KEYS = ("nx", "ny", "nz", "xsiz", "ysiz", "zsiz", "xmn", "ymn", "zmn")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""): h.update(block)
    return h.hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)), encoding="utf-8")


def load_config(path):
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    grid = cfg["preparation"]["grid"]
    cfg["preparation"]["grid"] = {key: grid[key] for key in GRID_KEYS}
    for key, value in dict(radius_xy=[300,300],tolerance=20,outlier_tolerance=5,trans_method="Johnson").items():
        cfg["preparation"].setdefault(key,value)
    if cfg["ensemble"] != {"parts": 20, "per_part": 5, "seed_base": 69069, "seed_step": 11}:
        raise ValueError("Compatibility profile requires the legacy 20 x 5 SISIM ensemble and seeds")
    return cfg


def parser(description, executables=False):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--config", type=Path, default=ROOT / "configs" / "example_1km.json")
    if executables:
        p.add_argument("--gslib-dir", type=Path, default=os.environ.get("GSLIB_DIR"))
    return p


@contextlib.contextmanager
def stage(args, name, needs=()):
    run = args.run_dir.resolve()
    cfg = load_config(args.config.resolve())
    for required in needs:
        manifest_path = run / required / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["status"] != "complete": raise ValueError(f"Incomplete stage: {required}")
        recorded = manifest.get("config_digest")
        if recorded and recorded != digest(args.config.resolve()):
            raise ValueError(f"Config differs from completed stage {required}; use a new run")
    out = run / name
    out.mkdir(parents=True, exist_ok=False)
    old = Path.cwd()
    gslib = getattr(args, "gslib_dir", None)
    if gslib is not None: args.gslib_dir = gslib.resolve()
    manifest = dict(status="started", config=cfg, config_digest=digest(args.config.resolve()),
                    python=sys.version, core_hashes={f: digest(ROOT / "scripts" / f) for f in ["UTILITIES_Wang.py", "GSLIB_parfl_011.py"]})
    try:
        os.chdir(run)
        with (out / "stage.log").open("w", encoding="utf-8") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            yield cfg, Path(name), manifest
        manifest["status"] = "complete"
    except Exception as exc:
        manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        os.chdir(old)
        dump(out / "run_manifest.json", manifest)
    print(f"{name}: complete ({out})")


def executable(directory, name):
    if directory is None: raise ValueError("Set --gslib-dir or GSLIB_DIR to your local GSLIB90 directory")
    path = Path(directory) / (name + (".exe" if os.name == "nt" else ""))
    if not path.is_file(): raise FileNotFoundError(path)
    return path


def run_gslib(directory, name, parfile):
    # Original RUN_GSLIB call, with Windows command-line quoting at the boundary.
    # All file names INSIDE the parameter file are short relative ASCII paths.
    path = executable(directory, name)
    if os.name != "nt":
        raise RuntimeError("Validated runner is Windows/GSLIB90; no silent backend substitution")
    U.RUN_GSLIB('"' + str(path) + '"', '"' + str(parfile) + '"')


def gamv_params(categorical=False, vertical=False):
    directions = [[90, 30, 30, 0, 2, 1], [0, 30, 30, 0, 2, 1]] if categorical else [[90, 45, 100, 0, 0, 2], [0, 45, 100, 0, 0, 2]]
    if vertical: directions = [[0, 30, 30, 90, 2, 5]]
    return dict(icolx=1, icoly=2, icolz=3 if categorical else 0, nvar=1,
                ivar1=4 if categorical else 3, tmin=-1.e21, tmax=1.e21,
                nlag=30 if categorical else 100, xlag=.5 if vertical else 5,
                lagtol=.25 if vertical else 2.5, ndir=len(directions), xdir=directions,
                standardize=1 if vertical else 0, nvarg=1, ivtail=1, ivhead=1, ivtype=1)


def varmap_params(categorical=False):
    if categorical:
        dims, sizes, lags, distances = (50,50,30), (10,10,.5), (20,20,20), (5.,5.,.5)
    else:
        dims, sizes, lags, distances = (20,20,1), (5.,5.,1), (50,50,0), (5,5,1)
    return dict(zip(["nx","ny","nz","xsiz","ysiz","zsiz","icolx","icoly","icolz","nvar","ivar1","tmin","tmax","igrid","nxlag","nylag","nzlag","dxlag","dylag","dzlag","minpairs","standardize","nvarg","ivtail","ivhead","ivtype","cut"],
        [*dims,*sizes,1,2,3 if categorical else 0,1,4 if categorical else 3,-1.e21,1.e21,0,*lags,*distances,2,0 if categorical else 1,1,1,1,1,.5]))


def sgsim_params(cfg, model, out):
    grid = copy.deepcopy(cfg["preparation"]["grid"])
    grid["nz"] = 1
    return dict(icolx=1,icoly=2,icolz=0,icolvr=3,icolwt=0,icolsec=0,tmin=-3,tmax=3,
        itrans=0,transfl="none.trn",ismooth=0,smthfl="none.dat",icolvrsmthfl=3,icolwtsmthfl=0,
        zmin=0,zmax=40,ltail=1,ltpar=0,utail=1,utpar=1,idbg=1,dbgfl=str(out)+".dbg",
        nsim=100,grid_params=grid,seed=12345,ndmin=1,ndmax=12,ncnode=12,sstrat=1,
        multgrid=1,nmult=3,noct=0,radiushmax=300,radiushmin=200,radiusvert=10,
        sang1=0,sang2=0,sang3=0,covtab1=60,covtab2=60,covtab3=60,ktype=0,rho=.6,varred=1,
        secfl="none.dat",icolsecfl=0,nst=model["nst"],c0=model["c0"],variogram=model["nest"])


def sisim_params(cfg, model):
    # W lines 1063-1074, before the experimental MIK parameter overwrite.
    return dict(vartype=0,ncat=2,catx=[0,1,"","",""],gpdfx=[.5,.5,"","",""],
        icolx=1,icoly=2,icolz=3,icolvr=4,directik="direct.ik",icolsx=1,icolsy=2,icolsz=3,
        icol=[4,5,6,7],imbsim=0,bz=[1,2,3,4,5],tmin=0,tmax=1.e21,zmin=0,zmax=1,
        ltail=1,ltbar=0,middle=1,midpar=1,utail=1,utpar=1,tabfl="cluster.dat",icolvrt=3,
        icolwtt=0,idbg=1,nsim=5,grid_params=copy.deepcopy(cfg["preparation"]["grid"]),
        seed=12345,ndmax=12,ncnode=12,maxsec=0,sstrat=1,multgrid=0,nmult=0,noct=5,
        radiushmax=100,radiushmin=100,radiusvert=10,sang1=0,sang2=0,sang3=0,
        covtab1=50,covtab2=50,covtab3=50,mik=0,mikcat=1,ktype=1,
        variogram=[[model["nst"],model["c0"],model["nest"]]]*2)


def simulate_gaussian(cfg, model, data, stem, gslib):
    stem = str(stem)
    params = sgsim_params(cfg, model, stem)
    G.GSLIB_SGSIM(str(executable(gslib,"sgsim")), str(data), stem+".par", stem+".out", params)
    run_gslib(gslib,"sgsim",stem+".par")
    return params


def simulate_indicators(cfg, model, data, folder, gslib):
    """Same 20 parts x 5, seeds and four batches of five as legacy runner.

    Parameter writer and executable unchanged. Added checked process completion
    and safe command quoting. Read each original output in numeric part order.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    params = sisim_params(cfg, model)
    for part in range(20):
        params["seed"] = 69069 + 11 * part
        stem = str(folder / f"sisim_{part}")
        G.GSLIB_SISIM(str(executable(gslib,"sisim")),stem+".par",str(data),stem+".dbg",stem+".out",params)
    for batch in range(4):
        with ThreadPoolExecutor(max_workers=5) as pool:
            jobs = [pool.submit(run_gslib,gslib,"sisim",folder/f"sisim_{part}.par") for part in range(batch*5,batch*5+5)]
            for job in jobs: job.result()
    grid = cfg["preparation"]["grid"]
    ncell = grid["nx"]*grid["ny"]*grid["nz"]
    simulations = np.lib.format.open_memmap(folder/"realizations.npy", mode="w+", dtype=np.uint8, shape=(100,ncell))
    for part in range(20):
        values = U.readsisim(str(folder/f"sisim_{part}.out"))
        if values.size != 5*ncell or not np.isin(values,[0,1]).all():
            raise ValueError(f"Invalid SISIM output in part {part}: shape/categories")
        simulations[part*5:part*5+5] = values.reshape(5,ncell)
    simulations.flush()
    return simulations


def read_variogram(gslib, data, stem, params):
    stem = str(stem)
    G.GSLIB_GAMV(str(executable(gslib,"gamv")),stem+".par",str(data),stem+".out",params)
    run_gslib(gslib,"gamv",stem+".par")
    tables = U.read_GAMV_GAM(stem+".out",params)
    for i, table in enumerate(tables): table.to_csv(stem+f"_direction_{i}.csv",index=False)
    return tables
