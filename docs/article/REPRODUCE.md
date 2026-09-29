# Reproduce the Munich reference example

The manuscript figures use `configs/exponential_20260916.json`. The other profiles remain useful for baseline/smoke regression but do not reproduce the article's selected covariance models. Numerical implementation is unchanged; release localization modifies comments only (see docs/core_source_localization.json). This guide describes the prepared reproduction route; the 29 September review checked stored evidence and the existing tests, without repeating the expensive 100-member production ensemble or five-fold simulations.

## Environment and supplied inputs

Use Windows and Python 3.12. Install the requirements and supply an existing GSLIB90 directory with `varmap.exe`, `gamv.exe`, `gam.exe`, `sgsim.exe` and `sisim.exe`. No GSLIB executable is redistributed. The article calculation used the executable hashes recorded in `recorded_run/03/run_manifest.json`; matching seeds alone does not ensure bitwise equivalence with different builds.

Run the following from the repository root. Choose a short ASCII work directory for the GSLIB parameter-file paths. The working directory must not already exist.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-figures.txt
$env:GSLIB_DIR = 'C:\tools\gslib90'
.\.venv\Scripts\python.exe scripts/verify_article_evidence.py
.\.venv\Scripts\python.exe scripts/prepare_article_workdir.py --work-dir runs/article
.\.venv\Scripts\python.exe run_workflow.py --config runs/article/exponential.json --run-dir runs/article/run --to-stage 4 --verify-preparation
```

This computes 100 realizations on 2,000,000 voxels each. It requires substantial local disk space and runtime. Raw GSLIB output is much larger than the compact byte-valued realization arrays. `runs/` is excluded from Git. Stage 01 uses the bundled `data/real_munich` records by default. `--wd` may select the complete source archive, but must point to the parent containing `data/bohrungen_epsg25832_shp/`.

`prepare_article_workdir.py` supplies the exact sampled GOK and boolean mask used for the manuscript, plus the saved directional fit record. Full 1 m terrain TIFFs are external; their original names and hashes are in `data/article/terrain.json`. The mask keeps cells whose **top face** is at or below the sampled terrain, leaving 857,561 cells. Stages 01–04 themselves retain the full elevation slab.

## Article-specific grouped five-fold validation

```powershell
.\.venv\Scripts\python.exe pics/indicator_cv.py --work-dir runs/article --repo . --gslib $env:GSLIB_DIR
.\.venv\Scripts\python.exe pics/make_revision_figures.py --work-dir runs/article --repo . --out runs/article/figures --part cv
```

This is the source of the manuscript's hierarchical three-class CV. It groups entire borehole IDs and connected original/composited XY grid columns; seed 12345; 100 realizations per fold. Variograms remain fixed. Fold trends and Quantile transforms use training observations. Production preprocessing uses Johnson instead, an inherited methodological difference. `run_workflow.py` stage 05 uses separate legacy submodel CV and is **not** the source of the manuscript's Figure 11.

## Regenerate the selected result figures

```powershell
.\.venv\Scripts\python.exe pics/make_revision_figures.py --work-dir runs/article --repo . --out runs/article/figures --part fits
.\.venv\Scripts\python.exe pics/make_uq_figures.py --work-dir runs/article --repo . --out runs/article/figures --part stability
.\.venv\Scripts\python.exe pics/make_uq_figures.py --work-dir runs/article --repo . --out runs/article/figures --part sections --x-cut 691000 --y-cut 5334800
```

These commands reproduce the numerical content of individual figures; the manuscript includes manually assembled panels and inherited illustrations. See `FIGURES.md`. The optional `--part map` command downloads OpenStreetMap tiles and requires network access. Existing map images are already supplied. The saved `vtk/model_gok_exponential.vtk` is delivered as a ready-to-view, hash-checked artifact, without requiring a rerun. The `--part vtk` script uses this delivered file as a geometry template.

The fit plot uses the recorded directional fits of the selected article run. It is distinct from the explicit, consolidated simulation parameters. Do not silently substitute automatically generated fitting candidates into the simulation profile.

## Additional model checks and restart

```powershell
.\.venv\Scripts\python.exe run_workflow.py --config runs/article/exponential.json --run-dir runs/article/run --from-stage 5 --skip-cv
.\.venv\Scripts\python.exe -B -m unittest discover -s tests -v
```

Stage 05 adds hard-data and variogram-reproduction diagnostics. These were not part of the saved September exponential article run and must not be advertised as already passed. Numerical consistency does not imply geological acceptance. Existing numbered stage folders are not overwritten; archive a failed stage outside its numbered path before retrying. If a CV fold fails, inspect its log and archive the incomplete `run/cv` directory before starting a fresh CV run; the inherited script waits on existing incomplete fold directories.

## Fixed settings that are not JSON options

| Setting | Article implementation |
| --- | --- |
| Grid | EPSG:25832; angle 0; outer X 690400–691400, Y 5334500–5335500, elevation 475–575 m |
| Cells | 100 × 100 × 200; 10 × 10 × 0.5 m; first centre (690405, 5334505, 475.25) |
| SGSIM | SK (`ktype=0`), 100 members, seed 12345; `ndmin=1`, `ndmax=12`, `ncnode=12`, `sstrat=1`, `multgrid=1`, `nmult=3`; input trim [-3,3] |
| SISIM | Full IK (`mik=0`), OK (`ktype=1`), 20 parts × 5; seed 69069+11×part; `ndmax=12`, `ncnode=12`, `sstrat=1`, `multgrid=0`, `noct=5` |
| SISIM priors | Fixed 0.5/0.5; code 0=clay group, 1=sand |
| Final classes | 1=clay group, 2=sand, 3=Quaternary/gravel |
| Pairing | Member i with member i; X fastest, then Y, then Z |
| Production trend | Radius 300 m; at most 100 neighbours; exp(-distance/100 m); residual exclusion >5 m |
| Entropy | Natural logarithm, unnormalized; zero terms replaced with 0.0001 only inside entropy calculation |

| Simulation covariance | Nugget | Structure contribution | Major N–S range (m) | Minor E–W range (m) | Vertical range (m) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Boundary: exponential | 0.136955506 | 0.763377625 | 63.942492 | 29.514768 | 6 (2D interface parameter) |
| Indicator: exponential | 0.025146356 | 0.117679560 | 64.879576 | 52.981169 | 7.120614 |

The total sill is nugget plus structure contribution. JSON stores the full precision. Indicator consolidation uses the original 90% horizontal/10% vertical sill/nugget weighting after converting standardized vertical semivariance back to raw units. The clay group combines several DIN fine-grained classes. The original compositing routine also retains a small XY offset for nominally vertical boreholes; the CV grouping accounts for both original and composited columns.
