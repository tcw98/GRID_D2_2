# GRID D2.2 — Geological modelling and uncertainty quantification

This repository provides a staged workflow for probabilistic subsurface modelling and uncertainty analysis from borehole data. It connects data quality assessment, geological interpretation, variography, conditional simulation, ensemble summaries and validation. The methodological background is provided in [Deliverable D2.1](docs/reports/D2_1_reference.docx).

Continuous geological boundaries are represented with sequential Gaussian simulation (SGSIM), and categorical lithology with sequential indicator simulation (SISIM). A hierarchical combination produces geological realizations from which class probabilities, the most likely class, entropy and variability are calculated.

The workflow is illustrated with a Munich borehole dataset. This reference case demonstrates the methods; its domain, geological interpretation and covariance parameters are not universal defaults for other sites.

## Workflow

```text
Borehole records
  -> quality assessment and geological preparation
  -> statistical preparation and variography
  -> continuous and categorical conditional simulation
  -> hierarchical geological realizations
  -> uncertainty summaries and validation
```

| Stage | Script | Main products |
| --- | --- | --- |
| 01 | `scripts/01_Data_Quality_and_Preparation.py` | Quality flags, composites, boundary observations, trends and transformed residuals |
| 02 | `scripts/02_Variography.py` | Directional experimental variograms, candidate fits and adopted parameters |
| 03 | `scripts/03_Simulation.py` | Continuous surfaces, categorical fields and paired geological realizations |
| 04 | `scripts/04_Uncertainty_Quantification.py` | Class probabilities, mode, entropy, variability, boundary statistics and stability diagnostics |
| 05 | `scripts/05_Validation.py` | Hard-data checks, spatial continuity diagnostics and submodel cross-validation |

See the [function mapping](docs/MAPPING.md) for the relationship between the staged scripts and the original numerical routines.

## Installation

The supported execution environment is Windows with Python 3.12. Run these commands from the repository root:

```powershell
py -3.12 -m venv .venv
..venvScriptspython.exe -m pip install -r requirements.txt
..venvScriptspython.exe -B -m unittest discover -s tests -v
```

Simulation requires an external GSLIB90 installation containing `varmap.exe`, `gamv.exe`, `gam.exe`, `sgsim.exe` and `sisim.exe`. These executables are not redistributed here.

## Run the workflow

A small-grid smoke profile exercises the staged implementation:

```powershell
$env:GSLIB_DIR = 'C:	oolsgslib90'
..venvScriptspython.exe run_workflow.py --config configs/smoke.json --run-dir runs/smoke --verify-preparation --sweep
..venvScriptspython.exe scripts/verify_run.py --config configs/smoke.json --run-dir runs/smoke
```

The smoke profile uses the included real borehole subset, a 24 × 24 × 40 grid and 100 realizations. Cross-validation runs by default; add `--skip-cv` to omit it. The `--sweep` flag evaluates the supplied continuous lag/tolerance combinations.

For a larger model, select an explicit configuration and a fresh run directory:

```powershell
..venvScriptspython.exe run_workflow.py --config configs/example_1km.json --run-dir runs/example_1km --verify-preparation
```

Use `--wd` to select an input-data directory containing `data/bohrungen_epsg25832_shp/`. Completed stages save logs and manifests. To resume, use `--from-stage`; failed stage directories must first be archived outside their numbered location. Use a new run directory when changing configuration.

## Adapting the method to another site

Set the coordinate system, domain, discretization and selected variogram parameters in a configuration file. Review the geological interpretation and category mapping separately: the supplied DIN grouping and Quaternary-base extraction are specific to the Munich example. Some search controls, kriging options and realization settings remain implementation constants in `scripts/workflow_support.py`.

Quality flags identify records for review rather than automatically repairing them. Existing preparation exclusions are recorded separately. Conditional realizations express uncertainty under the adopted data interpretation and spatial model; they do not cover every source of geological or engineering uncertainty.

The implementation reports the original **unnormalized natural-log entropy**, with its legacy handling of zero probabilities. This differs from the normalized entropy formulation discussed in D2.1. Variability is `1 - max(p)`. These indicators describe ensemble uncertainty and are not engineering failure probabilities or risk measures.

## Munich reference example

The selected example figures use `configs/exponential_20260916.json`. The baseline and smoke profiles have different covariance settings and are not substitutes for this configuration.

- [Example reproduction guide and parameter tables](docs/article/REPRODUCE.md)
- [Example validation evidence and interpretation limits](docs/article/VALIDATION.md)
- [Available figures](docs/article/FIGURES.md)
- [Borehole provenance and selection](data/real_munich/README.md)
- [Terrain-masked 3D model](vtk/model_gok_exponential.vtk)

Verify the packaged example evidence without running GSLIB:

```powershell
..venvScriptspython.exe scripts/verify_article_evidence.py
```

The example includes its selected figures, statistical tables, sampled terrain and VTK output. Full realization arrays and fold simulations are generated locally. Its grouped hierarchical five-fold validation uses `pics/indicator_cv.py`; it is distinct from stage 05's legacy submodel cross-validation.

## Repository and provenance

`scripts/` contains the staged workflow and numerical routines; `configs/` contains example profiles; `data/` contains the input subset and reference terrain; `tests/` contains regression checks. The interactive `legacy/Workflow_U9_Wang.py` is retained as a reference and must not be executed as the production runner.

English translations affect comment text only in two original source files. Executable tokens and numerical parameters are unchanged. [Source localization records](docs/core_source_localization.json) retain both original and release hashes. Historical run manifests retain the hashes of the versions that produced their results. Machine-local workspace prefixes have been replaced with a labelled placeholder in selected provenance records.

The existing repository [LICENSE](LICENSE), GNU Affero General Public License v3.0, is retained. See [THIRD_PARTY.md](THIRD_PARTY.md) for author notices and source/dependency provenance. Windows CI runs the regression tests and reference-evidence checks; GSLIB integration requires the external executables. [Release verification](docs/RELEASE_VERIFICATION.md) records the checks for this English-language package.
