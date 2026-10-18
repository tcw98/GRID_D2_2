# Reproduce the current article result - 18 October 2026

The current article uses `configs/exponential_sisim_final_20261018.json`. Its zero-nugget four-structure input and the corresponding new SISIM simulations replace the September article presentation. The numerical kernels and continuous SGSIM input remain unchanged.

## Redraw the supplied evidence

Install the repository requirements (Windows, Python 3.12). No GSLIB or new CV is needed to redraw the bundled tables/arrays:

```powershell
python scripts/plot_sgsim_reproduction.py
python scripts/plot_final_variograms.py
python scripts/refresh_final_model_figures.py
python scripts/plot_final_calibration.py
python scripts/publish_article_figures.py
python scripts/verify_final_validation.py
```

The first four commands write canonical outputs in `pics/final_20261018/`. The fifth replaces the current root `pics/` article assets and updates `pics/current_article.json`. See [FIGURES.md](FIGURES.md).

## Run the current simulation input

Supply an external GSLIB90 directory and use a fresh output directory:

```powershell
$env:GSLIB_DIR = 'C:\tools\gslib90'
python run_workflow.py --config configs/exponential_sisim_final_20261018.json --run-dir runs/article_final --verify-preparation --skip-cv
```

This uses the existing five-stage workflow. Original terrain masking and article postprocessing are separate steps. Stage 05's legacy CV is not the grouped validation behind the article calibration curves; the flag prevents silently substituting that design.

The [detailed replay guide](../final_20261018/REPRODUCE.md) explains the selected-array postprocessing and the independent adapter for the same saved grouped five-fold split. Replaying the exact CV requires the external original fold training/SGSIM/transform records and an external native SISIM executable. They are not all reconstructed by the stage-05 runner, and matching seeds with another executable build is not a bitwise-equivalence guarantee.

The supplied current VTK is `vtk/model_gok_sisim_final_20261018.vtk.gz`. Decompress it with a gzip-capable tool and open the resulting VTK in a visualization application. Its decompressed SHA-256 equals the selected model record. No numerical conversion or resampling is involved.

## Historical records

The September config, original validation tables, audit and original VTK are historical evidence retained for traceability; they are not current-model outputs. The earlier reproduction documentation remains available in [Git history](https://github.com/tcw98/GRID_D2_2/blob/d5d2df25032e991b6f1c4ac566734c0a087f2eaf/docs/article/REPRODUCE.md). The old consistency check remains available as `python scripts/verify_article_evidence.py`.

## SGSIM spatial reproduction

The added SGSIM plot reuses the 100 original Gaussian residual fields in EW/NS, keeps the non-zero nugget and scales all series by one input sill. See the [full scope and lag definitions](../final_20261018/README.md#continuous-sgsim-reproduction).

Default redraw uses the bundled compact arrays. To recalculate them from the same external production run without new simulation:

```powershell
python scripts/plot_sgsim_reproduction.py --source-run <original-run>
python scripts/publish_article_figures.py
```
