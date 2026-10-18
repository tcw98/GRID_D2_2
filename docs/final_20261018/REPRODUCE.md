# Reproduce the calibrated-input supplement

Run these from the repository root after installing the existing requirements. Plot redraw needs NumPy, pandas, SciPy and Matplotlib and does not invoke GSLIB:

```powershell
python scripts/plot_sgsim_reproduction.py
python scripts/plot_final_variograms.py
python scripts/refresh_final_model_figures.py
python scripts/plot_final_calibration.py
python scripts/publish_article_figures.py
python scripts/verify_final_validation.py
```

Canonical output stays under `pics/final_20261018/`; `publish_article_figures.py` then replaces the current root article figures. Historical numerical evidence and its checks remain separate. The plotting scripts use the small bundled figure arrays/tables, not a variogram of an ensemble mean. The final 07 and single 08 divide each direction/series by its own plotted maximum; 08 uses only the full binary SISIM grid. These are shape displays and remove amplitude differences. Raw full-grid common-reference RMSE is retained numerically in `pics/final_20261018/data/variogram_peak_normalization.json`. No alternative 08 plots are generated.

For a new production run, explicitly select the new input and a fresh output directory:

```powershell
python run_workflow.py --config configs/exponential_sisim_final_20261018.json --run-dir runs/sisim_final --verify-preparation --skip-cv
```

Supply GSLIB_DIR or --gslib-dir as in the root guide. This runs the existing five-stage implementation. Its stage-05 legacy submodel CV is not the grouped CV used by these figures; --skip-cv avoids silently substituting that design. Terrain masking and publication postprocessing remain separate, as in the original article route.

To regenerate model-dependent plot data from external selected simulation outputs, provide an original run with `03/qbasis_realizations.npy`, the selected 100-member binary SISIM array, the matching terrain folder (`gok.npy`, `below_gok.npy`), and a fresh external cache:

```powershell
python scripts/refresh_final_model_figures.py --source-run <original-run> --sisim-realizations <selected-realizations.npy> --terrain-dir <matching-terrain-folder> --cache-dir <new-cache> --vtk
```

This reuses the exact original numerical combination/UQ functions and does not simulate. The optional VTK is generated in the external cache. No old model VTK is overwritten.

The independent CV adapter is provided to replay the selected categorical input against the **same saved original grouped folds**. It reuses each fold's training-only categorical data and saved SGSIM/transform/trend, with 20 parts x 5 and seeds 69069+11*part. Existing fold outputs are external prerequisites; do not substitute stage-05 random submodel folds. See the baseline article guide if reconstructing that initial grouped run from scratch.

```powershell
python scripts/prepare_final_calibration.py --init --run <original-run> --cv-root <new-cv-output>/response_step_50 --sgsim-validation-dir <original-continuous-validation> --audit-dir <new-audit-folder> --exe <sisim.exe>
python scripts/calibrate_sisim.py cv --run <original-run> --repo . --exe <sisim.exe> --output <new-cv-output> --models docs/final_20261018/selected_variogram.json --nsim 100 --workers 5
python scripts/prepare_final_calibration.py --run <original-run> --cv-root <new-cv-output>/response_step_50 --sgsim-validation-dir <original-continuous-validation> --audit-dir <new-audit-folder> --exe <sisim.exe>
python scripts/plot_final_calibration.py
```

The helper is specific to this frozen model/grid/holdout design. It refuses incomplete batches or changes to other native simulation settings. The published protocol normalizes local path labels while retaining input hashes; `source_record_sha256` identifies the original local frozen record, not the normalized public JSON bytes. The table summary records both native and public protocol hashes.

## SGSIM spatial reproduction

The added SGSIM plot reuses the 100 original Gaussian residual fields in EW/NS, keeps the non-zero nugget and scales all series by one input sill. See the [full scope and lag definitions](README.md#continuous-sgsim-reproduction).

Default redraw uses the bundled compact arrays. To recalculate them from the same external production run without new simulation:

```powershell
python scripts/plot_sgsim_reproduction.py --source-run <original-run>
python scripts/publish_article_figures.py
```
