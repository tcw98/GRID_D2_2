> Historical record: source hashes below refer to the original calculation versions. Comment-only English localization and current release hashes are documented in `core_source_localization.json`.

# Historical synthetic verification

This file records the 12 September synthetic verification only. The default example was replaced by real boreholes on 13 September. Current evidence is in [real_example/VALIDATION.md](real_example/VALIDATION.md).

# Numerical compatibility and verification

## Preserved baseline

`core_hashes.json` records SHA-256 hashes of both original numerical modules and the unedited workflow. Tests enforce these hashes. The workflow reference contains interactive experiments and incomplete statements, so it is deliberately not imported or used as an executable historical baseline.

The profiles select its SGSIM plus full-IK SISIM branch, before later MIK/alternative experiments. The original reference reassigns parameter variables across experiments; it does not define a unique final production profile. The selected categorical model has `nst=1`, `c0=0.068`, and nest `[2, 0.125, 0, 0, 0, 50, 50, 6]`. The four continuous nests and all profile values are explicit in the JSON configuration. The full-IK setup corresponds to the pre-MIK parameter block around original lines 1063–1074. Fit candidates do not replace this baseline automatically.

SGSIM retains 100 realizations and seed 12345. SISIM retains 20 parts of 5 realizations, seed `69069 + 11 * part`, in four batches of five parts. Parameter writers and `RUN_GSLIB` are unchanged; wrappers quote executable paths and use relative ASCII paths inside parameter files. Numeric part order is retained. Streaming changes memory/storage use, with original combination and UQ functions still called directly.

Entropy retains natural logarithms, zero probability replacement by 0.0001, and original class tie rules. It can exceed 1 and is slightly positive even for a single certain class. No normalized entropy is introduced. New continuous statistics use NumPy linear percentiles and population variance.

## Executed checks (12 September 2026)

The complete five-stage smoke run used bundled synthetic boreholes, 24 × 24 × 40 voxels (23,040 per realization), 100 realizations, actual local GSLIB90 executables, the 25-case sensitivity sweep and two-fold validation for each submodel. Execution used Windows and Python 3.12.3 with the pinned requirements. Saved evidence is in `validation/`.

| Check | Result |
| --- | --- |
| Python QC/numerical contract tests | 10 passed |
| Synthetic stage 01 vs direct original calls | 14 passed; numeric differences zero |
| Saved GSLIB outputs, backtransform, combination and seven UQ fields vs direct originals | 11 passed; differences zero |
| Output shape, classes, probabilities, UQ and quantile checks | 8 passed |
| Continuous hard conditioning | 9 in-grid observations reproduced within 0.001 in Gaussian space |
| Indicator hard conditioning | 84 in-grid observations reproduced exactly |
| Final hierarchy hard conditioning | 360 in-grid observations reproduced exactly |
| Continuous fold validation | 9 evaluated; legacy ensemble MAE 0.5463778370640049 m |
| Categorical fold validation | 84 evaluated; accuracy 0.6547619047619048 |

The categorical test folds contain 74 and 10 in-grid points, with zero shared training/test XY locations and voxels. Outside-grid conditioning points are explicitly reported and are not counted as reproduction failures. Horizontal/vertical variogram reproduction figures are diagnostics; the synthetic input is not calibrated to the retained Munich variograms. In particular, the saved vertical reproduction figure shows discrepancies, so no blanket variogram-reproduction pass is claimed.

`real_data_stage01_comparison.json` additionally records the earlier real-borehole 1 km preprocessing comparison: qbasis (2087 × 3), tertiary_points (94498 × 4), coords_d (1340 × 3), and other intermediates/transform parameters matched direct original calls. Full 1 km stages 02–05 have not been executed. The smoke run verifies numerical orchestration and does not establish historical full-model equivalence or geological accuracy.

## Interpretation limits retained explicitly

- Cross-validation evaluates continuous boundary and binary indicator submodels separately, rather than a jointly held-out hierarchical model. Variograms are fixed, not refitted inside each fold.
- The old `save_train` actually performs a Quantile transform in CV; this behaviour is preserved even though production preprocessing uses Johnson. Its extra residual filtering also remains. Continuous MAE averages realization-wise absolute errors, not only ensemble-mean errors.
- Legacy fold splitting had no fixed workflow seed. A saved split seed of 12345 makes this new orchestration reproducible; it cannot recreate unknown historical folds.
- The volume is an absolute elevation slab in rotated model coordinates. No terrain mask/DEM is bundled, so above-ground voxels remain in volume summaries.
- Prefix stability compares 20/50/100 members from the same ensemble. It is a convergence diagnostic, not evidence of model correctness or a calibrated risk threshold.
- Local QC moving averages flag residuals only. Legacy model filtering remains active and is recorded separately. Insufficient neighbourhood support is reported instead of fabricating a reference.

## Repeat the regression

After the README smoke command, run:

```powershell
python scripts/verify_run.py --config configs/smoke.json --run-dir runs/smoke
python -m unittest discover -s tests -v
```

`verify_run.py` reads the same actual GSLIB files through the original all-at-once functions and compares against saved staged outputs. It intentionally limits grid size because those reference calls require more memory. Executable hashes are recorded in each relevant run manifest and in the bundled validation summary. Identical seeds alone do not guarantee identical results across different executable builds.
