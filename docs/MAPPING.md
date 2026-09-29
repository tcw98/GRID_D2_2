# Function mapping

| Current function/block | Proposed location | Treatment |
| --- | --- | --- |
| `LFU_2_YOU` | 01 geological input loading | Keep as-is; raw audit runs first |
| Missing IDs/records, duplicate rows, coordinate/CRS audit | 01 `audit_raw` | New QC; report only, stop on structurally unusable input |
| Collar elevation neighbourhood residual | 01 `local_plausibility` | New QC; excludes the same borehole from neighbours; no correction |
| Depth ordering, gaps/overlap/thickness, Endteufe | 01 `audit_intervals` | New QC; legacy sorting is recorded, not replaced |
| Unknown/ambiguous DIN and missing descriptions | 01 `audit_intervals` | New flags against existing classification rules |
| `qnd_compositing`, `len_wei_compositing` | 01 `geological_preprocessing` | Wrapper; retain Munich interpretation |
| Base-vs-collar/terminal depth and local base plausibility | 01 `audit_base_geometry`, `local_plausibility` | New QC flags |
| `moving_average_2d_kdtree`, `detrend_2D`, `transform` | 01 `statistical_preprocessing` | Wrapper; retain numerical operations and exclusions |
| `GSLIB_VARMAP`, `GSLIB_GAMV`, readers | 02 + shared helpers | Original writers/readers; relative parameter-file paths |
| `fit_2_curve`, `MAKENESTEDVARGM_Major/Minor` | 02 | Keep as-is; fit candidates separated from adopted parameters |
| Continuous lag/tolerance loop | 02 `--sweep` | Preserve original 25 combinations and bandwidth relation |
| `GSLIB_SGSIM`, `RUN_GSLIB`, `readsgsim`, `post_processing_2D` | 03 + shared helpers | Wrapper; original seed, writers and backtransform/trend order |
| `GSLIB_SISIM`, `readsisim`, batch scheduling | 03 + shared helpers | Original writers/readers; retain 20 parts × 5 and seed progression |
| `to_3D`, `combine_sim_ohnelm` | 03 | Call originals one realization at a time; preserve pairing/order |
| `evaluate_nmodel` | 04 `legacy_uq` | Original function on voxel chunks; original entropy unchanged |
| Continuous mean/variance/percentiles and prefix stability | 04 | New summaries; realizations unchanged |
| Uniform-grid export | 04 | Geometry-only PyVista `ImageData` compatibility wrapper |
| Hard-data and array consistency checks | 05 | New explicit reports; no input/output correction |
| `GSLIB_GAM` and readers | 05 variogram reproduction | Keep as-is; horizontal and vertical diagnostics |
| `split_array`, `split_array_3D`, `save_train`, `save_train_3D` | 05 cross-validation | Existing split/preprocessing; explicit saved split seed |
| `k_fold_cross_validation`, categorical CV and metrics | 05 | Existing evaluators; saved fold inputs and overlap diagnostics |
| `Workflow_U9_Wang.py` | legacy/ | Byte-identical reference; experimental alternatives retained |

The legacy already removes NaN/zero collar elevations, selects a spatial domain, sorts intervals, groups DIN codes, composites and trims detrending residuals. These operations did not provide a unified pre-deletion quality audit. Stage 01 adds that audit, including duplicate identity/row issues, geometry and CRS, spatial collar plausibility, interval integrity, unknown classifications and continuous-base flags. It reports swallowed legacy compositing diagnostics without repairing the affected data. Detrending remains modelling preprocessing, separate from QC.

General workflow: QA/QC → compositing → variography → simulation → UQ → validation. Case-specific interpretation: Quaternary-base extraction, DIN grouping, collar plausibility range and selected geological variograms. Other regions require an explicitly reviewed interpretation and configuration.
