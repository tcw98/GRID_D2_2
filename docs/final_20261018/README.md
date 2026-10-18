# Current calibrated SISIM article record - 18 October 2026

This record documents the current selected `response_step_50` input and its actually generated outputs. It retains the original September article example, source localization records, tests, license, terrain and historical VTK. Current article figures replace their root `pics/` paths; their canonical outputs and data are under `pics/final_20261018/`. `pics/current_article.json` identifies the current figure set.

## Selected input and actual calculations

`configs/exponential_sisim_final_20261018.json` changes only the profile name and categorical variogram relative to `configs/exponential_20260916.json`. The new categorical input has zero nugget and four positive exponential contributions 0.05/0.04/0.04/0.12. The calibration search held that four-structure layout and contributions fixed and adjusted EW/NS/Z ranges. The archived September profile instead has one categorical structure and a nugget of 0.025146355593272343; it is preserved as published. Explicitly select the new configuration to use the calibrated input. The continuous SGSIM input, conditioning procedure, search settings, ensemble ordering and original numerical simulation functions were retained.

The selected input received 100 spatial realizations. It subsequently received a new whole-borehole/group five-fold SISIM CV with **100 realizations per fold: 500 new realizations, 100 native batches**. Existing fold assignments and training conditions were reused; the native signatures, output dimensions, completion logs and protected input hashes were checked. `native_batches_verified.csv` records all batch hashes. `holdout_audit.json` checks the retained holdouts, including zero shared training/test voxels in every fold. Full raw simulation arrays/logs are external generated outputs, not bundled GitHub inputs.

## Spatial continuity

The Figure 8 style panels compare fixed experimental points, the selected theoretical input, one reproduced mean and one 5th-95th percentile envelope per direction. Horizontal nominal lag is 10 m, vertical 0.5 m. Experimental zero/near-zero records remain in the CSVs but are not displayed as regular lags. The reproduced mean/envelope includes the identity gamma(0)=0; this is not an added measured or fitted datum. Axes retain metres horizontally and normalized semivariance 0-1.2 vertically.

Shape panels use each reliable variance/theoretical sill as a unit reference. Formal panels divide every curve by the same reference variance 0.24815277295195295. No endpoint/maximum rescaling is used. Common-reference normalized RMSE is EW 0.13585777, NS 0.06213126 and Z 0.07241337; frozen J=0.10966702. Compared with the preceding candidate the objective improved modestly, but directional residuals remain. These results do not establish satisfactory reproduction of every spatial feature or a universal acceptance threshold.

- [Shape comparison, actual Neogene region](../../pics/final_20261018/08_variogram_neogene_shape.png)
- [Common-reference spatial validation](../../pics/final_20261018/08_variogram_neogene_common_variance.png)
- [Updated entropy sections](../../pics/final_20261018/13_sections_entropy.png)
- [Updated variability sections](../../pics/final_20261018/14_sections_variability.png)

## Probability calibration

[Wiegel et al. (2024)](https://doi.org/10.1002/gete.202400011), Figures 9 and 7a, inform the two displays. This project retains five folds for both submodels.

- **Binary SISIM:** 22,313 withheld observed Neogene clay/sand composites, 694 boreholes, 623 groups. Sand probability is the raw 0/1 SISIM frequency over 100 realizations, not hierarchical p_sand. The 50 fixed 2% probability bins show actual sand frequency as a staircase, with one bisector and no older-model overlay. Bins are (lower, upper], including p=0 in the first bin. All bins are populated; the smallest has 179 samples. There is no smoothing or artificial data insertion. ECE is 0.04422624, unweighted bin MAE 0.05229931, Brier 0.13174859 and binary classification accuracy 0.81172411 (ties classified as clay).
- **SGSIM:** the same saved 664 withheld boundary observations and original 100-member fold ensembles are reused without new SGSIM. The ECDF of their unchanged mid-PIT ranks is compared with the bisector, following the cumulative-quantile presentation. The existing nominal-quantile table is preserved, with mean absolute calibration gap 0.01085923. This ECDF display differs slightly from the previous 5%-quantile-point display because of finite-ensemble rank discretization; it is not a change to predictions.

SISIM still has local miscalibration: in the 0-10% prediction region, mean predicted sand probability is 0.0430 versus observed frequency 0.1427 (3441 samples); in the 40-80% region it is 0.5976 versus 0.5269 (4333 samples). These deviations are not removed or explained away as universally sparse bins. The CV was not used to tune the selected ranges further.

- [Binary-only SISIM calibration](../../pics/final_20261018/figure9_sisim_calibration.png)
- [Unchanged SGSIM data, cumulative calibration](../../pics/final_20261018/figure9_sgsim_calibration.png)
- [Two-panel calibration figure](../../pics/final_20261018/figure9_final_calibration.png)

This is fixed-full-data-variogram, whole-borehole/group conditional CV, not nested variogram-selection validation. Composite weighting is retained; samples are not independent equal-weight boreholes. The original production/CV preprocessing distinction remains. The original model VTK is retained as baseline evidence and does not represent this new categorical input. The selected-input VTK is supplied losslessly compressed as `vtk/model_gok_sisim_final_20261018.vtk.gz`; its decompressed hash is in the new model provenance.

All new simulation/UQ panels use the selected 100 SISIM fields paired with unchanged original SGSIM surfaces. Original entropy, class tie-breaking, variability and threshold rules are retained. The numerical-run utility hash and current English-release hash differ only in comments; this equivalence was independently checked and is documented in the existing source-localization records.

See [REPRODUCE.md](REPRODUCE.md) for commands. The English release, license, source data and original evidence remain in place. No new license rights are inferred for source data, figures or third-party materials.
