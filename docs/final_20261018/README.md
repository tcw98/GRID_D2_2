# Model parameters and calculation record

This record describes the `response_step_50` SISIM model, the SGSIM boundary model and their article results. Figures are in `pics/`; redraw inputs are in `pics/final_20261018/` and indexed by `pics/current_article.json`.

## Selected input and calculations

The article configuration is `configs/exponential_sisim_final_20261018.json`. SISIM uses zero nugget and four exponential structures with contributions 0.05/0.04/0.04/0.12. Calibration held these contributions fixed and adjusted the EW/NS/Z ranges. SGSIM uses the continuous model specified in the same configuration.

The spatial ensemble has 100 realizations. Grouped five-fold SISIM CV uses 100 realizations per fold: 500 realizations in 100 native batches. [Batch records](native_batches_verified.csv) contain output hashes; the [holdout audit](holdout_audit.json) records group assignments and zero shared training/test voxels in each fold. Full native realization arrays and logs are generated externally.

## Spatial continuity

### SISIM

Figures 07-08 show the experimental points and theoretical input; Figure 08 adds the mean reproduced variogram and 5th-95th percentile envelope over the complete binary SISIM grid. Horizontal nominal lag is 10 m and vertical lag is 0.5 m. Positive nominal-lag experimental bins are plotted at their actual distances. The reproduced origin is gamma(0)=0. Distances are in metres and the normalized y-axis is 0-1.2.

Within each direction and plotted lag window, experimental, theoretical and reproduced series are divided by their own maxima to compare shape. All realizations use the reproduced-mean peak as the envelope denominator. Scaling factors are in [variogram_peak_normalization.json](../../pics/final_20261018/data/variogram_peak_normalization.json).

Raw full-grid errors use common reference variance 0.24815277295195295:

| Direction | Common-reference normalized RMSE |
| --- | ---: |
| EW | 0.19886292 |
| NS | 0.13333520 |
| Z | 0.09661573 |

- [Experimental and theoretical input](../../pics/07_variogram_indicator_fit.png)
- [SISIM reproduced variogram](../../pics/08_variogram_reproduction.png)
- [Entropy sections](../../pics/13_sections_entropy.png)
- [Variability sections](../../pics/14_sections_variability.png)

### Continuous SGSIM reproduction

SGSIM reproduction uses 100 complete 100 x 100 Gaussian residual fields before trend restoration and back-transform. The two-dimensional interface has EW and NS spatial directions. All series share total input sill 0.9003331304926536, including nugget 0.13695550591961783 (nugget/sill=0.15211648). The theoretical curve shows the positive-distance nugget jump, with gamma(0)=0 displayed separately.

The original experimental GAMV uses nominal lag 5 m and tolerance 2.5 m. Reproduced curves use all axial grid pairs at 10 m increments. Comparisons use actual lag distances and pair counts; the display spans 0-155 m and normalized semivariance 0-1.2.

| Direction | Reproduced/input NRMSE | Reproduced/experimental NRMSE |
| --- | ---: | ---: |
| EW | 0.00873236 | 0.08626311 |
| NS | 0.00839186 | 0.06261752 |

These errors use the common input sill and overlapping distances. The [source record](../../pics/final_20261018/data/sgsim_reproduction_provenance.json) contains hashes, pair-estimator checks and calculation settings; per-realization lag curves are bundled with the plot data.

- [SGSIM reproduction](../../pics/06_sgsim_variogram_reproduction.png)

## Probability calibration

The displays follow [Wiegel et al. (2024)](https://doi.org/10.1002/gete.202400011), Figures 9 and 7a. Both submodels use five folds and 100 realizations per fold.

- **Binary SISIM:** 22,313 held-out observed Neogene clay/sand composites from 694 boreholes and 623 groups. Sand probability is the binary SISIM frequency over 100 realizations. Observed frequencies are displayed in 50 fixed 2% bins, (lower, upper], with p=0 in the first bin. Every bin is populated; the smallest has 179 samples. ECE is 0.04422624, unweighted bin MAE 0.05229931, Brier score 0.13174859 and classification accuracy 0.81172411; ties are classified as clay.
- **SGSIM:** 664 held-out boundary observations and their 100-member ensembles. The plot is the ECDF of mid-PIT quantile ranks. The saved nominal-quantile table has mean absolute calibration gap 0.01085923.

For predictions in the 0-10% region, mean sand probability is 0.0430 versus observed frequency 0.1427 (3441 samples). In the 40-80% region, the values are 0.5976 and 0.5269 (4333 samples).

The [validation summary](../article/VALIDATION.md#probability-calibration) describes holdout construction, weighting and preprocessing. The [protocol](protocol.json) records the model-selection and fold settings.

- [Binary SISIM calibration](../../pics/figure9_sisim_calibration.png)
- [SGSIM calibration](../../pics/figure9_sgsim_calibration.png)
- [Combined calibration panels](../../pics/figure9_final_calibration.png)

## Calculation records

The current model is delivered as `vtk/model_gok_sisim_final_20261018.vtk.gz`. Figure arrays, section sampling and model hashes are in [final_model_provenance.json](../../pics/final_20261018/data/final_model_provenance.json). Ensemble summaries pair SISIM fields with their SGSIM surfaces and use the recorded class, entropy and variability definitions.

The September example is archived in `configs/exponential_20260916.json`, `docs/article/recorded_run/`, `docs/article/source_audit.json`, the original `pics/data/` tables and `vtk/model_gok_exponential.vtk`. Its categorical model has one structure and nugget 0.025146355593272343. `python scripts/verify_article_evidence.py` checks that reference snapshot.

Implementation version checks are in [core_source_localization.json](../core_source_localization.json). See [REPRODUCE.md](REPRODUCE.md) for commands and [THIRD_PARTY.md](../../THIRD_PARTY.md) for source notices.
