# Current article validation - 18 October 2026

The current model uses `response_step_50`, zero nugget and four exponential structures. It has 100 actually generated selected-input spatial realizations and 500 new SISIM cross-validation realizations: five saved whole-borehole/group folds, 100 members per fold. Continuous SGSIM validation reuses its original saved predictions.

| Assessment | Current result |
| --- | --- |
| Spatial comparison | Fixed experimental indicator variograms versus the actual 100-member reproduced ensemble, EW/NS/Z |
| Common-reference normalized RMSE | EW 0.13585777, NS 0.06213126, Z 0.07241337 |
| Binary SISIM calibration sample | 22,313 observed Neogene composites; 694 boreholes; 623 groups |
| Binary Brier score | 0.13174859 |
| Binary 2% bin ECE | 0.04422624 |
| Binary classification accuracy | 0.81172411; clay tie rule retained |
| Continuous SGSIM calibration | 664 saved held-out boundary observations; five folds, 100 members per fold |
| Saved nominal-quantile calibration MAE | 0.01085923 |
| Native CV completion evidence | 100 completed batches x 5 realizations, with hashes in `../final_20261018/native_batches_verified.csv` |
| Holdout audit | Boreholes/groups exclusive; no shared training/test voxels in any fold |

Use [current spatial panels](../../pics/08_variogram_neogene_common_variance.png) and [probability calibration panels](../../pics/figure9_final_calibration.png). The two validation axes follow the spatial-continuity and probability-calibration approach discussed by [Wiegel et al. (2024)](https://doi.org/10.1002/gete.202400011); no universal numerical acceptance threshold is inferred from that paper.

Spatial calibration adjusts the input and actually reruns SISIM; it is not a rescaling of a previous reproduced curve. The selected objective improved modestly and residuals remain. In binary predictions below 10%, mean predicted P(Sand)=0.0430 versus observed frequency 0.1427; at 40-80%, 0.5976 versus 0.5269. The staircase is unsmoothed, with all 50 bins populated. These local biases must not be hidden by line thickness.

The current CV holds out the saved complete borehole/group folds but fixes variograms selected using the full data. It is conditional validation, not nested variogram selection. Observed composites are the weighting unit. Production uses the preserved Johnson preprocessing whereas CV retains the original training-only Quantile preprocessing. Legacy entropy is unnormalized natural-log entropy with its original zero replacement; it can exceed 1. The adopted geological interpretation, search settings, grid/support and finite ensemble remain additional limitations.

Run `python scripts/verify_final_validation.py` for consistency of current packaged predictions, calibration bins, spatial arrays and figure/VTK hashes. That is not a fresh simulation or independent validation. The historical `source_audit.json`, `recorded_run/`, original `pics/data/` and original VTK retain the September case; `verify_article_evidence.py` checks those old records only. They must not be quoted as current-model statistics.

See the [full calculation record](../final_20261018/README.md), [frozen protocol](../final_20261018/protocol.json), [holdout audit](../final_20261018/holdout_audit.json) and [reproduction guide](REPRODUCE.md).
