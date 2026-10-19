# Article validation

The categorical model uses `response_step_50`, with zero nugget and four exponential structures. Validation evaluates spatial continuity and predictive probability calibration using 100 spatial realizations and grouped five-fold SISIM CV with 100 realizations per fold. SGSIM results use the existing continuous ensembles.

| Assessment | Result |
| --- | --- |
| SISIM spatial comparison | Experimental indicator variograms versus 100 full-grid binary realization variograms, EW/NS/Z |
| Raw full-grid common-reference normalized RMSE | EW 0.19886292, NS 0.13333520, Z 0.09661573 |
| Binary SISIM calibration sample | 22,313 observed Neogene composites; 694 boreholes; 623 groups |
| Binary Brier score | 0.13174859 |
| Binary 2% bin ECE | 0.04422624 |
| Binary classification accuracy | 0.81172411; clay tie rule |
| SGSIM spatial reproduction | 100 Gaussian residual fields; common-sill mean/input NRMSE EW 0.00873236, NS 0.00839186 |
| Continuous SGSIM calibration | 664 held-out boundary observations; five folds, 100 members per fold |
| Saved nominal-quantile calibration MAE | 0.01085923 |
| Native CV batches | 100 batches x 5 realizations; [batch records](../final_20261018/native_batches_verified.csv) |
| Holdout grouping | Complete boreholes/groups; zero shared training/test voxels in each fold |

## Spatial continuity

The SISIM directional ranges were calibrated by adjusting the input variogram and rerunning simulations. Figures 07-08 compare directional shapes after each series is divided by its maximum over the plotted lag window. Raw amplitude errors are reported separately in the table above. The displayed ensemble envelope uses the 5th-95th percentiles.

SGSIM reproduction is evaluated in Gaussian residual space for the two-dimensional boundary. Experimental, theoretical and reproduced curves share the total input sill 0.9003331304926536, including nugget 0.13695550591961783. The reproduced mean closely follows the input in EW/NS, with experimental scatter visible in the plot.

- [SGSIM spatial reproduction](../../pics/06_sgsim_variogram_reproduction.png)
- [SISIM spatial reproduction](../../pics/08_variogram_reproduction.png)
- [Spatial calculation details and scaling factors](../final_20261018/README.md#spatial-continuity)

## Probability calibration

The displays follow the spatial-continuity and probability-calibration approach discussed by [Wiegel et al. (2024)](https://doi.org/10.1002/gete.202400011). Binary SISIM uses observed sand frequencies in 50 fixed 2% probability bins; SGSIM uses the empirical distribution of held-out quantile ranks.

Regional SISIM calibration differences are visible: below 10% predicted probability, mean P(Sand)=0.0430 versus observed frequency 0.1427; in the 40-80% region, the corresponding values are 0.5976 and 0.5269. All 50 bins are populated.

Five-fold CV holds out complete borehole groups while keeping variogram parameters estimated from the full dataset fixed. Scores are weighted by observed composites. Production uses the Johnson transformation; CV uses a Quantile transformation fitted to each training fold.

- [Calibration panels](../../pics/figure9_final_calibration.png)
- [Calibration definitions and sample counts](../final_20261018/README.md#probability-calibration)
- [Protocol](../final_20261018/protocol.json) and [holdout audit](../final_20261018/holdout_audit.json)

## Calculation checks

`python scripts/verify_final_validation.py` checks the bundled tables, arrays, figures and provenance. Uncertainty summaries use natural-log entropy with the original zero-probability replacement; entropy values can exceed 1. See the [reproduction guide](REPRODUCE.md) and [calculation records](../final_20261018/README.md#calculation-records).
