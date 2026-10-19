# Current article figures - 18 October 2026

The current categorical model is `configs/exponential_sisim_final_20261018.json`. Figures affected by its revised input and simulations have been replaced at the root `pics/` paths. Matching PDF versions are included. Redraw commands are in [REPRODUCE.md](REPRODUCE.md).

| File under `pics/` (PNG and PDF) | Content |
| --- | --- |
| `06_sgsim_variogram_reproduction` | SGSIM Gaussian-residual experimental/input/mean reproduced plus 5%-95% envelope, EW/NS, common input sill and non-zero nugget |
| `07_variogram_indicator_fit` | Experimental and four-structure zero-nugget input, each direction/series peak-scaled to 1 |
| `08_variogram_reproduction` | Complete binary SISIM grid, 100 realizations; peak-scaled mean and 5th-95th percentile envelope |
| `10_sisim_sections` | New raw SISIM realizations 1, 50 and 100 |
| `11_hierarchical_combination` | SGSIM boundary combined with SISIM |
| `12_uq_plan` | New class probabilities, mode, entropy and variability on a plan |
| `13_sections_entropy` | New mode, entropy and entropy tiers on two sections |
| `14_sections_variability` | New mode, variability and variability tiers on two sections |
| `16_realization_stability` | New ensemble-prefix mean entropy on the fixed terrain mask |
| `figure9_sisim_calibration` | Binary SISIM, 50 probability bins of width 2%, staircase plus bisector |
| `figure9_sgsim_calibration` | SGSIM held-out quantile ranks, ECDF plus bisector |
| `figure9_final_calibration` | Combined final calibration panels |

Figures `01`-`06` and `09` cover data preparation, continuous variography and SGSIM realizations. The `18_*` accuracy figures belong to the September reference case; its records are listed in the [calculation record](../final_20261018/README.md#calculation-records).

SISIM Figure 08 shows experimental points, theoretical input, reproduced mean and a 5th-95th percentile envelope. Figures 07-08 use the same direction-wise peak factors. Distances are in metres and normalized semivariance spans 0-1.2. [Spatial definitions and scaling factors](../final_20261018/README.md#spatial-continuity) accompany the numerical results.

`pics/current_article.json` lists the current figure hashes and source directory. `pics/final_20261018/data/` provides the compact arrays and tables used to redraw these figures. The current model record is `pics/final_20261018/data/final_model_provenance.json`.

## SGSIM spatial scope and scale

The SGSIM panel uses 100 full Gaussian residual fields in EW/NS and a common input sill, retaining its non-zero nugget. GAMV 5 m bins (2.5 m tolerance) and reproduced 10 m axial lags are displayed at their actual distances. See the [full spatial record](../final_20261018/README.md#continuous-sgsim-reproduction).
