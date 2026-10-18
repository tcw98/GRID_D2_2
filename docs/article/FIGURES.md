# Current article figures - 18 October 2026

The current categorical model is `configs/exponential_sisim_final_20261018.json`. Figures affected by its revised input and simulations have been replaced at the root `pics/` paths. Matching PDF versions are included. Reproduction uses the independent [current plotting commands](REPRODUCE.md); simulation kernels are unchanged.

| File under `pics/` (PNG and PDF) | Content |
| --- | --- |
| `07_variogram_indicator_fit` | Experimental and four-structure zero-nugget input, each direction/series peak-scaled to 1 |
| `08_variogram_reproduction` | Complete binary SISIM grid, 100 realizations; each direction/series peak-scaled to 1; thick reproduced mean and one envelope |
| `10_sisim_sections` | New raw SISIM realizations 1, 50 and 100 |
| `11_hierarchical_combination` | Unchanged SGSIM boundary combined with new SISIM |
| `12_uq_plan` | New class probabilities, mode, entropy and variability on a plan |
| `13_sections_entropy` | New mode, entropy and entropy tiers on two sections |
| `14_sections_variability` | New mode, variability and variability tiers on two sections |
| `16_realization_stability` | New ensemble-prefix mean entropy on the fixed terrain mask |
| `figure9_sisim_calibration` | Binary SISIM, 50 probability bins of width 2%, staircase plus bisector |
| `figure9_sgsim_calibration` | Reused SGSIM held-out quantile ranks, ECDF plus bisector |
| `figure9_final_calibration` | Combined final calibration panels |

The unchanged data/preparation, continuous variography and SGSIM illustrations (`01`-`06`, `09`) remain relevant. The original `18_*` accuracy figures and original `pics/data/` confusion matrices are September historical evidence and are not used as the new binary calibration result.

Each spatial panel contains one experimental point set, one theoretical line, one reproduced mean and one 5th-95th percentile envelope. Experimental zero/near-zero records remain in the source tables but are excluded from regular positive lags. The drawn reproduced origin is the identity gamma(0)=0. Distances are in metres and the normalized y-axis is 0-1.2. Only one 08 is delivered. For each direction, experimental points, theoretical input and reproduced mean are divided by their own maxima over the displayed lag window. The 07 theoretical and experimental curves use exactly the same factors as 08. This removes amplitude differences and is a shape display, not a common-variance validation. Raw full-grid RMSE and all scaling factors are retained in `pics/final_20261018/data/variogram_peak_normalization.json`.

`pics/current_article.json` lists the current figure hashes and source directory. `pics/final_20261018/data/` provides the compact arrays and tables used to redraw these figures. The old `pics/provenance.json` is explicitly historical; the current model record is `pics/final_20261018/data/final_model_provenance.json`.
