# Munich reference-example figures

These assets illustrate the supplied reference case. The general modelling workflow is described in the repository-root README and the D2.1 reference. Example-specific parameters and plotting commands are given in [REPRODUCE.md](REPRODUCE.md).

| File under `pics/` | Content |
| --- | --- |
| `01_grid_boreholes.png` | Modelling grid, borehole locations and attributed map background |
| `02_length_weighted_compositing.png` | Compositing example |
| `03_statistical_preparation.png` | Boundary statistical preparation |
| `04_varmap_boundary.png` | Continuous boundary variogram map |
| `05_varmap_indicator.png` | Indicator variogram maps |
| `06_variogram_boundary_fit.png` | Directional exponential fits for the boundary |
| `07_variogram_indicator_fit.png` | Horizontal and vertical indicator fits |
| `09_sgsim_realizations.png` | Selected continuous simulations |
| `13_sections_entropy.png` | Most likely class, entropy and entropy tiers on two sections |
| `14_sections_variability.png` | Most likely class, variability and variability tiers on two sections |
| `16_realization_stability.png` | Mean entropy over a fixed below-ground mask for ensemble prefixes |
| `18_five_fold_accuracy_table.png` | Overall and one-vs-rest accuracy by fold |
| `18_indicator_prediction_accuracy.png` | One-vs-rest class accuracies |
| `19_five_fold_plan.png` | Spatial fold assignments |

Numerical values, grouping definitions and plotting provenance are retained in `pics/data/` and `pics/provenance.json`. Entropy is unnormalized. The sections requested at X=691000 and Y=5334800 m use the positive-side cells centred at X=691005 and Y=5334805 m. Plotting uses the supplied terrain mask. Directional fitted curves are distinct from the consolidated simulation covariance parameters.
