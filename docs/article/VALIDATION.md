# Munich reference-example evidence — 29 September 2026

The available local production run was inspected and its saved numerical results rechecked. No new production or CV simulations were run for this review. Evidence from earlier synthetic and rotated-grid tests is retained as history and is not substituted for the article case.

| Check | Result | Evidence |
| --- | --- | --- |
| Numerical/reference source equivalence | Original numerical tokens retained; comment translations documented | `../core_source_localization.json` |
| Production config | Equal to checked-in exponential config | `source_audit.json` |
| Saved production ensemble | 100 × 2,000,000 | `source_audit.json` |
| Terrain mask | 857,561 below-ground cells | `../../data/article/below_gok.npy` |
| Recomputed legacy entropy vs saved output below ground | Maximum absolute difference 0 | `source_audit.json` |
| Entropy greater than 1 | 3,061 below-ground cells | `source_audit.json` |
| Maximum entropy | 1.0985126170507196 | `source_audit.json` |
| CV confusion matrices | 24,563 correct / 29,961 composites; 81.983244885% | `../../pics/data/fold_*_confusion.csv` |
| VTK hash | Matches original figure provenance | `source_audit.json` |

## Ensemble-prefix diagnostic

| Members | Mean legacy entropy | Fixed evaluation cells |
| ---: | ---: | ---: |
| 20 | 0.425537315621 | 857561 |
| 40 | 0.439058097830 | 857561 |
| 60 | 0.444032462378 | 857561 |
| 80 | 0.446831507415 | 857561 |
| 100 | 0.449302112721 | 857561 |

The last 40 members increase the mean by 0.005269650343, or about 1.17% of the final value. A flattening global mean does not demonstrate voxel-wise convergence, independent-seed stability or calibration. The ensemble size is prescribed by the implementation.

## Interpretation boundaries

- D2.1 presents normalized Shannon entropy; the preserved D2.2 implementation reports unnormalized natural-log entropy with its legacy zero-probability replacement. The entropy classification must use `H > 0.7` for the high class, including values above 1.
- Article CV holds out groups of complete boreholes and grid columns but keeps full-data variograms fixed. It is conditional CV, not nested refitting or a completely independent model-selection assessment.
- The score is weighted by observed composites, not independent voxels or equal-weight boreholes. Class curves are one-vs-rest accuracies, not recalls. Fold metrics also contain recall, precision and borehole-macro accuracy.
- CV preprocessing uses the original `save_train` Quantile transformation, whereas production uses Johnson. This difference remains unchanged and must be reported.
- The saved exponential run has stages 01–04 and separate hierarchical CV. It has no completed stage-05 hard-data/variogram-reproduction evidence. Full claims that all original observations and spatial structure are reproduced are unsupported by this run.
- The uncertainty output is conditional on fixed interpretation, class grouping and covariance assumptions. It is not an engineering failure probability or a quantified risk assessment.


English release checks are recorded in `../RELEASE_VERIFICATION.md`. The evidence above refers to the saved numerical example, not to a new simulation run.
