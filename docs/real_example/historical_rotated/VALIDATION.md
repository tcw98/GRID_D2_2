# Real-data verification

The relocated 1 km grid retains 100 × 100 × 200 cells with 10 × 10 × 0.5 m spacing. It uses the real CSV/shapefile subset selected with the Excel survey. The buffer contains 2,243 original location/collar rows and 23,804 layer records; the main square contains 883 borehole positions.

Stage 01 passed all 14 direct-original-function comparisons. There are 2,149 selected survey rows, 2,149 boundary records (including their missing-value cases), 97,581 lower-unit composites and 1,461 retained detrended conditioning points. The separate full-archive comparison is exact for qbasis, tertiary_points, coords_d, coords_t, grid_data and points_rotated. The stored qbasis array has four columns including borehole ID; its numeric coordinate/elevation component has three.

The smaller 24 × 24 × 40 real-data profile was run through stages 01–05 using actual GSLIB90 executables and 100 realizations, including two-fold submodel cross-validation. Eleven direct-original-function comparisons and all eight basic output checks passed. The ten Python tests also passed. The optional 25-case parameter sweep was not repeated. The full 1 km stages 02–05 have not been executed.

## Scientific diagnostics remain unresolved

Hard-data reproduction did not pass completely. The continuous field has 85 in-grid observations and 23 conflicting points; the indicator field has 2,141 in-grid observations and 36 conflicting points; the combined hierarchy has 2,641 in-grid observations and 79 conflicting points. Each check also reports unambiguous observations that were not reproduced. No observations or algorithms were altered to hide these failures.

The continuous legacy ensemble MAE is 0.9798037065 m over 85 evaluated observations. Binary categorical accuracy is 0.8589444185 over 2,141 observations. Although the two folds share no XY locations, each contains 133 training/test voxel overlaps after grid assignment. These are legacy diagnostic scores, not independent voxel-level validation or joint three-class hierarchical validation. Fold preprocessing retains the original Quantile transform and fixed variograms.

Mean entropy changes are 0.053757 for 20→50 realizations and 0.031344 for 50→100. The larger ensemble is more stable by this diagnostic, but neither comparison meets an illustrative 0.01 entropy-change criterion. Entropy remains the original, unnormalised natural-log function.

The JSON/CSV files in this directory contain the detailed evidence and the new 1 km quality report. Historical synthetic evidence remains in ../validation and is explicitly distinguished from this real-data run.
