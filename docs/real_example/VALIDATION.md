> Scope: this file records the earlier baseline/preparation profile. The article uses the later exponential configuration and terrain-masked postprocessing. See `docs/article/REPRODUCE.md` and `docs/article/VALIDATION.md` from the repository root. Statements below about no terrain mask or simulations not rerun refer to that earlier profile.

# Unrotated real-data example

The delivered example uses EPSG:25832 coordinates directly, with preparation.angle = 0 in both example_1km.json and smoke.json. No rotation or inverse rotation is required. The unchanged legacy core still receives its existing rotation argument, set to zero.

The main domain is X 690400–691400 m, Y 5334500–5335500 m, elevation 475–575 m. It has 100 × 100 × 200 cells at 10 × 10 × 0.5 m spacing. It contains 774 real borehole positions. The exact legacy preparation extent retains 2,265 location/collar rows and 24,110 original layer records.

All 14 stage-01 direct-original-function comparisons passed for this configuration. Counts: qbasis 2179 records, tertiary_points 97818 records, coords_d 1363 points. QA remains flag/report only; existing function exclusions remain unchanged. Ten Python regression tests passed.

Stages 02–05 have not been rerun for this unrotated location. Results in historical_rotated/ belong to the previous 25.67-degree configuration, including the full-source/subset comparison and small-grid simulations. They are not verification of the current grid. Changing the location and orientation changes the conditioning input, so output equality with that previous configuration is not expected.

Core functions, seeds, realization schedule and the original entropy implementation remain unchanged. Variogram parameters are retained example parameters and have not been refitted to the new location or orientation.
