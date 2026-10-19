> Historical preparation profile - September 2026. [Current article validation](../article/VALIDATION.md).

# Historical unrotated preparation - September 2026

The delivered example uses EPSG:25832 coordinates directly, with preparation.angle = 0 in both example_1km.json and smoke.json. No rotation or inverse rotation is required. The unchanged legacy core still receives its existing rotation argument, set to zero.

The main domain is X 690400–691400 m, Y 5334500–5335500 m, elevation 475–575 m. It has 100 × 100 × 200 cells at 10 × 10 × 0.5 m spacing. It contains 774 real borehole positions. The exact legacy preparation extent retains 2,265 location/collar rows and 24,110 original layer records.

All 14 stage-01 direct-original-function comparisons passed for this configuration. Counts: qbasis 2179 records, tertiary_points 97818 records, coords_d 1363 points. QA remains flag/report only; existing function exclusions remain unchanged. Ten Python regression tests passed.

This record covers stage 01 for the unrotated configuration. Earlier simulation diagnostics for the 25.67-degree configuration are in historical_rotated/.

This preparation profile uses the existing covariance parameters, seeds, realization schedule and entropy definition.
