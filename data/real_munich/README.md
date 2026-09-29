> Scope: this file records the earlier baseline/preparation profile. The article uses the later exponential configuration and terrain-masked postprocessing. See `docs/article/REPRODUCE.md` and `docs/article/VALIDATION.md` from the repository root. Statements below about no terrain mask or simulations not rerun refer to that earlier profile.

# Real Munich borehole example

This is a spatial subset of the user-supplied public Bavarian borehole archive. Its layer descriptions are real observations from `bohrungen_schichten.csv`; no lithology, depth, coordinate or elevation is simulated or repaired.

The selected 1 km square contains 774 borehole locations. A 500 m preparation buffer retains 2,265 shapefile/collar rows and 24,110 layer records. Records with missing collar information remain present for the existing QA and legacy exclusion rules. The shapefile and collar table retain the same paired source-row order because the legacy loader joins them by row index. `source_row_mapping.csv` records those original row positions.

`Schichtdaten_U9.xlsx` was used to select the location, not to overwrite CSV values. Its Survey sheet contains 5,080 records and its layer sheet contains 52,158 records. Matching IDs support interpretation of its coordinates as EPSG:31468: transformation to EPSG:25832 differs from the same-ID shapefile coordinates by a median of 0.282 m and a 95th percentile of 0.295 m. Modelling uses the supplied EPSG:25832 shapefile coordinates directly.

Selection ranks 1 km candidate squares on a 100 m lattice by the number of occupied 100 m cells, then capped per-cell support, then eligible hole count. Siting eligibility uses Excel-linked holes with collar elevation 400–700 m and terminal depth at least 30 m. These are selection criteria only; they do not filter the packaged modelling input. The selected square has eligible observations in 54 of its 100 cells and 249 eligible holes.

See `docs/real_example/domain_selection.json` for source SHA-256 hashes, original-row issues, coordinate bounds and the complete criteria; the companion CSV and PNG show candidate rankings and the selected footprint. The source archive was supplied at `GRID_D2.2/Example Results/WD/data/bohrungen_epsg25832_shp/`. The accompanying original quick guide is retained here for provenance. No new license is assigned to the source data.

The main grid is unchanged in size: 100 × 100 × 200 cells, spacing 10 × 10 × 0.5 m. In EPSG:25832 its outer boundaries are X 690400–691400 m, Y 5334500–5335500 m, and elevation 475–575 m. First cell centres are (690405, 5334505, 475.25). The example explicitly sets angle = 0; no horizontal rotation is applied. Ground masking is still not applied.

The exact legacy preparation bounds are X 689905–691905 m and Y 5334005–5336005 m: the old function measures the buffer from the first cell centre. This 5 m offset is preserved.
