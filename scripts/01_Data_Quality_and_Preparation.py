"""GRID D2.2 stage 01: audit -> Munich compositing -> statistical preparation.

Run beside the unmodified UTILITIES_Wang.py and GSLIB_parfl_011.py.
Default: original U9 grid, rotation, buffer, compositing and statistical settings.
--example: the real-data 1 km square in EPSG:25832, without rotation.
QC is observational: no QC flag changes inputs. Invalid inputs stop preparation.
Legacy LFU/detrend exclusions are preserved and exported separately from QC.
--verify compares against direct legacy function calls, NOT a historical run.
Never import Workflow_U9_Wang.py: it contains interactive experimental branches.
"""
from __future__ import annotations

import argparse
import ast
import contextlib
import hashlib
import importlib.metadata
import json
import io
import pickle
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import affinity


REPORT_COLUMNS = ["stage", "source_file", "source_row", "borehole_id", "field",
                  "issue_code", "severity", "observed_value", "reference_value",
                  "residual", "threshold", "neighbour_count", "action", "review_status"]


def defaults():
    # Order is part of the legacy interface: functions unpack grid.values().
    return {
        "grid": dict(nx=123, ny=692, nz=210, xsiz=10, ysiz=10, zsiz=0.5,
                     xmn=688440.0, ymn=5332620.0, zmn=475.25),
        "origin": [688000, 5332000], "angle": 25.67, "buffer": 500,
        "crs": "EPSG:25832", "interval": 0.5, "qbasistol": 3,
        "radius_xy": [300, 300], "tolerance": 20, "outlier_tolerance": 5,
        "trans_method": "Johnson", "vtk": True,
        # Review thresholds only; they do not trim observations.
        "collar_range": [400, 700], "collar_residual_limit": 10,
        "qc_radius": 300, "qc_min_neighbours": 3, "qb_residual_limit": 5,
        "interval_tolerance": 1e-6, "end_depth_review_limit": 500,
    }


class Report:
    """A small report accumulator; no model classes or algorithm abstractions."""
    def __init__(self):
        self.rows = []

    def add(self, stage, source, row, bh, field, code, value=None,
            reference=None, residual=None, threshold=None, neighbours=None,
            severity="warning", action="flag_for_review"):
        self.rows.append(dict(zip(REPORT_COLUMNS, [stage, str(source), row, bh,
            field, code, severity, value, reference, residual, threshold,
            neighbours, action, "pending"])))

    def save(self, path):
        pd.DataFrame(self.rows, columns=REPORT_COLUMNS).to_csv(path, index=False,
                                                               encoding="utf-8-sig")


def local_plausibility(frame, value_col, radius, min_neighbours, limit,
                       report, source, stage):
    """Independent QC at observations, leave-one-borehole-out, never filtering.

    Arithmetic moving average of OTHER borehole IDs within radius. Unsupported
    observations have NaN trend/residual, not a fabricated reference value.
    Coordinates are metres in the same planar frame. Not a modelling detrend.
    """
    result = frame.copy()
    result["local_reference"] = np.nan
    result["local_residual"] = np.nan
    result["neighbour_count"] = 0
    values = pd.to_numeric(frame[value_col], errors="coerce").to_numpy(float)
    xy = frame[["xcoord", "ycoord"]].to_numpy(float)
    valid = np.isfinite(values) & np.isfinite(xy).all(axis=1)
    positions = np.flatnonzero(valid)
    if len(positions):
        tree = cKDTree(xy[valid])
        ids = frame["borehole_id"].astype(str).to_numpy()
        for pos, indices in zip(positions, tree.query_ball_point(xy[valid], radius)):
            neighbours = positions[np.asarray(indices, dtype=int)]
            neighbours = neighbours[ids[neighbours] != ids[pos]]
            # Duplicate rows for one borehole must not multiply its QC weight.
            by_id = pd.Series(values[neighbours], index=ids[neighbours]).groupby(level=0).mean()
            n = len(by_id)
            result.iloc[pos, result.columns.get_loc("neighbour_count")] = n
            code = "insufficient_local_support"
            ref = residual = None
            if n >= min_neighbours:
                ref = float(by_id.mean())
                residual = float(values[pos] - ref)
                result.iloc[pos, result.columns.get_loc("local_reference")] = ref
                result.iloc[pos, result.columns.get_loc("local_residual")] = residual
                code = "spatial_residual_exceeds_threshold" if abs(residual) > limit else None
            if code:
                report.add(stage, source, frame.iloc[pos].get("source_row", pos),
                           ids[pos], value_col, code, values[pos], ref, residual,
                           limit, n)
    return result


def audit_raw(wd, cfg, u, report, out):
    """Read raw collar elevation BEFORE calling LFU_2_YOU (which drops NaN/0)."""
    folder = wd / "data" / "bohrungen_epsg25832_shp"
    shp = folder / "bohrungen.shp"
    collar_path = folder / "bohrungen_stammdaten.csv"
    layers_path = folder / "bohrungen_schichten.csv"
    geo = gpd.read_file(shp)
    collars = pd.read_csv(collar_path, sep=";", decimal=",")
    layers = pd.read_csv(layers_path, sep=";", decimal=",")
    required = [(geo, ["ObjektID", "geometry"], shp),
                (collars, ["ObjektID", "Ansatzhoeh", "Endteufe"], collar_path),
                (layers, ["ObjektID", "Obergrenz", "Untergrenz", "DIN", "PetBez"], layers_path)]
    for table, columns, source in required:
        for column in columns:
            if column not in table:
                report.add("raw", source, None, None, column, "missing_required_column",
                           severity="error", action="stop_preparation")
    if any(r["issue_code"] == "missing_required_column" for r in report.rows):
        raise ValueError("Missing required raw columns; see data_quality_report.csv")
    if geo.crs is None or geo.crs != cfg["crs"]:
        report.add("raw", shp, None, None, "CRS", "crs_mismatch_or_missing",
                   str(geo.crs), cfg["crs"], severity="error", action="stop_preparation")
        raise ValueError("CRS mismatch: no automatic reprojection in compatibility mode")
    # Global ID / coordinate audit, domain-limited geological checks below.
    for table, _, source in required:
        missing = table.ObjektID.isna() | table.ObjektID.astype(str).str.strip().eq("")
        duplicate = table.duplicated(keep=False)
        if source != layers_path:
            duplicate_id = table.ObjektID.duplicated(keep=False) & ~missing
            for idx in table.index[duplicate_id]:
                report.add("raw", source, int(idx) + 2, table.at[idx, "ObjektID"],
                           "ObjektID", "duplicate_borehole_id")
        for mask, code in [(missing, "missing_borehole_id"), (duplicate, "duplicate_record")]:
            for idx in table.index[mask]:
                report.add("raw", source, int(idx) + 2, table.at[idx, "ObjektID"], "record", code)
    usable = geo.geometry.notna() & ~geo.geometry.is_empty & geo.geometry.geom_type.eq("Point")
    for idx in geo.index[~usable]:
        report.add("raw", shp, int(idx), geo.at[idx, "ObjektID"], "geometry", "missing_or_nonpoint_geometry")
    xy = pd.DataFrame(index=geo.index, columns=["xcoord", "ycoord"], dtype=float)
    rotated = geo.loc[usable].geometry.apply(lambda p: affinity.rotate(p, cfg["angle"], origin=cfg["origin"]))
    xy.loc[usable, "xcoord"] = rotated.x
    xy.loc[usable, "ycoord"] = rotated.y
    for idx in xy.index[usable & ~np.isfinite(xy).all(axis=1)]:
        report.add("raw", shp, int(idx), geo.at[idx, "ObjektID"], "coordinates", "nonfinite_coordinates")
    extent = u.domain_plus_buffer(cfg["grid"], cfg["buffer"])
    xmin, xmax, ymin, ymax = extent
    selected = xy.xcoord.gt(xmin) & xy.xcoord.lt(xmax) & xy.ycoord.gt(ymin) & xy.ycoord.lt(ymax)
    ids = set(geo.loc[selected, "ObjektID"])
    # Audit precisely the legacy INDEX join, not a replacement ID merge.
    joined = geo.loc[selected, ["ObjektID"]].join(collars, lsuffix="ObjektID", rsuffix="2", how="inner")
    bad_join = joined.ObjektIDObjektID.ne(joined.ObjektID2)
    legacy_excluded = joined.Ansatzhoeh.isna() | joined.Ansatzhoeh.eq(0)
    for idx in joined.index[bad_join]:
        report.add("raw", collar_path, int(idx) + 2, joined.at[idx, "ObjektIDObjektID"],
                   "ObjektID", "legacy_index_join_id_mismatch", joined.at[idx, "ObjektID2"],
                   severity="warning" if legacy_excluded.loc[idx] else "error",
                   action="legacy_excluded" if legacy_excluded.loc[idx] else "stop_preparation")
    for bh in ids - set(collars.ObjektID):
        report.add("raw", collar_path, None, bh, "ObjektID", "missing_collar_record")
    for bh in ids - set(layers.ObjektID):
        report.add("raw", layers_path, None, bh, "ObjektID", "missing_layer_records")
    orphan = ~layers.ObjektID.isin(set(collars.ObjektID))
    for bh, count in layers.loc[orphan].groupby("ObjektID", dropna=False).size().items():
        report.add("raw", layers_path, None, bh, "ObjektID", "orphan_layer_records", int(count))
    raw = joined.copy()
    raw["borehole_id"] = raw.ObjektIDObjektID
    raw["source_row"] = raw.index + 2
    raw[["xcoord", "ycoord"]] = xy.loc[raw.index]
    raw["original_x"] = geo.loc[raw.index].geometry.x
    raw["original_y"] = geo.loc[raw.index].geometry.y
    elev = pd.to_numeric(raw.Ansatzhoeh, errors="coerce")
    raw["legacy_collar_excluded"] = raw.Ansatzhoeh.isna() | raw.Ansatzhoeh.eq(0)
    for idx in raw.index:
        val = elev.loc[idx]
        code = None
        if pd.isna(val): code = "missing_or_nonnumeric_collar_elevation"
        elif not np.isfinite(val): code = "nonfinite_collar_elevation"
        elif val == 0: code = "zero_collar_elevation"
        elif not cfg["collar_range"][0] <= val <= cfg["collar_range"][1]: code = "collar_elevation_outside_review_range"
        if code:
            action = "legacy_excluded" if raw.at[idx, "legacy_collar_excluded"] else "flag_for_review"
            report.add("raw", collar_path, int(idx) + 2, raw.at[idx, "borehole_id"],
                       "Ansatzhoeh", code, raw.at[idx, "Ansatzhoeh"], threshold=str(cfg["collar_range"]), action=action)
    # NaN/0 legacy-excluded collars remain in the audit, but cannot define a QC reference.
    spatial_input = raw.copy()
    spatial_input.loc[raw.legacy_collar_excluded, "Ansatzhoeh"] = np.nan
    spatial = local_plausibility(spatial_input, "Ansatzhoeh", cfg["qc_radius"],
                                cfg["qc_min_neighbours"], cfg["collar_residual_limit"],
                                report, collar_path, "raw")
    raw[spatial.columns.difference(raw.columns)] = spatial[spatial.columns.difference(raw.columns)]
    grid = cfg["grid"]
    raw["in_modelling_domain"] = (raw.xcoord.ge(grid["xmn"] - grid["xsiz"] / 2) &
        raw.xcoord.le(grid["xmn"] + (grid["nx"] - .5) * grid["xsiz"]) &
        raw.ycoord.ge(grid["ymn"] - grid["ysiz"] / 2) &
        raw.ycoord.le(grid["ymn"] + (grid["ny"] - .5) * grid["ysiz"]))
    for idx in raw.index[~raw.in_modelling_domain]:
        report.add("raw", shp, int(idx), raw.at[idx, "borehole_id"], "coordinates",
                   "outside_model_inside_legacy_buffer", severity="info", action="retain_for_conditioning")
    raw.to_csv(out / "raw_collar_audit.csv", index=False)
    relevant = layers.loc[layers.ObjektID.isin(ids)].copy()
    relevant["source_row"] = relevant.index + 2
    relevant.to_csv(out / "raw_layers_audited.csv", index=False)
    audit_intervals(relevant, raw, cfg, report, layers_path)
    if any(r["action"] == "stop_preparation" for r in report.rows):
        raise ValueError("Raw data linkage/preconditions failed; inputs have not been repaired")
    return (shp, collar_path, layers_path), extent, raw


def audit_intervals(layers, collars, cfg, report, source):
    tol = cfg["interval_tolerance"]
    end_map = collars.set_index("borehole_id").Endteufe
    for bh, group in layers.groupby("ObjektID", sort=False):
        top = pd.to_numeric(group.Obergrenz, errors="coerce")
        bottom = pd.to_numeric(group.Untergrenz, errors="coerce")
        numeric = np.isfinite(top) & np.isfinite(bottom)
        for idx in group.index[~numeric]:
            report.add("raw", source, int(idx) + 2, bh, "interval", "nonnumeric_or_missing_depth",
                       str((group.at[idx, "Obergrenz"], group.at[idx, "Untergrenz"])))
        for idx in group.index[numeric]:
            thickness = bottom.loc[idx] - top.loc[idx]
            if thickness <= tol or top.loc[idx] < 0:
                report.add("raw", source, int(idx) + 2, bh, "interval", "negative_or_zero_thickness_or_depth",
                           thickness, threshold=tol)
        if not bottom[numeric].is_monotonic_increasing:
            report.add("raw", source, None, bh, "Untergrenz", "original_interval_order", action="legacy_sorts_by_lower_depth")
        # Sorting is for the audit only; the input table is never changed.
        ordered = pd.DataFrame({"top": top[numeric], "bottom": bottom[numeric]}).sort_values("top", kind="stable")
        previous_bottom = 0.0
        for idx, interval in ordered.iterrows():
            delta = interval.top - previous_bottom
            if abs(delta) > tol:
                report.add("raw", source, int(idx) + 2, bh, "interval", "gap" if delta > 0 else "overlap",
                           interval.top, previous_bottom, delta, tol)
            previous_bottom = max(previous_bottom, interval.bottom)
        end = pd.to_numeric(end_map.get(bh, np.nan), errors="coerce")
        if not np.isscalar(end):
            report.add("raw", source, None, bh, "Endteufe", "ambiguous_duplicate_collar_end_depth")
        elif not np.isfinite(end) or end <= 0 or end > cfg["end_depth_review_limit"]:
            report.add("raw", source, None, bh, "Endteufe", "missing_or_implausible_end_depth", end)
        elif len(ordered) and abs(ordered.bottom.max() - end) > tol:
            report.add("raw", source, None, bh, "Endteufe", "end_depth_mismatch", end, ordered.bottom.max(), threshold=tol)
        for idx, row in group.iterrows():
            din = row.DIN
            if pd.isna(din) or not str(din).strip():
                code = "missing_DIN"
            else:
                main = str(din).replace("/", "*").split(",")[0]
                matches = [c for c in "GASUTFZ" if c in main]
                code = "unknown_DIN_legacy_default_10" if not matches else "ambiguous_DIN_legacy_precedence" if len(matches) > 1 else None
            if code:
                report.add("raw", source, int(idx) + 2, bh, "DIN", code, din)
            if pd.isna(row.PetBez) or not str(row.PetBez).strip():
                report.add("raw", source, int(idx) + 2, bh, "PetBez", "missing_lithological_description")


def geological_preprocessing(survey, layers, cfg, u, report=None):
    # CASE SPECIFIC: preserve all Munich DIN grouping and base identification.
    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured):
            comp = u.qnd_compositing(survey, layers, compositing_interval=cfg["interval"],
                                    vtkoption=cfg["vtk"], verticaloption=True,
                                    qbasisoption=True, qbasistol=cfg["qbasistol"])
    finally:
        print(captured.getvalue(), end="")
    # Legacy catches some per-hole exceptions internally and only prints IDs.
    # Copy that diagnostic into the unified report; do not rerun/repair the hole.
    if report is not None and "VTK - Error in:" in captured.getvalue():
        tail = captured.getvalue().split("VTK - Error in:", 1)[1].strip()
        faulty_ids = ast.literal_eval(tail.splitlines()[0])
        for bh in faulty_ids:
            report.add("geological", "qnd_compositing", None, bh, "compositing",
                       "legacy_compositing_exception", action="legacy_compositing_skipped")
    tertiary = u.len_wei_compositing(survey, layers, compositing_interval=cfg["interval"], verticaloption=True)
    return comp, tertiary


def audit_base_geometry(qb, survey, report):
    """Compare extracted base to the actual collar without changing extraction."""
    lookup = survey.set_index("ObjektIDObjektID")
    for idx, row in qb.iterrows():
        if row.borehole_id not in lookup.index:
            report.add("geological", "qnd_compositing", int(idx), row.borehole_id,
                       "ObjektID", "base_id_not_in_survey")
            continue
        collar = lookup.loc[row.borehole_id]
        if isinstance(collar, pd.DataFrame):
            continue  # Duplicate identity already reported in raw audit.
        if not np.array_equal([row.xcoord, row.ycoord], [collar.xcoord, collar.ycoord]):
            report.add("geological", "qnd_compositing", int(idx), row.borehole_id,
                       "coordinates", "base_collar_coordinate_mismatch", str([row.xcoord, row.ycoord]),
                       str([collar.xcoord, collar.ycoord]))
        if np.isfinite(row.qbasis):
            depth = float(collar.Ansatzhoeh) - row.qbasis
            end = pd.to_numeric(collar.Endteufe, errors="coerce")
            if depth < 0 or (np.isfinite(end) and depth > end):
                report.add("geological", "qnd_compositing", int(idx), row.borehole_id,
                           "qbasis", "base_outside_logged_depth", depth, end)


def statistical_preprocessing(qbasis, cfg, u):
    # GENERAL modelling preprocessing; these legacy functions DO filter data.
    # Diagnostic call is separate; detrend_2D makes its own identical trend.
    buffered = u.buffer_grid_params(cfg["grid"], max(cfg["radius_xy"]))
    data = np.asarray(qbasis)[:, 1:].astype(float)
    data = data[~np.isnan(data[:, 2])]
    data = data[data[:, 2] > np.mean(data[:, 2]) - cfg["tolerance"]]
    moving = u.moving_average_2d_kdtree(np.c_[np.ones(len(data)), data], buffered,
                                      *cfg["radius_xy"], tolerance=cfg["tolerance"])
    coords, residual, trend = u.detrend_2D(qbasis, cfg["grid"], moving_window=True,
        radius_xy=cfg["radius_xy"], tolerance=cfg["tolerance"], outlier_tolerance=cfg["outlier_tolerance"])
    mask = ~np.isnan(residual[:, 2].astype(float))
    coords, residual = coords[mask], residual[mask]
    transformed, method, parameters = u.transform(residual, cfg["trans_method"])
    return dict(coords=coords, coords_d=residual, grid_data=trend, coords_t=transformed,
                trans_method=method, trans_params=parameters, moving=moving,
                buffered_grid=buffered)


def selection_audit(qbasis, stats, cfg, u):
    result = pd.DataFrame(qbasis, columns=["borehole_id", "xcoord", "ycoord", "qbasis"])
    xyz = result[["xcoord", "ycoord", "qbasis"]].to_numpy(float)
    z = xyz[:, 2]
    finite_z = z[~np.isnan(z)]
    cut = np.mean(finite_z) - cfg["tolerance"]
    loc, inside = u.checkgrid(*xyz.T, stats["buffered_grid"])
    trend = np.full(len(z), np.nan)
    trend[inside] = stats["moving"][0][loc[inside]]
    result["legacy_local_trend"] = trend
    result["legacy_residual_before_selection"] = z - trend
    result["legacy_global_lower_cutoff"] = cut
    # Actual membership of returned data, rather than inferred acceptance.
    keys = set(map(tuple, stats["coords"].tolist()))
    result["legacy_selected"] = [tuple(row) in keys for row in xyz]
    result["selection_reason"] = np.select(
        [np.isnan(z), z <= cut, ~inside, np.isnan(trend), np.abs(z - trend) > cfg["outlier_tolerance"]],
        ["missing_base", "legacy_global_lower_cutoff", "outside_buffered_grid",
         "unsupported_legacy_trend", "legacy_residual_threshold"], default="retained")
    result.loc[~result.legacy_selected & result.selection_reason.eq("retained"), "selection_reason"] = "other_legacy_exclusion_review"
    return result


def compare_array(name, actual, expected):
    a, b = np.asarray(actual), np.asarray(expected)
    if a.shape != b.shape:
        return dict(name=name, passed=False, actual_shape=list(a.shape), expected_shape=list(b.shape))
    try:
        a, b = a.astype(float), b.astype(float)
        passed = np.array_equal(a, b, equal_nan=True)
        finite = np.isfinite(a) & np.isfinite(b)
        difference = float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.0
    except (ValueError, TypeError):
        passed = np.array_equal(a.astype(str), b.astype(str))
        difference = None
    return dict(name=name, passed=bool(passed), shape=list(a.shape), max_absolute_difference=difference)


def verify_direct_legacy(paths, extent, cfg, u, actual_comp, actual_tertiary, stats):
    """Independent direct-call reference from the applicable W blocks.

    W's stale checkpoint reloads/undefined interactive variables are intentionally
    not executed. This verifies refactoring against the current utility functions,
    not the provenance or geological correctness of historical saved results.
    """
    survey, layers = u.LFU_2_YOU(*map(str, paths), extent, False, (cfg["angle"], cfg["origin"]))
    comp = u.qnd_compositing(survey, layers, cfg["interval"], cfg["vtk"], True, True, cfg["qbasistol"])
    tertiary = u.len_wei_compositing(survey, layers, cfg["interval"], True)
    coords, coords_d, trend = u.detrend_2D(comp[3], cfg["grid"], True,
        cfg["radius_xy"], cfg["tolerance"], cfg["outlier_tolerance"])
    mask = ~np.isnan(coords_d[:, 2].astype(float))
    coords, coords_d = coords[mask], coords_d[mask]
    coords_t, method, parameters = u.transform(coords_d, cfg["trans_method"])
    results = [compare_array("qbasis_ids", np.asarray(actual_comp[3])[:, 0], np.asarray(comp[3])[:, 0]),
        compare_array("qbasis", np.asarray(actual_comp[3])[:, 1:], np.asarray(comp[3])[:, 1:]),
        compare_array("points", actual_comp[0], comp[0]),
        compare_array("tertiary_points", actual_tertiary, tertiary),
        compare_array("coords", stats["coords"], coords),
        compare_array("coords_d", stats["coords_d"], coords_d),
        compare_array("grid_data", stats["grid_data"], trend),
        compare_array("coords_t", stats["coords_t"], coords_t)]
    if method == "Johnson":
        results.extend(compare_array(f"trans_params_{i}", a, b)
                       for i, (a, b) in enumerate(zip(stats["trans_params"], parameters)))
    return {"reference": "direct legacy calls on the same raw inputs and settings",
            "historical_baseline": False, "passed": all(r["passed"] for r in results), "checks": results}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wd", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "real_munich")
    parser.add_argument("--output", type=Path, required=True, help="New empty run directory; existing results never overwritten")
    parser.add_argument("--utilities-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--example", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--qc-only", action="store_true")
    parser.add_argument("--config", type=Path, help="Repository profile JSON; modelling changes are explicit")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    cfg = defaults()
    if args.example:
        cfg["grid"].update(nx=100, ny=100, nz=200, xmn=690405.0, ymn=5334505.0)
        cfg["angle"] = 0
    if args.config:
        overrides = json.loads(args.config.read_text(encoding="utf-8"))["preparation"]
        for key, value in overrides.items():
            if key == "grid": cfg["grid"].update(value)
            elif key in cfg: cfg[key] = value
            else: raise ValueError(f"Unknown preparation setting: {key}")
    sys.path.insert(0, str(args.utilities_dir.resolve()))
    import UTILITIES_Wang as u
    report = Report()
    manifest = {"status": "started", "config": cfg, "wd": str(args.wd.resolve()),
                "qc_scope": "global IDs/geometry; collar and interval checks within legacy buffer",
                "vertical_coordinates": "absolute elevations; example is a 100 m thick slab, not terrain-following depth",
                "python": sys.version, "dependencies": {}, "source_hashes": {}}
    for package in ["numpy", "pandas", "scipy", "scikit-learn", "geopandas", "pyvista", "rasterio", "pykrige"]:
        manifest["dependencies"][package] = importlib.metadata.version(package)
    for name in ["UTILITIES_Wang.py", "GSLIB_parfl_011.py"]:
        manifest["source_hashes"][name] = sha256(args.utilities_dir / name)
    if args.config: manifest["config_digest"] = sha256(args.config)
    manifest["stage01_sha256"] = sha256(__file__)
    try:
        with (out / "preparation.log").open("w", encoding="utf-8") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            paths, extent, raw = audit_raw(args.wd, cfg, u, report, out)
            manifest["input_hashes"] = {str(p): sha256(p) for p in [*paths, paths[0].with_suffix(".dbf"), paths[0].with_suffix(".prj")]}
            report.save(out / "data_quality_report.csv")
            if not args.qc_only:
                survey, layers = u.LFU_2_YOU(*map(str, paths), extent, False, (cfg["angle"], cfg["origin"]))
                survey.to_csv(out / "survey_legacy_selected.csv", index=True)
                # Preconditions stop unsafe calls, never sort/repair/filter here.
                relevant = layers[layers.ObjektID.isin(survey.ObjektIDObjektID)]
                numeric = pd.to_numeric(relevant.Untergrenz, errors="coerce")
                elevations = pd.to_numeric(survey.Ansatzhoeh, errors="coerce")
                if len(survey) == 0 or not np.isfinite(numeric).all() or not np.isfinite(elevations).all():
                    raise ValueError("Empty/nonfinite modelling input; see raw audit")
                comp, tertiary = geological_preprocessing(survey, layers, cfg, u, report)
                points, cloud, line, qbasis = comp
                for name, values in [("points_rotated", points), ("tertiary_points", tertiary), ("qbasis", qbasis)]:
                    np.save(out / f"{name}.npy", values)
                if cfg["vtk"]:
                    cloud.save(out / "composited_points.vtp")
                    line.save(out / "borehole_lines.vtp")
                qb = pd.DataFrame(qbasis, columns=["borehole_id", "xcoord", "ycoord", "qbasis"])
                for column in ["xcoord", "ycoord", "qbasis"]:
                    qb[column] = pd.to_numeric(qb[column], errors="coerce")
                audit_base_geometry(qb, survey, report)
                qbqc = local_plausibility(qb, "qbasis", cfg["qc_radius"], cfg["qc_min_neighbours"],
                                           cfg["qb_residual_limit"], report, "qnd_compositing", "post_compositing_QC")
                qbqc.to_csv(out / "qbasis_plausibility.csv", index=False)
                for idx in qb.index[qb.qbasis.isna()]:
                    report.add("geological", "qnd_compositing", int(idx), qb.at[idx, "borehole_id"],
                               "qbasis", "base_not_identified", action="legacy_excluded")
                stats = statistical_preprocessing(qbasis, cfg, u)
                for name in ["coords", "coords_d", "grid_data", "coords_t"]:
                    np.save(out / f"{name}.npy", stats[name])
                np.save(out / "moving_average_buffer.npy", stats["moving"][0])
                np.save(out / "moving_average_neighbour_count.npy", stats["moving"][1])
                with (out / "transformation_parameters.pkl").open("wb") as stream:
                    pickle.dump((stats["trans_method"], stats["trans_params"]), stream)
                selection = selection_audit(qbasis, stats, cfg, u)
                selection.to_csv(out / "legacy_statistical_selection.csv", index=False)
                for _, row in selection.loc[~selection.legacy_selected].iterrows():
                    report.add("modelling_preprocessing", "detrend_2D", None, row.borehole_id,
                               "qbasis", row.selection_reason, row.qbasis, row.legacy_local_trend,
                               row.legacy_residual_before_selection, cfg["outlier_tolerance"], action="legacy_excluded")
                if not np.isfinite(stats["coords_t"]).all():
                    raise ValueError("Legacy transformation produced nonfinite values; outputs retained for review")
                u.PANDASDF2GSLIBGeoEAS(pd.DataFrame(stats["coords_t"], columns=["X", "Y", "Z"]), str(out / "qbasis_no_check_300.gslib"))
                cs = u.save_Cl_Sa_points(tertiary, str(out / "tertiary.gslib"))
                np.save(out / "Cl_Sa_points.npy", cs)
                manifest["counts"] = dict(survey=len(survey), qbasis=len(qbasis), tertiary_points=len(tertiary), coords_d=len(stats["coords_d"]))
                if args.verify:
                    verification = verify_direct_legacy(paths, extent, cfg, u, comp, tertiary, stats)
                    (out / "numerical_comparison.json").write_text(json.dumps(verification, indent=2), encoding="utf-8")
                    if not verification["passed"]:
                        raise AssertionError("Legacy numerical comparison failed")
            manifest["status"] = "qc_only_complete" if args.qc_only else "complete"
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        report.add("execution", "stage01", None, None, "preparation", "preparation_stopped",
                   manifest["error"], severity="error", action="stop_preparation")
        raise
    finally:
        manifest["core_sources_unchanged"] = all(sha256(args.utilities_dir / name) == digest
                                                for name, digest in manifest["source_hashes"].items())
        report.save(out / "data_quality_report.csv")
        (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Stage 01: {manifest['status']}. Results: {out}")


if __name__ == "__main__":
    main()
