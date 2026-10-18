"""Independent diagnostics for the saved GRID D2.2 simulations and borehole CV.

Run from any directory (original source and simulation files are read only)::

    python validation_model.py --run-dir <run> --output-dir <new-directory>

Requires numpy, scipy, pandas, matplotlib and scikit-learn (use the version that
created cv/fold_*/transform.pkl, here 1.5.2). No GSLIB execution is required.
Only load trusted local model outputs: continuous CV uses their pickled transforms.

The regular-grid estimator is GAM's non-standardized semivariogram:
sum((z(u+h)-z(u))**2)/(2*N(h)), using every valid axial pair, per realization.
FFT accelerates the pair sums without periodic wraparound or pair subsampling.
SGSIM is checked in Gaussian residual space; SISIM is checked as its 0/1 indicator.
Never calculate a variogram of the ensemble mean, or of numeric final class labels.
Normalized panels divide EACH realization by its own population variance and the
input curve by its total sill. They diagnose shape; they do not fix raw mismatch.

Categorical calibration uses saved OUT-OF-FOLD final hierarchical probabilities.
Continuous calibration reconstructs held-out boundary ensembles using each fold's
own transform and trend. No refitting/recalibration is done. Fixed full-data
variograms make this conditional CV, not nested model-selection validation.
Bootstrap intervals resample original CV borehole groups within folds; they do not
account for spatial dependence BETWEEN groups or refitting variability.

Acceptance limits are optional project choices, NEVER universal defaults. Use
--ece-limit, --coverage-tolerance and --variogram-nrmse-limit only after agreeing
the intended use. Existing output directories are refused to prevent overwrites.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.fft import rfft, irfft, next_fast_len

CLASSES = ["Clay", "Sand", "Quaternary"]
PROB_COLS = ["p_clay", "p_sand", "p_gravel"]


def dump(path, value):
    def convert(x):
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, Path):
            return str(x)
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False), encoding="utf-8")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def axial_variogram(values, axis, nlag):
    """Exact all-pair, nonperiodic grid estimator for finite, complete fields."""
    a = np.asarray(values, dtype=np.float64)
    if not np.isfinite(a).all():
        raise ValueError("Nonfinite simulation values: no implicit missing-data repair")
    n = a.shape[axis]
    k = np.arange(1, min(nlag, n - 1) + 1)
    if not len(k):
        raise ValueError("Direction has no cell pairs")
    # Remove a constant to avoid cancellation for large absolute coordinates.
    a = a - a.mean()
    other = tuple(i for i in range(a.ndim) if i != axis)
    squares = np.sum(a * a, axis=other)
    prefix = np.r_[0.0, np.cumsum(squares)]
    nfft = next_fast_len(2 * n - 1)
    ft = rfft(a, n=nfft, axis=axis, workers=1)
    power = np.sum(ft.real**2 + ft.imag**2, axis=other)
    cross = irfft(power, n=nfft)[k]
    pairs = (n - k) * (a.size // n)
    gamma = (prefix[n-k] + prefix[n] - prefix[k] - 2 * cross) / (2 * pairs)
    return k, np.maximum(gamma, 0.0), pairs


def theoretical(h, model, direction):
    """GSLIB spherical/exponential/Gaussian, zero rotation only (fail explicitly)."""
    h = np.asarray(h, dtype=float)
    gamma = np.full(h.shape, float(model["c0"]))
    for row in model["nest"]:
        kind, contribution, azimuth, dip, rotation, major, minor, vertical = row
        if any(abs(v) > 1e-10 for v in (azimuth, dip, rotation)):
            raise ValueError("This validation supports zero-angle models only")
        # GSLIB azimuth 0: major axis North; minor axis East.
        distance_range = {"E-W": minor, "N-S": major, "Vertical": vertical}[direction]
        if distance_range <= 0:
            raise ValueError("Nonpositive variogram range")
        r = h / distance_range
        if kind == 1:
            term = np.where(r < 1, 1.5*r - 0.5*r**3, 1.0)
        elif kind == 2:
            term = -np.expm1(-3*r)
        elif kind == 3:
            term = -np.expm1(-3*r*r)
        else:
            raise ValueError(f"Unsupported structure type {kind}")
        gamma += contribution * term
    gamma[h == 0] = 0.0
    return gamma


def total_sill(model):
    return float(model["c0"] + sum(row[1] for row in model["nest"]))


def load_run(run):
    m = json.loads((run / "03/run_manifest.json").read_text(encoding="utf-8"))
    if m["status"] != "complete":
        raise ValueError("Stage 03 is not complete")
    cfg = m["config"]
    sg = json.loads((run / "03/sgsim_parameters.json").read_text())
    si = json.loads((run / "03/sisim_parameters_template.json").read_text())
    # Actual simulation parameters take precedence over example configs or fits.
    continuous = dict(c0=sg["c0"], nest=sg["variogram"])
    _, c0, nest = si["variogram"][0]
    categorical = dict(c0=c0, nest=nest)
    if si["variogram"][0] != si["variogram"][1]:
        raise ValueError("Binary categories have different covariance models")
    for model in [continuous, categorical]:
        if total_sill(model) <= 0:
            raise ValueError("Nonpositive sill")
    return cfg, dict(continuous=continuous, categorical=categorical), si


def reproduction(run, out, grid, models):
    nx, ny, nz = (grid[k] for k in ["nx", "ny", "nz"])
    rows, summaries, histogram = [], [], []
    indicator_variance = float(np.var(np.load(run / "01/Cl_Sa_points.npy")[:, 3]))
    for name, filename, shape, directions in [
        ("continuous", "03/sgsim_gaussian.npy", (1, ny, nx),
         [("E-W", 2, grid["xsiz"], 50), ("N-S", 1, grid["ysiz"], 50)]),
        ("categorical", "03/sisim/realizations.npy", (nz, ny, nx),
         [("E-W", 2, grid["xsiz"], 30), ("N-S", 1, grid["ysiz"], 30),
          ("Vertical", 0, grid["zsiz"], 30)]),
    ]:
        sim = np.load(run / filename, mmap_mode="r")
        if sim.ndim != 2 or sim.shape[1] != int(np.prod(shape)):
            raise ValueError(f"Wrong shape: {filename}: {sim.shape}")
        model = models[name]
        sill = total_sill(model)
        cubes = {d[0]: [] for d in directions}
        variances = []
        for realization, field in enumerate(sim, 1):
            field = np.asarray(field).reshape(shape)
            if name == "categorical" and not np.isin(field, [0, 1]).all():
                raise ValueError("SISIM must be the binary 0/1 intermediate, not final labels")
            variance = float(np.var(field, dtype=np.float64))
            if variance <= 0:
                raise ValueError("Constant realization: normalized variogram undefined")
            variances.append(variance)
            histogram.append(dict(component=name, realization=realization,
                                  mean=float(field.mean()), variance=variance, input_sill=sill))
            for direction, axis, step, nl in directions:
                lag, gamma, pairs = axial_variogram(field, axis, nl)
                h = lag * step
                expected = theoretical(h, model, direction)
                cubes[direction].append(gamma)
                for j in range(len(lag)):
                    rows.append(dict(component=name, direction=direction, realization=realization,
                        lag=int(lag[j]), distance_m=h[j], pairs=int(pairs[j]),
                        semivariance=gamma[j], normalized_semivariance=gamma[j]/variance,
                        input_semivariance=expected[j], input_normalized=expected[j]/sill))
            if realization % 10 == 0:
                print(f"{name}: {realization}/{len(sim)} realizations", flush=True)
        fig, axes = plt.subplots(2, len(directions), figsize=(5*len(directions), 8), squeeze=False)
        for column, (direction, axis, step, nl) in enumerate(directions):
            curves = np.asarray(cubes[direction])
            h = np.arange(1, curves.shape[1]+1)*step
            target = theoretical(h, model, direction)
            mean = curves.mean(axis=0)
            lo, hi = np.quantile(curves, [.025, .975], axis=0)
            normcurves = curves / np.asarray(variances)[:, None]
            normmean = normcurves.mean(axis=0)
            ranges = {"E-W": 6, "N-S": 5, "Vertical": 7}
            relevant = h <= max(row[ranges[direction]] for row in model["nest"])
            summaries.append(dict(component=name, direction=direction,
                realizations=len(sim), input_sill=sill, ensemble_mean_variance=float(np.mean(variances)),
                raw_nrmse=float(np.sqrt(np.mean((mean-target)**2))/sill),
                raw_nrmse_within_range=float(np.sqrt(np.mean((mean[relevant]-target[relevant])**2))/sill) if relevant.any() else None,
                normalized_rmse=float(np.sqrt(np.mean((normmean-target/sill)**2))),
                input_inside_pointwise_95pct_fraction=float(np.mean((target>=lo)&(target<=hi))),
                mean_last_lag=float(mean[-1]), input_last_lag=float(target[-1]), max_distance_m=float(h[-1])))
            pd.DataFrame(dict(distance_m=h, input_semivariance=target, ensemble_mean=mean,
                realization_p025=lo, realization_p975=hi, normalized_mean=normmean,
                input_normalized=target/sill)).to_csv(out/f"variogram_{name}_{direction}.csv", index=False)
            expfile = (f"02/continuous_direction_{column}.csv" if name == "continuous" else
                       "02/categorical_vertical_direction_0.csv" if direction == "Vertical" else
                       f"02/categorical_horizontal_direction_{column}.csv")
            exp = pd.read_csv(run / expfile)
            exp = exp[(exp.npairs > 0) & (exp.avsepdist > 0) & (exp.avsepdist <= h[-1])]
            expraw = exp.svargval.to_numpy() * (indicator_variance if direction == "Vertical" else 1)
            for r, (values, theory) in enumerate([(curves, target), (normcurves, target/sill)]):
                ax = axes[r, column]
                low, high = np.quantile(values, [.025, .975], axis=0)
                ax.plot(h, values.T, color="#97afbd", alpha=.15, lw=.5)
                ax.fill_between(h, low, high, color="#6fa0b8", alpha=.25, label="95% realization envelope")
                ax.plot(h, values.mean(axis=0), color="#006b8f", lw=2, label="Ensemble mean variogram")
                ax.plot(h, theory, color="#c43831", lw=2, ls="--", label="Actual input model")
                if r == 0:
                    ax.scatter(exp.avsepdist, expraw, s=15, color="#353535", label="Input GAMV (angular tolerance)")
                ax.set(title=f"{direction} | {'raw' if r == 0 else 'normalized shape'}",
                       xlabel="Lag distance [m]", ylabel="Semivariance" if r == 0 else "Semivariance / variance")
                ax.grid(alpha=.2)
                ax.legend(fontsize=7)
        fig.suptitle(f"{name.upper()}: {len(sim)} complete-grid realizations", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, .96))
        fig.savefig(out/f"variogram_{name}.png", dpi=180)
        plt.close(fig)
    pd.DataFrame(rows).to_csv(out/"variogram_all_realizations.csv", index=False)
    pd.DataFrame(histogram).to_csv(out/"realization_moments.csv", index=False)
    pd.DataFrame(summaries).to_csv(out/"variogram_metrics.csv", index=False)
    return dict(metrics=summaries, conditioning_indicator_variance=indicator_variance)


def cluster_bootstrap_weights(groups, folds, repeats, seed):
    """Bootstrap CV groups within folds, retaining all composites in a group."""
    unique, idx = np.unique(groups, return_inverse=True)
    gf = np.zeros(len(unique), int)
    for j in range(len(unique)):
        f = np.unique(folds[idx == j])
        if len(f) != 1:
            raise ValueError("A validation group appears in multiple folds")
        gf[j] = f[0]
    rng = np.random.default_rng(seed)
    weights = np.zeros((repeats, len(unique)), dtype=float)
    for fold in np.unique(gf):
        take = np.flatnonzero(gf == fold)
        weights[:, take] = rng.multinomial(len(take), np.full(len(take), 1/len(take)), size=repeats)
    return idx, weights


def weighted_ratio_ci(numerator, denominator, idx, weights):
    n = np.bincount(idx, weights=np.asarray(numerator, float), minlength=weights.shape[1])
    d = np.bincount(idx, weights=np.asarray(denominator, float), minlength=weights.shape[1])
    den = weights @ d
    value = (weights @ n)[den > 0] / den[den > 0]
    return np.quantile(value, [.025, .975]) if len(value) else np.array([np.nan, np.nan])


def reliability(prob, outcome, groups, folds, bins, repeats, seed):
    edges = np.linspace(0, 1, bins+1)
    # Match sklearn uniform bins: (lower, upper], with p=0 in first bin.
    which = np.searchsorted(edges[1:-1], prob)
    idx, weights = cluster_bootstrap_weights(groups, folds, repeats, seed)
    rows = []
    for b in range(bins):
        mask = which == b
        n = int(mask.sum())
        if n:
            low, high = weighted_ratio_ci(mask*outcome, mask, idx, weights)
            mp, obs = float(prob[mask].mean()), float(outcome[mask].mean())
            gaplow, gaphigh = weighted_ratio_ci(mask*(outcome-prob), mask, idx, weights)
        else:
            mp = obs = low = high = gaplow = gaphigh = np.nan
        rows.append(dict(bin=b, lower=edges[b], upper=edges[b+1], count=n,
            groups=int(len(np.unique(groups[mask]))), mean_probability=mp, observed_frequency=obs,
            ci_low=low, ci_high=high, gap_ci_low=gaplow, gap_ci_high=gaphigh))
    table = pd.DataFrame(rows)
    good = table[table["count"] > 0]
    ece = float(np.sum(good["count"] * abs(good.mean_probability-good.observed_frequency))/len(prob))
    return table, ece


def categorical_calibration(run, out, repeats, seed):
    assignments = pd.read_csv(run/"cv/borehole_folds.csv", dtype={"borehole_id": str})
    if assignments.groupby("borehole_id").fold.nunique().max() != 1:
        raise ValueError("Borehole IDs cross folds")
    mapping = assignments.drop_duplicates("borehole_id").set_index("borehole_id")
    frames = []
    for fold in sorted(assignments.fold.unique()):
        f = pd.read_csv(run/f"cv/fold_{fold}/predictions.csv", dtype={"borehole_id": str})
        f["fold"] = fold
        matched = mapping.reindex(f.borehole_id)
        if matched["group"].isna().any() or not (matched.fold.to_numpy() == fold).all():
            raise ValueError("Predictions do not match withheld boreholes")
        metrics = pd.read_csv(run/f"cv/fold_{fold}/metrics.csv")
        if int(metrics.shared_train_test_voxels.iloc[0]) != 0:
            raise ValueError("Saved CV reports shared train/test voxels")
        f["group"] = matched["group"].to_numpy()
        frames.append(f)
    data = pd.concat(frames, ignore_index=True)
    p = data[PROB_COLS].to_numpy()
    truth = data.observed.to_numpy(dtype=int)
    if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any() or not np.allclose(p.sum(axis=1), 1, atol=1e-8):
        raise ValueError("Invalid held-out probabilities")
    if not np.isin(truth, [1, 2, 3]).all():
        raise ValueError("Unexpected observed category")
    target = truth[:, None] == np.arange(1, 4)
    pred = data.predicted.to_numpy(dtype=int)
    if not np.isin(pred, [1, 2, 3]).all() or not np.allclose(p[np.arange(len(p)), pred-1], p.max(axis=1)):
        raise ValueError("Saved predicted classes are not probability maximizers")
    group, fold = data["group"].to_numpy(), data.fold.to_numpy()
    correct = (pred == truth)
    evaluations = [(CLASSES[k], p[:, k], target[:, k]) for k in range(3)]
    evaluations += [("Top-label confidence", p.max(axis=1), correct)]
    scores = []
    fig, axes = plt.subplots(2, 4, figsize=(17, 7), gridspec_kw={"height_ratios": [3, 1]})
    for k, (name, prob, outcome) in enumerate(evaluations):
        table, ece = reliability(prob, outcome, group, fold, 10, repeats, seed)
        table.to_csv(out/f"calibration_{k}_bins.csv", index=False)
        other_bins = {}
        for nbin in [5, 20]:
            extra, score = reliability(prob, outcome, group, fold, nbin, repeats, seed)
            extra.to_csv(out/f"calibration_{k}_{nbin}_bins.csv", index=False)
            other_bins[str(nbin)] = score
        scores.append(dict(target=name, ece_10_bins=ece, ece_other_bins=other_bins,
            mean_probability=float(prob.mean()), observed_frequency=float(np.mean(outcome)),
            binary_brier=float(np.mean((prob-outcome)**2))))
        ax = axes[0, k]
        ax.plot([0, 1], [0, 1], ls="--", color="black", lw=1, label="Ideal")
        # NaNs preserve gaps for empty bins.
        ax.plot(table.mean_probability, table.observed_frequency, "o-", color="#006b8f", ms=4)
        ax.vlines(table.mean_probability, table.ci_low, table.ci_high, color="#006b8f", lw=1.5)
        ax.set(xlim=(-.02, 1.02), ylim=(-.02, 1.02), title=f"{name}\nECE = {ece:.3f}",
               xlabel="Mean predicted probability", ylabel="Observed frequency")
        ax.grid(alpha=.2)
        axes[1, k].bar((table.lower+table.upper)/2, table["count"], width=.085, color="#6fa0b8")
        axes[1, k].set(xlim=(0, 1), xlabel="Probability bin", ylabel="Composites")
    fig.suptitle("Out-of-fold hierarchical model | group bootstrap 95% intervals", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(out/"calibration_categorical.png", dpi=180)
    plt.close(fig)
    brier = np.sum((p-target)**2, axis=1)
    baseline = np.zeros_like(p)
    for f in np.unique(fold):
        # Baseline uses other folds' outcomes only, never this fold's outcomes.
        baseline[fold == f] = target[fold != f].mean(axis=0)
    base_brier = np.sum((baseline-target)**2, axis=1)
    foldrows = []
    for f in np.unique(fold):
        mask = fold == f
        foldrows.append(dict(fold=int(f), n=int(mask.sum()), accuracy=float(correct[mask].mean()),
            brier=float(brier[mask].mean()), baseline_brier=float(base_brier[mask].mean()),
            mean_confidence=float(p[mask].max(axis=1).mean())))
    pd.DataFrame(foldrows).to_csv(out/"calibration_fold_metrics.csv", index=False)
    data.assign(correct=correct, confidence=p.max(axis=1), brier=brier).to_csv(out/"heldout_predictions.csv", index=False)
    chosen = p[np.arange(len(p)), truth-1]
    idx, weights = cluster_bootstrap_weights(group, fold, repeats, seed)
    accuracy_ci = weighted_ratio_ci(correct, np.ones(len(data)), idx, weights)
    gap_ci = weighted_ratio_ci(p.max(axis=1)-correct, np.ones(len(data)), idx, weights)
    return dict(n=len(data), boreholes=data.borehole_id.nunique(), groups=len(np.unique(group)),
        accuracy=float(correct.mean()), accuracy_cluster_ci=accuracy_ci,
        borehole_macro_accuracy=float(data.assign(correct=correct).groupby("borehole_id").correct.mean().mean()),
        mean_confidence=float(p.max(axis=1).mean()), confidence_minus_accuracy_ci=gap_ci,
        multiclass_brier_sum=float(brier.mean()), baseline_brier_sum=float(base_brier.mean()),
        brier_skill_vs_other_fold_prevalence=float(1-brier.mean()/base_brier.mean()),
        baseline_scope="Prevalence from other held-out folds; descriptive baseline, not raw training-data prevalence",
        zero_probability_observations=int((chosen == 0).sum()),
        raw_log_loss="infinite" if (chosen == 0).any() else float(-np.log(chosen).mean()),
        classwise=scores, folds=foldrows)


def grid_xy(xy, grid):
    """Exact legacy checkgrid right=True boundary rule and x-fastest indexing."""
    indices = []
    valid = np.ones(len(xy), dtype=bool)
    for j, axis in enumerate(["x", "y"]):
        n, step, centre = grid["n"+axis], grid[axis+"siz"], grid[axis+"mn"]
        edges = np.arange(n+1)*step + centre - .5*step
        i = np.digitize(xy[:, j], edges, right=True)
        valid &= (i > 0) & (i <= n)
        indices.append(i-1)
    return indices[1]*grid["nx"]+indices[0], valid


def read_sg(path, cells):
    with path.open() as f:
        f.readline()
        # GSLIB SGSIM appends grid metadata after the variable count.
        nvar = int(f.readline().split()[0])
        if nvar != 1:
            raise ValueError("Expected single-variable SGSIM GeoEAS output")
        f.readline()
        values = np.loadtxt(f)
    if values.size % cells:
        raise ValueError("Truncated SGSIM output")
    return values.reshape(-1, cells)


def continuous_calibration(run, out, grid, repeats, seed):
    source = np.load(run/"cv/inputs.npz")
    qb, qfold = source["qb"], source["qfold"]
    assignments = pd.read_csv(run/"cv/borehole_folds.csv")
    mapping = assignments.drop_duplicates(["x", "y"]).set_index(["x", "y"])
    levels = np.arange(.05, 1, .05)
    coverage_levels = np.array([.1, .2, .3, .4, .5, .6, .7, .8, .9, .95])
    frames, all_quantiles, all_cover, all_width, all_samples = [], [], [], [], []
    excluded = 0
    cells = grid["nx"]*grid["ny"]
    for fold in sorted(np.unique(qfold)):
        test = qb[qfold == fold]
        loc, inside = grid_xy(test[:, :2], grid)
        excluded += int((~inside).sum())
        test, loc = test[inside], loc[inside]
        matched = mapping.reindex(pd.MultiIndex.from_arrays([test[:, 0], test[:, 1]], names=["x", "y"]))
        if matched["group"].isna().any() or not (matched.fold.to_numpy() == fold).all():
            raise ValueError("Boundary validation points do not map to withheld groups")
        folder = run/f"cv/fold_{fold}"
        # These files are trusted saved outputs from the user's existing workflow.
        with (folder/"transform.pkl").open("rb") as f:
            method, transform, trend = pickle.load(f)
        if method != "Quantile":
            raise ValueError(f"Continuous fold transform {method!r} unsupported; no silent substitution")
        gaussian = read_sg(folder/"sgsim.out", cells)
        samples = transform.inverse_transform(gaussian[:, loc].reshape(-1, 1)).reshape(len(gaussian), -1)
        samples += np.asarray(trend)[loc]
        if not np.isfinite(samples).all():
            raise ValueError("Nonfinite held-out boundary predictions")
        y = test[:, 2]
        quantiles = np.quantile(samples, levels, axis=0)
        intervals = np.quantile(samples, np.r_[(1-coverage_levels)/2, (1+coverage_levels)/2], axis=0)
        lower, upper = np.split(intervals, 2)
        cover = (y >= lower) & (y <= upper)
        pit = ((samples < y).sum(axis=0) + .5*(samples == y).sum(axis=0))/len(samples)
        frames.append(pd.DataFrame(dict(fold=fold, group=matched["group"].to_numpy(),
            x=test[:, 0], y=test[:, 1], observed=y, mean=samples.mean(axis=0),
            median=np.median(samples, axis=0), pit=pit, ensemble_size=len(samples))))
        all_quantiles.append(quantiles)
        all_cover.append(cover)
        all_width.append(upper-lower)
        all_samples.append(samples)
        print(f"Boundary calibration fold {fold+1}: {len(y)} held-out observations", flush=True)
    data = pd.concat(frames, ignore_index=True)
    samples = np.concatenate(all_samples, axis=1)
    q = np.concatenate(all_quantiles, axis=1)
    cover = np.concatenate(all_cover, axis=1)
    widths = np.concatenate(all_width, axis=1)
    y = data.observed.to_numpy()
    idx, weights = cluster_bootstrap_weights(data["group"].to_numpy(), data.fold.to_numpy(), repeats, seed)
    quantilerows, coverrows = [], []
    for j, level in enumerate(levels):
        observed = y <= q[j]
        low, high = weighted_ratio_ci(observed, np.ones(len(y)), idx, weights)
        quantilerows.append(dict(nominal=level, observed=float(observed.mean()), ci_low=low, ci_high=high))
    for j, level in enumerate(coverage_levels):
        low, high = weighted_ratio_ci(cover[j], np.ones(len(y)), idx, weights)
        coverrows.append(dict(nominal=level, observed=float(cover[j].mean()), ci_low=low, ci_high=high,
                              mean_width_m=float(widths[j].mean())))
    qt, ct = pd.DataFrame(quantilerows), pd.DataFrame(coverrows)
    qt.to_csv(out/"continuous_quantile_calibration.csv", index=False)
    ct.to_csv(out/"continuous_interval_coverage.csv", index=False)
    data.to_csv(out/"continuous_heldout_predictions.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    for ax, table, title in zip(axes[:2], [qt, ct], ["Quantile calibration", "Central interval coverage"]):
        ax.plot([0, 1], [0, 1], "k--", lw=1)
        ax.plot(table.nominal, table.observed, "o-", color="#006b8f", ms=4)
        ax.fill_between(table.nominal, table.ci_low, table.ci_high, alpha=.2, color="#6fa0b8")
        ax.set(title=title, xlabel="Nominal probability", ylabel="Observed fraction", xlim=(0, 1), ylim=(0, 1))
        ax.grid(alpha=.2)
    axes[2].hist(data.pit, bins=np.linspace(0, 1, 11), color="#6fa0b8", edgecolor="white")
    axes[2].axhline(len(data)/10, color="black", ls="--", lw=1)
    axes[2].set(title="Held-out boundary ranks (mid-PIT)", xlabel="Empirical rank in ensemble", ylabel="Boundary observations")
    fig.suptitle(f"Boundary cross-validation | {len(data)} observations | 95% group bootstrap intervals")
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(out/"calibration_continuous.png", dpi=180)
    plt.close(fig)
    # Empirical-ensemble CRPS: E|X-y| - 0.5 E|X-X'|, no pairwise matrix.
    n = len(samples)
    coefficient = (2*np.arange(1, n+1)-n-1)[:, None]
    crps = abs(samples-y).mean(axis=0) - np.sum(coefficient*np.sort(samples, axis=0), axis=0)/n**2
    return dict(n=len(data), groups=int(data["group"].nunique()), excluded_outside_grid=excluded,
        mae_ensemble_mean_m=float(abs(data["mean"]-data.observed).mean()),
        mae_ensemble_median_m=float(abs(data["median"]-data.observed).mean()),
        legacy_ensemble_mae_m=float(abs(samples-y).mean()), mean_crps_m=float(crps.mean()),
        coverage=coverrows, quantile_mean_absolute_gap=float(abs(qt.nominal-qt.observed).mean()),
        scope="All in-grid withheld qbasis, without truth-based residual filtering; training-only Quantile transform and trend")


def self_test():
    """Independent brute-force estimator and calibration edge-case checks."""
    rng = np.random.default_rng(21)
    for shape in [(4, 5, 7), (1, 4, 5)]:
        for a in [rng.normal(size=shape), rng.integers(0, 2, size=shape)]:
            for axis in range(3):
                if shape[axis] < 2:
                    continue
                lags, actual, pairs = axial_variogram(a, axis, 10)
                for j, lag in enumerate(lags):
                    left = [slice(None)]*3
                    right = [slice(None)]*3
                    left[axis] = slice(None, -lag)
                    right[axis] = slice(lag, None)
                    delta = a[tuple(left)].astype(float)-a[tuple(right)].astype(float)
                    np.testing.assert_allclose(actual[j], .5*np.mean(delta**2), atol=1e-12)
                    assert pairs[j] == delta.size
    p = np.array([0., 0., .5, .5, 1., 1.])
    outcome = np.array([0, 0, 0, 1, 1, 1])
    table, ece = reliability(p, outcome, np.arange(6), np.zeros(6), 10, 30, 2)
    assert table["count"].sum() == 6 and ece == 0
    assert table["count"].iloc[0] == 2 and table["count"].iloc[-1] == 2
    model = dict(c0=.1, nest=[[2, .9, 0, 0, 0, 100, 50, 5]])
    for d, r in [("N-S", 100), ("E-W", 50), ("Vertical", 5)]:
        np.testing.assert_allclose(theoretical([0, r], model, d), [0, .1+.9*(1-np.exp(-3))])
    print("PASS: all-pair FFT vs independent direct estimator, axes/pair counts, bins 0/1, anisotropic input curves")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261015)
    parser.add_argument("--ece-limit", type=float)
    parser.add_argument("--coverage-tolerance", type=float)
    parser.add_argument("--variogram-nrmse-limit", type=float)
    parser.add_argument("--self-test", action="store_true")
    a = parser.parse_args()
    if a.self_test:
        self_test()
        return
    if a.run_dir is None or a.output_dir is None:
        parser.error("--run-dir and --output-dir are required")
    if a.bootstrap < 100:
        parser.error("Use at least 100 bootstrap replicates")
    for name in ["ece_limit", "coverage_tolerance", "variogram_nrmse_limit"]:
        if getattr(a, name) is not None and getattr(a, name) < 0:
            parser.error("Acceptance limits must be nonnegative")
    run, out = a.run_dir.resolve(), a.output_dir.resolve()
    cfg, models, si = load_run(run)
    out.mkdir(parents=True, exist_ok=False)
    import sklearn
    manifest = dict(status="started", source_run=run, script_sha256=digest(Path(__file__)),
        python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
        sklearn=sklearn.__version__, bootstrap=a.bootstrap, seed=a.seed, grid=cfg["preparation"]["grid"],
        models=models, indicator_global_pdf=si["gpdfx"][:2],
        cv_scope=json.loads((run/"cv/scope.json").read_text()),
        caveats=["Conditional CV with fixed full-data variograms; not nested model selection",
                 "Composite-weighted calibration; spatial dependence between borehole groups remains",
                 "Realization envelope is pointwise ensemble spread, not a simultaneous test or mean confidence interval",
                 "Grid axial GAM and input angular-tolerance GAMV use different pair geometries",
                 "Complete rectangular component fields, without terrain or hierarchical masking",
                 "No input covariance is specified for the combined three-class hierarchy"],
        sources={})
    sources = [run/"03/run_manifest.json", run/"03/sgsim_parameters.json", run/"03/sisim_parameters_template.json",
               run/"03/sgsim_gaussian.npy", run/"03/sisim/realizations.npy", run/"cv/inputs.npz",
               run/"cv/borehole_folds.csv", run/"cv/scope.json"]
    sources += sorted((run/"cv").glob("fold_*/predictions.csv"))
    sources += sorted((run/"cv").glob("fold_*/sgsim.out"))
    sources += sorted((run/"cv").glob("fold_*/transform.pkl"))
    sources += sorted((run/"cv").glob("fold_*/metrics.csv"))
    sources += sorted((run/"02").glob("*_direction_*.csv"))
    sources += [run/"01/Cl_Sa_points.npy"]
    manifest["sources"] = {str(p.relative_to(run)): digest(p) for p in sources}
    dump(out/"manifest.json", manifest)
    try:
        v = reproduction(run, out, cfg["preparation"]["grid"], models)
        cat = categorical_calibration(run, out, a.bootstrap, a.seed)
        cont = continuous_calibration(run, out, cfg["preparation"]["grid"], a.bootstrap, a.seed)
        dump(out/"summary.json", dict(variogram=v, categorical=cat, continuous=cont))
        acceptance = dict(status="not_assessed_without_project_limits", checks=[],
            note="These diagnostics do not establish universal geological/engineering acceptance criteria.")
        for metric, limit, values in [
            ("classwise_and_confidence_ECE_10", a.ece_limit, [s["ece_10_bins"] for s in cat["classwise"]]),
            ("coverage_absolute_gap_at_80_90_95pct", a.coverage_tolerance,
             [abs(c["observed"]-c["nominal"]) for c in cont["coverage"] if c["nominal"] >= .8]),
            ("raw_variogram_NRMSE_all_directions", a.variogram_nrmse_limit, [m["raw_nrmse"] for m in v["metrics"]]),
        ]:
            acceptance["checks"].append(dict(metric=metric, limit=limit, observed_max=max(values),
                passed=None if limit is None else bool(max(values) <= limit)))
        if all(c["limit"] is not None for c in acceptance["checks"]):
            acceptance["status"] = "passed_user_limits" if all(c["passed"] for c in acceptance["checks"]) else "failed_user_limits"
        dump(out/"acceptance.json", acceptance)
        manifest["status"] = "complete"
    except Exception as e:
        manifest.update(status="failed", error=f"{type(e).__name__}: {e}")
        raise
    finally:
        dump(out/"manifest.json", manifest)
    print(f"Complete: {out}", flush=True)


if __name__ == "__main__":
    main()
