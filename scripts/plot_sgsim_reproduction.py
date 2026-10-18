"""SGSIM experimental/input/reproduced variograms in Gaussian residual space.

Default: redraw from compact bundled arrays, without simulation or CV.
--source-run <run> recalculates GAM-equivalent pair sums from saved Gaussian
realizations. It does not simulate or write to the source run.
All curves share the actual input total sill; the non-zero nugget is retained.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validation_model import axial_variogram, theoretical

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'pics/final_20261018/data'
OUT = ROOT / 'pics/final_20261018'
ORIGINAL = ROOT / 'pics/data'
RECORD = DATA / 'sgsim_reproduction_provenance.json'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(value):
    RECORD.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def extract(run):
    cfg = json.loads((DATA / 'exponential.json').read_text())
    model = cfg['variograms']['continuous']
    params = json.loads((run / '03/sgsim_parameters.json').read_text())
    assert params['c0'] == model['c0'] and params['variogram'] == model['nest']
    assert params['nsim'] == 100 and params['grid_params']['nz'] == 1
    assert (params['grid_params']['nx'], params['grid_params']['ny']) == (100, 100)
    paths = ['03/sgsim_gaussian.npy', '03/sgsim_parameters.json', '01/qbasis_no_check_300.gslib',
             '02/continuous_direction_0.csv', '02/continuous_direction_1.csv']
    before = {name: sha(run / name) for name in paths}
    if RECORD.exists():
        previous = json.loads(RECORD.read_text())
        assert previous['source_hashes'] == before, 'Source differs from bundled production SGSIM'
    values = np.load(run / '03/sgsim_gaussian.npy', mmap_mode='r')
    assert values.shape == (100, 10000) and np.isfinite(values).all()
    curves = {'EW': [], 'NS': []}
    moments = []
    worst = 0.0
    for i, row in enumerate(values):
        field = np.asarray(row).reshape(100, 100)
        moments.append(dict(realization=i+1, mean=float(field.mean()), variance=float(field.var())))
        for direction, axis in [('EW', 1), ('NS', 0)]:
            lag, gamma, pairs = axial_variogram(field, axis, 50)
            np.testing.assert_array_equal(pairs, (100-lag)*100)
            if i in [0, 49, 99]:
                for k in [1, 5, 15, 50]:
                    left, right = [slice(None)]*2, [slice(None)]*2
                    left[axis], right[axis] = slice(None, -k), slice(k, None)
                    direct = .5*np.mean((field[tuple(left)]-field[tuple(right)])**2)
                    error = abs(float(gamma[k-1])-float(direct))
                    worst = max(worst, error)
                    np.testing.assert_allclose(gamma[k-1], direct, rtol=0, atol=2e-12)
            curves[direction].append(gamma)
    for j, direction in enumerate(['EW', 'NS']):
        assert sha(run / f'02/continuous_direction_{j}.csv') == sha(ORIGINAL / f'continuous_direction_{j}.csv')
        np.save(DATA / f'sgsim_{direction}_curves.npy', np.asarray(curves[direction]))
    pd.DataFrame(moments).to_csv(DATA / 'sgsim_realization_moments.csv', index=False)
    assert before == {name: sha(run / name) for name in paths}
    save(dict(component='SGSIM', space='Gaussian residuals before trend/back-transform',
              realizations=100, grid_shape_yx=[100, 100], directions=['EW', 'NS'],
              model=model, source_hashes=before, source_files_unchanged=True,
              estimator='GAM-equivalent .5 * mean squared differences over all nonperiodic axial grid pairs, per realization',
              grid_lag_m=10, nlag=50, direct_pair_max_abs_error=worst,
              experimental='Unchanged original GAMV: nominal lag 5 m, tolerance 2.5 m; zero and near-zero bins omitted; use actual distances and pair counts',
              rerun_sgsim=False, rerun_cv=False))


def plot():
    record = json.loads(RECORD.read_text())
    model = record['model']
    cfg = json.loads((DATA / 'exponential.json').read_text())
    assert model == cfg['variograms']['continuous']
    sill = model['c0'] + sum(n[1] for n in model['nest'])
    assert model['c0'] > 0 and sill > 0
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.55), layout='constrained')
    metrics, exports = {}, []
    input_hashes = {}
    limit = 155
    for j, (direction, ax) in enumerate(zip(['EW', 'NS'], axes)):
        file = DATA / f'sgsim_{direction}_curves.npy'
        full = np.load(file)
        assert full.shape == (100, 50) and np.isfinite(full).all()
        all_h = np.arange(1, 51)*10
        h, curves = all_h[all_h <= limit], full[:, all_h <= limit]
        exp_file = ORIGINAL / f'continuous_direction_{j}.csv'
        exp = pd.read_csv(exp_file)
        # Original GAMV outputs zero, near-zero, then nominal 5 m, 10 m, ...
        exp = exp[(exp.lag > 2) & (exp.npairs > 0) & (exp.avsepdist > 0) & (exp.avsepdist <= limit)]
        mean = curves.mean(axis=0)
        lower, upper = np.quantile(curves, [.05, .95], axis=0)
        # Show the positive-distance right limit, without erasing the nugget.
        x = np.linspace(1e-8, limit, 600)
        label = 'E-W' if direction == 'EW' else 'N-S'
        expected = theoretical(x, model, label)
        ax.fill_between(h, lower/sill, upper/sill, color='#007f86', alpha=.18, linewidth=0,
                        label='5th-95th percentile', zorder=1)
        ax.plot(x, expected/sill, '--', color='#ce8a1a', lw=1.7, label='Theoretical input', zorder=2)
        ax.plot(h, mean/sill, color='#007f86', lw=3.2, label='Mean reproduced', zorder=3)
        ax.scatter(exp.avsepdist, exp.svargval/sill, s=24, color='#222222', label='Experimental', zorder=5)
        # gamma(0)=0 is an isolated identity, not a zero-nugget assumption.
        ax.scatter([0], [0], s=24, facecolors='none', edgecolors='#333333', clip_on=False, zorder=6)
        ax.set(xlim=(0, limit), ylim=(0, 1.2), title=direction,
               xlabel='Separation distance [m]', ylabel='Semivariance / common input sill')
        ax.grid(alpha=.16)
        target_at_grid = theoretical(h, model, label)
        overlap = exp[exp.avsepdist.between(all_h[0], all_h[-1])]
        metrics[direction] = dict(
            mean_vs_input_NRMSE=float(np.sqrt(np.mean((mean-target_at_grid)**2))/sill),
            mean_vs_experimental_NRMSE=float(np.sqrt(np.mean((np.interp(overlap.avsepdist, all_h, full.mean(axis=0))-overlap.svargval)**2))/sill),
            experimental_count=len(exp), reproduced_count=len(h),
            theoretical_inside_pointwise_90pct_fraction=float(np.mean((target_at_grid >= lower) & (target_at_grid <= upper))))
        for row in exp.itertuples():
            exports.append(dict(direction=direction, series='experimental', distance=row.avsepdist,
                                raw_gamma=row.svargval, normalized_gamma=row.svargval/sill, pairs=int(row.npairs)))
        for distance, raw, lo, hi in zip(h, mean, lower, upper):
            exports.append(dict(direction=direction, series='mean_reproduced', distance=distance,
                                raw_gamma=raw, normalized_gamma=raw/sill, lower=lo/sill, upper=hi/sill,
                                pairs=int((100-distance/10)*100)))
        for distance, raw in zip(x, expected):
            exports.append(dict(direction=direction, series='theoretical_input', distance=distance,
                                raw_gamma=raw, normalized_gamma=raw/sill))
        for p in [file, exp_file]:
            input_hashes[p.relative_to(ROOT).as_posix()] = sha(p)
        assert len(ax.lines) == 2
    handles, labels = axes[0].get_legend_handles_labels()
    lookup = dict(zip(labels, handles))
    order = ['Experimental', 'Theoretical input', 'Mean reproduced', '5th-95th percentile']
    axes[0].legend([lookup[n] for n in order], order, fontsize=8, loc='lower right')
    fig.suptitle('SGSIM spatial reproduction | 100 Gaussian residual realizations\n'
                 f'Common input sill = 1; non-zero nugget/sill = {model["c0"]/sill:.3f} retained', fontsize=11)
    for ext in ['png', 'pdf']:
        fig.savefig(OUT / f'06_sgsim_variogram_reproduction.{ext}', dpi=240)
    plt.close(fig)
    pd.DataFrame(exports).to_csv(DATA / 'sgsim_reproduction_plot_values.csv', index=False)
    record.update(common_input_sill=sill, nugget_over_sill=model['c0']/sill,
                  display='All three series share the same actual input total sill; no independent peak normalization',
                  origin='gamma(0)=0 shown separately; positive-distance input retains the nugget jump',
                  plot_input_hashes=input_hashes, display_xlim_m=limit, display_ylim=[0, 1.2], metrics=metrics)
    save(record)
    print(json.dumps(dict(input_sill=sill, nugget=model['c0'], metrics=metrics), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', type=Path)
    args = parser.parse_args()
    if args.source_run:
        extract(args.source_run.resolve())
    plot()
