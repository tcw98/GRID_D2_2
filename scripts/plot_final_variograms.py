"""Final Figures 07/08: direction-wise peak normalization, complete binary SISIM grid.

Each series is divided by its plotted maximum within each direction.
Raw amplitude errors are recorded alongside the shape comparison.
The script redraws figures and plot metadata from bundled realization curves.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'pics/final_20261018'
DATA = OUT / 'data'


def theoretical(h, model, direction):
    index = {'EW': 6, 'NS': 5, 'Z': 7}[direction]
    assert model['c0'] == 0 and model['nst'] == 4
    assert all(n[0] == 2 for n in model['nest'])
    return sum(n[1] * (-np.expm1(-3 * h / n[index])) for n in model['nest'])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inputs = [DATA / 'exponential.json']
    inputs += [DATA / f'experimental_{d}.csv' for d in ['EW', 'NS', 'Z']]
    inputs += [DATA / f'full_{d}_curves.npy' for d in ['EW', 'NS', 'Z']]
    hashes = {p.name: digest(p) for p in inputs}
    cfg = json.loads((DATA / 'exponential.json').read_text())
    model = cfg['variograms']['categorical']
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.55), layout='constrained')
    fig7, axes7 = plt.subplots(1, 3, figsize=(13.3, 4.55), layout='constrained')
    factors = []
    exports = []
    for direction, ax, ax7 in zip(['EW', 'NS', 'Z'], axes, axes7):
        limit = 15 if direction == 'Z' else 155
        step = .5 if direction == 'Z' else 10
        target = pd.read_csv(DATA / f'experimental_{direction}.csv')
        valid = ((target.pairs > 0) & (target.distance > 0) &
                 (target.nominal_distance > 0) & (target.distance <= limit))
        target = target[valid].copy()
        curves = np.load(DATA / f'full_{direction}_curves.npy')
        assert curves.shape[0] == 100
        h = np.arange(1, curves.shape[1] + 1) * step
        displayed = h <= limit
        h, curves = h[displayed], curves[:, displayed]
        mean = curves.mean(axis=0)
        experimental_peak = float(target.gamma.max())
        reproduced_peak = float(mean.max())
        # A unit-total-sill model is computed first, then divided by its own
        # plotted peak: total sill alone does not give a plotted maximum of one.
        total_sill = sum(n[1] for n in model['nest'])
        x = np.linspace(0, max(float(h.max()), float(target.distance.max())), 500)
        input_raw = theoretical(x, model, direction)
        input_unit_sill = input_raw / total_sill
        input_display_peak = float(input_unit_sill.max())
        input_normalized = input_unit_sill / input_display_peak
        experimental_normalized = target.gamma.to_numpy() / experimental_peak
        reproduced_normalized = mean / reproduced_peak
        # One mean-derived denominator is applied to all members, preserving
        # their spread. Individual member peaks are not forced to one.
        normalized_members = curves / reproduced_peak
        lower, upper = np.quantile(normalized_members, [.05, .95], axis=0)
        hp = np.r_[0, h]
        ax.fill_between(hp, np.r_[0, lower], np.r_[0, upper],
                        color='#007f86', alpha=.18, linewidth=0,
                        label='5th-95th percentile', zorder=1)
        ax.plot(x, input_normalized, '--', color='#ce8a1a', lw=1.7,
                label='Theoretical input', zorder=2)
        ax.plot(hp, np.r_[0, reproduced_normalized], color='#007f86', lw=3.2,
                solid_capstyle='round', label='Mean reproduced', zorder=3)
        ax.scatter(target.distance, experimental_normalized, s=25,
                   color='#222222', label='Experimental', zorder=5)
        ax.set(xlim=(0, limit), ylim=(0, 1.2), title=direction,
               xlabel='Separation distance [m]', ylabel='Peak-normalized semivariance')
        ax.grid(alpha=.16)
        ax7.plot(x, input_normalized, '--', color='#ce8a1a', lw=1.7,
                 label='Theoretical input', zorder=2)
        ax7.scatter(target.distance, experimental_normalized, s=25, color='#222222',
                    label='Experimental', zorder=5)
        ax7.set(xlim=(0, limit), ylim=(0, 1.2), title=direction,
                xlabel='Separation distance [m]', ylabel='Peak-normalized semivariance')
        ax7.grid(alpha=.16)
        assert len(ax7.lines) == 1 and len(ax.lines) == 2
        for value in [experimental_normalized.max(), input_normalized.max(), reproduced_normalized.max()]:
            np.testing.assert_allclose(value, 1, rtol=0, atol=1e-14)
        assert hp[0] == 0 and ax.lines[1].get_ydata()[0] == 0
        factors.append(dict(direction=direction, experimental_peak=experimental_peak,
                            input_total_sill=total_sill,
                            input_plotted_peak=float(input_raw.max()),
                            reproduced_mean_peak=reproduced_peak,
                            experimental_N=len(target), reproduced_N=len(h),
                            experimental_last_distance=float(target.distance.max()),
                            reproduced_last_distance=float(h.max()), linewidth=3.2))
        for j, row in target.iterrows():
            exports.append(dict(direction=direction, series='experimental', distance=row.distance,
                                raw_gamma=row.gamma, normalized_gamma=row.gamma / experimental_peak,
                                pairs=int(row.pairs)))
        for distance, raw, value, lo, hi in zip(h, mean, reproduced_normalized, lower, upper):
            exports.append(dict(direction=direction, series='mean_reproduced', distance=distance,
                                raw_gamma=raw, normalized_gamma=value, lower=lo, upper=hi))
        for distance, raw, value in zip(x, input_raw, input_normalized):
            exports.append(dict(direction=direction, series='theoretical_input', distance=distance,
                                raw_gamma=raw, normalized_gamma=value))
    handles, labels = axes[0].get_legend_handles_labels()
    lookup = dict(zip(labels, handles))
    order = ['Experimental', 'Theoretical input', 'Mean reproduced', '5th-95th percentile']
    axes[0].legend([lookup[n] for n in order], order, fontsize=8, loc='lower right')
    fig.suptitle('Figure 8 | Complete binary SISIM, 100 realizations\n'
                 'Each direction and series scaled to its own plotted maximum = 1 (shape only)', fontsize=11)
    handles7, labels7 = axes7[0].get_legend_handles_labels()
    lookup7 = dict(zip(labels7, handles7))
    order7 = ['Experimental', 'Theoretical input']
    axes7[0].legend([lookup7[n] for n in order7], order7, fontsize=8, loc='lower right')
    fig7.suptitle('Figure 7 | Experimental and theoretical input\n'
                  'Each direction and series scaled to its own plotted maximum = 1 (shape only)', fontsize=11)
    for ext in ['png', 'pdf']:
        fig.savefig(OUT / f'08_variogram_reproduction.{ext}', dpi=240)
        fig7.savefig(OUT / f'07_variogram_indicator_fit.{ext}', dpi=240)
    plt.close(fig)
    plt.close(fig7)
    pd.DataFrame(exports).to_csv(DATA / 'variogram_plot_values.csv', index=False)
    assert hashes == {p.name: digest(p) for p in inputs}
    record = dict(scope='Complete binary SISIM simulation grid, not hierarchical class labels or ensemble-mean variogram',
                  normalization='For each direction, experimental and mean reproduced divide by their own plotted maxima; input divides by total sill and then by its plotted maximum.',
                  envelope='100 members share the mean reproduced peak denominator in each direction.',
                  origin='Only gamma(0)=0 identity for reproduced; no experimental zero/near-zero bin plotted.',
                  interpretation='Direction-wise peak-normalized shapes; raw amplitude errors are recorded separately.',
                  input_hashes=hashes, factors=factors, source_files_unchanged=True)
    reference = json.loads((DATA / 'variogram_normalization.json').read_text())['common_reference_variance']
    rmse = {}
    for direction in ['EW', 'NS', 'Z']:
        t = pd.read_csv(DATA / f'experimental_{direction}.csv')
        t = t[(t.pairs > 0) & (t.distance > 0) & (t.nominal_distance > 0)]
        a = np.load(DATA / f'full_{direction}_curves.npy')
        h = np.arange(1, a.shape[1] + 1) * (.5 if direction == 'Z' else 10)
        t = t[t.distance.between(h[0], h[-1])]
        delta = np.interp(t.distance, h, a.mean(axis=0)) - t.gamma
        rmse[direction] = float(np.sqrt(np.mean(delta ** 2)) / reference)
    record['common_reference_variance'] = reference
    record['raw_full_grid_common_reference_RMSE'] = rmse
    (DATA / 'variogram_peak_normalization.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(factors, indent=2))


if __name__ == '__main__':
    main()
