"""Check the bundled calibrated-input evidence without external simulation files.

Checks cover the packaged tables, spatial arrays, figures and provenance.
Stored source hashes identify the external production records.
Run the original verify_article_evidence.py as well for the September baseline.
"""
from pathlib import Path
import hashlib
import gzip
import json
import numpy as np
import pandas as pd


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(root=None):
    root = Path(root) if root else Path(__file__).resolve().parents[1]
    data = root / 'pics/final_20261018/data'
    docs = root / 'docs/final_20261018'
    read = lambda p: json.loads(p.read_text(encoding='utf-8'))
    cfg = read(root / 'configs/exponential_sisim_final_20261018.json')
    assert cfg == read(data / 'exponential.json')
    original = read(root / 'configs/exponential_20260916.json')
    assert cfg['variograms']['continuous'] == original['variograms']['continuous']
    selected = read(docs / 'selected_variogram.json')['response_step_50']
    assert selected == cfg['variograms']['categorical']
    assert selected['nst'] == 4 and selected['c0'] == 0
    np.testing.assert_allclose([n[1] for n in selected['nest']], [.05, .04, .04, .12])
    assert all(n[0] == 2 and min(n[5:]) > 0 for n in selected['nest'])
    summary = read(data / 'final_cv_summary.json')
    assert summary['status'] == 'complete'
    assert summary['selected_variogram_sha256'] == digest(docs / 'selected_variogram.json')
    # The distributed snapshot includes a path-normalized copy of the native protocol.
    if 'public_protocol_sha256' in summary:
        assert summary['public_protocol_sha256'] == digest(docs / 'protocol.json')
        assert summary['native_local_protocol_sha256'] == read(docs / 'protocol.json')['source_record_sha256']

    pred = pd.read_csv(data / 'binary_sisim_cv_predictions.csv', dtype={'borehole_id': str})
    assert len(pred) == summary['binary_N'] == 22313
    assert pred.observed.value_counts().to_dict() == {2: 12876, 1: 9437}
    assert pred.realizations.eq(100).all() and pred.sand_count.between(0, 100).all()
    assert pred.sand_count.eq(pred.sand_count.astype(int)).all()
    np.testing.assert_allclose(pred.p_sand, pred.sand_count / 100, atol=1e-15)
    np.testing.assert_array_equal(pred.is_sand, (pred.observed == 2).astype(int))
    np.testing.assert_array_equal(pred.predicted, np.where(pred.p_sand > .5, 2, 1))
    assert pred.borehole_id.nunique() == summary['binary_boreholes'] == 694
    assert pred.group.nunique() == summary['binary_groups'] == 623
    assert pred.groupby('borehole_id').fold.nunique().max() == 1
    assert pred.groupby('group').fold.nunique().max() == 1

    table = pd.read_csv(data / 'binary_sisim_calibration_2pct.csv')
    edges = np.linspace(0, 1, 51)
    assert len(table) == 50
    np.testing.assert_allclose(table.lower, edges[:-1], atol=1e-15)
    np.testing.assert_allclose(table.upper, edges[1:], atol=1e-15)
    bins = np.searchsorted(edges[1:-1], pred.p_sand.to_numpy())
    for b, row in table.iterrows():
        sample = pred[bins == b]
        assert len(sample) == row['count'] > 0
        assert sample.is_sand.sum() == row['sand_count']
        np.testing.assert_allclose([sample.p_sand.mean(), sample.is_sand.mean()],
                                   [row.mean_probability, row.observed_frequency], atol=1e-14)
    assert table['count'].sum() == len(pred)
    gap = abs(table.mean_probability - table.observed_frequency)
    actual = {
        'binary_accuracy': float((pred.predicted == pred.observed).mean()),
        'binary_brier': float(((pred.p_sand - pred.is_sand) ** 2).mean()),
        'binary_ECE_2pct': float((table['count'] * gap).sum() / len(pred)),
        'binary_unweighted_bin_MAE': float(gap.mean()),
    }
    for key, value in actual.items():
        np.testing.assert_allclose(value, summary[key], rtol=0, atol=1e-12)

    folds = pd.read_csv(data / 'final_cv_fold_metrics.csv')
    assert folds.fold.tolist() == [1, 2, 3, 4, 5]
    assert folds.realizations.eq(100).all()
    assert folds.test_N.tolist() == [6115, 5401, 6999, 5412, 6034]
    assert folds.binary_N.tolist() == pred.groupby('fold').size().tolist()
    assert folds.test_N.sum() == summary['total_heldout_N'] == 29961
    assert summary['realizations_per_fold'] == [100] * 5
    batches = pd.read_csv(docs / 'native_batches_verified.csv')
    assert len(batches) == summary['actual_native_batches'] == 100
    assert batches.realizations.eq(5).all() and batches.seconds.gt(0).all()
    for fold in range(1, 6):
        assert sorted(batches.loc[batches.fold == fold, 'part']) == list(range(20))
    assert batches.realizations.sum() == summary['actual_new_sisim_realizations'] == 500
    assert batches.binary_sha256.str.fullmatch(r'[0-9a-f]{64}').all()
    audit = read(docs / 'holdout_audit.json')
    assert audit['whole_borehole_ids_exclusive'] and audit['fold_groups_exclusive']
    assert all(r['same_train_test_voxels'] == 0 and
               r['native_training_data_match_saved_holdout_partition'] for r in audit['rows'])

    continuous = pd.read_csv(data / 'continuous_heldout_predictions.csv')
    assert len(continuous) == summary['continuous_N'] == 664
    assert continuous.fold.nunique() == 5 and continuous.ensemble_size.eq(100).all()
    assert continuous.pit.between(0, 1).all() and summary['sgsim_rerun'] is False
    qt = pd.read_csv(data / 'continuous_quantile_calibration.csv')
    np.testing.assert_allclose(abs(qt.nominal - qt.observed).mean(),
                               summary['continuous_saved_quantile_MAE'], atol=1e-14)

    norm = read(data / 'variogram_normalization.json')
    rmse = {}
    for direction in ['EW', 'NS', 'Z']:
        exp = pd.read_csv(data / f'experimental_{direction}.csv')
        exp = exp[(exp.pairs > 0) & (exp.distance > 0) & (exp.nominal_distance > 0)]
        step = .5 if direction == 'Z' else 10
        count = 30 if direction == 'Z' else 16
        assert len(exp) == (30 if direction == 'Z' else 15)
        for scope in ['neogene', 'full']:
            curves = np.load(data / f'{scope}_{direction}_curves.npy')
            assert curves.shape == (100, count)
            assert np.isfinite(curves).all() and (curves >= 0).all()
            if scope == 'full':
                h = np.arange(1, count + 1) * step
                valid = exp.distance.between(h[0], h[-1])
                target = exp[valid]
                residual = np.interp(target.distance, h, curves.mean(axis=0)) - target.gamma
                rmse[direction] = float(np.sqrt(np.mean(residual ** 2)) / norm['common_reference_variance'])
    provenance = read(data / 'final_model_provenance.json')
    record = read(root / 'docs/core_source_localization.json')['files']['scripts/UTILITIES_Wang.py']
    assert provenance['numerical_run_core_sha256'] == record['original_sha256']
    assert provenance['repository_release_core_sha256'] == record['release_sha256']
    assert digest(root / 'scripts/UTILITIES_Wang.py') == record['release_sha256']
    assert provenance['realizations'] == 100 and provenance['mask_count'] == 857561
    with np.load(data / 'final_model_sections.npz') as section:
        for axis in ['x', 'y', 'plan']:
            p = sum(section[f'p_{kind}_{axis}'] for kind in ['clay', 'sand', 'gravel'])
            mask = section[f'mask_{axis}'].astype(bool)
            np.testing.assert_allclose(p[mask], 1, atol=1e-12)
    # The current article aliases must match the canonical figure outputs.
    current = read(root / 'pics/current_article.json')
    assert current['candidate'] == 'response_step_50'
    assert current['config'] == 'configs/exponential_sisim_final_20261018.json'
    assert len(current['figure_sha256']) == 24
    for name, expected in current['figure_sha256'].items():
        assert digest(root / 'pics' / name) == digest(root / 'pics/final_20261018' / name) == expected
    for folder in [root / 'pics', root / 'pics/final_20261018']:
        assert {p.name for p in folder.glob('08*.png')} == {'08_variogram_reproduction.png'}
        assert {p.name for p in folder.glob('08*.pdf')} == {'08_variogram_reproduction.pdf'}
    peak_record = read(data / 'variogram_peak_normalization.json')
    for name, expected in peak_record['input_hashes'].items():
        assert digest(data / name) == expected
    for direction, value in rmse.items():
        np.testing.assert_allclose(value, peak_record['raw_full_grid_common_reference_RMSE'][direction], atol=1e-14)
    plot_values = pd.read_csv(data / 'variogram_plot_values.csv')
    np.testing.assert_allclose(plot_values.groupby(['direction', 'series']).normalized_gamma.max(), 1, atol=1e-14)
    sg = read(data / 'sgsim_reproduction_provenance.json')
    assert sg['realizations'] == 100 and sg['grid_shape_yx'] == [100, 100]
    assert sg['model'] == cfg['variograms']['continuous']
    assert sg['source_hashes']['03/sgsim_gaussian.npy'] == '209c9316e67fd079e2f12d23e793e7941ba061aafab6492607ea6dde676cd5a8'
    assert sg['rerun_sgsim'] is False and sg['rerun_cv'] is False
    assert sg['source_files_unchanged'] and sg['direct_pair_max_abs_error'] < 2e-12
    sill = sg['model']['c0'] + sum(n[1] for n in sg['model']['nest'])
    np.testing.assert_allclose(sill, sg['common_input_sill'], atol=1e-14)
    moments = pd.read_csv(data / 'sgsim_realization_moments.csv')
    assert len(moments) == 100 and moments.variance.gt(0).all()
    table = pd.read_csv(data / 'sgsim_reproduction_plot_values.csv')
    np.testing.assert_allclose(table.normalized_gamma, table.raw_gamma/sill, atol=1e-14)
    for name, expected in sg['plot_input_hashes'].items():
        assert digest(root / name) == expected
    for direction in ['EW', 'NS']:
        curves = np.load(data / f'sgsim_{direction}_curves.npy')
        assert curves.shape == (100, 50) and np.isfinite(curves).all() and (curves >= 0).all()
        rows = table[(table.direction == direction) & (table.series == 'mean_reproduced')]
        np.testing.assert_allclose(rows.raw_gamma, curves[:, :15].mean(axis=0), atol=1e-14)
        for column, q in [('lower', .05), ('upper', .95)]:
            np.testing.assert_allclose(rows[column], np.quantile(curves[:, :15], q, axis=0)/sill, atol=1e-14)
    vtk_hash = hashlib.sha256()
    with gzip.open(root / current['vtk'], 'rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            vtk_hash.update(chunk)
    assert vtk_hash.hexdigest() == provenance['vtk_sha256']
    return dict(binary_N=len(pred), native_realizations=500, continuous_N=len(continuous),
                metrics=actual, full_grid_common_reference_spatial_RMSE=rmse)


if __name__ == '__main__':
    result = verify()
    print('PASS: calibrated-input tables, spatial arrays, figures and provenance')
    print(json.dumps(result, indent=2))
