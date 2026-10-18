"""Publish the checked current article figures at the root pics/ paths.

Run after the plotting commands. --check verifies without overwriting.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
STEMS = ['06_sgsim_variogram_reproduction', '07_variogram_indicator_fit', '08_variogram_reproduction', '10_sisim_sections', '11_hierarchical_combination', '12_uq_plan', '13_sections_entropy', '14_sections_variability', '16_realization_stability', 'figure9_final_calibration', 'figure9_sgsim_calibration', 'figure9_sisim_calibration']


def publish(check=False):
    records = {}
    for stem in STEMS:
        for ext in ['png', 'pdf']:
            name = f'{stem}.{ext}'
            source = ROOT / 'pics/final_20261018' / name
            target = ROOT / 'pics' / name
            if check:
                assert source.read_bytes() == target.read_bytes(), name
            else:
                shutil.copy2(source, target)
            records[name] = hashlib.sha256(source.read_bytes()).hexdigest()
    record = dict(version='20261018', candidate='response_step_50',
                  config='configs/exponential_sisim_final_20261018.json',
                  canonical_outputs='pics/final_20261018',
                  data='pics/final_20261018/data',
                  model_provenance='pics/final_20261018/data/final_model_provenance.json',
                  sgsim_reproduction_provenance='pics/final_20261018/data/sgsim_reproduction_provenance.json',
                  vtk='vtk/model_gok_sisim_final_20261018.vtk.gz',
                  figure_sha256=records)
    path = ROOT / 'pics/current_article.json'
    if check:
        assert json.loads(path.read_text(encoding='utf-8')) == record
    else:
        path.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', action='store_true')
    publish(p.parse_args().check)
    print('PASS: 12 current article figure pairs')
