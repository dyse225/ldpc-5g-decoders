"""Validate and combine Violeta's saved MATLAB caches; no simulation.

Run with Python, NumPy and SciPy. Original submissions are never overwritten.
"""
from pathlib import Path
import csv
import math

import numpy as np
from scipy.io import loadmat, savemat

HERE = Path(__file__).resolve().parent
OUT = HERE / 'analysis'


def vector(value):
    return np.atleast_1d(value).astype(float)


def threshold(res, metric, target):
    x, y, n = map(vector, (res['EbN0dB'], res[metric], res['numFrames']))
    for i in range(len(x) - 1):
        if (n[i:i+2] > 0).all() and (y[i:i+2] > 0).all() and y[i] >= target >= y[i+1]:
            if y[i] == y[i+1]:
                return float(x[i])
            f = (math.log10(target) - math.log10(y[i])) / (math.log10(y[i+1]) - math.log10(y[i]))
            return float(x[i] + f * (x[i+1] - x[i]))
    return float('nan')


def write_csv(name, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    runs, raw, summary, inventory = [], [], [], []
    low_count = zero_count = 0
    keys = set()
    for path in sorted((HERE / 'results' / 'cache').glob('*.mat')):
        item = loadmat(path, simplify_cells=True)['item']
        meta, res = item['meta'], item['res']
        key = (int(meta['bgn']), meta['algorithm'], int(meta['maxIter']))
        assert key not in keys, f'Duplicate configuration: {path}'
        keys.add(key)
        runs.append(item)
        n, be, fe = map(vector, (res['numFrames'], res['numBitErrors'], res['numFrameErrors']))
        valid = n > 0
        assert valid.any(), path
        np.testing.assert_allclose(vector(res['fer'])[valid], fe[valid] / n[valid], rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(vector(res['ber'])[valid], be[valid] / (float(res['K']) * n[valid]), rtol=1e-12, atol=1e-15)
        assert (fe[valid] <= n[valid]).all() and (be[valid] <= n[valid] * res['K']).all()
        assert res['params']['modulation'] == 'qpsk'
        inventory.append(dict(BG=key[0], algorithm=key[1], maxIter=key[2],
                              points=int(valid.sum()), skipped=int((~valid).sum()), source=path.relative_to(HERE).as_posix()))
        for i in np.flatnonzero(valid):
            row = dict(BG=key[0], K=int(res['K']), E=int(res['E']), algorithm=key[1], maxIter=key[2],
                       alpha=float(meta['scaling']), beta=float(meta['offset']),
                       EbN0dB=float(vector(res['EbN0dB'])[i]), BER=float(vector(res['ber'])[i]),
                       FER=float(vector(res['fer'])[i]), avgIter=float(vector(res['avgIter'])[i]),
                       numFrames=int(n[i]), numBitErrors=int(be[i]), numFrameErrors=int(fe[i]),
                       FER_upper95_zero=float(1 - .05 ** (1 / n[i])) if fe[i] == 0 else float('nan'))
            raw.append(row)
            low_count += int(0 < fe[i] < 100)
            zero_count += int(fe[i] == 0)
        if key[1] != 'uncoded':
            summary.append(dict(BG=key[0], K=int(res['K']), algorithm=key[1], maxIter=key[2],
                                alpha=float(meta['scaling']), beta=float(meta['offset']),
                                EbN0dB_at_FER1e_2=threshold(res, 'fer', 1e-2),
                                EbN0dB_at_FER1e_3=threshold(res, 'fer', 1e-3),
                                EbN0dB_at_BER1e_4=threshold(res, 'ber', 1e-4)))
    for bg in (1, 2):
        for alg in ('bp', 'offset-min-sum'):
            assert (bg, alg, 50) in keys
        for alg in ('layered-bp', 'norm-min-sum'):
            for it in (5, 10, 20, 50):
                assert (bg, alg, it) in keys
    cells = np.empty((1, len(runs)), dtype=object)
    for i, run in enumerate(runs):
        cells[0, i] = run
    savemat(OUT / 'combined_results.mat', {'runs': cells, 'source': 'Violeta: results/cache; no simulation'}, do_compression=True)
    restored = loadmat(OUT / 'combined_results.mat', simplify_cells=True)['runs']
    assert len(restored) == len(runs)
    write_csv('raw_results_both_bg.csv', raw)
    write_csv('summary_both_bg.csv', summary)
    write_csv('source_inventory.csv', inventory)
    print(f'Validated {len(runs)} runs, {len(raw)} measured points; {low_count} points with 1..99 frame errors; {zero_count} zero-error points.')
    for row in summary:
        if row['maxIter'] == 50:
            print(row)


if __name__ == '__main__':
    main()
