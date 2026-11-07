#!/usr/bin/env python3
"""Open N0 — check the map built by build_n0.py. MIT license.

1. Against the ITU-R P.452 N050 map, if you have your own copy (download it from the ITU and pass its path with
   --itu path/to/N050.TXT). Only difference statistics are printed: the ITU file is not part of this repository and
   must not be added to it (ITU: no reproduction without written permission).
2. Zonal values for February and August, to read against Figs 1 and 2 of ITU-R P.453.
3. Sample points in the ITM climate classes (NTIA Report 82-100, Table of suggested N_s values).
Usage: python3 tools/compare_n0.py [npz] [--itu path/to/N050.TXT]
"""
import os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
args = [a for a in sys.argv[1:] if not a.startswith('--itu')]
ITU = next((sys.argv[i + 1] for i, a in enumerate(sys.argv[:-1]) if a == '--itu'), '')
NPZ = args[0] if args and args[0] != ITU else os.path.join(ROOT, 'data', 'n0_r2_1991-2020.npz')
d = np.load(NPZ); n0, lat, lon = d['n0'], d['lat'], d['lon']


def at(layer, la, lo):   # bilinear on the 1.5 degree grid
    lo %= 360.0
    y = (90.0 - la) / 1.5; x = lo / 1.5
    i, j = int(min(np.floor(y), len(lat) - 2)), int(min(np.floor(x), len(lon) - 2)); fy, fx = y - i, x - j
    a = n0[layer]
    return (a[i, j] * (1 - fx) + a[i, j + 1] * fx) * (1 - fy) + (a[i + 1, j] * (1 - fx) + a[i + 1, j + 1] * fx) * fy


itu = ITU
if itu and os.path.exists(itu):
    ref = np.loadtxt(itu)
    assert ref.shape == n0[0].shape, (ref.shape, n0[0].shape)
    diff = n0[0] - ref
    w = np.cos(np.radians(lat))[:, None] * np.ones_like(diff)          # area weights
    print(f'vs ITU N050 (annual): mean {np.average(diff, weights=w):+.1f}, rms {np.sqrt(np.average(diff**2, weights=w)):.1f}, '
          f'|diff| p50 {np.percentile(abs(diff), 50):.1f} p90 {np.percentile(abs(diff), 90):.1f} max {abs(diff).max():.1f} N-units')
    for lo_, hi_ in [(60, 90), (30, 60), (0, 30), (-30, 0), (-60, -30), (-90, -60)]:
        r = (lat >= lo_) & (lat <= hi_)
        print(f'   lat {lo_:+3d}..{hi_:+3d}: ours {np.average(n0[0][r], weights=w[r]):6.1f}  ITU {np.average(ref[r], weights=w[r]):6.1f}  '
              f'mean diff {np.average(diff[r], weights=w[r]):+5.1f}')
    for name, la, lo_ in [('Messina', 38.19, 15.55), ('Roma', 41.9, 12.5), ('London', 51.5, 0.0), ('Kinshasa', -4.3, 15.3),
                          ('Sahara (In Salah)', 27.2, 2.5), ('Singapore', 1.35, 103.8), ('Tromso', 69.6, 18.9), ('Denver', 39.7, -105.0)]:
        y = (90.0 - la) / 1.5; x = (lo_ % 360) / 1.5; i, j = int(round(y)), int(round(x))
        print(f'   {name:18s} ours {at(0, la, lo_):6.1f}  ITU(nearest node) {ref[i, j]:6.1f}')
else:
    print('no ITU N050 copy given (--itu): comparison skipped')

print('\nzonal medians (P.453 Figs 1-2 are monthly means at sea level):')
for la in (80, 60, 50, 40, 30, 20, 10, 0, -10, -20, -30, -40, -50, -60, -70, -80):
    r = np.argmin(abs(lat - la))
    print(f'   {la:+3d}°  Feb {np.median(n0[2][r]):6.1f}  Aug {np.median(n0[8][r]):6.1f}  annual {np.median(n0[0][r]):6.1f}')

print('\nITM climate classes (suggested N_s, NTIA 82-100) vs this map, annual:')
for name, sug, la, lo_ in [('Equatorial (Congo)', 360, -1.0, 23.0), ('Continental subtropical (Sudan)', 320, 13.0, 30.0),
                           ('Maritime subtropical (W Africa coast)', 370, 9.0, -14.0), ('Desert (Sahara)', 280, 23.0, 10.0),
                           ('Continental temperate (US Midwest)', 301, 41.0, -93.0), ('Maritime temperate over land (UK)', 320, 52.5, -1.5),
                           ('Maritime temperate over sea (North Sea)', 350, 56.0, 3.0)]:
    print(f'   {name:42s} ITM {sug}  map {at(0, la, lo_):6.1f}')
