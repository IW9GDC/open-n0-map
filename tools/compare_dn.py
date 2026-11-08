#!/usr/bin/env python3
"""Open N0 — compare the open DN map (build_dn.py) with the ITU-R P.452/P.1812 map DN50 (private copy, never published:
only the differences are printed). Overall, by latitude band and by surface height, plus a few places.
Usage: python3 tools/compare_dn.py [npz] --itu path/to/your/own/DN50.TXT. MIT license (https://github.com/IW9GDC/open-n0-map)."""
import os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
args = [a for i, a in enumerate(sys.argv[1:]) if a != '--itu' and sys.argv[i] != '--itu']
ITU = next((sys.argv[i + 1] for i, a in enumerate(sys.argv[:-1]) if a == '--itu'), '')
if not ITU or not os.path.exists(ITU):
    sys.exit('the ITU map is not part of this repository: pass your own copy with --itu path/to/DN50.TXT')
f = np.load(args[0] if args else os.path.join(ROOT, 'data', 'dn_r2_1991-2020.npz'))
itu = np.loadtxt(ITU)
lat, lon = f['lat'], f['lon']
LA = np.repeat(lat[:, None], len(lon), 1)
w = np.cos(np.radians(LA))                                     # area weights
zs = np.array([np.interp(lon, np.append(f['lon25'], 360.0), np.append(r, r[0])) for r in f['zs25'][::-1]])
zs = np.array([np.interp(lat, f['lat25'][::-1], zs[:, j]) for j in range(len(lon))]).T
for key in ('dn50', 'dnmean'):
    d = f[key] - itu
    st = lambda m: '%6d  mean %+5.1f  rms %4.1f' % (m.sum(), np.average(d[m], weights=w[m]), np.sqrt(np.average(d[m] ** 2, weights=w[m])))
    print('==', key, '- ITU DN50 (N-units/km); this map ranges %.1f..%.1f' % (f[key].min(), f[key].max()))
    print('  all cells           ', st(np.ones_like(d, bool)))
    for a, b in [(60, 90), (30, 60), (0, 30), (-30, 0), (-60, -30), (-90, -60)]:
        print('  lat %+3d..%+3d         ' % (a, b), st((LA >= a) & (LA < b) if b < 90 else (LA >= a)))
    for a, b in [(-1e9, 500), (500, 1500), (1500, 3000), (3000, 1e9)]:
        print('  surface %5s-%5s m  ' % (max(a, 0) if a > -1e9 else '<', b if b < 1e9 else ''), st((zs >= a) & (zs < b)))
P = {'Messina': (38.19, 15.55), 'London': (51.5, -0.1), 'Sahara (Tamanrasset)': (22.8, 5.5), 'Singapore': (1.3, 103.8), 'Denver': (39.7, -105.0), 'Lhasa': (29.65, 91.1)}
def at(a, la, lo):
    lo = lo % 360; i, j = (90 - la) / 1.5, lo / 1.5; i0, j0 = int(i), int(j); fi, fj = i - i0, j - j0
    return (a[i0, j0] * (1 - fi) * (1 - fj) + a[i0, j0 + 1] * (1 - fi) * fj + a[i0 + 1, j0] * fi * (1 - fj) + a[i0 + 1, j0 + 1] * fi * fj)
for n, (la, lo) in P.items():
    print('  %-22s ours median %5.1f  mean %5.1f   minus ITU %+5.1f' % (n, at(f['dn50'], la, lo), at(f['dnmean'], la, lo), at(f['dn50'], la, lo) - at(itu, la, lo)))
