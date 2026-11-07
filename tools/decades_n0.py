#!/usr/bin/env python3
"""Open N0 — change of the sea-level refractivity N0 between 1991-2000 and 2011-2020 (NCEP-DOE Reanalysis 2), with the
part due to the dry term (pressure, temperature) and to the wet term (water vapour). Same physics as build_n0.py.
MIT license. Usage: python3 tools/decades_n0.py [raw dir]   (default raw/; writes data/n0_decades.npz)"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from build_n0 import field, H0
RAW = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'raw')
hz, _, lat, lon = field(os.path.join(RAW, 'hgt.sfc.gauss.nc'), 'hgt'); red = np.exp(np.maximum(hz[0], 0) / H0)
res = {}
for name, years in (('1991-2000', range(1991, 2001)), ('2011-2020', range(2011, 2021))):
    n0s, dry, wet, tt, ee = [], 0, 0, 0, 0; cnt = 0
    for y in years:
        T, *_ = field(os.path.join(RAW, f'air.2m.gauss.{y}.nc'), 'air')
        q, *_ = field(os.path.join(RAW, f'shum.2m.gauss.{y}.nc'), 'shum')
        P, *_ = field(os.path.join(RAW, f'pres.sfc.gauss.{y}.nc'), 'pres')
        k = min(len(T), len(q), len(P)); T, q, P = T[:k], q[:k], P[:k] / 100.0
        e = q * P / (0.622 + 0.378 * q)
        d = 77.6 * (P - e) / T * red; w = (72.0 * e / T + 3.75e5 * e / (T * T)) * red
        n0s.append(np.round((d + w) * 10).astype(np.uint16))
        dry = dry + d.sum(0); wet = wet + w.sum(0); tt = tt + T.sum(0); ee = ee + e.sum(0); cnt += k
        print(name, y, flush=True)
    allv = np.concatenate(n0s)
    med = np.empty(allv.shape[1:], np.float32)
    for i in range(allv.shape[1]): med[i] = np.median(allv[:, i, :], axis=0) / 10.0
    res[name] = dict(med=med, dry=dry / cnt, wet=wet / cnt, T=tt / cnt, e=ee / cnt)
a, b = res['1991-2000'], res['2011-2020']
w = np.cos(np.radians(lat))[:, None] * np.ones_like(a['med'])
avg = lambda x, m=None: float(np.average(x if m is None else x[m], weights=w if m is None else w[m]))
dN, dD, dW, dT, dE = b['med'] - a['med'], b['dry'] - a['dry'], b['wet'] - a['wet'], b['T'] - a['T'], b['e'] - a['e']
print(f'\nGLOBAL (area-weighted): median N0 {avg(dN):+.2f} N-units | mean dry term {avg(dD):+.2f} | mean wet term {avg(dW):+.2f} | '
      f'temperature {avg(dT):+.2f} K | vapour pressure {avg(dE):+.2f} hPa ({100 * avg(dE) / avg(a["e"]):+.1f} %)')
print('\nband          dN0(median)  dry    wet    dT(K)  de(hPa)')
for lo_, hi_ in [(60, 90), (30, 60), (0, 30), (-30, 0), (-60, -30), (-90, -60)]:
    m = (lat >= lo_) & (lat <= hi_); m2 = np.broadcast_to(m[:, None], dN.shape)
    print(f'{lo_:+3d}..{hi_:+3d}    {avg(dN, m2):+6.2f}     {avg(dD, m2):+5.2f}  {avg(dW, m2):+5.2f}  {avg(dT, m2):+5.2f}  {avg(dE, m2):+5.2f}')
def at(x, la, lo):
    i = np.argmin(abs(lat - la)); j = np.argmin(abs(lon - lo % 360)); return float(x[i, j])
print('\nplace             1991-2000  2011-2020  change   dT(K)')
for name, la, lo in [('Messina', 38.13, 15.56), ('London', 51.5, -0.1), ('Kinshasa', -4.3, 15.3), ('Singapore', 1.35, 103.8),
                     ('In Salah (Sahara)', 27.2, 2.5), ('Tromso', 69.6, 18.9), ('Denver', 39.7, -105.0), ('Delhi', 28.6, 77.2)]:
    print(f'{name:18s} {at(a["med"], la, lo):8.1f}  {at(b["med"], la, lo):8.1f}  {at(dN, la, lo):+6.1f}  {at(dT, la, lo):+5.2f}')
p = np.percentile(dN, [5, 50, 95]); print(f'\ngrid cells: 5th / 50th / 95th percentile of the change {p[0]:+.1f} / {p[1]:+.1f} / {p[2]:+.1f}; '
      f'cells increased {100 * (dN > 0).mean():.0f} %')
np.savez_compressed(os.path.join(os.path.dirname(__file__), '..', 'data', 'n0_decades.npz'), lat=lat, lon=lon, **{f'{k}_{n}': v for n, r in (('a', a), ('b', b)) for k, v in r.items()})
