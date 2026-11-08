#!/usr/bin/env python3
"""Open N0 — map of the radio-refractivity lapse rate through the lowest 1 km of the atmosphere, DN (the quantity of the
ITU-R P.452 / P.1812 map DN50: "the average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere",
used as a positive number; P.2001 calls it the median value), built from public-domain data. MIT license (https://github.com/IW9GDC/open-n0-map).

What it computes
  Every 6 hours for 30 years (1991-2020): DN = N(surface) - N(surface + 1 km), N-units/km; then its median (the value
  exceeded for 50 % of the time) and its mean per grid cell, on the 1.5 degree layout of the ITU map (latitude +90 to
  -90, longitude 0 to 360).

Physics (see README.md for the references)
  N   = 77.6 (P - e)/T + 72 e/T + 3.75e5 e/T^2               ITU-R P.453-14 eq. (6), the same formula as build_n0.py
  surface: 2 m air temperature, 2 m specific humidity q, surface pressure; e = q P / (0.622 + 0.378 q)
  levels 1000..300 hPa: temperature and relative humidity; e = RH/100 * EF * 6.1121 exp((18.678 - t/234.5) t/(t + 257.14))
            (P.453-14 eqs. (9)-(10), saturation over water, t in deg C, EF = 1 + 1e-4 (7.2 + P (0.0320 + 5.9e-6 t^2)))
  heights of the levels: hypsometric equation upward from the surface, z2 = z1 + (Rd/g) Tv_mean ln(p1/p2),
            Tv = T (1 + 0.608 q), Rd = 287.05 J/(kg K), g = 9.80665 m/s2; only levels above the ground are used
  N(surface + 1 km): ln N interpolated linearly in height between the two points around it (surface or levels)

Data
  NCEP-DOE Reanalysis 2 (Kanamitsu et al. 2002, BAMS 83, 1631-1643), NOAA PSL, public domain: the surface fields of
  tools/fetch_r2.sh (T62 Gaussian grid, interpolated bilinearly to the 2.5 degree grid of the levels) and the pressure-level
  air temperature and relative humidity of fetch_r2_levels.sh. The model orography (hgt.sfc) gives the surface height.

Usage
  tools/fetch_r2.sh raw && tools/fetch_r2_levels.sh raw-levels      (about 5.3 + 15 GB)
  python3 tools/build_dn.py [surface raw dir] [levels raw dir] [out dir]      (defaults: raw, raw-levels, data)
  -> data/dn.bin (uint16 DN x 10, the median, 121 x 241, row = latitude from +90 down, column = longitude from 0 east),
     data/DN_annual.txt (the same grid as text), data/dn.json (layout), data/dn_r2_1991-2020.npz (median, mean and
     standard deviation, also on the 2.5 degree grid of the levels)
  python3 tools/build_dn.py --files     only rewrites dn.bin, DN_annual.txt and dn.json from the existing .npz
"""
import json, os, sys, datetime as dt
import numpy as np
import netCDF4

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
WEB_ONLY = '--files' in sys.argv
if WEB_ONLY:
    sys.argv.remove('--files')
SFC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'raw')
LEV = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'raw-levels')
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, 'data')
YEARS = range(int(os.environ.get('DN_FIRST', 1991)), int(os.environ.get('DN_LAST', 2020)) + 1)   # env: a shorter run to test
NLEV = 8                                        # 1000, 925, 850, 700, 600, 500, 400, 300 hPa
RD, G = 287.05, 9.80665
B0, BW, NB = -150.0, 0.05, 12000                # histogram of DN: -150 .. +450 N-units/km in 0.05 steps


def read(path, name, levels=False):
    n = netCDF4.Dataset(path)
    v = n.variables[name]
    a = v[:, :NLEV] if levels else v[:]
    t = n.variables['time']
    times = np.asarray(t[:], dtype=np.float64)
    units = t.units
    out = dict(v=np.asarray(np.ma.filled(a.astype(np.float32), np.nan)), t=times, units=units,
               lat=n.variables['lat'][:].astype(np.float64), lon=n.variables['lon'][:].astype(np.float64))
    if levels:
        out['lev'] = n.variables['level'][:NLEV].astype(np.float64)
    n.close()
    return out


def refractivity(T, e, P):
    """N-units. T kelvin, e and P hPa."""
    return 77.6 * (P - e) / T + 72.0 * e / T + 3.75e5 * e / (T * T)


def weights(glat, glon, lat, lon):
    """Bilinear weights from the Gaussian grid (glat north to south) to the regular grid (lat, lon): 4 indices + 4 weights
    per target cell, flattened row-major."""
    asc = glat[::-1]
    yi = np.clip(np.searchsorted(asc, lat) - 1, 0, len(asc) - 2)
    fy = np.clip((lat - asc[yi]) / (asc[yi + 1] - asc[yi]), 0, 1)
    rows0, rows1 = len(glat) - 1 - yi, len(glat) - 2 - yi                 # back to north-to-south indices
    dlon = glon[1] - glon[0]
    xi = np.floor((lon - glon[0]) / dlon).astype(int) % len(glon)
    fx = ((lon - glon[0]) / dlon) - np.floor((lon - glon[0]) / dlon)
    xj = (xi + 1) % len(glon)
    nx = len(glon)
    R0, R1 = np.repeat(rows0, len(lon)), np.repeat(rows1, len(lon))
    FY = np.repeat(fy, len(lon))
    X0, X1, FX = np.tile(xi, len(lat)), np.tile(xj, len(lat)), np.tile(fx, len(lat))
    idx = np.stack([R0 * nx + X0, R0 * nx + X1, R1 * nx + X0, R1 * nx + X1])
    w = np.stack([(1 - FY) * (1 - FX), (1 - FY) * FX, FY * (1 - FX), FY * FX])
    return idx, w.astype(np.float32)


def regrid(a, idx, w):
    """a (..., ny*nx) on the Gaussian grid -> (..., cells) on the regular grid."""
    return (a[..., idx[0]] * w[0] + a[..., idx[1]] * w[1] + a[..., idx[2]] * w[2] + a[..., idx[3]] * w[3])


def write_web(npz):
    """dn.bin, DN_annual.txt and dn.json (in OUT) from the median layer of the .npz"""
    dn = np.load(npz)['dn50']
    web = OUT; os.makedirs(web, exist_ok=True)
    np.round(dn * 10).astype('<u2').tofile(os.path.join(web, 'dn.bin'))
    np.savetxt(os.path.join(web, 'DN_annual.txt'), dn, fmt='%.1f')
    meta = {
        'quantity': 'radio-refractivity lapse rate through the lowest 1 km of the atmosphere, DN = N(surface) - N(surface + 1 km), '
                    'N-units/km, median (exceeded for 50 % of the time): the quantity of the ITU-R P.452/P.1812 map DN50',
        'layers': ['annual'], 'rows': int(dn.shape[0]), 'cols': int(dn.shape[1]), 'lat0': 90.0, 'lon0': 0.0, 'step': 1.5, 'scale': 0.1,
        'type': 'uint16 little-endian',
        'source': 'NCEP-DOE Reanalysis 2 (NOAA PSL, public domain): 2 m T and q, surface pressure, and pressure-level T and RH, 6-hourly 1991-2020',
        'method': 'ITU-R P.453 eq. (6) at the surface and at 1 km above it (hypsometric heights, ln N interpolated); median per grid cell; bilinear to 1.5 degrees',
        'built': dt.date.today().isoformat(), 'license': 'MIT (Open N0); input data public domain',
        'files': {'dn.bin': '1 layer uint16 little-endian (DN x 10), 121 rows (lat +90 to -90) x 241 columns (lon 0 to 360)',
                  'DN_annual.txt': 'the same grid as text, one row per latitude, N-units/km with one decimal'},
    }
    json.dump(meta, open(os.path.join(web, 'dn.json'), 'w'), indent=1)
    print('dn.bin, DN_annual.txt, dn.json written, range %.1f..%.1f' % (dn.min(), dn.max()))


def main():
    if WEB_ONLY:
        return write_web(os.path.join(OUT, 'dn_r2_1991-2020.npz'))
    hs = read(os.path.join(SFC, 'hgt.sfc.gauss.nc'), 'hgt')
    first = read(os.path.join(LEV, 'air.1991.nc'), 'air', levels=True)
    lat, lon, lev = first['lat'], first['lon'], first['lev']
    assert np.allclose(lev, [1000, 925, 850, 700, 600, 500, 400, 300]), lev
    del first
    idx, w = weights(hs['lat'], hs['lon'], lat, lon)
    ncell = len(lat) * len(lon)
    zs = regrid(hs['v'][0].reshape(-1), idx, w)                            # model surface height on the 2.5 degree grid
    hist = np.zeros(ncell * NB, np.int32)
    s, s2, cnt = np.zeros(ncell), np.zeros(ncell), np.zeros(ncell)
    for y in YEARS:
        T2 = read(os.path.join(SFC, f'air.2m.gauss.{y}.nc'), 'air')
        q2 = read(os.path.join(SFC, f'shum.2m.gauss.{y}.nc'), 'shum')
        ps = read(os.path.join(SFC, f'pres.sfc.gauss.{y}.nc'), 'pres')
        TL = read(os.path.join(LEV, f'air.{y}.nc'), 'air', levels=True)
        RL = read(os.path.join(LEV, f'rhum.{y}.nc'), 'rhum', levels=True)
        # same instants in all five files (all 6-hourly; hours since 1800-01-01)
        k = min(len(T2['t']), len(TL['t']))
        for d in (q2, ps, RL):
            k = min(k, len(d['t']))
        tt = [netCDF4.num2date(d['t'][:k], d['units'], only_use_cftime_datetimes=False) for d in (T2, q2, ps, TL, RL)]
        assert all((np.array(x) == np.array(tt[0])).all() for x in tt), y
        Ts = regrid(T2['v'][:k].reshape(k, -1), idx, w)
        qs = regrid(q2['v'][:k].reshape(k, -1), idx, w)
        Ps = regrid(ps['v'][:k].reshape(k, -1), idx, w) / 100.0             # hPa
        es = qs * Ps / (0.622 + 0.378 * qs)
        Ns = refractivity(Ts, es, Ps)
        Tvs = Ts * (1 + 0.608 * qs)
        T = TL['v'][:k].reshape(k, NLEV, -1)
        RH = np.clip(RL['v'][:k].reshape(k, NLEV, -1), 0, 100)
        del T2, q2, ps, TL, RL
        z1 = zs[None, :] + 1000.0
        pz, pN, pTv, pp = np.broadcast_to(zs, Ts.shape).copy(), Ns.copy(), Tvs.copy(), Ps.copy()
        N1 = np.full(Ts.shape, np.nan, np.float32)
        for j in range(NLEV):
            p = lev[j]
            t = T[:, j] - 273.15
            EF = 1 + 1e-4 * (7.2 + p * (0.0320 + 5.9e-6 * t * t))
            e = RH[:, j] / 100.0 * EF * 6.1121 * np.exp((18.678 - t / 234.5) * t / (t + 257.14))
            N = refractivity(T[:, j], e, p)
            q = 0.622 * e / (p - 0.378 * e)
            Tv = T[:, j] * (1 + 0.608 * q)
            above = p < pp                                                 # level above the previous point (and the ground)
            z = pz + RD / G * 0.5 * (pTv + Tv) * np.log(pp / p)
            hit = above & np.isnan(N1) & (z >= z1)
            f = np.where(hit, (z1 - pz) / np.where(z > pz, z - pz, 1), 0)
            N1 = np.where(hit, np.exp(np.log(pN) + f * (np.log(N) - np.log(pN))), N1)
            pz, pN, pTv, pp = (np.where(above, z, pz), np.where(above, N, pN), np.where(above, Tv, pTv), np.where(above, p, pp))
        DN = (Ns - N1).astype(np.float32)                                  # N-units per km (1 km layer)
        ok = np.isfinite(DN)
        b = np.clip(((DN - B0) / BW).astype(np.int64), 0, NB - 1)
        cells = np.broadcast_to(np.arange(ncell), DN.shape)
        hist += np.bincount((cells[ok] * NB + b[ok]), minlength=ncell * NB).astype(np.int32)
        s += np.where(ok, DN, 0).sum(0); s2 += np.where(ok, DN * DN, 0).sum(0); cnt += ok.sum(0)
        print(y, k, 'steps, DN median of the year (all cells)', round(float(np.nanmedian(DN)), 2),
              'missing', int((~ok).sum()), flush=True)
    H = hist.reshape(ncell, NB)
    c = np.cumsum(H, axis=1)
    half = c[:, -1] / 2.0
    j = np.argmax(c >= half[:, None], axis=1)                              # first bin reaching half of the samples
    below = np.where(j > 0, c[np.arange(ncell), j - 1], 0)
    frac = (half - below) / np.maximum(H[np.arange(ncell), j], 1)
    med = (B0 + (j + frac) * BW).reshape(len(lat), len(lon)).astype(np.float32)
    mean = (s / np.maximum(cnt, 1)).reshape(len(lat), len(lon)).astype(np.float32)
    sd = np.sqrt(np.maximum(s2 / np.maximum(cnt, 1) - (s / np.maximum(cnt, 1)) ** 2, 0)).reshape(len(lat), len(lon)).astype(np.float32)

    # bilinear regrid to the 1.5 degree ITU layout (lat +90..-90, lon 0..360 inclusive); the 2.5 degree grid already
    # runs from +90 to -90 with both poles
    glat = np.arange(90.0, -90.0 - 1e-9, -1.5)
    glon = np.arange(0.0, 360.0 + 1e-9, 1.5)
    lonw = np.concatenate([lon, [lon[0] + 360.0]])
    asc = lat[::-1]

    def itu(a):
        a = np.concatenate([a, a[:, :1]], axis=1)[::-1]
        tmp = np.array([np.interp(glon, lonw, row) for row in a])
        return np.array([np.interp(glat, asc, tmp[:, i]) for i in range(len(glon))]).T.astype(np.float32)

    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(os.path.join(OUT, 'dn_r2_1991-2020.npz'), dn50=itu(med), dnmean=itu(mean), dnsd=itu(sd), lat=glat, lon=glon,
                        dn50_25=med, dnmean_25=mean, lat25=lat, lon25=lon, zs25=zs.reshape(len(lat), len(lon)))
    print('written: DN50 range', float(med.min()), float(med.max()), 'mean of cells', float(med.mean()))
    write_web(os.path.join(OUT, 'dn_r2_1991-2020.npz'))


if __name__ == '__main__':
    main()
