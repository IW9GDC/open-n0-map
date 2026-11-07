#!/usr/bin/env python3
"""Open N0 — sea-level surface refractivity map, built from public-domain data. MIT license (https://github.com/IW9GDC/open-n0-map).

What it computes
  The radio refractivity of the air at the surface, every 6 hours for 30 years (1991-2020), reduced to sea level, and
  its median (the value exceeded for 50 % of the time) for the average year and for each month, on a 1.5 degree grid
  (the layout of the ITU-R P.452 map: latitude +90 to -90, longitude 0 to 360).

Physics (see README.md for the references)
  N   = 77.6 (P - e)/T + 72 e/T + 3.75e5 e/T^2              ITU-R P.453-14 eq. (6); Smith & Weintraub 1953, Thayer 1974
  e   = q P / (0.622 + 0.378 q)                              water vapour pressure from specific humidity q (kg/kg)
  N0  = N exp(z / 9460 m)                                    exponential atmosphere, Bean & Thayer 1959 / P.453 eq. (12);
                                                             9.46 km is the scale height ITM itself uses to bring N0
                                                             back to the height of the path (TN101 eq. 4.3), so a value
                                                             read from this map goes through the model unchanged
  P hPa, e hPa, T kelvin, z = height of the surface above sea level.

Data
  NCEP-DOE Reanalysis 2 (Kanamitsu et al. 2002, BAMS 83, 1631-1643), NOAA PSL, public domain: 2 m air temperature,
  2 m specific humidity, surface pressure, 6-hourly, T62 Gaussian grid; model surface height. tools/fetch_r2.sh

Usage
  tools/fetch_r2.sh raw          (about 5.3 GB)
  python3 tools/build_n0.py [raw dir] [out dir]          (defaults: raw, data)
  -> data/n0.bin (uint16 N0*10, 13 layers: annual, Jan..Dec, each 121 x 241, row = latitude from +90 down, column =
     longitude from 0 east), data/N0_annual.txt and N0_01..12.txt (the same grids as text), data/lat.txt, data/lon.txt,
     data/n0.json (layout), data/n0_r2_1991-2020.npz (all fields, also on the native Gaussian grid).
"""
import json, os, sys, glob, datetime as dt
import numpy as np
import netCDF4

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RAW = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'raw')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'data')
YEARS = range(1991, 2021)
H0 = 9460.0                                    # m, ITM's scale height


def field(path, name):
    n = netCDF4.Dataset(path)
    v = n.variables[name][:]                   # netCDF4 applies scale_factor / add_offset and masks fill values
    t = n.variables['time']
    times = netCDF4.num2date(t[:], t.units, only_use_cftime_datetimes=False)
    lat, lon = n.variables['lat'][:].astype(np.float64), n.variables['lon'][:].astype(np.float64)
    n.close()
    v = np.asarray(np.ma.filled(v.astype(np.float32), np.nan)).reshape(len(times), len(lat), len(lon))
    return v, np.array([d.month for d in times]), lat, lon


def refractivity(T, q, P):
    """N at the surface, N-units. T kelvin, q kg/kg, P pascal."""
    P = P / 100.0
    e = q * P / (0.622 + 0.378 * q)
    return 77.6 * (P - e) / T + 72.0 * e / T + 3.75e5 * e / (T * T)


def main():
    hz, _, lat, lon = field(os.path.join(RAW, 'hgt.sfc.gauss.nc'), 'hgt')
    red = np.exp(np.maximum(hz[0], 0.0) / H0).astype(np.float32)     # sea and depressions: no reduction below 0 m
    ny, nx = len(lat), len(lon)
    per_year, months = [], []
    for y in YEARS:
        T, m, la, lo = field(os.path.join(RAW, f'air.2m.gauss.{y}.nc'), 'air')
        q, m2, _, _ = field(os.path.join(RAW, f'shum.2m.gauss.{y}.nc'), 'shum')
        P, m3, _, _ = field(os.path.join(RAW, f'pres.sfc.gauss.{y}.nc'), 'pres')
        k = min(len(m), len(m2), len(m3))
        assert (m[:k] == m2[:k]).all() and (m[:k] == m3[:k]).all(), y
        assert np.allclose(la, lat) and np.allclose(lo, lon), y
        N0 = refractivity(T[:k], q[:k], P[:k]) * red
        per_year.append(np.round(N0 * 10).astype(np.uint16).reshape(k, ny * nx))   # 0.1 N-unit steps
        months.append(m[:k])
        print(y, k, 'steps, N0 median', float(np.nanmedian(N0)), flush=True)
    allm = np.concatenate(months)
    layers = []
    for sel in [None] + list(range(1, 13)):
        med = np.empty(ny * nx, np.float32)
        for c0 in range(0, ny * nx, 2048):                                 # chunks of cells: bounded memory
            block = np.concatenate([a[(mm == sel) if sel else slice(None), c0:c0 + 2048] for a, mm in zip(per_year, months)])
            med[c0:c0 + 2048] = np.median(block, axis=0) / 10.0
        layers.append(med.reshape(ny, nx))
        print('layer', 'annual' if sel is None else sel, 'done', flush=True)
    layers = np.stack(layers)                                              # (13, 94, 192) on the Gaussian grid

    # bilinear regrid to the 1.5 degree layout of ITU-R P.452 (lat +90..-90, lon 0..360 inclusive)
    glat = np.arange(90.0, -90.0 - 1e-9, -1.5)
    glon = np.arange(0.0, 360.0 + 1e-9, 1.5)
    asc = lat[::-1]                                                        # Gaussian latitudes run north to south
    out = np.empty((13, len(glat), len(glon)), np.float32)
    lonw = np.concatenate([lon, [lon[0] + 360.0]])
    for i in range(13):
        a = layers[i][::-1]                                                # south to north
        a = np.concatenate([a, a[:, :1]], axis=1)                          # wrap in longitude
        tmp = np.array([np.interp(glon, lonw, row) for row in a])          # along longitude
        out[i] = np.array([np.interp(glat, asc, tmp[:, j]) for j in range(len(glon))]).T   # along latitude (edge values at the poles)
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(os.path.join(OUT, 'n0_r2_1991-2020.npz'), n0=out, lat=glat, lon=glon, gauss=layers, glat=lat, glon=lon)
    web = OUT
    np.round(out * 10).astype('<u2').tofile(os.path.join(web, 'n0.bin'))
    for i, nm in enumerate(['annual'] + [f'{m:02d}' for m in range(1, 13)]):
        np.savetxt(os.path.join(web, f'N0_{nm}.txt'), out[i], fmt='%.1f')
    np.savetxt(os.path.join(web, 'lat.txt'), glat, fmt='%.1f'); np.savetxt(os.path.join(web, 'lon.txt'), glon, fmt='%.1f')
    meta = {
        'quantity': 'sea-level surface refractivity N0, N-units, median (exceeded for 50 % of the time)',
        'layers': ['annual'] + [dt.date(2000, i, 1).strftime('%b') for i in range(1, 13)],
        'rows': len(glat), 'cols': len(glon), 'lat0': 90.0, 'lon0': 0.0, 'step': 1.5, 'scale': 0.1, 'type': 'uint16 little-endian',
        'source': 'NCEP-DOE Reanalysis 2 (NOAA PSL, public domain), 2 m T, 2 m q, surface pressure, 6-hourly 1991-2020',
        'method': 'ITU-R P.453 eq. (6); N0 = Ns exp(h / 9460 m); median per grid cell; bilinear to 1.5 degrees',
        'built': dt.date.today().isoformat(), 'license': 'MIT; input data public domain',
    }
    json.dump(meta, open(os.path.join(web, 'n0.json'), 'w'), indent=1)
    print('written', out.shape, 'annual range', float(out[0].min()), float(out[0].max()))


if __name__ == '__main__':
    main()
