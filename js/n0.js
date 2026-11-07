// Open N0 — read the sea-level surface refractivity map (data/n0.bin, tools/build_n0.py) in a browser or Node: median of
// the refractivity computed every 6 hours in 1991-2020 from NCEP-DOE Reanalysis 2 (public domain) with ITU-R P.453
// eq. (6), reduced to sea level with the scale height ITM uses (9.46 km). 13 layers (whole year, January..December) of
// 121 x 241 points, 1.5 degrees, rows from +90 down, columns from 0 east; uint16 = N0 x 10, little-endian.
// MIT license.
const ROWS = 121, COLS = 241, STEP = 1.5, LAYER = ROWS * COLS;
let grid = null;

export function loadN0(url) {
  if (!grid) grid = fetch(url).then(r => {
    if (!r.ok) throw new Error(`N0 map: HTTP ${r.status}`);
    return r.arrayBuffer();
  }).then(b => {
    const a = new Uint16Array(b);
    if (a.length !== 13 * LAYER) throw new Error('N0 map: unexpected size');
    return a;
  }).catch(e => { grid = null; throw e; });
  return grid;
}

// N0 at a point, bilinear between the four nearest grid points (as ITU-R P.1144 describes). layer 0 = whole year,
// 1..12 = month.
export async function n0At(url, lat, lon, layer = 0) {
  const g = await loadN0(url), off = layer * LAYER;
  const y = Math.min(Math.max((90 - lat) / STEP, 0), ROWS - 1.000001);
  const x = (((lon % 360) + 360) % 360) / STEP;
  const i = Math.floor(y), j = Math.min(Math.floor(x), COLS - 2), fy = y - i, fx = x - j;
  const v = (r, c) => g[off + r * COLS + c] / 10;
  return (v(i, j) * (1 - fx) + v(i, j + 1) * fx) * (1 - fy) + (v(i + 1, j) * (1 - fx) + v(i + 1, j + 1) * fx) * fy;
}
