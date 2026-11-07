# Open N0: a free world map of sea-level surface refractivity

A world map of **N0**, the radio refractivity of the air at the surface reduced to sea level, for the whole year and
for each month, built only from public-domain data with the physics of the open literature, and released under the
**MIT licence**: data and code.

| | |
|---|---|
| Quantity | sea-level surface refractivity N0, N-units, **median** (exceeded for 50 % of the time) |
| Period | 1991-2020, every 6 hours (about 43 800 samples per point) |
| Layers | whole year, January ... December (13) |
| Grid | 1.5°, latitude +90 to −90 (121 rows), longitude 0 to 360 (241 columns): the layout of the ITU-R P.452 maps |
| Source data | NCEP-DOE Reanalysis 2, NOAA PSL (public domain) |
| Licence | MIT |

It is used by [RFscape](https://rfscape.com/), a free radio coverage planner, to suggest the value of this parameter
for the place and the month of a link or a coverage map.

## What it is for

Radio propagation models need the refractivity of the air near the ground:

- **ITM / Longley-Rice** takes it as an input (`N_0`) and uses it for the effective earth radius and for the
  diffraction and troposcatter loss. Most programs leave it at **301**, the value for a continental temperate climate
  in the United States.
- **ITU-R P.452 and P.1812** use N0 in their troposcatter models.
- The **effective earth factor k** and the **radio horizon** follow from it: from a 1000 m site the radio horizon is
  130 km with N = 301 and 150 km with N = 400.

The real value is far from constant. It depends above all on the water vapour: about 380-390 at the equator, 280-300
over deserts, 310-370 at mid latitudes (higher in summer and near the sea), about 310-320 near the poles. Around the
Mediterranean it goes from about 325-335 in winter to 350-375 in summer. With this map a program can use the right
value for the place and the month instead of a fixed one.

## Why this exists

The ITU publishes an equivalent map, `N050.TXT`, as an *integral digital product* of Recommendation
[ITU-R P.452](https://www.itu.int/rec/R-REC-P.452) (annual, 1.5°, derived from radiosonde ascents of 1983-1992).
Every ITU Recommendation, with its digital products, carries this notice:

> © ITU. All rights reserved. No part of this publication may be reproduced, by any means whatsoever, without written
> permission of ITU.

So open-source software, websites and browser applications cannot ship that map, and a web page cannot read it from
the ITU site either. The physics behind it, however, is public and decades old, and the meteorological data needed to
apply it are in the public domain. This repository rebuilds an equivalent map from those, and publishes data and code
under the MIT licence, so anyone can use, modify and redistribute it. The ITU map is not included here and was used
only privately, to compare statistics (see [Validation](#validation)).

## Physics

1. **Radio refractivity** of moist air (ITU-R P.453-14, eq. 6; constants after Smith & Weintraub 1953 and Thayer 1974):

   N = 77.6 (P − e) / T + 72 e / T + 3.75·10⁵ e / T²

   P total pressure (hPa), e water vapour pressure (hPa), T temperature (K).

2. **Water vapour pressure** from the specific humidity q (kg/kg): e = q P / (0.622 + 0.378 q).

3. **Reduction to sea level** with the exponential reference atmosphere (Bean & Thayer 1959; ITU-R P.453 eq. 12):

   N0 = Ns · exp(h / h0), with h the height of the ground and **h0 = 9.46 km**.

   9.46 km is the scale height that ITM itself uses to bring N0 back to the mean height of a path (NBS Technical
   Note 101, eq. 4.3; `InitializePointToPoint.cpp` in the NTIA code). With the same value on both sides, a number read
   from this map goes through ITM and gives back the refractivity of the ground of that path.

4. **Statistics**: for every grid point, the median of all the 6-hourly values of the period (whole year, or the
   values of one calendar month), then bilinear interpolation from the reanalysis grid to the 1.5° grid. Between grid
   points, interpolate bilinearly (as ITU-R P.1144 describes).

## Data

NCEP-DOE Reanalysis 2 (Kanamitsu et al. 2002): 2 m air temperature, 2 m specific humidity and surface pressure every
6 hours, 1991-2020, T62 Gaussian grid (about 1.9°), and the model ground height.

*NCEP-DOE Reanalysis 2 data provided by the NOAA PSL, Boulder, Colorado, USA, from their website at
https://psl.noaa.gov.* The data are a work of the U.S. Government and in the public domain.

## Files

| File | Content |
|---|---|
| `data/n0.bin` | 13 layers (whole year, then January ... December), each 121 × 241, `uint16` little-endian = N0 × 10. Row 0 is latitude +90, column 0 is longitude 0 (east). 758 KB |
| `data/N0_annual.txt`, `data/N0_01.txt` ... `N0_12.txt` | the same grids as text, one row per latitude, N-units with one decimal |
| `data/lat.txt`, `data/lon.txt` | grid coordinates (degrees) |
| `data/n0.json` | layout and provenance |
| `data/n0_r2_1991-2020.npz` | all layers as float32, also on the native Gaussian grid |
| `tools/fetch_r2.sh` | downloads the Reanalysis 2 files (about 5.3 GB, resumable) |
| `tools/build_n0.py` | builds the map |
| `tools/compare_n0.py` | prints the checks of [Validation](#validation) |
| `js/n0.js` | reader for browsers and Node (bilinear) |

## Use

Python:

```python
import numpy as np
n0 = np.fromfile('data/n0.bin', '<u2').reshape(13, 121, 241) / 10.0   # [layer, row, col]

def n0_at(lat, lon, layer=0):              # layer 0 = whole year, 1..12 = month
    y, x = (90 - lat) / 1.5, (lon % 360) / 1.5
    i, j = min(int(y), 119), min(int(x), 239)
    fy, fx = y - i, x - j
    a = n0[layer]
    return ((a[i, j] * (1 - fx) + a[i, j + 1] * fx) * (1 - fy)
            + (a[i + 1, j] * (1 - fx) + a[i + 1, j + 1] * fx) * fy)

print(round(n0_at(38.13, 15.56)), round(n0_at(38.13, 15.56, 8)))   # Messina: 344 (year), 373 (August)
```

JavaScript:

```js
import { n0At } from './js/n0.js';
const year = await n0At('data/n0.bin', 51.5, -0.1);        // London, whole year
const july = await n0At('data/n0.bin', 51.5, -0.1, 7);     // London, July
```

The text grids have the same layout as the ITU-R P.452 map files, so code written for those can read them directly.

## Validation

Details and tables in [validation/VALIDATION.md](validation/VALIDATION.md).

- **Against the ITU-R P.452 `N050` map** (annual): area-weighted mean difference **+2.0 N-units**, rms 9.7; half of
  the grid points within 4.9, 90 % within 14.8. Below 500 m of ground height (80 % of the points, oceans included) the
  mean difference is +2.6. The two maps differ in period (1991-2020 against 1983-1992) and in source (reanalysis against
  radiosondes); a slightly wetter recent climate is consistent with the small positive bias.
- **High ground**: the difference grows with the height of the ground (about −14 between 1500 and 3000 m, −19 above
  3000 m). It depends on the convention of the reduction to sea level: with h0 = 7.35 km (the global mean of P.453)
  instead of 9.46 km it reverses sign (+15 above 3000 m). This map keeps 9.46 km on purpose, to match ITM. The coarse
  model orography of the reanalysis also smooths valleys and plateaus.
- **Against the monthly maps of ITU-R P.453** (Figs 1 and 2, February and August): the same latitude pattern and
  seasonal swing, within about 10-20 N-units of the values read on the charts.
- **ITM climate classes**: the map is higher than the values suggested by the ITM guide in most classes (for example
  389 in the Congo basin against 360, 320 in the US Midwest against 301). Those are surface values (not reduced to sea
  level) of a 1960s climatology.

## Limitations

- About 1.9° source resolution: coasts, islands, narrow valleys and mountain ranges are smoothed. Use the value as
  climatology for a region, not for a single valley.
- A median: refractivity changes from day to day and from hour to hour (and ducting is a different matter, described
  by the refractivity gradient, not by N0).
- Not a substitute for local measurements, nor for the ITU maps where a Recommendation requires them in a formal
  procedure.

## Rebuild

```sh
tools/fetch_r2.sh raw                      # about 5.3 GB from psl.noaa.gov
python3 tools/build_n0.py raw data         # about 10 minutes, 3 GB of memory
python3 tools/compare_n0.py data/n0_r2_1991-2020.npz [--itu path/to/your/own/N050.TXT]
```

Python 3 with `numpy` and `netCDF4`.

## References

- Recommendation ITU-R P.453-14 (08/2019), *The radio refractive index: its formula and refractivity data*.
- Recommendation ITU-R P.452-18 (10/2023), *Prediction procedure for the evaluation of interference between stations on
  the surface of the Earth at frequencies above about 0.1 GHz* (compared only).
- E. K. Smith, S. Weintraub, *The constants in the equation for atmospheric refractive index at radio frequencies*,
  Proc. IRE 41 (8), 1953. doi:10.1109/JRPROC.1953.274297
- G. D. Thayer, *An improved equation for the radio refractive index of air*, Radio Science 9 (10), 803-807, 1974.
  doi:10.1029/RS009i010p00803
- B. R. Bean, G. D. Thayer, *Models of the atmospheric radio refractive index*, Proc. IRE 47 (5), 740-755, 1959.
  doi:10.1109/JRPROC.1959.287422
- B. R. Bean, J. D. Horn, A. M. Ozanich, *Climatic charts and data of the radio refractive index for the United States
  and the world*, NBS Monograph 22, 1960.
- B. R. Bean, E. J. Dutton, *Radio Meteorology*, NBS Monograph 92, 1966.
- P. L. Rice, A. G. Longley, K. A. Norton, A. P. Barsis, *Transmission loss predictions for tropospheric communication
  circuits*, NBS Technical Note 101, revised 1967 (eq. 4.3, the 9.46 km scale height used by ITM).
- G. A. Hufford, A. G. Longley, W. A. Kissick, *A guide to the use of the ITS Irregular Terrain Model in the area
  prediction mode*, NTIA Report 82-100, 1982 (suggested values per climate).
- M. Kanamitsu, W. Ebisuzaki, J. Woollen, S.-K. Yang, J. J. Hnilo, M. Fiorino, G. L. Potter, *NCEP-DOE AMIP-II
  Reanalysis (R-2)*, Bull. Amer. Meteor. Soc. 83 (11), 1631-1643, 2002. doi:10.1175/BAMS-83-11-1631

## Licence and citation

MIT, see [LICENSE](LICENSE). Data derived from public-domain NOAA data; if you use the map, please cite this
repository (see [CITATION.cff](CITATION.cff)) and acknowledge NOAA PSL as above.

Emanuele IW9GDC, for [RFscape](https://rfscape.com/).
