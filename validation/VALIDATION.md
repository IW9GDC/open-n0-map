# Validation

All numbers from `tools/compare_n0.py` and the build of 07/11/2026. For the comparison with the ITU-R P.452 `N050` map only **differences** (this map minus ITU) are given: the ITU values themselves are not reproduced.

## Against ITU-R P.452 N050 (annual)

- area-weighted mean difference **+2.0** N-units, rms **9.7**
- absolute difference: median 4.9, 90th percentile 14.8, maximum 50.4

| Latitude band | mean difference |
|---|---|
| +60° to +90° | +0.5 |
| +30° to +60° | +0.9 |
| +0° to +30° | +3.3 |
| -30° to +0° | +4.8 |
| -60° to -30° | +0.9 |
| -90° to -60° | -6.3 |

| Ground height | grid points | mean difference | rms |
|---|---|---|---|
| 0-500 m | 23385 | +2.6 | 7.6 |
| 500-1500 m | 2928 | -5.4 | 11.2 |
| 1500-3000 m | 2028 | -14.2 | 16.3 |
| 3000-7000 m | 820 | -18.9 | 22.3 |

The dependence on height is the convention of the reduction to sea level. Re-reduced with h0 = 7.35 km (the global mean of ITU-R P.453) instead of 9.46 km, the mean differences per height band become: 0-500 m +3.2, 500-1500 m +3.4, 1500-3000 m +7.8, 3000-7000 m +15.4. The map keeps 9.46 km, the scale height ITM uses, so that ITM gives back the refractivity of the ground of each path.

The largest differences (about −45 to −50) are on the Himalayan foothills and the Gangetic plain, where the coarse model orography of the reanalysis puts the plain several hundred metres too high, on the Antarctic plateau (height convention), and over the Sahel (16.5° N), where the two sources disagree on the moisture of the monsoon season.

## Zonal medians (N-units)

To read against Figs 1 and 2 of ITU-R P.453 (monthly means of N0, February and August).

| Latitude | February | August | Whole year |
|---|---|---|---|
| +80° | 323 | 315 | 316 |
| +60° | 316 | 332 | 318 |
| +50° | 316 | 337 | 319 |
| +40° | 317 | 357 | 328 |
| +30° | 335 | 371 | 352 |
| +20° | 358 | 383 | 372 |
| +10° | 378 | 388 | 384 |
| +0° | 385 | 381 | 383 |
| -10° | 383 | 371 | 380 |
| -20° | 376 | 356 | 366 |
| -30° | 365 | 339 | 351 |
| -40° | 344 | 328 | 334 |
| -50° | 324 | 317 | 319 |
| -60° | 311 | 309 | 309 |
| -70° | 303 | 309 | 306 |
| -80° | 304 | 316 | 310 |

## ITM climate classes

Values suggested by the ITM guide (NTIA Report 82-100; surface values of a 1960s climatology, not reduced to sea level) and this map (whole year) at a representative point.

| Climate | ITM guide | this map |
|---|---|---|
| Equatorial (Congo) | 360 | 389 |
| Continental subtropical (Sudan) | 320 | 332 |
| Maritime subtropical (West Africa coast) | 370 | 383 |
| Desert (Sahara) | 280 | 299 |
| Continental temperate (US Midwest) | 301 | 320 |
| Maritime temperate over land (UK) | 320 | 326 |
| Maritime temperate over sea (North Sea) | 350 | 326 |

## Sample points (this map, N-units)

| Place | Whole year | January | July |
|---|---|---|---|
| Messina | 344 | 330 | 369 |
| Rome | 337 | 325 | 357 |
| London | 328 | 322 | 339 |
| Kinshasa | 381 | 387 | 357 |
| In Salah (Sahara) | 301 | 310 | 287 |
| Singapore | 391 | 387 | 390 |
| Tromsø | 316 | 311 | 329 |
| Denver | 314 | 311 | 334 |

## Change between 1991-2000 and 2011-2020

From `tools/decades_n0.py`: annual median N0 of each decade (same data and physics as the map), and the change of the
mean dry term (77.6 (P − e)/T) and wet term (72 e/T + 3.75·10⁵ e/T²), all reduced to sea level. The two decades are
20 years apart.

| | Change |
|---|---|
| 2 m temperature | +0.47 K |
| Water vapour pressure | +0.44 hPa (+2.6 %) |
| Dry term | −0.58 N-units |
| Wet term | +1.91 N-units |
| **N0 (median)** | **+1.26 N-units** (area-weighted; 62 % of the grid points increased) |

Warming lowers the dry term (about −1 N-unit per kelvin) and raises the wet term through the extra water vapour;
the vapour wins, so N0 increases. At constant relative humidity +0.47 K would give about +3.3 % of vapour
(Clausius-Clapeyron); +2.6 % is consistent with a slight drop of relative humidity over land.

| Latitude band | ΔN0 | dry | wet | ΔT (K) |
|---|---|---|---|---|
| +60° to +90° | +0.0 | −2.3 | +2.2 | +1.86 |
| +30° to +60° | +1.3 | −0.5 | +2.1 | +0.48 |
| 0° to +30° | +2.6 | −0.5 | +3.1 | +0.34 |
| −30° to 0° | +1.2 | −0.3 | +1.6 | +0.27 |
| −60° to −30° | +0.7 | −0.3 | +1.1 | +0.31 |
| −90° to −60° | −0.8 | −1.2 | +0.3 | +0.88 |

The Arctic warmed four times the global mean, but with so little vapour the two terms cancel. In dry regions the
change can be negative (Sahara −0.3, Denver −2.1); in humid ones it is largest (Delhi +5.0, Messina +2.8).

The trend, about +0.06 N-units per year worldwide and +0.13 between the equator and 30° N, accounts for roughly half of
the difference from the ITU-R P.452 `N050` map, whose radiosonde data (1983-1992) are about 17 years older than the
centre of this map's period: about +1.1 worldwide and +2.3 in the northern tropics, against +2.0 and +3.3 observed.
The rest comes from the instruments (radiosonde humidity sensors against a reanalysis), the observing times and the
station coverage. Trends in a reanalysis also carry changes of the observing system: take the order of magnitude,
not the decimals.

## ΔN against ITU-R P.452 / P.1812 DN50

From `tools/compare_dn.py` and the build of 08/11/2026 (median map; only differences, this map minus ITU, are given).

- area-weighted mean difference **+2.7** N-units/km, rms **7.1**

| Latitude band | mean difference | rms |
|---|---|---|
| +60° to +90° | +0.8 | 2.4 |
| +30° to +60° | +2.7 | 5.6 |
| +0° to +30° | +5.3 | 9.9 |
| -30° to +0° | +4.5 | 8.0 |
| -60° to -30° | -0.9 | 3.7 |
| -90° to -60° | -1.5 | 3.1 |

| Ground height | grid points | mean difference | rms |
|---|---|---|---|
| below 500 m | 23379 | +3.1 | 7.0 |
| 500-1500 m | 2947 | +1.8 | 7.8 |
| 1500-3000 m | 2059 | -2.7 | 5.4 |
| above 3000 m | 776 | -5.1 | 6.9 |

The mean of ΔN (instead of the median) is about 1 N-unit/km higher everywhere. A single year (1991) already gives
+2.9 / rms 7.3, so the differences are not a matter of sample size: the likely causes are the different period, the
different source (reanalysis against radiosondes, which see thin low-level inversions) and the coarse orography.

### In use: ITU-R P.1812 on the ITU-R SG3 databank

The reference implementation of Recommendation ITU-R P.1812-8 was run on 53,467 paths of the ITU-R Study Group 3
databank of VHF/UHF measurements, once with the ITU maps (DN50, N050) and once with this repository's maps (ΔN, and the
annual N0). Error = predicted minus measured field strength:

| Case | ITU maps: mean / sd | Open maps: mean / sd |
|---|---|---|
| transmitter at 50 m or more, receiver at 8 m or more (35,650) | +5.7 / 8.5 dB | +5.7 / 8.5 dB |
| low transmitting antennas (10,154) | +12.8 / 13.9 dB | +13.1 / 13.9 dB |
| all (53,352) | +7.3 / 10.3 dB | +7.4 / 10.3 dB |

For a propagation prediction the two maps are interchangeable.

