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
