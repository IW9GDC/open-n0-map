#!/bin/bash
# Open N0 — download the NCEP-DOE Reanalysis 2 fields used to build the surface refractivity map (public domain, NOAA PSL).
# 2 m air temperature, 2 m specific humidity, surface pressure, 6-hourly, 1991-2020, plus the model surface height.
# Usage: tools/fetch_r2.sh [dest]   (default raw/). Resumable. MIT license.
set -u
DEST="${1:-$(dirname "$0")/../raw}"; mkdir -p "$DEST"; cd "$DEST" || exit 1
BASE=https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/gaussian_grid
get() { [ -s "$1" ] && return 0; curl -sf -L --retry 4 --retry-delay 10 -o "$1.part" "$BASE/$1" && mv "$1.part" "$1" && echo "ok $1" || echo "FAILED $1"; }
export -f get; export BASE
get hgt.sfc.gauss.nc
# the server limits each connection (~0.4 MB/s), not the total: six files at a time
for y in $(seq 1991 2020); do for v in air.2m shum.2m pres.sfc; do echo "$v.gauss.$y.nc"; done; done | xargs -P 6 -I{} bash -c 'get "$@"' _ {}
echo "done $(ls *.nc | wc -l) files"
