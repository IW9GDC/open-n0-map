#!/bin/bash
# Open N0 — download the NCEP-DOE Reanalysis 2 pressure-level fields used to build the refractivity lapse-rate map DN
# (public domain, NOAA PSL): air temperature and relative humidity, 6-hourly, 1991-2020, 17 levels,
# 2.5 degree grid (~495 MB per year, ~15 GB; the heights of the levels are computed hydrostatically by build_dn.py). The surface fields come from tools/fetch_r2.sh.
# Usage: tools/fetch_r2_levels.sh [dest]   (default raw-levels). Resumable. MIT license (https://github.com/IW9GDC/open-n0-map).
set -u
DEST="${1:-$(dirname "$0")/../raw-levels}"; mkdir -p "$DEST"; cd "$DEST" || exit 1
BASE=https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/pressure
get() { [ -s "$1" ] && return 0; curl -sf -L -C - --retry 6 --retry-delay 20 -o "$1.part" "$BASE/$1" && mv "$1.part" "$1" && echo "ok $1" || echo "FAILED $1"; }
export -f get; export BASE
# the server limits each connection (~0.4 MB/s), not the total: sixteen files at a time
for y in $(seq 1991 2020); do for v in air rhum; do echo "$v.$y.nc"; done; done | xargs -P 16 -I{} bash -c 'get "$@"' _ {}
echo "done $(ls *.nc | wc -l) files"
