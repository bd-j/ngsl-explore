#!/bin/bash
# Fetch NGSL v2 documentation, star listing, and spectra from MAST.
# Source: https://archive.stsci.edu/prepds/stisngsl/
set -euo pipefail
cd "$(dirname "$0")/.."
BASE=https://archive.stsci.edu
mkdir -p data docs/ngsl_delivery

for f in aaareadme File_format_contents_V2 Table_V2 Lindler_AAS_Jan2010 v2_versus_v1; do
    curl -sSL -o "docs/ngsl_delivery/${f}.pdf" "${BASE}/prepds/stisngsl/docs/${f}.pdf"
done

curl -sSL -o data/stis_ngsl_v2.zip "${BASE}/pub/hlsp/stisngsl/v2/stis_ngsl_v2.zip"

# Full 374-star listing as CSV from the MAST search interface
curl -sSL -o data/ngsl_listing.csv \
  "${BASE}/prepds/stisngsl/search.php?ordercolumn1=targname&max_records=500&max_rpp=500&action=Search&outputformat=CSV&nondefault=&selectedColumnsCsv=targname,ra,dec,bmag,vmag,b-v,spectral_type"

unzip -oq data/stis_ngsl_v2.zip -d data/spectra

# STIS model line spread functions, for the true (not sampling-limited) resolution
# https://www.stsci.edu/hst/instrumentation/stis/performance/spectral-resolution
LSF="https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/stis/performance/spectral-resolution/_documents/LSF"
mkdir -p data/stis_lsf
for f in LSF_G230L_1700 LSF_G230L_2400 LSF_G430L_3200 LSF_G430L_5500 LSF_G750L_7000; do
    curl -sSL -A "Mozilla/5.0" -o "data/stis_lsf/${f}.txt" "${LSF}/${f}.txt"
done

# Pickles (1998) stellar spectral flux atlas, UVILIB
# https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/pickles-atlas
PK="https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/pickles"
mkdir -p data/pickles
curl -sSL -A "Mozilla/5.0" -o data/pickles/AA_README "${PK}/AA_README"
curl -sSL -A "Mozilla/5.0" -o data/pickles/pickles_index.fits "${PK}/dat_uvi/pickles.fits"
for i in $(seq 1 131); do
    curl -sSL -A "Mozilla/5.0" -o "data/pickles/pickles_${i}.fits" \
        "${PK}/dat_uvi/pickles_${i}.fits"
done
