#!/bin/bash
# Fetch NGSL v2 documentation, star listing, and spectra from MAST.
# Source: https://archive.stsci.edu/prepds/stisngsl/
set -euo pipefail
cd "$(dirname "$0")/.."
BASE=https://archive.stsci.edu
mkdir -p data docs

for f in aaareadme File_format_contents_V2 Table_V2 Lindler_AAS_Jan2010 v2_versus_v1; do
    curl -sSL -o "docs/${f}.pdf" "${BASE}/prepds/stisngsl/docs/${f}.pdf"
done

curl -sSL -o data/stis_ngsl_v2.zip "${BASE}/pub/hlsp/stisngsl/v2/stis_ngsl_v2.zip"

# Full 374-star listing as CSV from the MAST search interface
curl -sSL -o data/ngsl_listing.csv \
  "${BASE}/prepds/stisngsl/search.php?ordercolumn1=targname&max_records=500&max_rpp=500&action=Search&outputformat=CSV&nondefault=&selectedColumnsCsv=targname,ra,dec,bmag,vmag,b-v,spectral_type"

unzip -oq data/stis_ngsl_v2.zip -d data/spectra
