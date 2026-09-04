#!/bin/bash
cd /Users/bjohnson/Projects/ngsl
while IFS=, read -r star teff logg feh; do
  python3 scripts/make_atlas_model.py --star "$star" --teff "$teff" --logg "$logg" --feh "$feh"
done <<'ROWS'
HD040573,10200,4.2,-0.4
HD147550,10074,3.9,-0.0
HD128801,10123,3.7,-1.9
HD143459,9878,3.6,-0.6
ROWS
