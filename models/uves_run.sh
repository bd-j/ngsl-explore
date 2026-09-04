#!/bin/bash
cd /Users/bjohnson/Projects/ngsl
python3 scripts/make_atlas_model.py --star HD162393 --teff 9955 --logg 4.05 --feh -0.55
python3 scripts/make_atlas_model.py --star HD162678 --teff 9908 --logg 3.53 --feh 0.03
python3 scripts/make_atlas_model.py --star HD188294 --teff 11016 --logg 4.04 --feh 0.05
