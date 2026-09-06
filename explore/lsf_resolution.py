"""NGSL's actual spectral resolution, from the STIS model line spread functions.

Source: https://www.stsci.edu/hst/instrumentation/stis/performance/spectral-resolution
Column 52x0.2 is NGSL's aperture. FWHM is measured in pixels by interpolating
the half-maximum crossings, then converted to Angstroms using the dispersion
measured from the delivered NGSL wavelength grid.

Caveat carried in the output: that page tabulates the MAMA G230L, not the CCD
G230LB that NGSL actually used. G230L FWHM is therefore reported in pixels
only - converting it would require the MAMA dispersion, a different detector.

Writes data/stis_lsf_resolution.csv
"""
import csv
import numpy as np

# Dispersion measured from the NGSL v2 wavelength grid (grating_summary.txt).
DISP = {'G430L': 2.744, 'G750L': 4.877}
SEG = {'G430L': (3058, 5647), 'G750L': (5647, 10198)}

FILES = [('LSF_G230L_1700', 'G230L', 1700), ('LSF_G230L_2400', 'G230L', 2400),
         ('LSF_G430L_3200', 'G430L', 3200), ('LSF_G430L_5500', 'G430L', 5500),
         ('LSF_G750L_7000', 'G750L', 7000)]
APERTURE = '52x0.2'
BALMER = 3646.0


def fwhm_pixels(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    col = lines[1].split().index(APERTURE)
    d = np.array([[float(v) for v in ln.split()] for ln in lines[2:] if ln.strip()])
    x, y = d[:, 0], d[:, col] / d[:, col].max()
    pk = int(np.argmax(y))
    left = np.interp(0.5, y[:pk + 1], x[:pk + 1])
    right = np.interp(0.5, y[pk:][::-1], x[pk:][::-1])
    return right - left


rows = []
for stem, grating, wave in FILES:
    fw_px = fwhm_pixels(f'data/stis_lsf/{stem}.txt')
    r = {'grating': grating, 'wavelength_A': wave, 'aperture': APERTURE,
         'fwhm_pixels': round(float(fw_px), 3), 'ngsl_mode': 'yes'}
    if grating in DISP:
        fw_a = fw_px * DISP[grating]
        r.update(dispersion_A_per_pix=DISP[grating],
                 fwhm_angstrom=round(float(fw_a), 2),
                 R_actual=round(float(wave / fw_a)),
                 R_if_2pixel=round(float(wave / (2 * DISP[grating]))))
    else:
        # MAMA G230L, not NGSL's CCD G230LB: pixels only, no conversion.
        r.update(dispersion_A_per_pix='', fwhm_angstrom='', R_actual='',
                 R_if_2pixel='', ngsl_mode='no (MAMA G230L, not CCD G230LB)')
    rows.append(r)

# G430L FWHM varies slowly and near-linearly with wavelength; interpolate
# between the two tabulated points to reach the segment ends and the Balmer limit.
g = [r for r in rows if r['grating'] == 'G430L']
xw = [g[0]['wavelength_A'], g[1]['wavelength_A']]
xf = [g[0]['fwhm_pixels'], g[1]['fwhm_pixels']]
for wave, label in [(SEG['G430L'][0], 'G430L blue end'), (BALMER, 'Balmer limit'),
                    (SEG['G430L'][1], 'G430L red end')]:
    fw_px = float(np.interp(wave, xw, xf))
    fw_a = fw_px * DISP['G430L']
    rows.append({'grating': 'G430L', 'wavelength_A': round(wave), 'aperture': APERTURE,
                 'fwhm_pixels': round(fw_px, 3), 'ngsl_mode': f'yes (interpolated: {label})',
                 'dispersion_A_per_pix': DISP['G430L'], 'fwhm_angstrom': round(fw_a, 2),
                 'R_actual': round(wave / fw_a),
                 'R_if_2pixel': round(wave / (2 * DISP['G430L']))})

with open('data/stis_lsf_resolution.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

hdr = ['grating', 'wavelength_A', 'fwhm_pixels', 'fwhm_angstrom', 'R_actual',
       'R_if_2pixel', 'ngsl_mode']
print(''.join(f'{h:>16}' for h in hdr[:-1]) + '  ' + hdr[-1])
for r in rows:
    print(''.join(f'{str(r[h]):>16}' for h in hdr[:-1]) + '  ' + str(r['ngsl_mode']))
print('\n-> data/stis_lsf_resolution.csv')
