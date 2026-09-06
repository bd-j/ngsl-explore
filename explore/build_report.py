"""Build report.html from report.src.html, inlining the figures as data URIs.

Artifacts cannot load external images, so each {{FIG_*}} placeholder is
replaced with a base64 data URI of the corresponding PNG.
"""
import base64
import pathlib

FIGURES = {
    '{{FIG_SNR}}': 'figures/snr_vs_wavelength.png',
    '{{FIG_COV}}': 'figures/parameter_coverage.png',
    '{{FIG_BALMER}}': 'figures/balmer_break_candidates.png',
}

src = pathlib.Path('report.src.html').read_text()
for token, path in FIGURES.items():
    if token not in src:
        raise SystemExit(f'placeholder {token} not found in report.src.html')
    b64 = base64.b64encode(pathlib.Path(path).read_bytes()).decode()
    src = src.replace(token, f'data:image/png;base64,{b64}')

if '{{' in src:
    raise SystemExit('unsubstituted placeholder remains')

pathlib.Path('report.html').write_text(src)
print(f'report.html: {len(src)/1024/1024:.2f} MB')
