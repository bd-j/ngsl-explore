"""Extract text from the NGSL documentation PDFs.

build_catalog.py parses the stellar-parameter table out of aaareadme.txt, so
this must run after fetch_ngsl.sh and before build_catalog.py.
v2_versus_v1.pdf is skipped: 374 pages of plots with no extractable text.
"""
import pathlib
import pypdf

SKIP = {'v2_versus_v1'}

for p in sorted(pathlib.Path('docs').glob('*.pdf')):
    if p.stem in SKIP:
        continue
    reader = pypdf.PdfReader(str(p))
    txt = "\n".join((pg.extract_text() or "") for pg in reader.pages)
    out = p.with_suffix('.txt')
    out.write_text(txt)
    print(f'{p.stem}: {len(reader.pages)} pages, {len(txt)} chars -> {out}')
