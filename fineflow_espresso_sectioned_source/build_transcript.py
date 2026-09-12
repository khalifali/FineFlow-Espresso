"""Build the side-by-side presenter PDF from the sectioned deck and notes.

Usage (from this directory):
  python3 -m pip install pymupdf
  pdflatex fineflow_espresso_sectioned.tex  # run twice after section edits
  python3 build_transcript.py
Edit transcript.json to change notes. Page keys must match the compiled deck.
"""
from pathlib import Path
import json
import fitz

root = Path(__file__).resolve().parent
deck = fitz.open(root / 'fineflow_espresso_sectioned.pdf')
notes = json.loads((root / 'transcript.json').read_text())
assert set(notes) == {str(i + 1) for i in range(len(deck))}
out = fitz.open()
for i, slide in enumerate(deck):
    p = out.new_page(width=1200, height=675)
    p.insert_text((28, 40), 'FineFlow-Espresso: Presenter transcript', fontsize=20, color=(0, .396, .741))
    p.insert_text((28, 65), f'Slide {i+1} of {len(deck)}', fontsize=12)
    p.show_pdf_page(fitz.Rect(20, 110, 720, 504), deck, i)
    p.draw_line((740, 96), (740, 615), color=(.8, .8, .8), width=.6)
    p.insert_text((770, 118), 'Presenter transcript', fontsize=17, color=(0, .396, .741))
    spare = p.insert_textbox(fitz.Rect(770, 146, 1170, 630), notes[str(i+1)], fontsize=16, lineheight=1.4, fontname='helv')
    assert spare >= 0, (i+1, spare)
out.save(root / 'fineflow_espresso_sectioned_with_transcript.pdf', garbage=4, deflate=True)
print(f'Built {len(out)} presenter pages')
