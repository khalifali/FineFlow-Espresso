Compile fineflow_espresso_sectioned.tex twice with pdflatex, then run python3 build_transcript.py (requires pymupdf).

49 slides with seven section dividers matching PoreAccess-CO2. Edit transcript.json for presenter notes.
Run ../fineflow_espresso_solver/fineflow_espresso/compare_pressure_modes.py to regenerate the pressure comparison.
Baseline field figures use constant 9 bar and epsilon0=0.55; pressure and closure comparisons are regenerated with the 0.55 baseline.
See Ali_Khalifa_default_deck_sectioning.md for the approved section style.

The results section includes an axial porosity-profile figure generated from the saved solver porosity history.
