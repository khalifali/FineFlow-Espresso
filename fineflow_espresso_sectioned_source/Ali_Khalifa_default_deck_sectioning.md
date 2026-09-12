# Ali Khalifa — default presentation sectioning style

Explicit user preference, 12 September 2026: use the PoreAccess-CO2 Beamer sectioning style for all future decks unless the user requests another style. Preserve the relevant deck branding and content. This is a reusable style reference, not an automatically installed application setting.

Reference: https://github.com/khalifali/PoreAccess-CO2/blob/main/beamer/PoreAccess_CO2_presentation/poreaccess_co2_results.tex

- Start with the title slide, then an automatic Beamer table of contents.
- Use a full section-divider slide before every section, including the first.
- Divider: white background, gray large “Section X of Y”, 0.4 cm gap, large bold blue section title, 0.5 cm gap, one short explanatory subtitle. Vertically centre the group using vfill above and below.
- Preserve the TUM logo at top right when using this TUM theme.
- Footer: authors and project on the left, short section name in the middle, slide number/total at right.
- Declare sections with short and long names: \section[Methods]{Numerical methods and verification}.
- Use normal numbered divider frames, so PDF page counts and transcript keys include the dividers.
- Adapt the section count and names to the project. Do not force seven sections for every deck.
- Compile twice and synchronize the side-by-side presenter transcript after reordering or adding dividers.

The preferred treatment is full section transitions, not only a small section label in the footer.
