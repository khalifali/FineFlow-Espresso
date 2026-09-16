#!/usr/bin/env python3
"""Add the axial porosity-profile output and matching results slide.

The script is intentionally idempotent so the associated GitHub Actions workflow
can be re-run without duplicating code, slides, or presenter notes.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def update_solver() -> None:
    path = ROOT / "fineflow_espresso_solver/fineflow_espresso/fineflow_espresso.py"
    text = path.read_text(encoding="utf-8")
    if 'f"{prefix}_porosity_profiles"' in text:
        return

    marker = '    _save_publication_figure(fig, output_dir / f"{prefix}_relative_permeability")\n'
    if marker not in text:
        raise RuntimeError("Could not locate publication-figure insertion point")

    addition = r'''

    # Axial porosity profiles at representative extraction times. These are
    # the same 1D observable that can be compared with interrupted CT profiles.
    target_times_s = (0.0, min(10.0, time[-1]), min(30.0, time[-1]), time[-1])
    profile_indices = []
    for target_time in target_times_s:
        index = int(np.argmin(np.abs(time - target_time)))
        if index not in profile_indices:
            profile_indices.append(index)

    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    colors = ("#9A9A9A", "#66B5D8", "#005293", "#6F4E37")
    for color, index in zip(colors, profile_indices):
        ax.plot(
            z_mm,
            result.porosity[index],
            color=color,
            lw=2.1,
            label=f"{time[index]:g} s",
        )
    selected_porosity = result.porosity[profile_indices]
    ymin = float(np.min(selected_porosity))
    ymax = float(np.max(selected_porosity))
    margin = max(0.0015, 0.08 * max(ymax - ymin, 1.0e-6))
    ax.set_ylim(max(0.0, ymin - margin), min(1.0, ymax + margin))
    ax.set_xlabel("Depth from inlet / top of puck [mm]")
    ax.set_ylabel("Porosity [-]")
    ax.grid(alpha=0.22)
    ax.legend(title="Extraction time", loc="best", frameon=False)
    _save_publication_figure(fig, output_dir / f"{prefix}_porosity_profiles")
'''
    path.write_text(text.replace(marker, marker + addition, 1), encoding="utf-8")


def update_tests() -> None:
    path = ROOT / "fineflow_espresso_solver/fineflow_espresso/test_fineflow.py"
    text = path.read_text(encoding="utf-8")
    if "named_test_porosity_profiles.png" in text:
        return

    marker = '            self.assertIn("named_test_relative_permeability.png", names)\n'
    if marker not in text:
        raise RuntimeError("Could not locate publication-output test insertion point")
    replacement = (
        marker
        + '            self.assertIn("named_test_porosity_profiles.png", names)\n'
        + '            self.assertIn("named_test_porosity_profiles.pdf", names)\n'
    )
    path.write_text(text.replace(marker, replacement, 1), encoding="utf-8")


def update_deck() -> None:
    path = ROOT / "fineflow_espresso_sectioned_source/fineflow_espresso_sectioned.tex"
    text = path.read_text(encoding="utf-8")
    title = "Porosity profile evolves as fines redistribute"
    if title in text:
        return

    anchor = r"\begin{frame}{Fine puck: constant pressure versus PI}"
    if anchor not in text:
        raise RuntimeError("Could not locate results-slide insertion point")

    slide = r'''
\begin{frame}{Porosity profile evolves as fines redistribute}
{\small\color{tumblue}Constant 9 bar across puck + basket. Exponential permeability. Illustrative, uncalibrated case.}\par\medskip
\begin{columns}[T,onlytextwidth]
 \column{.66\textwidth}\centering
  \safeimage{porosity_profiles.pdf}{\linewidth}{4.05cm}
 \column{.30\textwidth}\small
  \hi{Predicted structure}
  \begin{itemize}
   \item Direct model observable: $\varepsilon(z,t)$.
   \item Initial profile: $\varepsilon_0=0.36$.
   \item Deposition lowers porosity, most strongly near the outlet/bottom.
   \item At 60 s: $\varepsilon\approx0.338$--$0.355$.
  \end{itemize}
\end{columns}
\vspace{.08cm}
\takeaway{The predicted axial porosity evolution can be compared directly with interrupted CT profiles at matched extraction times.}
\slidesource{Model $t=0$ is the start of extraction in the saturated 1D model; dry and preinfusion stages are not represented separately.}
\end{frame}

'''
    path.write_text(text.replace(anchor, slide + anchor, 1), encoding="utf-8")


def update_slide_map() -> None:
    path = ROOT / "fineflow_espresso_sectioned_source/slide_map.json"
    slides = json.loads(path.read_text(encoding="utf-8"))
    title = "Porosity profile evolves as fines redistribute"
    if any(item.get("title") == title for item in slides):
        return

    for item in slides:
        if int(item["slide"]) >= 36:
            item["slide"] = int(item["slide"]) + 1
    insert_at = next(i for i, item in enumerate(slides) if int(item["slide"]) == 35) + 1
    slides.insert(
        insert_at,
        {
            "slide": 36,
            "section": "Illustrative results and model sensitivity",
            "title": title,
        },
    )
    path.write_text(json.dumps(slides, indent=2) + "\n", encoding="utf-8")


def update_transcript() -> None:
    path = ROOT / "fineflow_espresso_sectioned_source/transcript.json"
    transcript = json.loads(path.read_text(encoding="utf-8"))
    note = (
        "This plot uses the axial porosity history already returned by the solver. "
        "The initial illustrative bed is uniform at epsilon 0.36. As fines are released, "
        "transported and recaptured, deposition lowers porosity throughout the puck, with "
        "the strongest decrease near the outlet because the current capture closure is weighted there. "
        "At sixty seconds the profile spans about 0.338 to 0.355. This is not a dry-to-wet comparison: "
        "model time zero is the start of full extraction in an already saturated one-dimensional bed. "
        "The direct validation target is therefore the CT porosity profile at matched extraction times."
    )
    if "36" in transcript and "axial porosity history" in transcript["36"]:
        return

    shifted: dict[str, str] = {}
    for key, value in transcript.items():
        number = int(key)
        shifted[str(number + 1 if number >= 36 else number)] = value
    shifted["36"] = note
    last = max(map(int, shifted))
    ordered = {str(i): shifted[str(i)] for i in range(1, last + 1)}
    path.write_text(
        json.dumps(ordered, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def update_readme() -> None:
    path = ROOT / "fineflow_espresso_sectioned_source/README_BUILD.txt"
    text = path.read_text(encoding="utf-8")
    text = text.replace(
        "48 slides with seven section dividers",
        "49 slides with seven section dividers",
    )
    if "axial porosity-profile figure" not in text:
        text += (
            "\nThe results section includes an axial porosity-profile figure generated "
            "from the saved solver porosity history.\n"
        )
    path.write_text(text, encoding="utf-8")


def main() -> None:
    update_solver()
    update_tests()
    update_deck()
    update_slide_map()
    update_transcript()
    update_readme()


if __name__ == "__main__":
    main()
