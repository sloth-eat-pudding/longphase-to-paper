# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is **not a software project** — it is the LaTeX manuscript, supplementary material, and
figure-generation scripts for the paper *"Joint Somatic Phasing and Purity Estimation
Incorporating Loss of Heterozygosity in Tumor-Only Long-Read Sequencing"* (the LongPhase-TO
method). The paper is written for a Springer Nature journal and uses the official
`sn-jnl` document class with the `sn-nature` option.

**`paper.tex` is the target deliverable** (the journal submission). All manuscript work happens
here plus its `supplementary.tex`. `thesis.tex` is a legacy backup — do not edit it or keep it in
sync with `paper.tex`.

Work here is prose editing (`.tex`), reference management (`.bib`), and regenerating benchmark
figures (`.py` → `.png`). There is no application to run or test.

## Build / compile

There is no Makefile; compile with the TeX Live toolchain directly. The manuscript cross-references
the supplementary via `\externaldocument{Supplementary}` (the `xr-hyper` package), so
**`Supplementary` must be compiled first** to produce `Supplementary.aux` before `paper.tex` can
resolve `\ref` calls into the supplement. Correct order:

```sh
# 1. Build the supplement first (generates Supplementary.aux for cross-refs)
pdflatex Supplementary && bibtex Supplementary && pdflatex Supplementary && pdflatex Supplementary

# 2. Then the main manuscript
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```

`latexmk -pdf paper` also works and handles the rerun loop, but still requires `Supplementary.aux`
to exist for the supplementary cross-references to resolve.

Note: macOS is case-insensitive, so `Supplementary.*` (capital, the built artifacts) and
`supplementary.tex` (lowercase, as git tracks it) refer to the same files.

## Manuscript files

- **`paper.tex`** — the journal manuscript and sole editing target (Introduction / Results / Discussion / Methods).
- **`supplementary.tex`** — Supplementary Results/Methods/Figures/Tables plus command-line arguments.
- **`thesis.tex`** — legacy backup only; ignore it (see "Legacy / ignore" below).

## Figures

Final figures are `Figure_1.pdf` … `Figure_4.pdf`, included in `paper.tex` via
`\includegraphics`. Each has a matching **`Figure_N.pptx`** — the PowerPoint is the editable
source; the PDF is exported from it. To change a main figure, edit the `.pptx` and re-export
the `.pdf` (do not hand-edit the PDF).

Some figure panels are matplotlib charts produced by standalone Python scripts. Each script has
its **benchmark data hardcoded as inline arrays** (there is no external data file) and writes one
PNG next to itself:

| Script | Output | Content |
|---|---|---|
| `Phasing_N50.py` | `phasing_n50_grouped_bar_kbp.png` | Block N50 (Kbp) across phasing tools |
| `Block_sum.py` | `block_sum_grouped_bar.png` | Total phased block sum (bp) |
| `LOH_N50.py` | `loh_n50_grouped.png` | LOH segment N50 vs ASCAT |
| `LOH_burden.py` | `loh_burden_grouped_bar.png` | Genome-wide LOH burden |
| `Q1Q3.py` | (purity Q1–Q3 spread plot) | Purity estimates vs expected purity |

Run a script with `python3 Phasing_N50.py`. They require `numpy`, `pandas`, and `matplotlib`.
To update a benchmark number, edit the array literal in the corresponding script and rerun it;
the styling is intentionally "Nature-style" grouped bars.

The eight benchmark datasets are consistent across scripts: `COLO829_ONT`, `COLO829_NYGC`,
`H1437_UCSC`, `H2009_UCSC`, `HCC1395_HKU`, `HCC1395_NYGC`, `HCC1937_UCSC`, `HCC1954_UCSC`.

## Editing conventions (paper.tex & Supplementary.tex)

- **Cross-references live in the Supplement.** `paper.tex` contains only Figures 1–4
  (`Figure_1.pdf`…`Figure_4.pdf`); **every table and every Methods figure is defined in
  `Supplementary.tex`** and reached via `\externaldocument`. Cite them as
  **`Supplementary Fig.~\ref{…}`** (plural **`Supplementary Figs.~\ref{…}`**) /
  **`Supplementary Table~\ref{…}`** (not bare "Figure"/"Table"). Springer Nature abbreviates
  "Fig." for figures but spells out "Table"; the supplement's caption name
  (`\captionsetup[figure]{name=Supplementary Fig.}`) must match the in-text callout. Numbering is
  plain arabic (no "S" prefix). Note this differs from the main-figure style below, which keeps the
  full word "Figure" by deliberate author choice.
  Supplement labels are pandoc-derived: figures `fig:met-page-N-cropped-jpg` /
  `fig:app-page-N-cropped-jpg`; tables `tab:*` (e.g. `tab:cell-lines`, `tab:callers`,
  `tab:loh_summary`, `tab:phasing-performance`, `tab:somatic_purity`, `tab:tumor-purity`).
  **Confirm a label exists in `Supplementary.tex` before adding a `\ref`** — mismatches
  (`tab:loh-performance` vs `tab:loh_summary`) and missing figures render silently as `??`.
- **Main-figure style.** Main figures use the filename as the label (`\label{Figure_1.pdf}`) and are
  cited `Figure~\ref{Figure_1.pdf}(a)` — full word "Figure", panel letter in parentheses. Keep this.
- **Terminology (keep identical in both files).** The core model is the **"triplet graph"** (never
  "Tri-nodal-edge" or "triplate"), contrasted with the **"pairwise graph"**. Somatic callers always
  carry the **`-TO`** suffix (ClairS-TO, DeepSomatic-TO, Severus-TO). Use American spelling
  (**"tumor"**, not "tumour").
- **The estimated quantity is "tumor DNA fraction", never "tumor purity".** The method infers the
  fraction of tumor-derived reads directly from germline haplotype imbalance; this is deliberately
  distinguished from *tumor cellular purity* (which requires aneuploidy/CNV modeling). "cellular
  purity" appears only in the **Discussion**, as the explicit contrast term cited against the
  purity/ploidy callers ASCAT, PURPLE, and Battenberg (`vanloo2010`, `priestley2019`, `nikzainal2012`). Note: the LaTeX
  **label identifiers** `tab:somatic_purity`, `tab:somatic_indel_purity`, `tab:tumor-purity`, and
  `app:purity-aware-workflow` still contain "purity" on purpose — they are internal keys referenced
  across both files; **do not rename them**. Only display text says "DNA fraction".
- **Math symbols:** prefer LaTeX macros (`$\ge$`, `$\le$`, `$\times$`) over literal Unicode — a
  literal `≥` in the Supplement broke the `pdflatex` build.
- **No em dashes in prose.** Do not use em dashes (`—` or LaTeX `---`) in manuscript text. Recast
  with commas, parentheses, or a colon instead (e.g. a dash-delimited list becomes a parenthetical).
  This is a fixed author preference for both `paper.tex` and `Supplementary.tex`.
- **After any edit, rebuild (Supplement then paper) and grep the `.log` for `undefined` and
  `Missing character`** — the cross-file `\ref`s make broken references easy to introduce and easy
  to miss.

## Do not edit

- **`sn-jnl.cls`** and **`sn-nature.bst`** — the unmodified Springer Nature template class and
  BibTeX style. Treat as read-only vendor files.
- **`references.bib.*.bak`** — backups of the bibliography; edit `references.bib` only.

## Legacy / ignore

**`thesis.tex`** is a legacy backup of the manuscript; it is not the deliverable and should not be
edited or kept in sync with `paper.tex`.

Many `*_chapter*.md`, `*_section*.txt`, and `*_citation_*.json` files are **deleted** in the
working tree (they appear in `git status` as `D`). These were intermediate drafts from an earlier
pandoc-based chapter→LaTeX workflow (which is also why some labels look auto-generated, e.g.
`\ref{fig:met-page-14-cropped-jpg}`). They are not part of the current workflow — do not
resurrect or rely on them. The `page_*_cropped.*` images are cropped source scans used for some
method figures.
