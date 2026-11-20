# NAR Results/Discussion restructure: three options to review

Manuscript: `NAR/paper.tex` (Results start at line 192, Discussion at line 237). The current Results subsections are:

1. Overview of the LongPhase-TO framework
2. LongPhase-TO improves haplotype contiguity and delineates chromosome-scale LOH (cell lines, no truth; headline "median 2.9-fold block N50 over germline phasers")
3. Haplotype-aware recalibration improves somatic SNV and indel calls (8 ONT cell lines, dilution series)
4. Phasing-informed copy-number fitting recovers tumor DNA fraction and ploidy (ONT, PacBio, HG008 with NIST CN/LOH truth)
5. Somatic recalibration transfers to PacBio HiFi data

Framing (author decision): one tumor-only graph gives somatic haplotypes, LOH, and phasing-informed DNA fraction/ploidy/CN, each better than the dedicated tool. A downstream use of somatic haplotypes will close Results later.

## New results to place

Data in this folder (copied from the student's commits):

- `issue2_hg008_runs.tsv`: switch error (SER), Hamming, N50, phase-correct N50 (blocks split at each switch error) on five HG008 ONT mixtures, against an assembly-derived HG008-N truth (GIAB v6.3, dipcall). Tools: LongPhase-TO, LongPhase 2.1.2, WhatsHap, HapCUT2.
- `issue2_celllines_concordance.tsv`: the same metrics on 8 ONT cell lines at TDF 1.0, against LongPhase phasing of the matched normal ("concordance", not error). The LongPhase row here is the Figure 2 run.
- `issue9_hg008_snv_indel_metrics.tsv`: somatic SNV/indel precision/recall/F1 on the five HG008 mixtures against the NIST HG008-T v0.4 somatic benchmark. Use only `ClairS-TO raw PASS`, `DeepSomatic-TO raw PASS`, `LongPhase-TO(ClairS-TO)`, `LongPhase-TO(DeepSomatic-TO)`; ignore the `+ClairS-TO PoN` rows.

Summary:

- SER: LongPhase-TO 0.05–0.15% vs WhatsHap/HapCUT2 0.16–0.48% on all datasets. Against LongPhase: HG008 0.2–0.8 slightly higher (0.058–0.083 vs 0.049–0.068), pure tumor lower; lower in 6 of 8 cell lines.
- Phase-correct N50 gain over LongPhase: HG008 1.2× at 0.2–0.8, 2.3× pure; cell lines 0.95–2.9×. Over WhatsHap/HapCUT2: 1.8–3.8×. The current headline uses raw N50 (2.9-fold median).
- Hamming 7–16% for all tools (cause unknown); author: do not report completeness; no LongPhase rerun.
- HG008 SNV F1 gain +0.15–0.21 at 0.2–0.8 for both callers; no gain in pure tumor (ClairS-TO 0.487→0.478, DeepSomatic-TO 0.351→0.377). Indel F1 ≤ 0.13 for every method. Recall drop at t40_n10 (0.85→0.73) under investigation.
- HG008 is the only genome not used in development.

## The three options

| | A. Add in place | B. HG008 in its own subsection | C. One subsection per task, truth first (proposed) |
|---|---|---|---|
| Results subsections | Keep the 5. Phasing accuracy goes into subsection 2. HG008 calls go into subsection 3. | 1 Overview, 2 Phasing and LOH (cell lines), 3 Recalibration (cell lines), 4 TDF and ploidy (cell lines), 5 PacBio, 6 HG008: phasing accuracy, somatic calls, CN and LOH against NIST truth | 1 Overview, 2 Phasing: HG008 SER and phase-correct N50 first, then cell-line contiguity, 3 LOH (separate), 4 Recalibration: ONT cell lines, HG008 and PacBio in one subsection, 5 TDF, ploidy and CN |
| Advantage | Smallest edit | Independent-genome result in one place | Each task starts with truth-based evidence; subsection 2 answers the editor's switch-error comment first |
| Disadvantage | Truth result sits in the middle of contiguity text | Evidence per task split in two places; subsection 4 loses its HG008 paragraph | Largest edit; PacBio no longer a separate subsection |
| Figures | Fig 2: add HG008 accuracy panel. Fig 3: add HG008 column. | Possibly new Fig 5, or the same panels as A | Same as A |

Discussion changes (all options): paragraph 1 adds truth-based phasing accuracy and states the N50 gain as phase-correct N50; near-pure cost paragraph adds HG008 pure tumor; limitations paragraph: HG008 is the independent genome for three analyses, still one cell line.

Editor's comment that triggered this: "no switch error rate is reported, so it is unclear whether the longer blocks are also correct"; earlier rejection reason: "advance too modest".
