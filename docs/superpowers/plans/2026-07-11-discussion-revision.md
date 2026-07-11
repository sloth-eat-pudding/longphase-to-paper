# Discussion Revision Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current seven-paragraph Discussion in `paper.tex` with the approved six-paragraph, evidence-led Nature-style argument.

**Architecture:** The revision treats LOH-aware somatic haplotyping as the central conceptual contribution and organizes all other findings around it. The six paragraphs move from advance, to methodological positioning, to bounded implications, to estimand definition, to mechanistic and validation limitations, and finally to a concise conclusion.

**Tech Stack:** LaTeX (`paper.tex`), BibTeX (`references.bib`), TeX Live (`pdflatex`, `bibtex`), repository-prescribed Supplementary-first build order.

## Global Constraints

- Edit only the Discussion section of `paper.tex`.
- Preserve **tumor DNA fraction** versus **cellular purity** terminology.
- Preserve reported results and do not introduce new analyses.
- Use only defensible claims and citations.
- Do not edit `thesis.tex`, vendor files, figures, or Supplementary text.
- Compile `Supplementary` before `paper` and inspect logs for undefined references and missing characters.

---

### Task 1: Replace the opening with the principal advance

**Files:**
- Modify: `paper.tex:124-126`

**Interfaces:**
- Consumes: Results reported in `paper.tex:89-119`.
- Produces: Discussion paragraph 1, which defines the central advance used by paragraphs 2-6.

- [ ] **Step 1: Replace the two generic summary paragraphs with one evidence-led paragraph**

Use a direct problem-advance-evidence sequence. Retain the distinction between biological concordance in HCC1395 and technical reproducibility in the other datasets.

- [ ] **Step 2: Check against the Results**

Run: `rg -n "longer and less fragmented|CytoScan HD|technical reproducibility|tumor DNA fraction" paper.tex`

Expected: every summarized claim has a corresponding Results statement.

### Task 2: Add explicit methodological positioning

**Files:**
- Modify: `paper.tex` Discussion paragraph 2.

**Interfaces:**
- Consumes: The central haplotype representation defined in paragraph 1.
- Produces: A comparison with germline phasers and tumor-only callers that motivates paragraph 3.

- [ ] **Step 1: Explain the distinct assumptions of the methods**

Contrast conventional diploid germline phasing with lineage-aware somatic phasing across LOH. Explain why local haplotype evidence complements ClairS-TO and DeepSomatic-TO rather than replacing their caller-specific models.

- [ ] **Step 2: Use existing citations**

Use `martin2016`, `edge2017`, and `lin2022` for germline phasers and `park2024` and `chen2025` for the tumor-only somatic callers.

- [ ] **Step 3: Confirm citation keys**

Run: `rg -n "^@article\{(martin2016|edge2017|lin2022|park2024|chen2025)," references.bib`

Expected: five matching entries.

### Task 3: Narrow the implications to demonstrated scope

**Files:**
- Modify: `paper.tex` Discussion paragraph 3.

**Interfaces:**
- Consumes: The methodological distinction established in paragraph 2.
- Produces: A bounded statement of biological and practical utility before limitations.

- [ ] **Step 1: State the demonstrated value**

Describe parental-lineage assignment, allelic context, preserved phase across LOH, and potential utility when a matched normal is unavailable.

- [ ] **Step 2: State what was not demonstrated**

Explicitly limit the evidence to cell-line dilution benchmarks and avoid claims of diagnostic, prognostic, primary-tumor, or clinical utility.

### Task 4: Define tumor DNA fraction versus cellular purity

**Files:**
- Modify: `paper.tex` Discussion paragraph 4.

**Interfaces:**
- Consumes: The bounded tumor-only application described in paragraph 3.
- Produces: The estimand definition needed to interpret HCC1954 in paragraph 5.

- [ ] **Step 1: Define the two quantities**

Identify LongPhase-TO's output as the proportion of sequenced DNA molecules of tumor origin. Contrast it with cellular purity estimated jointly with allele-specific copy number and ploidy.

- [ ] **Step 2: Retain appropriate citations**

Use `vanloo2010`, `priestley2019`, and `nikzainal2012` for ASCAT, PURPLE, and Battenberg.

- [ ] **Step 3: Remove the repeated future-work sentence**

Reserve the copy-number-integration recommendation for paragraph 5.

### Task 5: Integrate mechanistic and validation limitations

**Files:**
- Modify: `paper.tex` Discussion paragraph 5.

**Interfaces:**
- Consumes: The tumor-DNA-fraction definition in paragraph 4.
- Produces: A mechanistic boundary and validation boundary that support the final conclusion.

- [ ] **Step 1: Explain the HCC1954 boundary case**

Connect whole-genome duplication and low true-LOH burden to weak haplotype-imbalance signal and compression of estimates. Retain `ghandi2019` and `bielski2018`, but remove the unsupported assertion about the exact internal cause of ASCAT's fit.

- [ ] **Step 2: State the practical interpretation**

Recommend interpreting the DNA fraction estimate together with LOH burden and identify allele-specific copy-number modeling as the principled extension for genome-doubled, low-LOH samples.

- [ ] **Step 3: State the validation limitation**

Clarify that HCC1395 is the only orthogonally validated biological case, whereas other datasets establish contiguity and technical reproducibility. Identify primary tumors, orthogonal allele-specific references, and controlled simulations as immediate validation needs.

### Task 6: Replace the generic future-work paragraph with a focused conclusion

**Files:**
- Modify: `paper.tex` Discussion paragraph 6.

**Interfaces:**
- Consumes: The advance and boundaries established in paragraphs 1-5.
- Produces: A concise final statement of present contribution and immediate next step.

- [ ] **Step 1: Write the final paragraph**

Conclude that tumor-only long reads support LOH-aware somatic phasing, call refinement, and read-level tumor DNA fraction estimation in one framework. State that primary-tumor validation with orthogonal references is needed. Remove speculative structural-variant, multi-omics, and cancer-treatment claims.

### Task 7: Verify source, citations, and build

**Files:**
- Verify: `paper.tex`
- Verify: `references.bib`
- Generated: `Supplementary.log`, `paper.log`, `Supplementary.pdf`, `paper.pdf`

**Interfaces:**
- Consumes: Revised six-paragraph Discussion.
- Produces: Compile-backed evidence that the manuscript remains valid.

- [ ] **Step 1: Confirm the Discussion contains exactly six paragraphs**

Inspect `paper.tex` between `\\section{Discussion}` and `\\section{Methods}`.

- [ ] **Step 2: Check terminology and generic language**

Run: `sed -n '/\\section{Discussion}/,/\\section{Methods}/p' paper.tex`

Expected: no display-text use of "tumor purity" or "tumour" and no generic closing claims about transforming cancer treatment.

- [ ] **Step 3: Build Supplementary**

Run: `pdflatex Supplementary && bibtex Supplementary && pdflatex Supplementary && pdflatex Supplementary`

Expected: exit code 0.

- [ ] **Step 4: Build paper**

Run: `pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper`

Expected: exit code 0.

- [ ] **Step 5: Inspect logs**

Run: `rg -n -i "undefined|Missing character|Citation.*undefined|Reference.*undefined" Supplementary.log paper.log`

Expected: no matches.

- [ ] **Step 6: Confirm output files**

Run: `pdfinfo Supplementary.pdf | rg "Pages|File size"` and `pdfinfo paper.pdf | rg "Pages|File size"`

Expected: both PDFs exist, have non-zero size, and report valid page counts.
