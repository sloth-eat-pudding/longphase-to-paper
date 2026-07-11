# Discussion Structure Revision Design

## Objective

Revise the Discussion in `paper.tex` to achieve a concise, evidence-led Nature-style argument. Preserve technically defensible claims and appropriate citations, qualify or remove claims that exceed the evidence, and add citations only when they materially support interpretation or positioning.

## Scope

- Edit only the Discussion section of `paper.tex`.
- Preserve the manuscript's terminology, especially **tumor DNA fraction** versus **cellular purity**.
- Preserve the reported results and numerical findings; do not introduce new analyses.
- Do not edit `thesis.tex`, vendor files, figures, or the Supplementary text.
- Rebuild `Supplementary` before `paper` and inspect the logs for undefined references and missing characters.

## Approved Six-Paragraph Structure

### 1. Principal advance and key findings

Open with the unresolved problem and the central conceptual advance: LongPhase-TO reconstructs parental and somatic haplotypes from tumor-only long reads while remaining aware of chromosome-scale LOH. Present LOH reconstruction, somatic-call recalibration, and tumor DNA fraction estimation as consequences of this common haplotype-resolved representation rather than as a list of independent pipeline features. Retain the distinction between orthogonal biological concordance in HCC1395 and technical reproducibility elsewhere.

### 2. Positioning against existing approaches

Explain how LongPhase-TO differs from conventional germline phasers, tumor-only somatic callers, and allele-specific copy-number methods. Focus on the assumptions that distinguish the methods: diploid germline phasing, local haplotype evidence for somatic-call refinement, and copy-number-free estimation of read-level tumor content. Use existing bibliography entries where they directly support this positioning; add references only if an important comparison is otherwise unsupported.

### 3. Evidence-bounded implications

State the demonstrated biological and practical value of assigning somatic variants to parental haplotypes and maintaining phase across LOH. Describe tumor-only analysis as potentially useful when matched normal material is unavailable, while making clear that clinical, diagnostic, prognostic, clonal-evolution, and primary-tumor utility were not established in this study.

### 4. Tumor DNA fraction versus cellular purity

Define LongPhase-TO's output as the fraction of sequenced DNA molecules derived from tumor cells. Contrast it with ploidy-adjusted cellular purity inferred by ASCAT, PURPLE, and Battenberg. Explain concisely why the quantities can diverge under aneuploidy or whole-genome duplication. Reserve the copy-number-integration proposal for the following paragraph to avoid repetition.

### 5. Mechanistic and validation limitations

Use HCC1954 as the mechanistic boundary case: whole-genome duplication and low true-LOH burden weaken the haplotype-imbalance signal and compress estimates toward the middle of the dilution range. Retain the practical recommendation to interpret DNA fraction together with LOH burden, but avoid unsupported claims about the exact internal cause of ASCAT's fit. Then state the validation limitation: only HCC1395 has an orthogonal chromosome-scale LOH reference; other datasets demonstrate contiguity and reproducibility rather than biological accuracy. Identify primary tumors, additional allele-specific references, controlled simulations, and explicit copy-number modeling as the immediate validation and methodological priorities.

### 6. Focused conclusion

End with the present contribution: a single tumor-only long-read sample can support LOH-aware somatic phasing, haplotype-aware call refinement, and read-level tumor DNA fraction estimation within one framework. Avoid speculative lists involving multi-omics, additional structural-variant classes, or transformation of cancer treatment.

## Style Requirements

- Lead each paragraph with a substantive claim rather than a generic transition.
- Maintain one principal argumentative function per paragraph.
- Prefer mechanistic interpretation over adjectives such as "comprehensive", "significant", "promising", or "transformative".
- Avoid repeating Results-level numerical detail unless needed to establish a limitation.
- Use calibrated language that distinguishes demonstrated results from future applications.
- Preserve American spelling, LaTeX citation syntax, and repository cross-reference conventions.

## Verification

1. Compare the revised section against the Results to remove redundant restatement and prevent new claims.
2. Confirm all citation keys exist in `references.bib`.
3. Compile `Supplementary` first, then `paper`, using the repository-prescribed `pdflatex` and BibTeX sequence.
4. Search both logs for undefined references/citations and missing characters.
5. Re-read the compiled Discussion for paragraph-to-paragraph continuity and a concise final landing.
