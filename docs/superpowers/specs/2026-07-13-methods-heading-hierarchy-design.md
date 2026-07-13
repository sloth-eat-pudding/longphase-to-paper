# Methods heading specificity and consistency design

## Goal

Make the Methods headings in `paper.tex` and `supplementary.tex` immediately informative, specific to the operations described beneath them, and consistent with Nature-style scientific writing, while preserving the journal template and technical content.

## Design

Use concise, sentence-case noun phrases throughout the main and Supplementary Methods. Each heading must identify the methodological operation, evidence type, or model described beneath it rather than use generic labels such as "methodology," "principle," "scoring," or "predictive modeling." Keep terminology parallel across the two files when they describe the same procedure.

Within somatic variant recalibration, retain candidate selection and triplet-graph recalibration as distinct components. Rename the three triplet-graph run-in blocks to distinguish graph construction and pattern classification, regression-based scoring of individual triplets, and candidate-level integration of evidence. Use the same conceptual sequence in the Supplementary Methods.

For tumor DNA fraction estimation, replace generic model language with headings that name the actual predictors: haplotype imbalance and LOH. Rename other Methods headings where their current wording is overly broad, indirect, or inconsistent with the paragraph content, including LOH detection, haplotype reconstruction, datasets, and benchmark design.

In `supplementary.tex`, revise the headings within `Supplementary Methods` but leave `Supplementary Results`, figures, tables, and `Commandline arguments` outside this editorial scope. Correct only unambiguous content inconsistencies exposed by the heading review, including the duplicated "LOH regions" subject in the paragraph contrasting heterozygosity ratios in LOH and non-LOH regions. Do not change `sn-jnl.cls`, labels, cross-references, mathematical definitions, thresholds, or scientific claims.

## Verification

Compile `supplementary.tex` before `paper.tex`, reject undefined-reference or missing-character warnings, extract the rendered headings, and inspect the affected Methods pages for ambiguous spacing, orphaned headings, overlap, and clipping. Confirm that the source contains no title-case drift within the revised Methods headings and that parallel procedures use parallel terminology across both files.
