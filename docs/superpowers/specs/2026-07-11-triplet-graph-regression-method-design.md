# Triplet-Graph Regression Method Rewrite Design

## Objective

Revise the description of LongPhase-TO somatic-variant recalibration so that the Results, main Methods, Supplementary Methods, and retained supplementary figures describe the current released implementation consistently. The rewrite must explain the regression upgrade without incorrectly implying that regression replaces all site-level evidence aggregation.

## Authoritative interpretation

The current implementation uses a two-level classifier:

1. Each local triplet graph is evaluated against predefined high-confidence somatic patterns. For a matching pattern, three read-supported path counts are scored with a mutation-type-specific logistic model: a linear predictor for indels and a cubic polynomial predictor for SNVs. A high-confidence pattern is accepted when its sigmoid probability is at least 0.5.
2. If no high-confidence pattern is accepted, the graph is evaluated against refined low-confidence patterns. Accepted high- and low-confidence patterns are then aggregated across all candidate-containing triplets to make the site-level decision using the existing somatic-path and low-confidence evidence ratios.

This interpretation is supported jointly by the student slides, the accompanying description, and the released `PhasingGraph.cpp`. The released implementation takes precedence where the student prose loosely describes the path counts as pooled across all 12 high-confidence patterns or describes regression as the sole final site classifier.

## Manuscript changes

### `paper.tex`

- Revise the Results overview so it says that LongPhase-TO extracts triplet-graph path counts, applies separate SNV and indel logistic classifiers to high-confidence patterns, and aggregates accepted high- and low-confidence evidence across triplets. Avoid saying that a regression model directly makes the entire site-level decision.
- Replace the current pattern-mining and final-decision paragraphs with a reproducible description of:
  - the eight possible reference/alternate triplet paths;
  - the 12 high-confidence patterns spanning left, middle, and right candidate placements;
  - the three model inputs `x`, `y`, and `z`, defined as the read counts for the two germline paths and the somatic-descendant path in a candidate pattern;
  - minimum path-support screening before regression;
  - the fixed indel linear predictor and SNV cubic polynomial predictor;
  - sigmoid transformation and the `0.5` graph-level threshold;
  - the nine refined low-confidence patterns and rejection of third paths that traverse the candidate somatic allele;
  - aggregation over candidate-containing triplets and the retained site-level thresholds of `0.8` for mean somatic-path support and `0.2` for the low-confidence evidence fraction.
- Use “logistic regression” or “logistic classifier,” not “linear regression followed by a sigmoid.”
- Do not claim a training cohort, cross-validation design, held-out test scheme, or feature normalization procedure because these details are absent from the supplied materials.

### `supplementary.tex`

- Replace the obsolete pattern-filtering prose with a concise supplementary explanation of the upgraded graph-level regression and the retained site-level evidence aggregation.
- Retain the empirical supplementary figures supporting the `0.8` somatic-path and `0.2` low-confidence thresholds, because those thresholds remain part of the released implementation.
- Remove the two obsolete figure environments labeled `fig:met-page-22-cropped-jpg` and `fig:met-page-23-cropped-jpg`. Their graphics explicitly show the retired pre-regression pattern-artifact rule and the incomplete `V_H >= 1` formula and would contradict the revised text.
- Retain the triplet-graph and combinatorial-pattern figures labeled `fig:met-page-20-cropped-jpg` and `fig:met-page-21-cropped-jpg`, which remain conceptually correct.

## Scientific and editorial constraints

- Preserve “triplet graph,” “tumor DNA fraction,” and all label identifiers except the two obsolete supplementary-figure labels explicitly removed above.
- Distinguish graph-level pattern acceptance from candidate-site classification throughout.
- Use `V_H`, `V_L`, and `V_N` only for high-confidence, low-confidence, and non-somatic evidence categories; do not call these probabilities.
- Describe the coefficients as fixed model parameters in the released implementation. Avoid unsupported statements about generalization, optimality, class-imbalance correction, or statistical significance.
- Keep the prose compact and procedural, consistent with a Nature Methods Methods section.

## Verification

1. Search `paper.tex` and `supplementary.tex` for the retired formula and obsolete descriptions.
2. Compile `supplementary.tex` first, then `paper.tex`, using the repository-prescribed BibTeX and repeated-pdflatex order.
3. Confirm both builds exit successfully.
4. Search both logs for `undefined` and `Missing character`.
5. Inspect the final diff to confirm that `thesis.tex`, vendor files, figure sources, and unrelated user changes were untouched.
