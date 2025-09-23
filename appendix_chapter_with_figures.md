# Appendix

## Supplementary Methods for Somatic Variant Analysis

This appendix provides supplementary details on the computational methodologies developed and employed for the detection of somatic variants and structural alterations from long-read sequencing data. The sections below describe the algorithms for Loss of Heterozygosity (LOH) detection, the rationale and implementation of pattern-based filtering for somatic variant calls, the integrated workflow for variant analysis, and a case study demonstrating the application of these methods.

### Methodology for Loss of Heterozygosity (LOH) Detection

Loss of Heterozygosity (LOH) is a common genomic event in cancer wherein one parental allele of a chromosome segment is lost. To detect LOH regions, we developed a two-step classification method based on the local density of heterozygous variants. The approach first assigns a genotype to each variant based on its allele frequency and then calculates a "Heterozygosity Ratio" for genomic windows to classify them as LOH or Non-LOH (@fig:app-page-37-cropped-jpg).

![A two-panel plot illustrating the methodology for Loss of Heterozygosity (LOH) detection. The left violin plot shows the distribution of Variant Allele Frequencies (VAF) in Non-LOH (blue) and LOH (green) regions, demonstrating that LOH regions lack the characteristic heterozygous peak at VAF=0.5. The right box plot shows the distribution of the calculated Heterozygosity Ratio for the two region types, revealing a clear separation that allows for robust classification using a threshold (Hr = 0.09).](page_37_cropped.jpg){#fig:app-page-37-cropped-jpg}


***
**Figure: Variant Allele Frequency Distributions in LOH and Non-LOH Regions.** The violin plot illustrates the distribution of Variant Allele Frequencies (VAF) for genetic variants categorized by their location in pre-defined LOH (green) or Non-LOH (blue) regions. In Non-LOH regions, the distribution is bimodal, with a prominent peak centered at a VAF of 0.5, characteristic of heterozygous variants, and a smaller peak near 1.0 for homozygous variants. In contrast, LOH regions exhibit a unimodal distribution strongly skewed towards a VAF of 1.0, reflecting the loss of one allele and the consequent absence of heterozygous variants. The dashed line at a VAF of 0.8 indicates the threshold used to classify variants as homozygous in the subsequent analysis.
***
**Figure: LOH Classification by Heterozygosity Ratio.** This box plot displays the distribution of the Heterozygosity Ratio for regions classified as LOH versus Non-LOH. The ratio, which quantifies the proportion of heterozygous variants in a genomic window, shows a clear and distinct separation between the two classes. Non-LOH regions consistently exhibit a high ratio with a median around 0.55, whereas LOH regions have a ratio approaching zero. The dashed line at a ratio of 0.09 represents the classification threshold, **H<sub>r</sub>**, below which a region is classified as LOH.
***

The classification logic is formalized by the following equations. The process begins by categorizing variants within a given genomic window as either homozygous (Hom) or heterozygous (Het). This categorization is determined by the observed Variant Allele Frequency (VAF) and any pre-existing genotype information. A variant is classified as homozygous if its VAF is greater than or equal to 0.8 or if it was previously genotyped as such.

$$
\begin{cases}
  \text{Hom count ++,} & \text{if (VAF} \ge \text{0.8) or (Genotype == Hom)} \\
  \text{Het count ++,} & \text{otherwise}
\end{cases}
$$

Next, the Heterozygosity Ratio ($H_r$) is calculated as the proportion of heterozygous variants relative to the total number of informative (heterozygous and homozygous) variants in the window.

$$
H_r (\text{Heterozygosity ratio}) = \frac{\text{Het Count}}{\text{Het Count} + \text{Hom Count}}
$$

Finally, a genomic window is classified as LOH if its $H_r$ falls below a data-driven threshold of 0.09.

$$
\text{LOH type} = \begin{cases}
  \text{LOH,} & \text{if } H_r < 0.09 \\
  \text{Non-LOH,} & \text{otherwise}
\end{cases}
$$

This quantitative approach provides a robust and automated method for identifying LOH events across the genome by leveraging the clear statistical signature of heterozygote depletion.

### Pattern-Based Filtering for Somatic Variant Calling

Alongside the detection of large-scale genomic events, the accurate identification of individual somatic variants is a critical objective that is frequently challenged by sequencing artifacts and other sources of noise. To address this challenge, we developed a filtering strategy that evaluates the specific patterns of evidence supporting each candidate variant. This approach is based on the rationale that true somatic mutations exhibit characteristic, high-confidence evidence patterns that differ from those of artifacts.

#### Rationale for Pattern-Based Evidence Selection

The fundamental premise of this approach is that the evidence supporting a variant call is not of uniform quality. By categorizing candidate variants based on predefined patterns derived from sequencing read characteristics, such as read alignment features, it becomes possible to distinguish high-confidence calls from likely false positives. The performance of these patterns was evaluated by comparing their relative abundance in sets of known true positive (TP) and false positive (FP) variant calls (@fig:app-page-38-cropped-jpg).

![A bar chart comparing the percentage of True Positive (TP, blue) and False Positive (FP, red) variant calls attributed to different evidence patterns. The patterns are grouped into high-confidence (V_H), low-confidence (V_L), and non-somatic (V_N) categories, highlighting that V_H patterns are highly enriched for TPs with minimal FPs, while the 'DISAGREE' pattern is a major source of FPs, justifying a pattern-based filtering strategy.](page_38_cropped.jpg){#fig:app-page-38-cropped-jpg}


***
**Figure: Comparison of True Positive and False Positive Rates Across Variant Evidence Patterns.** The bar chart displays the percentage of all true positive (TP, blue) and false positive (FP, red) variant calls that are associated with various predefined evidence patterns (x-axis). The analysis reveals three distinct categories: **V<sub>H</sub> (High)** patterns, which are highly enriched for TPs with a minimal contribution to FPs; **V<sub>L</sub> (Low)** patterns, which are less frequent but still reliable indicators of true variants; and **V<sub>N</sub> (NonSomatic)** patterns, dominated by the "DISAGREE" category, which is a major source of both TPs and FPs and is thus highly ambiguous and unreliable. This analysis provides a strong rationale for a filtering strategy that prioritizes variants supported by V<sub>H</sub> and V<sub>L</sub> patterns.
***

#### Pattern Artifact Filtering using a Somatic Path Ratio

To formalize the pattern-based filtering, we introduced the "Somatic path ratio" (also referred to as VOTE_RATIO), a metric calculated from a graph-based representation of sequencing reads. This ratio quantifies the proportion of evidence supporting the putative somatic variant path relative to all other paths at that genomic locus (@fig:app-page-39-cropped-jpg).

![Distribution and robustness of the Somatic Path Ratio (VOTE_RATIO) for artifact filtering. The histograms show that true positives (TP, blue) consistently have a high ratio near 1.0, whereas false positives (FP, red) have a lower, more dispersed distribution, enabling effective filtering with a threshold (τ ≈ 0.8). Additional plots demonstrate this separation is maintained across a wide range of tumor purities (1.0 to 0.2).](page_39_cropped.jpg){#fig:app-page-39-cropped-jpg}


$$
\text{Somatic path ratio} > \tau
$$

For instance, if the evidence supporting the somatic path consists of two reads, while the total evidence for all paths at the locus is three (e.g., comprising one reference read and two somatic reads), the ratio is calculated as follows:

$$
\frac{2}{1 + 0 + 0 + 2} = \frac{2}{3} \approx 0.67
$$

If this value falls below the empirically determined threshold $\tau$ (e.g., 0.8), the variant is filtered out as a likely artifact.

***
**Figure: Distribution and Robustness of the Somatic Path Ratio.** **(a)** A conceptual diagram illustrates how a genomic locus with a C->T somatic mutation is represented as a graph. Evidence is partitioned into high-confidence somatic paths (V<sub>H</sub>) and non-somatic paths (V<sub>N</sub>). **(b)** Histograms show the distribution of the Somatic Path Ratio (VOTE_RATIO) for true positive (TP, blue) and false positive (FP, red) variants. The vast majority of TPs exhibit a ratio close to 1.0, while FPs have a more dispersed distribution. A threshold, $\tau \approx 0.8$ (dashed line), effectively separates the two populations. **(c)** Paired histograms of the VOTE_RATIO for TPs and FPs are shown across a range of tumor purity levels (1.0 down to 0.2). The analysis demonstrates that while the TP distribution shifts slightly leftward with decreasing purity, it remains well-separated from the FP distribution, confirming the metric's robustness for filtering across diverse clinical samples.
***

#### Low-Confidence Variant Detection using Vlow Ratio

For variants supported by less definitive evidence patterns (categorized as $V_L$), a complementary filtering metric, the "Vlow Ratio," was developed. This ratio is equivalent to the variant allele frequency calculated from the evidence graph, representing the proportion of evidence supporting the low-confidence variant path ($V_L$) relative to the total evidence from both the variant and non-variant ($V_N$) paths.

$$
\frac{V_L}{V_L + V_N} \ge \theta
$$

A candidate variant is accepted if its Vlow Ratio exceeds a minimum threshold, $\theta$, empirically set to approximately 0.2.

***
**Figure: Distribution and Robustness of the Vlow Ratio.** **(a)** Histograms show the distribution of the Vlow Ratio for true positive (TP, blue) and false positive (FP, red) variants. The majority of FPs have a ratio near 0.0, while TPs are more broadly distributed. This allows for effective filtering using a threshold $\theta \approx 0.2$ (dashed line), which removes a large fraction of FPs while retaining most TPs. **(b)** Paired histograms of the Vlow Ratio are shown across varying tumor purity levels (1.0 down to 0.2). As purity decreases, the TP distribution shifts toward lower ratios. However, the FP distribution remains concentrated at zero, demonstrating that the threshold remains effective at filtering artifacts even in low-purity samples.
***

## Integrated Workflow and Application

The methodologies described above are integrated into a comprehensive computational workflow designed to leverage the unique advantages of long-read sequencing for somatic variant analysis. This pipeline enables not only accurate variant calling but also haplotype phasing, which is crucial for understanding the clonal architecture of tumors.

### Principle of Purity-Informed Haplotype Phasing

A key principle underlying our approach is the differential behavior of germline and somatic variant allele frequencies in mixed tumor-normal samples. The allele frequency of a heterozygous somatic variant is directly proportional to tumor purity, whereas the allele frequency of a heterozygous germline variant remains stable around 0.5 regardless of purity. This statistical difference allows for the deconvolution of somatic and germline signals, which in turn enables the phasing of somatic mutations onto their correct germline haplotype backbone (@fig:app-page-41-cropped-jpg).

![A diagram illustrating how tumor purity affects allele frequencies and enables haplotype phasing. As purity decreases, the allele frequency of somatic variants (orange curve) shifts toward zero, while that of germline variants (blue curve) remains at 0.5. This separation allows for the resolution of genomic patterns, showing how a somatic T allele can be linked to a specific germline G allele, a task impossible at 100% purity where the signals overlap.](page_41_cropped.jpg){#fig:app-page-41-cropped-jpg}


***
**Figure: The Effect of Tumor Purity on Allele Frequency Distributions and its Utility for Resolving Variant Haplotypes.** This diagram illustrates the relationship between tumor purity and variant allele frequencies (VAFs). **(Left)** Normal cells contain only heterozygous germline variants (blue), while tumor cells contain both germline and a new heterozygous somatic variant (orange). **(Center)** As tumor purity decreases from 1.0 to 0.2, the VAF of the somatic variant (orange curve) decreases proportionally, while the VAF of the germline variant (blue curve) remains stable at 0.5. **(Right)** At 100% purity, the VAF distributions overlap, obscuring the underlying genomic pattern. At lower purities (e.g., 0.6), the distinct VAFs allow for the statistical separation of reads, enabling the reconstruction of the haplotype and revealing that the somatic T allele is physically linked to the germline G allele.
***

### Integrated Iterative Refinement Workflow for Somatic Variant Calling

The complete computational pipeline integrates LOH detection, pattern-based filtering, and purity-informed phasing into a cohesive workflow. A key feature is an iterative refinement loop specifically designed to handle high-purity samples, where standard phasing methods may fail (@fig:app-page-42-cropped-jpg). In such cases, information from LOH detection is used to inform and improve the phasing of variants, leading to more accurate downstream analysis and tumor purity estimation.

![A schematic of the integrated computational workflow for somatic variant analysis using long-read sequencing. The pipeline shows the flow from candidate variant generation to phasing. It highlights a critical iterative refinement loop where Loss of Heterozygosity (LOH) detection is used to improve phasing accuracy in high-purity samples, ultimately leading to a more precise tumor purity estimation.](page_42_cropped.jpg){#fig:app-page-42-cropped-jpg}


***
**Figure: A Comprehensive Workflow for Somatic Variant Calling and Phasing.** This diagram outlines the integrated pipeline, starting from long-read sequencing data containing a mixture of reference (grey), germline (blue), somatic (orange), and non-somatic/error (yellow) alleles. The workflow proceeds through candidate generation and filtering to a central phasing step, which assigns variants to their respective haplotypes (HP1-1, HP2). A parallel path for high-purity samples performs LOH detection, and this information is fed back in an iterative refinement loop to improve the accuracy of phasing and subsequent tumor purity prediction.
***

### Exemplar Case Study: Haplotagging of Somatic Variants and LOH

The efficacy of the integrated workflow is demonstrated through its application to long-read sequencing data from a representative tumor sample. The process of assigning variants and reads to their specific haplotype of origin is termed "haplotagging." The visualization below shows the successful haplotagging of both small somatic variants and a large-scale LOH event (@fig:app-page-43-cropped-jpg).

![A genome browser visualization demonstrating the successful haplotagging of somatic variants and a large structural event. Long sequencing reads are separated by haplotype (HP1 and HP2). The image clearly shows a somatic variant present only on HP1 reads, followed by a large region of Loss of Heterozygosity (LOH) where HP1 reads are absent. Within the LOH region, a second somatic variant is identified on the remaining HP2 reads, showcasing the pipeline's ability to resolve complex cancer genomes.](page_43_cropped.jpg){#fig:app-page-43-cropped-jpg}


***
**Figure: Genome Browser Visualization Demonstrating Successful Haplotagging of Phased Somatic Variants and a Large-Scale Loss of Heterozygosity Event.** This image displays long sequencing reads aligned to a reference genome, with reads computationally separated and colored by their inferred haplotype of origin (e.g., HP1, HP2, and refined sub-haplotypes HP1-1, HP2-1). The visualization clearly shows: **(1)** A heterozygous somatic variant (left, red ticks) present exclusively on reads assigned to Haplotype 1 (HP1). **(2)** A large genomic region of Loss of Heterozygosity (LOH, right), identified by the near-complete absence of reads corresponding to Haplotype 1. **(3)** A second somatic variant (right, green ticks) that is present on the remaining Haplotype 2 reads within the LOH region. This result validates the ability of the pipeline to resolve complex genomic architectures by accurately phasing variants and identifying concurrent structural changes.
***