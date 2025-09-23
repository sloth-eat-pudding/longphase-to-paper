# Methods

## Introduction and Overall Objective

The primary objective of this research is the development and application of **LongPhase-TO**, a comprehensive bioinformatics pipeline designed for the analysis of tumor-only samples using long-read sequencing data (@fig:met-page-10-cropped-jpg). The absence of a matched normal sample in many clinical and research settings poses significant analytical challenges, including the confident identification of somatic mutations and the characterization of tumor-specific genomic aberrations. LongPhase-TO is engineered to overcome these challenges by integrating multiple, haplotype-aware analyses into a single, cohesive workflow.

![Flowchart illustrating the overall objective of the LongPhase-TO pipeline. The process starts with a heterogeneous tumor-only sample, which undergoes long-read sequencing. The resulting data is processed by the LongPhase-TO pipeline, which performs three core analyses: chromosomal-scale LOH detection, somatic phasing to reconstruct germline and somatic haplotypes, and tumor purity prediction.](page_10_cropped.jpg){#fig:met-page-10-cropped-jpg}


The pipeline is structured to process raw long-read sequencing data from a tumor sample and produce a detailed genomic profile of the cancer. This process encompasses three core analytical modules:

1.  **Chromosomal-Scale Loss of Heterozygosity (LOH) Detection:** This module identifies large genomic regions where one of the two parental chromosome copies has been lost, a frequent and significant event in tumorigenesis. The detection is performed at a chromosomal scale, enabling the characterization of major structural alterations.

2.  **Somatic Phasing and Variant Calling:** This central component leverages the length of the sequencing reads to perform phasing—the assignment of genetic variants to their parental chromosome of origin (haplotype). The pipeline first reconstructs the patient's two germline haplotypes. Subsequently, it identifies and reconstructs novel somatic haplotypes that have emerged during tumor evolution through the accumulation of somatic mutations, thereby providing a haplotype-resolved view of the tumor's genetic lineage.

3.  **Tumor Purity Prediction:** This module estimates the proportion of cancerous cells within the bulk tumor sample. Tumor purity is a critical covariate for interpreting somatic variant allele frequencies and understanding the tumor microenvironment. The prediction is derived from the relative abundance of the reconstructed germline and somatic haplotypes.

By combining these three analyses, LongPhase-TO aims to provide a detailed and accurate genomic portrait of a tumor from a single sample, offering insights into its clonal heterogeneity, evolutionary history, and structural landscape.

## Overview of the Analytical Workflow

To achieve its objectives, the LongPhase-TO pipeline employs a structured, multi-step methodology designed to progressively refine raw sequencing data into a comprehensive genomic characterization. The workflow integrates several distinct analytical stages, each building upon the output of the previous one (@fig:met-page-11-cropped-jpg). The six primary steps are outlined below.

![Overview of the six-step analytical method. The workflow starts with (1) CNV/BFB interval calling, followed by (2) chromosomal-scale LOH detection via heterozygosity ratio, (3) somatic variant candidate selection using a Panel of Normals, (4) graph-based somatic variant calling, (5) phasing of germline and somatic haplotypes, and concludes with (6) tumor purity prediction based on haplotype frequencies.](page_11_cropped.jpg){#fig:met-page-11-cropped-jpg}


1.  **$\frac{CNV}{BFB}$ Interval Calling:** The analysis begins with a genome-wide screen for regions indicative of large-scale structural variations. By analyzing the alignment patterns of long reads, specifically the locations of soft-clipped sequences, this step identifies candidate genomic intervals likely containing Copy Number Variations (CNVs) or complex rearrangements resulting from Breakage-Fusion-Bridge (BFB) cycles.

2.  **Chromosomal-Scale LOH Detection:** Using the information from the initial screen, this step focuses on detecting large regions of Loss of Heterozygosity (LOH). The method quantifies the ratio of heterozygous to homozygous variants across the genome. A significant and sustained drop in this heterozygosity ratio over a large chromosomal segment is identified as an LOH event, signifying the loss of one parental allele.

3.  **Somatic Variant Candidate Selection:** To confidently identify tumor-specific mutations in the absence of a matched normal sample, this step utilizes a rigorous filtering strategy. An initial set of all detected variants is filtered against a Panel of Normals (PON)—a comprehensive database of common germline variants and recurrent sequencing artifacts. This subtractive process enriches the candidate pool for high-confidence somatic mutations by removing known polymorphisms.

4.  **Somatic Variant Calling:** The refined list of candidates is subjected to a sophisticated variant calling process. This step utilizes a novel graph-based approach that models local haplotype structures. By analyzing the linkage of variants on the same long reads, the method distinguishes true, low-frequency somatic variants, which form consistent new haplotypes, from sporadic sequencing errors.

5.  **Phasing:** This core step reconstructs the full-length haplotypes present in the sample. It assigns the identified somatic and germline variants to their respective chromosome of origin, resolving the two primary germline haplotypes (HP1 and HP2). Crucially, it also identifies and reconstructs any novel somatic haplotypes that have evolved from the germline state (e.g., HP1-1 derived from HP1). This process integrates information from the LOH detection step to correctly model haplotype composition in regions of allelic loss.

6.  **Tumor Purity Prediction:** In the final step, the pipeline leverages the quantitative output of the phasing module to estimate tumor purity. The relative frequencies of the reconstructed germline and somatic haplotypes are used as input for a predictive model. As tumor purity increases, the abundance of somatic haplotypes relative to their germline counterparts increases in a predictable manner, allowing for the calculation of an accurate purity score.

This integrated workflow ensures that each analytical component informs the others, leading to a robust and detailed characterization of the tumor genome from a single long-read sequencing experiment.

## CNV/BFB Interval Calling

The identification of genomic intervals affected by large structural variants, such as CNVs and BFB-induced rearrangements, is a critical first step in characterizing the tumor genome. The method accomplishes this by systematically analyzing soft-clipping patterns in long-read sequencing alignments. Soft-clipping occurs when a portion of a sequencing read does not align to the reference genome, often because the read spans a structural variant breakpoint.

### Principle of Clipping Pattern Analysis

The alignment of long reads across structural variant breakpoints generates characteristic and non-random patterns of soft-clipping. These patterns, encoded in the CIGAR (Compact Idiosyncratic Gapped Alignment Report) string of an alignment file, serve as powerful signatures for specific SV types. A CIGAR string such as `10S5M2D1M` indicates that the first 10 bases of the read were soft-clipped (`10S`), followed by 5 matching bases (`5M`), a 2-base deletion (`2D`), and another matching base (`1M`). Clipping is categorized into two types:
*   **Up-clipping:** Soft-clipping at the beginning of a read's alignment (e.g., `10S...`).
*   **Down-clipping:** Soft-clipping at the end of a read's alignment (e.g., `...10S`).

Distinct SVs produce unique spatial arrangements of these clipping types. For instance, a tandem duplication often presents as a sharp peak of up-clipped reads immediately followed by a sharp peak of down-clipped reads at the duplication boundaries. A fold-back inversion, a hallmark of BFB cycles, can produce a more complex signature of alternating and overlapping up- and down-clipping signals. The algorithm is designed to computationally detect and delineate these signature-rich intervals.

### Algorithmic Detection of Clipping Signatures

The formal identification of $\frac{CNV}{BFB}$ intervals is performed using a four-step signal processing algorithm that transforms raw alignment data into defined genomic regions.

#### Parsing and Quantifying Clipping Events

The algorithm first parses the alignment file ($e$.$g$., BAM format) for every soft-clipped read. For each clipped read, the genomic coordinate of the clipping breakpoint is recorded. To create a directional signal, up-clipping events are assigned a positive value ($e$.$g$., +1) and down-clipping events are assigned a negative value ($e$.$g$., -1). This process converts the spatial distribution of clipped reads into a quantitative scatter plot of "Clipping Count" versus genomic position, where each point represents the net clipping signal at a single base.

#### Signal Smoothing via Forward Pileup

The raw clipping count signal is often sparse and noisy. To enhance true signals and reduce noise, a smoothing step is applied using a sliding window summation technique termed "Forward Pileup." A window of a predefined size `w` slides across the genome, and at each position, the algorithm sums the clipping count values of all points within the window. For example, with a window size `$w=6$`, the pileup value is calculated as:

$$
\text{Pileup} = \sum_{i=1}^{w} \text{ClippingCount}_i
$$

A window containing `[1, 1, 1, 1, 1, 1]` would yield a strong positive pileup value of 6, while a window with mixed signals like `[-1, -1, -1, 1, -1, 5]` would yield a value of 0. This transformation converts the discrete scatter plot into a continuous signal profile, where sustained clipping events are amplified into distinct peaks and valleys (@fig:met-page-14-cropped-jpg).

![The "Forward Pileup" method for signal smoothing. The algorithm converts a sparse scatter plot of raw clipping counts into a continuous signal by summing counts within a sliding window. This process reduces noise and amplifies consistent clipping signals into identifiable peaks, which are then used as candidates for structural variant breakpoints.](page_14_cropped.jpg){#fig:met-page-14-cropped-jpg}


#### Advanced Signal Feature Detection

Simple peak detection on the smoothed profile is insufficient to capture all relevant SV signatures. Therefore, more sophisticated criteria are employed to identify high-confidence breakpoint candidates.

*   **Calling Gentle Clipping:** This criterion identifies regions characterized not by a single high-amplitude peak, but by a high *density* of consistent, low-level clipping events. It scans the genome for windows of size `w` containing a high number of clipping events of the same sign (e.g., all positive). Such regions are flagged as "Gentle Clipping Candidates," representing breakpoints that produce a sustained but less sharp signal.

*   **Amplify Gentle Clipping:** This method detects a powerful signature of many SV breakpoints: a rapid, localized transition from one clipping direction to the other (@fig:met-page-15-cropped-jpg). The algorithm searches the "Forward Pileup" profile for a local minimum (valley) immediately followed by a local maximum (peak). The magnitude of this "signal swing" is quantified by calculating the difference in their amplitudes. For a peak of -2 following a valley of -8, the amplitude is:

![Advanced methods for detecting structural variant breakpoints. The "Amplify Gentle Clipping" method identifies high-confidence breakpoints by searching for a rapid transition from a local minimum (valley) in the clipping pileup signal to a local maximum (peak). The amplitude of this signal swing is quantified to score the confidence of the breakpoint.](page_15_cropped.jpg){#fig:met-page-15-cropped-jpg}


$$
-2 - (-8) = 6
$$

A large amplitude signifies a strong, high-confidence transition indicative of a complex breakpoint.

#### Formal Interval Calling

The final step formalizes the detection of candidate regions using a set of quantitative rules. A significant event is identified if it meets the criteria for a "Signal Swing," which is defined by two conditions. A pair of signal points, `i` and `j`, with clipping pileup values `c_i` and `c_j` respectively, constitute a signal swing if their signs are opposite and they satisfy:

$$
(c_j / c_i \ge \alpha) \land (|j - i| \le \beta)
$$

Here, `α` is a threshold for the minimum signal amplitude ratio, ensuring a significant swing, and `β` is a threshold for the maximum genomic distance between the points, ensuring localization. An "Event Burst" can also be called based on a simpler amplitude threshold, `c > λ`, to identify isolated but strong peaks. The genomic coordinates of validated signal swings are used to define the boundaries of the final called CNV/BFB intervals (@fig:met-page-16-cropped-jpg).

![The final step of interval calling using "Signal Swings." The algorithm identifies pairs of opposing peaks and valleys in the processed clipping signal that meet criteria for amplitude and proximity. These signal swings pinpoint the precise locations of breakpoints, which are then used to define the start and end of the final called genomic intervals, as shown in the result schematic.](page_16_cropped.jpg){#fig:met-page-16-cropped-jpg}


## Chromosomal-Scale Loss of Heterozygosity Detection

Loss of Heterozygosity (LOH) is a hallmark of cancer genomes, reflecting the loss of one parental allele over a large genomic region. LOH events are detected at a chromosomal scale by integrating variant information with the previously identified CNV/BFB intervals.

### LOH Detection via Heterozygosity Ratio

The primary signal for LOH is a stark imbalance in allelic ratios at heterozygous variant sites. For each known heterozygous single nucleotide variant (SNV), a heterozygosity ratio is calculated based on the read counts supporting the heterozygous versus homozygous state. This is defined as:

$$
\text{Heterozygosity ratio} = \frac{\text{Het Count}}{\text{Het Count} + \text{Hom Count}}
$$

In a normal diploid region, this ratio is expected to be close to 1.0 (or a high value, accounting for mapping biases). In a region of LOH, where one allele is lost, the vast majority of reads will support the remaining allele, causing the heterozygosity ratio to drop to a value near zero. An LOH event is called if this ratio falls below a predefined threshold, `σ`:

$$
\text{Heterozygosity ratio} < \sigma
$$

Visually, this corresponds to the complete disappearance of long reads supporting one of the two haplotypes over a large genomic span.

### Merging LOH Segments for Chromosomal-Scale Analysis

A key challenge in LOH detection is that large LOH regions can be interrupted by small, localized genomic events such as tandem duplications (CNVs) or complex rearrangements from BFB cycles. These events can locally restore a heterozygous state or create complex copy number profiles, causing transient spikes in the heterozygosity ratio and fragmenting the LOH call.

To identify the true, contiguous chromosomal-scale event, the method incorporates a crucial filtering and merging step. The algorithm first makes an initial pass, labeling regions with a low heterozygosity ratio as "LOH" and identifying the interrupting CNV/BFB events based on the interval calls from the previous stage. It then systematically "skips" over these small, intervening CNV/BFB events, merging the flanking LOH segments (@fig:met-page-17-cropped-jpg). This logic allows the pipeline to reconstruct a single, continuous LOH event, providing a more biologically accurate representation of a major chromosomal alteration rather than a collection of disjointed, smaller events.

![Method for chromosomal-scale LOH detection. The workflow shows that an initial pass may identify fragmented LOH regions interrupted by small CNV or BFB events that locally restore heterozygosity. The final step involves a merging process that filters out these small intervening events to define a single, continuous, large-scale LOH region.](page_17_cropped.jpg){#fig:met-page-17-cropped-jpg}


## Somatic Variant Identification and Phasing

The accurate identification of somatic variants from tumor-only long-read data requires a multi-stage process that combines rigorous filtering, advanced haplotype-aware variant calling, and high-resolution phasing.

### Somatic Variant Candidate Selection

The initial list of variant candidates produced by a standard caller is a heterogeneous mixture of true somatic mutations, common germline variants, and technical artifacts. The first step is to enrich this list for true somatic events.

#### Filtering Against a Panel of Normals (PON)

In the absence of a matched normal sample, a Panel of Normals (PON) is used to filter out common germline variants. The PON is a comprehensive database constructed from a large cohort of healthy individuals (e.g., 1000 Genomes Project, gnomAD) and also includes recurrent, technology-specific artifacts. All variant candidates from the tumor sample that are present in the PON are removed in a subtractive process that effectively eliminates the majority of non-somatic variants (@fig:met-page-18-cropped-jpg).

![Filtering somatic variant candidates using a Panel of Normals (PON). This schematic shows that an initial set of variant candidates, containing a mix of somatic (brown), germline (blue), and other variants, is filtered by subtracting all variants found in the PON. This process enriches the final candidate list for true somatic mutations by removing common germline polymorphisms and known artifacts.](page_18_cropped.jpg){#fig:met-page-18-cropped-jpg}


#### Filtering of Clustered Variants via Allele Concordance

Systematic sequencing and alignment errors often manifest as clusters of false positive variants in close proximity. To identify and remove these artifacts, an "Allele Concordance" filter is employed. This method evaluates a candidate variant at position `i` by analyzing its co-occurrence with neighboring variants on the same long reads. For a candidate variant `i` and a nearby variant `j`, two concordance ratios are calculated:

$$
\frac{\text{AltAlt}_{ij}}{\text{AltAlt}_{ij} + \text{AltRef}_{ij}} \ge \mu
$$

$$
\frac{\text{AltAlt}_{ij}}{\text{AltAlt}_{ij} + \text{RefAlt}_{ij}} \ge \nu
$$

where `AltAlt_ij` is the number of reads with alternate alleles at both `i` and `j`, `AltRef_ij` is the count of reads with an alternate allele at `i` and reference at `j`, and `RefAlt_ij` is the count of reads with reference at `i` and alternate at `j`. The thresholds `μ` and `ν` control the stringency. If a candidate variant `i` shows high concordance with multiple neighbors (i.e., its aggregate `Allele Concordance` score exceeds a final threshold `ρ`), it is classified as a clustered artifact and removed from the candidate list (@fig:met-page-19-cropped-jpg).

![Filtering of clustered artifacts using Allele Concordance. The method analyzes the co-occurrence of a candidate variant (at position `i`) with its neighbors (at position `j`) on the same long reads. Variants that consistently appear together (high concordance) are likely systematic artifacts. The workflow illustrates how pairs of variants are tested and, if the aggregate concordance score is high, the candidate is filtered out.](page_19_cropped.jpg){#fig:met-page-19-cropped-jpg}


### Graph-Based Somatic Variant Calling

The filtered list of high-confidence candidates is then subjected to a novel, graph-based variant calling algorithm designed to leverage local haplotype information.

#### The Tri-nodal-edge Graph Model

The core of this method is the "Tri-nodal-edge Graph," a data structure that models the linkage of alleles across three adjacent variant positions. This provides a richer representation of the local haplotype context than simpler models that only consider pairs of variants (a "Bi-nodal-edge Graph"), which were found to be non-informative for resolving complex cases (@fig:met-page-20-cropped-jpg). In this graph, nodes represent the alleles (Reference or Alternate), and weighted edges represent the observed 3-variant haplotypes and their supporting read counts.

![Comparison of graph models for somatic variant calling. The slide contrasts the proposed "Tri-nodal-edge Graph," which models haplotypes across three adjacent variant sites, with a simpler "Bi-nodal-edge Graph" that only considers pairs. The tri-nodal graph captures more linkage information, making it more powerful for distinguishing true somatic variants from artifacts, while the bi-nodal graph is labeled "noninformative" for this complex task.](page_20_cropped.jpg){#fig:met-page-20-cropped-jpg}


#### Combinatorial Haplotype Pattern Classification

A variant is classified not in isolation, but by matching the combination of local haplotypes observed at its locus against a predefined catalog of patterns. These patterns fall into three categories:

*   **High-Confidence Somatic (V_H):** These are patterns derived from a "3-haplotype" scenario, where reads support a reference haplotype, a germline haplotype, and a somatic haplotype. There are 12 such specific combinations of the possible `C(8,3) = 56` 3-haplotype sets that are classified as `V_H` (@fig:met-page-21-cropped-jpg).
*   **Low-Confidence Somatic (V_L):** These patterns arise from simpler "2-haplotype" scenarios, typically involving a reference and a somatic haplotype. There are 3 specific combinations of the `C(8,2) = 28` 2-haplotype sets classified as `V_L`.
*   **Non-Somatic (V_N):** Any observed haplotype combination not matching a `V_H` or `V_L` pattern.

![Catalog of high-confidence somatic variant patterns. The figure displays the 12 specific tri-haplotype patterns that are classified as high-confidence somatic variants (V_H). A variant is called with high confidence if the combination of local haplotypes observed in the sequencing data matches one of these predefined graph structures, which typically represent a reference, a germline, and a somatic haplotype.](page_21_cropped.jpg){#fig:met-page-21-cropped-jpg}


#### Pattern Mining and Artifact Filtering

To ensure high specificity, the classification is refined through a two-stage filtering process. First, the algorithm identifies the most frequent local haplotypes (Top 3 for `V_H` candidates, Top 2 for `V_L` candidates). Second, it applies a ratio-based artifact filter.
*   For a `V_H` candidate, the ratio of read support for a minor path versus the main somatic path must not exceed a threshold `τ` (@fig:met-page-22-cropped-jpg). This prevents misclassification due to a noisy, low-support path.
*   For a `V_L` candidate, the ratio of read support for the 3rd most frequent path to the 2nd must be below a threshold `δ`. This ensures the site is truly a simple 2-haplotype case and not a misidentified complex site.

![Artifact filtering for high-confidence (V_H) variant patterns. After a candidate V_H pattern is identified from the top 3 most frequent haplotypes, a ratio-based filter is applied. This filter checks that the read support for minor, potentially erroneous paths is not significant compared to the main somatic haplotype path. If the ratio exceeds a threshold τ, the pattern is rejected as a likely artifact.](page_22_cropped.jpg){#fig:met-page-22-cropped-jpg}


#### Final Somatic Decision Formula

A final call for a variant site is made using a "voting" system that aggregates the evidence from all local patterns (@fig:met-page-23-cropped-jpg). A site is called as a somatic variant if it meets either of two conditions:

![The final decision-making process for somatic variant calling. This diagram illustrates how all local haplotype patterns at a site are classified as high-confidence (V_H), low-confidence (V_L), or non-somatic (V_N). These classifications are tallied in a "voting" system, and a final decision is made based on a formula that prioritizes V_H evidence or a significant proportion of V_L evidence.](page_23_cropped.jpg){#fig:met-page-23-cropped-jpg}


$$
V_H \ge 1 \quad \text{or} \quad \frac{V_L}{V_L + V_N} \ge \theta
$$

This hierarchical rule ensures high sensitivity by calling any site with strong (`V_H`) evidence, while maintaining specificity by requiring weaker (`V_L`) evidence to surpass a statistical threshold `θ`.

### Haplotype Phasing

The foundation for the graph-based variant calling is a robust phasing algorithm that reconstructs the full haplotypes present in the sample.

#### Graph-Based Haplotype Reconstruction

The method first constructs a `k-variant graph` (typically with `k=2`) where nodes are alleles at variant sites and directed edges connect alleles observed on the same read (@fig:met-page-24-cropped-jpg). Edge weights correspond to the number of supporting reads. The algorithm then finds the most heavily supported paths through this graph to reconstruct the underlying haplotypes. It starts by seeding the two most common alleles at the first position as the germline haplotypes (HP1 and HP2) and extends them greedily by following the edges with the highest weights.

![Graph-based haplotype phasing workflow. This diagram illustrates how sequencing read data is transformed into a k-variant graph, where nodes are alleles and weighted edges represent their co-occurrence on reads. The algorithm then traverses this graph to reconstruct the most likely haplotypes, including distinguishing between germline (HP1, HP2) and newly evolved somatic (HP2-1) haplotypes.](page_24_cropped.jpg){#fig:met-page-24-cropped-jpg}


#### Discovery of Somatic Haplotypes

Less frequent paths that diverge from the main germline tracks are identified as potential somatic haplotypes. For example, if HP2 splits into two paths, one with high read support and one with low support, the low-support path is reconstructed as a new somatic haplotype (e.g., HP2-1). This process effectively deconvolutes the mixture of reads into a clear set of distinct germline and somatic haplotypes, providing the essential input for the variant calling and purity prediction modules.

## Tumor Purity Prediction

Accurate estimation of tumor purity—the fraction of cancer cells in a sample—is essential for the correct interpretation of somatic variant data. Purity is predicted by modeling the genomic consequences of mixing tumor and normal cell populations.

### Principle of Haplotype Imbalance

The core principle of this method is that somatic copy number alterations in tumor cells create a measurable Haplotype Imbalance (HI) in the bulk sequencing data. At any heterozygous site, HI is defined as the fraction of reads supporting the more abundant of the two germline haplotypes:

$$
\text{HI} = \frac{\max(\text{Count}_{\text{HP1}}, \text{Count}_{\text{HP2}})}{\text{Count}_{\text{HP1}} + \text{Count}_{\text{HP2}}}
$$

In a pure normal sample, HI is expected to be 0.5 across the genome. In a tumor sample containing regions of LOH, the HI will shift towards 1.0. The genome-wide statistical distribution of HI values is therefore strongly correlated with tumor purity; as purity increases, the distribution systematically shifts towards higher values (@fig:met-page-25-cropped-jpg).

![Principle of tumor purity prediction using Haplotype Imbalance (HI). The figure shows that as tumor purity increases from 0.2 to 1.0, the proportion of somatic haplotypes increases. This causes the genome-wide distribution of HI, shown in histograms and box plots, to shift from a balanced state (centered at 0.5) to a highly imbalanced state with higher median values and a peak at 1.0 (LOH).](page_25_cropped.jpg){#fig:met-page-25-cropped-jpg}


### Predictive Modeling with Genomic Features

To formalize this relationship, a machine learning approach is used. Two key sets of features are extracted from the genomic data:

1.  **LOH Ratio:** The fraction of the genome affected by LOH, which serves as a measure of large-scale allelic loss. It is calculated as:
    $$
    \text{LOH Ratio} = \frac{L_{\text{LOH}}}{L_{\text{Genome}}}
    $$
    where `$L_{LOH}$` is the total length of all detected LOH regions and `$L_{Genome}$` is the total analyzable genome length.

2.  **Haplotype Imbalance Statistics:** Key statistical properties of the genome-wide HI distribution, specifically the first (Q1) and third (Q3) quartiles, which capture the central tendency and spread of the imbalance.

These features (LOH Ratio, Q1 of HI, Q3 of HI) are used as input for a **Polynomial Regression** model (@fig:met-page-26-cropped-jpg). The model is trained to learn the non-linear relationship between these genomic features and tumor purity. Once trained, the model can take these features from any new tumor-only sample and produce a continuous, quantitative prediction of its purity.

![Schematic of the tumor purity prediction model. The model takes two types of genomic features as input: the LOH Ratio (fraction of the genome with LOH) and statistics from the Haplotype Imbalance distribution (represented by the Q1 and Q3 quartiles of a box plot). These features are fed into a Polynomial Regression model, which outputs a continuous prediction of tumor purity.](page_26_cropped.jpg){#fig:met-page-26-cropped-jpg}


## Data and Computational Resources

### Data Sources and Benchmark Sets

This study utilizes publicly available long-read sequencing data from a panel of well-characterized cancer cell lines. The use of these standard reference materials ensures the reproducibility of our findings and allows for rigorous benchmarking against established ground-truth call sets. The cell lines and their respective data and benchmark sources are detailed below.

| Cell Line | Material Source | Benchmark Source |
|-----------|-----------------|------------------|
| COLO829   | ONT, NYGC       | NYGC             |
| H1437     | UCSC            | Google           |
| H2009     | UCSC            | Google           |
| HCC1395   | HKU, NYGC       | SEQC2            |
| HCC1937   | UCSC            | Google           |
| HCC1954   | UCSC            | Google           |

### Somatic Variant Calling Tools

To evaluate the performance of the pipeline, several state-of-the-art somatic variant callers designed for long-read or tumor-only analysis were used for comparison.

| Somatic Variant Callers           |
|-----------------------------------|
| ClairS-TO v0.3.0 (ssrs)             |
| ClairS-TO v0.3.0 (ss)               |
| DeepSomatic v1.8.0 (tumor-only)   |

### Generation of In-Silico Tumor-Normal Mixtures

To train and validate the tumor purity prediction model, a benchmark dataset with precisely known purity levels was required. This was generated using an *in-silico* titration approach. High-coverage long-read sequencing data were obtained from a pure tumor cell line and its matched normal cell line. Synthetic mixture samples were then created by downsampling and mixing the raw reads from the tumor and normal samples at predefined ratios (@fig:met-page-27-cropped-jpg). For a target total coverage of 50x, samples with purities ranging from 1.0 (50x tumor, 0x normal) to 0.2 (10x tumor, 40x normal) were generated. This procedure yielded a high-quality benchmark dataset with a ground-truth purity label for each sample, enabling robust training and unbiased evaluation of the regression model.

![Generation of in-silico tumor-normal mixtures for benchmarking. This diagram illustrates the method of creating synthetic samples with known tumor purity. Raw sequencing reads from a pure tumor sample (orange) and a pure normal sample (blue) are mixed in specific proportions to create samples with precise purity levels (e.g., 0.8 purity from 40x tumor and 10x normal reads), providing a ground-truth dataset for model training and validation.](page_27_cropped.jpg){#fig:met-page-27-cropped-jpg}
