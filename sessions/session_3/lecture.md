---
layout: page
permalink: sessions/session_3/lecture
menubar_toc: true
---

<script src="{{ site.baseurl }}/assets/js/vanilla-back-to-top.min.js"></script> <script>addBackToTop()</script>

## Introduction to NGS sequencing:

A typical Next Generation Sequencing (NGS) workflow can be divided into 5 steps:

Specimen collection - e.g. primary vs secondary tumors, tumors and/or normal tissues, liquid biopsy (blood, urine) or traditional biopsy, preservation techniques (FFPE, fresh frozen, circulating tumor cells, cell free DNA, etc.)
Sample extraction: e.g. genomic DNA or RNA
Library preparation: prepare samples to be compatible with the sequencer. Sequencing libraries are typically created by fragmenting DNA and adding specialized adapters to both ends. In the Illumina sequencing workflow, these adapters contain complementary sequences that allow the DNA fragments to bind to the flow cell. Fragments can then be amplified and purified.
Sequencing: libraries are loaded onto a flow cell and placed on the sequencer. Two common sequencing platforms are Illumina and Ion Torrent.
Data analysis & interpretation: After sequencing, the instrument software identifies nucleotides (a process called base calling) and the predicted accuracy of those base calls. The raw output from sequencers are often in FASTA format, which can be aligned to a reference genome in BAM format.


<img> From Qiagen


Sample collection, preservation, and manipulation are important pre-analytical factors to consider. See the overview of the most commonly used biopsy techniques, preservation methods, and genomic analytes below.

<img>

Traditional biopsy methods include fine- or core-needle biopsy or surgical resection. These biopsies typically only access the primary tumor site. From traditional tissue biopsy the most common pathological preservation path is through formalin fixation and paraffin embedding (FFPE), though fresh frozen tissue or disaggregated cells are sometimes also available. From each of these material types, both DNA and RNA can be extracted. Liquid biopsy usually involves blood draw, though some groups are now testing urine and other body fluids. Liquid biopsy can have representative somatic lesions from more than one tumor site. Circulating tumor cells (CTCs), cell-free DNA (cfDNA), and exosomes or extracellular vesicles (EVs) are the most common components of liquid biopsy that are isolated for somatic analysis. DNA and RNA can be isolated from CTCs, but only DNA is represented in the cfDNA extraction, and RNA is most commonly targeted from EVs

Common pre-analytical issues, impacts, and contingencies associated with different sample types are summarized in the Table below:

<img>
Lennon NJ, Adalsteinsson VA, Gabriel SB. Technological considerations for genome-guided diagnosis and management of cancer. Genome Med. 2016 Oct 26;8(1):112. doi: 10.1186/s13073-016-0370-4. PMID: 27784341; PMCID: PMC5080740.

All DNA extraction procedures aim to isolate the DNA from all other cell components. Steps for preprocessing, lysis will vary substantially depending on specimen type. Quality and quantity of extracted DNA will vary depending on specimen type, storage conditions, and extraction techniques.

—--

## Quality control - DNA extraction and library preparation

### Measuring DNA quality and quantity

After extraction, DNA quantity and quality should be assessed to decide if your sample is up to your specifications before you invest time and money on sequencing.
#### Quantitation Methods
Spectrophotometry-absorbance measurement:
Measures DNA concentration by absorbance of UV light
Accessible to many labs
Can overestimate quantity due to other molecules absorbing at the same wavelength as DNA

<img>
Fluorometry:
Measures DNA via binding dye which is fluorescent when bound to dsDNA
Not standard in all labs; but becoming more accessible
Binding dye is specific to dsDNA (this is important for many NGS applications)
Utilizes a standard curve
In general, more accurate

<img>
For both quality and quantity, requirements will vary depending on downstream sequencing platform and application.  Amplification-free approaches generally require more material.
#### Qualitative approaches
Electrophoresis
Slab Gel (older technology, most labs can do this)
Size distribution measured against a ladder with DNA fragments of known sizes to determine quality

<img>
BioAnalyzer, TapeStation, FemtoPulse
Automated electrophoresis, detection, analysis
Determine fragment size distribution, assign a qualitative score

<img>
Metrics like GQN, DIN, DV200 are based on the size distribution of DNA fragments
qPCR
Ratios of amplification for larger amplicons vs smaller amplicons used to determine quality

<img>
Both methods can be informative or predictive for performance in downstream NGS assays, or determine feasibility
For both quality and quantity, requirements will vary depending on downstream sequencing platform and application.  Long Read approaches generally require longer fragments; Short Read approaches can accommodate lower quality DNA to a point; even then, data quality may be impacted.
FFPE samples in particular are notorious for being degraded/lower quality
—
### Library preparation & QC:

<img>

All library preparation methods are converting the analyte (DNA, RNA) into a molecule that is usable by the sequencer. Fragment size, adapter structure, and whether the library is amplified will depend on the sequencing platform and application. Fragmentation can be either enzymatic or mechanical, and adapters can be y-shaped, linear, or hairpins. Adapters may incorporate binding sites (for flow cells or beads), sequencing primer sites, and barcodes for demultiplexing or unique molecular identification
After library preparation, obtaining the size distribution and molar concentration (molecule count) of the library is critical to ensuring that yields on expensive sequencing runs are maximized, and data quality is as good as it can be. qPCR for measuring concentration, BioAnalyzer/Tapestation/FemtoPulse are common methods for library QC before sequencing.


<img><img>
Loading too much  library on a sequencing flowcell/chip/SMRTcell can impact data quality/reduce yield; loading too little will reduce yield. Each sequencing platform has a range of acceptable library sizes and amounts required for a sequencing run. Presence of adapter/primer dimers can also reduce yield of usable data. Of note, achieving a library of the right size and purifying/excluding dimers from the final library is much more difficult when the input material is degraded (as in FFPE samples).
—
## Tailoring sequencing strategies for study’s design and purposes:
Major sequencing platforms and their use cases are summarized in the table and figure below:

<img>

<img>

### Illumina Sequencing

The process for Illumina sequencing is pictured below:

<img>

[Chaitankar et al. 2020]
(https://www.sciencedirect.com/science/article/pii/S1350946216300301?via%3Dihub)

Library molecules are loaded onto flow cells and individual library molecules are clonally amplified into clusters via bridge amplification. Clusters are sequenced using a Sequencing By Synthesis (SBS) method using reversible terminator nucleotides that are fluorescently labeled. Fluorescence is captured in images for each cycle of sequencing, and the signal for each cluster is interpreted. Poor fluorescence or mixed fluorescent signals within a cluster will reduce base calling certainty and therefore sequence quality.

### Ion Torrent
Library molecules undergo emulsion PCR on beads. After amplification, the emulsion is broken and the amplified beads are loaded into wells on the sequencing chip. Nucleotides are incorporated sequentially onto the chip; base incorporations for a given well/bead are interpreted as pH changes due to the release of Hydrogen.


<img>

### Pacbio Long Read Sequencing

DNA Polymerase and sequencing primer are bound to the hairpin adapters on library molecules and loaded into wells (ZMWs) on a SMRTcell. Once the polymerase is activated, it unzips the library molecule into a circle as it incorporates nucleotides, which are fluorescently labeled. Nucleotide incorporation is captured in real time as a movie of fluorescent signals. Each individual molecule can be read through many times; these are analytically combined into a single HiFi (or CCS) read that is highly accurate.


<img>
Long read sequencing is great for repetitive regions, STR  or AT/GC rich loci that are difficult to sequence or map with short reads, or where phasing across long regions is useful. Among long read platforms, Pacbio has the highest accuracy at 99.9%.

### Oxford Nanopore Sequencing
A motor protein unwinds dsDNA and drives single-stranded DNA (negatively charged) through the nanopore due to voltage applied across the membrane.
As nucleotides pass through the nanopore, the current change is measured and used to determine nucleotides as they pass. Nanopore sequencing has the opportunity for the longest reads of any platform, albeit at lower accuracy than Pacbio.

<img>
https://www.nature.com/articles/s41587-021-01108-x

## NGS strategies comparison

The three most typical NGS approaches are whole genome sequencing (WGS), whole exome sequencing (WES or WXS), and targeted sequencing, plus RNA sequencing for transcriptomic analysis. Each strategy has unique limitations that should be considered when selecting methodology for any study. See the figures below for a summary of each method’s capabilities and limitations:


<img>
[Hess et al. 2020](https://www.sciencedirect.com/science/article/pii/S0734975020300343
)

<img>

The general question when choosing sequencing methods is cost vs discovery power. WGS, and in particular long read sequencing, provides the most opportunity for discovery but is much more expensive than WES and targeted sequencing and has sequencing depth limitations due to cost. Additionally, WGS produces much more data than other methods which imposes additional considerations regarding data storage and computational analysis.

|                         |                                                         Targeted Panel (capture)                                                        |                                                   WES                                                  |                                    WGS                                    |
|-------------------------|:---------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------:|
| Cost (per sample)       | iSeq 100: Higher MiniSeq: Mid MiSeq: Mid NextSeq 1000/2000: Low NovaSeq 6000: Lowest                                                    | NextSeq 550: Highest NextSeq 1000/2000: Mid NovaSeq 6000: Lowest                                       | NextSeq 1000/2000: Higher  NovaSeq 6000: Low                              |
| Run Time (hours)        | iSeq 100: 9-19 MiniSeq: 2-24 MiSeq: 5-55 NextSeq 1000/2000: 11-48 NovaSeq 6000: 13-44                                                   | NextSeq 550: 11-29  NextSeq 1000/2000: 11-48  NovaSeq 6000: 13-44                                      | NextSeq 1000/2000: 11-48  NovaSeq 6000: 13-44                             |
| DNA Quantity Required   | 50-200ng                                                                                                                                | 50-200ng                                                                                               | 100-1000ng                                                                |
| Standard Coverage Depth | >500x                                                                                                                                   | >50-100x                                                                                               | >30x                                                                      |
| Samples (Per Run)       | iSeq 100: 1-48 MiniSeq: 1-96 MiSeq: 1-96 NextSeq 1000/2000: 1-28 NovaSeq 6000: 24-384                                                   | NextSeq 550: 8 NextSeq 1000/2000: 16-48 NovaSeq 6000: 24-384                                           | NextSeq 1000/2000: 1-3  NovaSeq 6000: 1-48                                |
| Reads Per Run           | iSeq 100: 4 million MiniSeq: 7-25 million MiSeq: 1-25 million NextSeq 1000/2000: 100 million-1.2 billion NovaSeq 6000: Up to 20 billion | NextSeq 550: 130-400 million NextSeq 1000/2000: 100 million-1.2 billion NovaSeq 6000: Up to 20 billion | NextSeq 1000/2000: 100 million-1.2 billion NovaSeq 6000: Up to 20 billion |
*Illumina statistics per sequencing strategy*

Note that these requirements and numbers are likely to change with time and between sequencing centers. Numbers for targeted panels in particular are highly variable based on the number of sequencing targets in your panel.

|                         | Targeted Panel or Long-Range Amplicon | WES |                       WGS                      |
|-------------------------|:-------------------------------------:|:---:|:----------------------------------------------:|
| Cost (per sample)       | $10-$200                              | N/A | ~$1500-5000                                    |
| Run Time (hours)        | 30hrs per SMRT cell                   |     | 30hrs per SMRT cell                            |
| DNA Quantity Required   | 50-500 ng                             |     | >3-5 ug*                                       |
| Standard Coverage Depth | >100x                                 |     | 10-30x                                         |
| Samples (Per Run)       | Up to 384 per SMRT Cell               |     | 1 sample per SMRT cell (10x coverage per cell) |
| Reads Per Run           | Up to 4,000,000                       |     | Up to 4,000,000                                |
*Pacbio long read statistics per sequencing strategy*


https://www.nature.com/articles/s41576-020-0236-x
<img>

### Discovery limitations

In addition to the practical limitations, each strategy has its own limitations for scientific discovery potential. See the tables below for details on each strategy.

|      Study Design     |     Sanger Sequencing     | SNP Array |          Targeted Panel         |               WES               |               WGS               |
|:---------------------:|:-------------------------:|:---------:|:-------------------------------:|:-------------------------------:|:-------------------------------:|
| GWAS                  | No                        | Yes       | Pre-defined genes/variants only | Yes- limited to exome           | Best                            |
| Driver genes          | Pre-defined genes only    | Yes       | Pre-defined genes only          | Yes- De novo discovery possible | Yes- De novo discovery possible |
| Non-coding            | Pre-defined regions only  | Limited   | Pre-defined regions only        | No                              | Best                            |
| Structural variant    | Pre-defined variants only | No        | No                              | Limited                         | Best                            |
| Tumor evolution       | No                        | Limited   | No                              | Limited                         | Best                            |
| Copy number analysis  | No                        | Limited   | No                              | Limited                         | Best                            |
| Mutational signatures | No                        | No        | No                              | Limited and potentially biased  | Best                            |
| Gene fusion           | Pre-defined fusions only  | No        | Limited                         | Limited                         | Best                            |

|                      Study Design                      |       Illumina Short Read WGS      | PacBio Long Reads |
|:------------------------------------------------------:|:----------------------------------:|:-----------------:|
| Structural Variants                                    | Limited to small SVs               | Best              |
| Genome Assembly and Resolution of Repeat-Heavy Regions | Cannot resolve long repeat regions | Best              |
| Copy Number Analysis                                   | Biased by GC content               | Unbiased          |
| Haplotype Phasing                                      | Limited                            | Best              |
| Pseudogene Detection                                   | Limited                            | Best              |

### Epigenomics

In addition to DNA-based and RNA-based bulk-sequencing, there are also epigenomic related bulk-sequencing techniques to explore the epigenome. Hi-C can be used to study the 3D genome organization and topology-associated domain. To examine chromatin accessibility (closed and open chromatin regions) a number of different approaches could be used such as MNase-seq, ATAC-seq, DNase-seq, FAIRE-seq. This is important to understand which genomic regions are transcriptionally active or repressed.

<img>
https://www.sciencedirect.com/science/article/pii/B9780128122150000042

For protein-DNA interactions (e.g. the binding of specific transcription factors to DNA), ChIP-seq is particularly helpful. To study small RNA, miRNA-seq can be used.
For DNA methylation, there are various approaches including Whole-Genome Bisulfite sequencing, Reduced Representation Bisulfite sequencing, and Methylated DNA immunoprecipitation.

### Single-cell methods

In recent years many new technologies have been developed to sequence at the individual cell level. We will not discuss these methods in this session, but see the figure below for a brief, non-comprehensive summary of single-cell methods:


<img>
Ren et al. 2018

## Targeted Sequencing Details

Targeted sequencing sequences key genes or regions of interest to high depth (500–1000× or higher), allowing identification of rare/established variants. This provides cost-effective findings for studies of disease-related genes and delivers accurate, easy-to-interpret results, identifying causative novel or inherited mutations at low allele frequencies (down to 5%). Two common methods for targeted sequencing include target enrichment (hybridization capture) and amplicon generation.

|                                    Targeted Enrichment                                    |                                Amplicon Sequencing                               |
|:-----------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------:|
| Captured by hybridization to biotinylated probes and then isolated by magnetic pulldown.  | Amplified and purified using highly multiplexed oligo pools.                     |
| Target enrichment captures 20b–62 Mb regions                                              | Allows sequencing a few genes to hundreds of genes in a single run               |
| Larger gene content, typically >50 genes                                                  | Smaller gene content, typically < 50 genes                                       |
| More comprehensive profiling for all variant types                                        | Ideal for analyzing single nucleotide variants and insertions/deletions (indels) |
| More comprehensive method, but with longer hands-on time and turnaround time              | More affordable, easier workflow                                                 |


<img>
*A) Target enrichment, B) Amplicon sequencing. [Source](https://www.thermofisher.com/us/en/home/life-science/sequencing/sequencing-learning-center/next-generation-sequencing-information/ngs-basics/targeted-sequencing-approaches.html)*

### Targeted Panel Types

There are targeted panels for many purposes, but some of the most common panels for cancer research are highlighted below.

|           Gene Panel          | Gene Count |   Sample Type   |             Variants            |                                                                     Notes                                                                    |
|:-----------------------------:|:----------:|:---------------:|:-------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------------------:|
| Oncomine Comprehensive Plus   | 500 +      | FFPE, DNA,  RNA | Hotspots, CNVs, fusions, splice | Broad, PanCancer assay for comprehensive genomic profiling of key biomarkers across targeted and immuno-oncology applications                |
| Oncomine v3                   | 161        | FFPE, DNA,  RNA | Hotspots, CNVs, fusions         | Multi-biomarker assay based on latest clinical oncology research for targeted solid tumor applications                                       |
| MSK-IMPACT                    | 505        | DNA             | Hotspots, CNVs, fusions         | MSK-IMPACT is the first laboratory-developed tumor-profiling test to receive authorization from the US Food and Drug Administration.         |
| TruSight RNA Pan-Cancer Panel | 1385       | RNA, FFPE       | Hotspots, Fusions, Expression   | Enables quantitative measurement of gene expression as well as the detection of gene fusions with both known and novel gene fusion partners. |

Information on other targeted panels can be found on the Illumina website (see examples below) or from other vendors.


<img>

Many instrument-specific softwares are available with various analysis modules for data analysis, or alternatively standard exome sequencing pipelines can be used.


<img>

## Whole-Exome Bulk Sequencing (WES), short reads

As the name implies, WES sequences DNA coding regions which are also known as exons. Many of the genetic mutations that lead to genetic disorders happen in the exome, which is why WES can still be very useful despite the fact that it doesn’t analyze the entire genome. Importantly it is much more cost efficient than WGS. The standard bioinformatics pipeline for WES is somewhat similar to WGS, but with the need to specify WES capture platforms.

<img>Choi et al. 2018
To perform exome capture, genomic DNA is sheared (commonly by sonication i.e. physical force such as ultrasonicators using burst of ultrasounds) to create evenly sized DNA fragments which are then followed by adaptor ligation. Size selection is performed on the library prior to capture. Size-selected libraries are then incubated with biotinylated probes complementary to the exome regions (like cDNA library or RNA baits). The bait-DNA hybrids are then “pulled” out of the complex mixture by incubation with the magnetic beads. RNA baits are digested such that the remaining nuclear acid is the targeted DNA of interest. Captured DNA is amplified and the targeted samples are sequenced.

### WES Capture Platforms
|               Platform               |      Target Capture Region Length      | Required input quantity |                                                     Bedfile links                                                    |
|:------------------------------------:|:--------------------------------------:|:-----------------------:|:--------------------------------------------------------------------------------------------------------------------:|
| Agilent SureSelect Human All Exon v8 | 35.1 Mb                                | 10-400ng of DNA         | https://kb.10xgenomics.com/hc/en-us/articles/115004150923-Where-can-I-find-the-Agilent-Target-BED-files-             |
| Roche KAPA HyperExome v3.0           | 43 Mb - targeting hg38 genome assembly | 100ng DNA               | Download HG38 Design Files for the KAPA HyperExome Probes Download hg19 Design Files for the KAPA HyperExome Probes  |
| Illumina TruSeq                      | 45 Mb                                  | 100ng of DNA            | https://emea.support.illumina.com/downloads/truseq-exome-product-files.html                                          |
Note that these details and common platforms may change with time.
## Whole-Genome Bulk Sequencing (WGS), short reads

WGS is the most comprehensive method for analyzing entire genomes. WGS provides a high-resolution, base-by-base view of the genome and captures both large and small variants that might be missed with targeted approaches. WGS also delivers large volumes of data in a short amount of time to support assembly of novel genomes.
### WGS Analysis
SNV and INDELs can be called using a variety of callers, such as Mutect2, Strelka2, Vardict, Muse, or LoFreq, or an ensemble caller such as [Somatic Combiner](https://github.com/mingyi-wang/somatic-combiner) which is used in the NCI-CGR pipeline. We will discuss the practical details of mutation calling in the next session.


<img>*SomaticCombiner: Wang, M.  et al. Scientific Reports,  10:12898 (2020)*

Alternatively there is the Sentieon variant discovery pipeline which was used for Sherlock Lung. We will discuss the Sentieon pipeline further in the next session on variant calling.

<img>
## Quality Control: DNA Sequencing Data

There are several metrics of quality control downstream of NGS sequencing at various steps in the process.
Sample quality control - DNA/RNA integrity and quantity
Library quality control
Check library size distribution
Minimize adapter-dimers
Molar concentration
Sequencing data quality control
FastQC/MultiQC
Adapter trimming
Illumina Sequence analysis viewer (SAV)
Flowcell level run metrics such as yield, error rate, %>=Q30,, Cluster PF(%)
Quality control after mapping
Deduplication, duplication rates
Insert size
Coverage
Contamination
Sex concordance
There are several tools available for QC of sequencing data, but two of the most common are FastQC/MultiQC. We will practice using both in the practical section for this session.

#### MultiQC
MultiQC is a general use tool, perfect for summarizing output from numerous bioinformatics tools (currently [114 supported](https://multiqc.info/)). MultiQC searches a given directory for analysis logs and compiles an HTML report encompassing all samples and tools.

<img>
There are many QC tools, but some of the more popular include:
- FastQC
- Samtools
  - Stats, flagstats, idxstats, rmdup
- Picard
  - MarkDuplicates, WgsMetrics, AlignmentSummaryMetrics, etc
- GATK
  - BaseRecalibrator, VariantEval
- Peddy
  - het_check, ped_check, sex_check, summary_table, background_pca
- VerifyBamID: ancestry-agnostic DNA contamination estimation method
- Somalier: checking relatedness between sample
### Check Sequencing Quality Using FastQC

<img>

### Checking Insert Size with

<img>
### Assessing Gender and Relatedness with Somalier
 <img>
