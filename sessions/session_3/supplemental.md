---
layout: page
permalink: sessions/session_3/supplemental
menubar_toc: true
---

### Illumina Sequencing

The process for Illumina sequencing is pictured below:

<img class="center" style="display: block;margin-left: auto; margin-right: auto;max-width:66%;" src="lecture_assets/illumina_sequencing_method.png">
<figcaption class="has-text-centered is-size-7 is-italic"><a href="https://www.sciencedirect.com/science/article/pii/S1350946216300301?via%3Dihub">Chaitankar et al. 2020</a></figcaption><br>

Library molecules are loaded onto flow cells and individual library molecules are clonally amplified into clusters via bridge amplification. Clusters are sequenced using a Sequencing By Synthesis (SBS) method using reversible terminator nucleotides that are fluorescently labeled. Fluorescence is captured in images for each cycle of sequencing, and the signal for each cluster is interpreted. Poor fluorescence or mixed fluorescent signals within a cluster will reduce base calling certainty and therefore sequence quality.

---

### Ion Torrent
Library molecules undergo emulsion PCR on beads. After amplification, the emulsion is broken and the amplified beads are loaded into wells on the sequencing chip. Nucleotides are incorporated sequentially onto the chip; base incorporations for a given well/bead are interpreted as pH changes due to the release of Hydrogen.

<img class="center" style="display: block;margin-left: auto; margin-right: auto;max-width:66%;" src="lecture_assets/ionTorrent_sequencing_method.jpg">
<figcaption class="has-text-centered is-size-7 is-italic"><a href="https://academic.oup.com/bioinformatics/article/29/13/i344/188472?login=false">Golan and Medvedev, 2013</a></figcaption><br>

---

### Pacbio Long Read Sequencing

DNA Polymerase and sequencing primer are bound to the hairpin adapters on library molecules and loaded into wells (ZMWs) on a SMRTcell. Once the polymerase is activated, it unzips the library molecule into a circle as it incorporates nucleotides, which are fluorescently labeled. Nucleotide incorporation is captured in real time as a movie of fluorescent signals. Each individual molecule can be read through many times; these are analytically combined into a single HiFi (or CCS) read that is highly accurate.

<img class="center" style="display: block;margin-left: auto; margin-right: auto;max-width:66%;" src="lecture_assets/Pacbio_sequencing_method.jpg">
<figcaption class="has-text-centered is-size-7 is-italic">From https://www.pacb.com/</figcaption><br>

Long read sequencing is great for repetitive regions, STR  or AT/GC rich loci that are difficult to sequence or map with short reads, or where phasing across long regions is useful. Among long read platforms, Pacbio has the highest accuracy at 99.9%.

---

### Oxford Nanopore Sequencing
A motor protein unwinds dsDNA and drives single-stranded DNA (negatively charged) through the nanopore due to voltage applied across the membrane.
As nucleotides pass through the nanopore, the current change is measured and used to determine nucleotides as they pass. Nanopore sequencing has the opportunity for the longest reads of any platform, albeit at lower accuracy than Pacbio.

<img class="center" style="display: block;margin-left: auto; margin-right: auto;max-width:66%;" src="lecture_assets/Nanopore_sequencing_method.jpg">
<figcaption class="has-text-centered is-size-7 is-italic">From <a href="https://www.nature.com/articles/s41587-021-01108-x">Wang et al., 2021</a></figcaption>

---

### Epigenomics

In addition to DNA-based and RNA-based bulk-sequencing, there are also epigenomic related bulk-sequencing techniques to explore the epigenome. Hi-C can be used to study the 3D genome organization and topology-associated domain. To examine chromatin accessibility (closed and open chromatin regions) a number of different approaches could be used such as MNase-seq, ATAC-seq, DNase-seq, FAIRE-seq. This is important to understand which genomic regions are transcriptionally active or repressed.

<img class="center" style="display: block; margin-left: auto; margin-right: auto;" src="lecture_assets/epigenomics_seq.png">
<figcaption class="has-text-centered is-size-7 is-italic"><a href="https://www.sciencedirect.com/science/article/pii/B9780128122150000042">Hsu et al., 2018</a></figcaption>

For protein-DNA interactions (e.g. the binding of specific transcription factors to DNA), ChIP-seq is particularly helpful. To study small RNA, miRNA-seq can be used.
For DNA methylation, there are various approaches including Whole-Genome Bisulfite sequencing, Reduced Representation Bisulfite sequencing, and Methylated DNA immunoprecipitation.

---

### Single-cell methods

In recent years many new technologies have been developed to sequence at the individual cell level. We will not discuss these methods in this session, but see the figure below for a brief, non-comprehensive summary of single-cell methods:

<img class="center" style="display: block; margin-left: auto; margin-right: auto;" src="lecture_assets/singleCell_seqs.png">
<figcaption class="has-text-centered is-size-7 is-italic"><a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1593-z">Ren et al., 2018</a></figcaption>

---
