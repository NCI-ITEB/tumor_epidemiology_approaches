---
layout: page
permalink: sessions/session_4/practical
menubar_toc: true
---

<!--script link="practical_assets/{{ site.baseurl }}/assets/js/vanilla-back-to-top.min.js"></script> <script>addBackToTop()</script-->

<script src="{{ site.baseurl }}/assets/js/copyCodeSnippet.js" defer></script>
<script src="{{ site.baseurl }}/assets/js/copyCodeBlock.js" defer></script>

<style>
pre {
  max-height: 500px;
  overflow-y: auto;
  max-width: 120%;
}
</style>

Before we begin, please login to Biowulf and request an interactive session:

For a reminder on how to log-in to Biowulf, we refer you to this Biowulf HPC guide. In short:

- (Windows) use PuTTy
- (Mac) enter <code>ssh USERNAME@biowulf.nih.gov</code>{% include code-snippet-copy.html %} to your terminal, then enter your password

Request an interactive session with the following command:

<code>sinteractive --mem=12g --cpus-per-task=2</code>{% include code-snippet-copy.html %}

---

This practical section will focus on somatic mutational calling and downstream analysis.
Today's tasks include (A) performing somatic mutational calling using the GATK pipeline, (B) annotating the mutation calls, (C) visualizing somatic mutations in IGV, and (D) downstream analysis using MAFtools.

---

## The GATK pipeline of somatic short variant (SNVs and indels) discovery

We will use the targeted sequencing data of paired tumor and normal samples from kidney for this practice. The scripts for this task will takes about 25 mins running time, so we will submit the scripts to Biowulf at the beginning of the session and then walk through the scripts together.

1. First let's create a new folder in your personal data directory, and copy all the sample data to this directory.

<code>cd /data/$USER/<br>
mkdir practical_session_4<br>
cd practical_session_4<br>
cp /data/classes/DCEG_Somatic_Workshop/Practical_session_4/* .</code>{% include code-snippet-copy.html %}

2. Next, we will submit two sets of scripts to Biowulf. Use the following commands to submit the scripts for step 1.

<code>cd shfiles<br>
jobid1=$(swarm -t 8 -g 24  --gres=lscratch:24 --job-name step1 --time 1:00:00 scripts_for_step1.swarm)<br>
echo $jobid1</code>{% include code-snippet-copy.html %}

<code>swarm</code> is a command to submit a group of command lines to biowulf and those command lines will run as parallel processes. The options <code>-g</code>, <code>-t</code> and <code>--time</code> specify the memory (in GB), threads and running time per process. We request 8CPUs and 24GB of memory for each process. The default running time is 2hrs. We assign a job name using the option <code>--job-name</code>. The option <code>--gres</code> is to allocate a local scratch disk space for each process. <code>lscratch:N</code> will require N GB for the scratch space. The directory of this scratch space can be accessed in the script as <code>/lscratch/$SLURM_JOB_ID/</code>. More information could be found at [https://hpc.nih.gov/apps/swarm.html](https://hpc.nih.gov/apps/swarm.html).

The <code>swarm</code> will return a job ID and it will be saved to <code>$jobid1</code> in this script. We then submit the swarm file for step 2 and let it start only after the completion of all commands in the first job. To set up the dependency, we use the option <code>--dependency=afterany:jobid</code> where the jobid is $jobid1. Step 1 will take about 9 mins.


<code>jobid2=$(swarm -t 8 -g 20 --gres=lscratch:20 --dependency=afterany:$jobid1 --job-name step2 scripts_for_step2.swarm)<br>
echo $jobid2</code>{% include code-snippet-copy.html %}

In a few seconds, let's check the job status using <code>sjobs<code>, the second job will sit in a pending state until all the processes in job1 are finished.

### Data pre-processing for somatic short variant discovery

While we are waiting for the two jobs, let's examine the script.

<code>vi scripts_for_step1.swarm</code>{% include code-snippet-copy.html %}
```bash
#!/bin/bash
sh Step1_preprocess_variant_discovery.sh GPK7017_2000 ../fastq ../Results
sh Step1_preprocess_variant_discovery.sh GPK4013_0401 ../fastq ../Results
```

The first step for somatic or germline variant discovery is to pre-process the raw sequence data to produce analysis-ready BAM files. As the following figure shows, we start from the raw FASTQ files (or BAM files), align them to the reference genome, and clean up the data to correct for technical biases (e.g. duplicate reads, strand bias, cross-sample contamination, etc).

{% include image-modal.html link="practical_assets/GATK_data_preprocess.png" max-width="30%" %}

**1\.** In the swarm file, each line corresponds to one sample, and will be submitted as one process in the swarm submission. They are all in the format:
<code>sh ./Step1_preprocess_variant_discovery.sh $SAMPLE $INPUT_DIR $OUTPUT_DIR</code>.
All the fields after the bash file (.sh) are arguments, and they will be passed to the bash (.sh) file. All the arguments in the swarm commands should be processed within the bash file. In this way, we can prepare the 'core' commands as bash files to process samples in batches.

**2\.** Now let's quit the swarm file by typing <code>:q</code> and hitting enter, and open the corresponding bash file.

<code>vi Step1_preprocess_variant_discovery.sh</code>{% include code-snippet-copy.html %}

```bash
#!/bin/bash

#### Data preprocess for somatic short variant discovery using GATK workflow ###

module load samtools
module load bwa
module load fastp
module load GATK/4.3.0.0
module load picard/2.27.3

GATK_Bundle=/fdb/GATK_resource_bundle/hg38-v0
GENOME=$GATK_Bundle/Homo_sapiens_assembly38.fasta

SAMPLE=$1
INDIR=$2
DIR=$3
logs=$DIR/logs
read1=$INDIR/${SAMPLE}_R1.fastq.gz
read2=$INDIR/${SAMPLE}_R2.fastq.gz
id=$SAMPLE
lb=$id
sm=$id

SECONDS=0

echo -e "sample:$SAMPLE\nindir:$INDIR\noutdir:$DIR"

if [ ! -d "$DIR" ]; then
        mkdir -p $DIR
fi
if [ ! -d "$logs" ]; then
        mkdir -p $logs
fi
### Perform adaptor trimming on fastq files ###
fastp -i $read1 -I $read2 \
      --stdout --thread 2 \
      -j ${logs}/fastp-${SAMPLE}.json \
      -h ${logs}/fastp-${SAMPLE}.html \
      2> ${logs}/fastp-${SAMPLE}.log | \
bwa mem -M -t 8 \
      -R "@RG\tID:$id\tPL:ILLUMINA\tLB:$lb\tSM:$sm" \
      $GENOME - 2> ${logs}/bwa-${SAMPLE}.log | \
samtools sort -T /lscratch/$SLURM_JOB_ID/ -m 2G -@ 4 -O BAM \
      -o $DIR/${SAMPLE}_sort.bam  2> ${logs}/samtools-${SAMPLE}.log

###/lscratch/$SLURM_JOBID
duration=$SECONDS
echo "Alignment completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

### Duplicate marking in coordinate-sorted raw BAM files ###
java -Xmx2g -jar $PICARDJARPATH/picard.jar MarkDuplicates \
     -I $DIR/${SAMPLE}_sort.bam \
    -O /dev/stdout \
    -M marked_dup_metrics.txt 2> ${logs}/markdup-${SAMPLE}.log \
|java -Xmx2g -jar $PICARDJARPATH/picard.jar SortSam \
      -I /dev/stdin -O $DIR/${SAMPLE}_markdup_sorted.bam \
      -SORT_ORDER coordinate 2>> ${logs}/markdup-${SAMPLE}.log

duration=$SECONDS
echo "Duplicate marking completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

DBSNP=/fdb/GATK_resource_bundle/hg38/dbsnp_138.hg38.vcf.gz
INDEL=/fdb/GATK_resource_bundle/hg38-v0/Homo_sapiens_assembly38.known_indels.vcf.gz
GOLD_INDEL=/fdb/GATK_resource_bundle/hg38-v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

### Recaliberating base quality score###
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms6G -Xmx6G -XX:ParallelGCThreads=2" BaseRecalibrator \
  -I $DIR/${SAMPLE}_markdup_sorted.bam \
  -R $GENOME \
  -O  $DIR/${SAMPLE}_markdup_bqsr.report \
  --known-sites $DBSNP \
  --known-sites $INDEL \
  --known-sites $GOLD_INDEL \
  2> ${logs}/BQSR-${SAMPLE}.log

gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms6G -Xmx6G -XX:ParallelGCThreads=2" ApplyBQSR \
  -I $DIR/${SAMPLE}_markdup_sorted.bam  \
  -R $GENOME \
  --bqsr-recal-file  $DIR/${SAMPLE}_markdup_bqsr.report \
  -O  $DIR/${SAMPLE}_markdup_bqsr.bam \
  2>> ${logs}/BQSR-${SAMPLE}.log

duration=$SECONDS
echo "Base recaliration completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
```

We first load several modules using <code>module load</code>. Then we specify the reference/source files to be used. Most of the reference files corresponding to a pre-installed application in Biowulf can be found in the folder /fdb/.
Note the way we pass the arguments to the bash file. Arguments passed to a script are processed in the same order in which they’re sent (in this case, the order in the swarm file). The indexing of the arguments starts at one, and the first argument can be accessed inside the script using $1. Similarly, the second argument can be accessed using $2, and so on. Here we assign the first argument to variable $SAMPLE, and so on.

<!--{% include image-modal.html link="practical_assets/Step1_1.PNG" %}-->
```bash
GATK_Bundle=/fdb/GATK_resource_bundle/hg38-v0
GENOME=$GATK_Bundle/Homo_sapiens_assembly38.fasta

SAMPLE=$1
INDIR=$2
DIR=$3
logs=$DIR/logs
read1=$INDIR/${SAMPLE}_R1.fastq.gz
read2=$INDIR/${SAMPLE}_R2.fastq.gz
id=$SAMPLE
lb=$id
sm=$id
```

**3\.** We use the [fastp](https://github.com/OpenGene/fastp) command to trim the adapters, [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) to perform raw alignment to the reference genome, and [samtools sort](http://www.htslib.org/doc/samtools-sort.html) to sort the files by genomic coordinates. The three commands are connected by 'pipes' <code>|</code> to bypass storage of large intermediate files. The adapter trimming step is optional. The backslash at the end of each line is used to skip the line breaker and the command continues in the next line.

<!--{% include image-modal.html link="practical_assets/Step1_2.PNG" %}-->
```bash
### Perform adaptor trimming on fastq files ###
fastp -i $read1 -I $read2 \
      --stdout --thread 2 \
      -j ${logs}/fastp-${SAMPLE}.json \
      -h ${logs}/fastp-${SAMPLE}.html \
      2> ${logs}/fastp-${SAMPLE}.log | \
bwa mem -M -t 8 \
      -R "@RG\tID:$id\tPL:ILLUMINA\tLB:$lb\tSM:$sm" \
      $GENOME - 2> ${logs}/bwa-${SAMPLE}.log | \
samtools sort -T /lscratch/$SLURM_JOB_ID/ -m 2G -@ 4 -O BAM \
      -o $DIR/${SAMPLE}_sort.bam  2> ${logs}/samtools-${SAMPLE}.log
```

Let's take a look at the three command separately:

In the last line of the <code>fastp</code> command lines, '2' is a file descriptor. It redirects(writes) any 'standard error' messages to the log file fastp-${SAMPLE}.log, which will be otherwise output to the screen.
In <code>bwa mem</code> options, <code>-t</code>specifies the number of threads. <code>-R</code> specifies the read header lines. <code>-</code> indicates taking the input from the standard input. In this example, we used the pipes to redirect the standard output of last command (aka, fastp) as the standard input into the next command (bwa). So we use this option to tell the script take the output from command as the input data.
In <code>samtools sort</code> options, <code>-T</code> is the prefix of temporary files. <code>-m</code> and <code>-@</code> specify the memory per thread and the number of threads, respectively. <code>-O</code> specifies the output format.   

**4\.** Now we have a sorted raw alignment BAM file. In the next section of the script, we mark read pairs that are likely to have originated from duplicates of the same DNA fragments uing [picard tools](https://broadinstitute.github.io/picard/). The picard tools package contains multiple commands, and the command lines look like this:
<code>java [java opts] -jar $PICARDJARPATH/picard.jar COMMAND [options] </code>

<!--{% include image-modal.html link="practical_assets/Step1_3.PNG" %}-->
```bash
### Duplicate marking in coordinate-sorted raw BAM files ###
java -Xmx2g -jar $PICARDJARPATH/picard.jar MarkDuplicates \
     -I $DIR/${SAMPLE}_sort.bam \
    -O /dev/stdout \
    -M marked_dup_metrics.txt 2> ${logs}/markdup-${SAMPLE}.log \
| java -Xmx2g -jar $PICARDJARPATH/picard.jar SortSam \
      -I /dev/stdin -O $DIR/${SAMPLE}_markdup_sorted.bam \
      -SORT_ORDER coordinate 2>> ${logs}/markdup-${SAMPLE}.log
```

We use a pipe to connect two picard commands <code>MarkDuplicates</code> and <code>SortSam</code>. Similar to the bwa example, in the second picard command, we use <code>-I /dev/stdin</code> to specify reading from standard input.

**5\.** In the last part of the script, we recalibrate base quality score to correct for patterns of systematic errors using the [GATK tools for BQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-). The GATK command lines look like this:
<code>gatk --java-options "[java opts]" ToolName [options] </code>

<!--{% include image-modal.html link="practical_assets/Step1_4.PNG" %}-->
```bash
### Recaliberating base quality score###
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms6G -Xmx6G -XX:ParallelGCThreads=2" BaseRecalibrator \
  -I $DIR/${SAMPLE}_markdup_sorted.bam \
  -R $GENOME \
  -O  $DIR/${SAMPLE}_markdup_bqsr.report \
  --known-sites $DBSNP \
  --known-sites $INDEL \
  --known-sites $GOLD_INDEL \
  2> ${logs}/BQSR-${SAMPLE}.log

gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms6G -Xmx6G -XX:ParallelGCThreads=2" ApplyBQSR \
  -I $DIR/${SAMPLE}_markdup_sorted.bam  \
  -R $GENOME \
  --bqsr-recal-file  $DIR/${SAMPLE}_markdup_bqsr.report \
  -O  $DIR/${SAMPLE}_markdup_bqsr.bam \
  2>> ${logs}/BQSR-${SAMPLE}.log
```

The first tool <code>BaseRecalibrator</code> builds the recalibration model. As we calculate the mismatched bases, we exclude the loci known to vary in the population, which requires the input of known variants resource. And this is specified by the option <code>--known-sites</code>.  The second tool, <code>ApplyBQSR</code> ,adjusts the score based on the model.

**6\.** By the end of step 1, you should have generated analysis-ready BAM files for each individual sample (tumor and normal). We prepared the sample output data in the /Results/ directory for you to check.

<code>ls ../Results</code>

---

### Somatic short variant discovery using MuTect2

Next we will use the analysis-ready BAM files to proceed the GATK pipeline for somatic short variant discovery. There are two main steps: generating candidate variants and filtering. We use GATK tools for both steps. For more details please refer to: [https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-).

{% include image-modal.html link="practical_assets/GATK_somatic_variant_calling.png" %}

Let's examine the swarm file.
<code>vi scripts_for_step2.swarm</code>{% include code-snippet-copy.html %}
```bash
#!/bin/bash
sh Step2_somatic_variant_discovery_MuTect2.sh  GPK7017_2000 GPK4013_0401 PAP1_1708_01_T01 ../Results ../Results
```

Simlar to the first step, we use one bash file to process tumor/normal sample pairs in batches. In this example, we have only one tumor/normal pair, but we can easily scale up in this way. Now let's quit the swarm file by with <code>:q</code> and open the script file:

<code>vi Step2_somatic_variant_discovery_MuTect2.sh</code>{% include code-snippet-copy.html %}

```bash
#!/bin/bash

#### Call somatic short variant discovery using MuTect2 and filtering ###

module load samtools
module load bwa
module load fastp
module load GATK/4.3.0.0
module load picard/2.27.3

GATK_Bundle=/fdb/GATK_resource_bundle/hg38-v0
GENOME=$GATK_Bundle/Homo_sapiens_assembly38.fasta
COMMONVAR=/data/classes/DCEG_Somatic_Workshop/Practical_session_4/Reference/small_exac_common_3.hg38.vcf.gz

NSAMPLE=$1
TSAMPLE=$2
PREFIX=$3
INDIR=$4
DIR=$5
logs=$DIR/logs

BAM_NORMAL=$INDIR/${NSAMPLE}_markdup_bqsr.bam
BAM_TUMOR=$INDIR/${TSAMPLE}_markdup_bqsr.bam
BASE_TUMOR=`samtools view -H $BAM_TUMOR |awk '$1~/^@RG/ {for (i=1;i<=NF;i++) {if ($i~/SM/) {split($i,aa,":"); print aa[2]}}}'|sort|uniq`
BASE_NORMAL=`samtools view -H $BAM_NORMAL | awk '$1~/^@RG/ {for (i=1;i<=NF;i++) {if ($i~/SM/) {split($i,aa,":"); print aa[2]}}}'|sort|uniq`

OUT_VCF=$DIR/${PREFIX}.vcf
OUT_FILTERED_VCF=$DIR/${PREFIX}_filtered.vcf
OUT_PASSED_VCF=$DIR/${PREFIX}_passed.vcf
OUT_STATS=$DIR/${PREFIX}.vcf.stats

SECONDS=0

if [ ! -d "$DIR" ]; then
        mkdir -p $DIR
fi
if [ ! -d "$logs" ]; then
        mkdir -p $logs
fi
### SNV and Indel calling  ###
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=1" Mutect2 \
  -R $GENOME \
  -I $BAM_NORMAL \
  -I $BAM_TUMOR \
  -normal $BASE_NORMAL \
  -tumor $BASE_TUMOR \
  -O $OUT_VCF \
  2> ${logs}/Mutect2-${SAMPLE}.log

duration=$SECONDS
echo "MuTect2 completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

### Filtering the MuTect2 variant calls ###
### Estimate cross-sample contamination  ###

gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=2" GetPileupSummaries \
   -I $BAM_TUMOR \
   -V $COMMONVAR \
   -L $COMMONVAR \
   -O $DIR/${TSAMPLE}_pileups.table \
   2> ${logs}/VarFilter-${SAMPLE}.log

gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=2" GetPileupSummaries \
   -I $BAM_NORMAL \
   -V $COMMONVAR \
   -L $COMMONVAR \
   -O $DIR/${NSAMPLE}_pileups.table \
   2>> ${logs}/VarFilter-${SAMPLE}.log

gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=2" CalculateContamination \
     -I $DIR/${TSAMPLE}_pileups.table \
     -matched $DIR/${NSAMPLE}_pileups.table \
     -tumor-segmentation $DIR/${TSAMPLE}_segments.table \
     -O $DIR/${TSAMPLE}_calculatecontamination.table \
     2>> ${logs}/VarFilter-${SAMPLE}.log

### Filter variants  ###

gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2" FilterMutectCalls \
  -R $GENOME \
  --contamination-table $DIR/${TSAMPLE}_calculatecontamination.table \
  --stats $OUT_STATS \
  --tumor-segmentation $DIR/${TSAMPLE}_segments.table \
  -O $OUT_FILTERED_VCF \
  -V $OUT_VCF \
  2>> ${logs}/VarFilter-${SAMPLE}.log
awk '($1 ~/^#/) || ($7 ~ /PASS/) {print}' $OUT_FILTERED_VCF >$OUT_PASSED_VCF

###/lscratch/$SLURM_JOBID
duration=$SECONDS
echo "Filtering completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
```

**1\.** In the first part, we load necessary modules, specify reference resources and arguments from input. We pass five arguments from the swarm commands in the following order: normal sample name, tumor sample name, prefix for output file name (e.g. patient or subject ID), input file directory, and output file directory.

<!--{% include image-modal.html link="practical_assets/Step2_1.PNG" %}-->
```bash
GATK_Bundle=/fdb/GATK_resource_bundle/hg38-v0
GENOME=$GATK_Bundle/Homo_sapiens_assembly38.fasta
COMMONVAR=/data/classes/DCEG_Somatic_Workshop/Practical_session_4/Reference/small_exac_common_3.hg38.vcf.gz

NSAMPLE=$1
TSAMPLE=$2
PREFIX=$3
INDIR=$4
DIR=$5
logs=$DIR/logs
```

**2\.** In the next part, we use the GATK tool [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) to call candidate variants.

<!--{% include image-modal.html link="practical_assets/Step2_2.PNG" %}-->

```bash
### SNV and Indel calling  ###
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=1" Mutect2 \
  -R $GENOME \
  -I $BAM_NORMAL \
  -I $BAM_TUMOR \
  -normal $BASE_NORMAL \
  -tumor $BASE_TUMOR \
  -O $OUT_VCF \
  2> ${logs}/Mutect2-${SAMPLE}.log
```

The options <code>-I</code> were used twice to specify the input BAM files for normal and tumor samples, respectively. The names of the tumor and normal samples are specified by the options <code>-tumor</code> and <code>-normal</code>, which have to match the sample names in the read group lines of the BAM file. Note that the read group lines take this format: <code>@RG\tID:$id\tPL:ILLUMINA\tLB:$lb\tSM:$sm</code>, and the sample names are in the last field.

If we are working with whole exome or genome sequencing data, then this step requires a large amount of memory. They can be changed in the java options <code>-Xms</code> and <code>-Xmx</code>. MuTect2 does not allow multithreading in its options. So to spped up this process, you may consider parallel computing using a 'scatter-gather' approach, that is, run separate GATK commands to process a proportion of the input data and collate all the results into the final output. For more information, please refer to [https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather).

**TODO:Confirm with Tongwu the memory requested for WGS/WES analysis**

**3\.** The next part of the script calculates the cross-sample contamination and estimates the allelic copy number segmentation for each tumor sample. It will generate a contamination-table file that will be used for filtering  This is an optional step.

<!--{% include image-modal.html link="practical_assets/Step2_3.PNG" %}-->
```bash
### Filtering the MuTect2 variant calls ###
### Estimate cross-sample contamination  ###

gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=2" GetPileupSummaries \
   -I $BAM_TUMOR \
   -V $COMMONVAR \
   -L $COMMONVAR \
   -O $DIR/${TSAMPLE}_pileups.table \
   2> ${logs}/VarFilter-${SAMPLE}.log

gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=2" GetPileupSummaries \
   -I $BAM_NORMAL \
   -V $COMMONVAR \
   -L $COMMONVAR \
   -O $DIR/${NSAMPLE}_pileups.table \
   2>> ${logs}/VarFilter-${SAMPLE}.log

gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=2" CalculateContamination \
     -I $DIR/${TSAMPLE}_pileups.table \
     -matched $DIR/${NSAMPLE}_pileups.table \
     -tumor-segmentation $DIR/${TSAMPLE}_segments.table \
     -O $DIR/${TSAMPLE}_calculatecontamination.table \
     2>> ${logs}/VarFilter-${SAMPLE}.log
```

The tool <code>GetPileupSummaries</code> summarizes read counts that support reference, alternative and other alleles for given sites. In its options <code>-V</code> and <code>-L</code>, it requires an input file for common germline variant sites, e.g. this file from the gnomAD resource. We provide a reference file small_exac_common_3.hg38.vcf.gz in our folder. Other reference files for this purpose (e.g. gnomAD) can be found at the GATK Google bucket. In an interactive session, you may use the following commands to check and download the reference files for common germline variant sites.

<code>module load google-cloud-sdk<br>
gsutil ls gs://gatk-best-practices/somatic-hg38/<br>
gsutil ls gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz* ./</code>

**4\.** Then we will filter the mutation candidates using the tool [FilterMutationCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls). The options <code>--contamination-table</code> and <code>--tumor-segmentation</code> use the output from the previous step and are optional.

<!--{% include image-modal.html link="practical_assets/Step2_4.PNG" %}-->
```bash
### Filter variants  ###

gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2" FilterMutectCalls \
  -R $GENOME \
  --contamination-table $DIR/${TSAMPLE}_calculatecontamination.table \
  --stats $OUT_STATS \
  --tumor-segmentation $DIR/${TSAMPLE}_segments.table \
  -O $OUT_FILTERED_VCF \
  -V $OUT_VCF \
  2>> ${logs}/VarFilter-${SAMPLE}.log
awk '($1 ~/^#/) || ($7 ~ /PASS/) {print}' $OUT_FILTERED_VCF >$OUT_PASSED_VCF
```

This process annotates the variants as 'PASS' in the 'FILTER' field in its output VCF file, but will keep all the candidate variants, including those that failed the filtering, in its output. So we use an <code>awk</code> command to extract all the variants that pass the filtering.

**5\.** At the end of step 2, we have VCF files for original candidates and candidates that pass the filtering. How do we check the number of candidate mutations in each file? We use these commands.

<code>grep -v "^#" PAP1_1401_01_T01.vcf | wc -l<br>
grep -v "^#" PAP1_1401_01_T01_passed.vcf | wc -l </code>{% include code-snippet-copy.html %}

<code>grep -v</code> will extract all lines that do NOT contain a specific pattern. The pattern we specify is "^#", which is lines starting with the "#" symbol. This will exclude all the header lines in the VCF file. Then we redirect the output to the command <code>wc -l</code> which calculates the number of lines.

---

### Annotation of somatic short variants

Now we have a list of variants ready for the downstream analysis. To prioritize the variants for downstream analysis. we will functionally annotate these variants. There are multiple applications for this analysis, and we are going to show two methods, [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) and [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/).

**1\.** We first annotate with Funcotator. Confirm that you are in the interactive session, and enter the following commands.

<code>cd /data/$USER/practical_session_4/shfile<br>
sh /data/classes/DCEG_Somatic_Workshop/Practical_session_4/shfiles/Step4_funcotator_annotation.sh  PAP1_1401_01_T01 ../Results ../Results</code>{% include code-snippet-copy.html %}

The script will take a few seconds. Let's first check the bash scripts.

<code>vi /data/classes/DCEG_Somatic_Workshop/Practical_session_4/shfiles/Step4_funcotator_annotation.sh</code>{% include code-snippet-copy.html %}

<!--{% include image-modal.html link="practical_assets/Step3_1.PNG" %}-->
```bash
#!/bin/bash

#### Annotate VCF file with Funcotator ###

module load samtools
module load GATK/4.3.0.0

GATK_Bundle=/fdb/GATK_resource_bundle/hg38-v0
GENOME=$GATK_Bundle/Homo_sapiens_assembly38.fasta
FUNCOTATORDB=/fdb/GATK_resource_bundle/funcotator/funcotator_dataSources.v1.7.20200521s

PREFIX=$1
INDIR=$2
DIR=$3
logs=$DIR/logs
IN_VCF=$INDIR/${PREFIX}_passed.vcf
OUT_ANNOT_MAF=$DIR/${PREFIX}_funcotator.maf

if [ ! -d "$DIR" ]; then
        mkdir -p $DIR
fi
if [ ! -d "$logs" ]; then
        mkdir -p $logs
fi

SECONDS=0

gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms10G -Xmx10G -XX:ParallelGCThreads=2" Funcotator -R $GENOME \
     -V $IN_VCF \
     -O $OUT_ANNOT_MAF \
     --output-file-format MAF \
     --data-sources-path $FUNCOTATORDB \
     --ref-version hg38 \
     2> ${logs}/funcotator-${PREFIX}.log

duration=$SECONDS
echo "Annotation completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
```

Funcotator is a tool in the GATK package. It requires five input files:
<code>-R</code> A reference genome sequence.
<code>-V</code> A VCF of variant calls to annotate.
<code>--data-sources-path</code> The path to data sources or required formats.
<code>--ref_version</code> The version of the reference genome sequence being used.
<code>--output-file-format</code> The desired output format for the annotated variants file (MAF or VCF).
We also specify the output file name with the <code>-O</code> option.

```bash
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms10G -Xmx10G -XX:ParallelGCThreads=2" Funcotator -R $GENOME \
     -V $IN_VCF \
     -O $OUT_ANNOT_MAF \
     --output-file-format MAF \
     --data-sources-path $FUNCOTATORDB \
     --ref-version hg38 \
     2> ${logs}/funcotator-${PREFIX}.log
```

All the data source files have already been downloaded in the Biowulf, so we specify the data source folder using this command: <code>FUNCOTATORDB=/fdb/GATK_resource_bundle/funcotator/funcotator_dataSources.v1.7.20200521s</code>
The script takes three arguments to generate the input/output paths and file names: the prefix of the input VCF file, the input and output directories.
The basic output files in VCF or MAF format contains all variants from the input file, with added annotation from the data sources.

We will run the second annotation package ANNOVAR and check the results together later.

**2\.** Next we annotate with ANNOVAR. The command lines are very similar:

<code>cd /data/$USER/practical_session_4/shfile<br>
sh /data/classes/DCEG_Somatic_Workshop/Practical_session_4/shfiles/Step4_annovar_annotation.sh  PAP1_1401_01_T01 ../Results ../Results</code>{% include code-snippet-copy.html %}

Next let's check the bash scripts.

<code>vi /data/classes/DCEG_Somatic_Workshop/Practical_session_4/shfiles/Step4_annovarannotation.sh</code>{% include code-snippet-copy.html %}

<!--{% include image-modal.html link="practical_assets/Step3_2.PNG" %}-->
```bash
#!/bin/bash

#### Annotate VCF file with ANNOVAR ###
module load annovar/2020-06-08

PREFIX=$1
INDIR=$2
DIR=$3
logs=$DIR/logs
IN_VCF=$INDIR/${PREFIX}_passed.vcf

if [ ! -d "$DIR" ]; then
        mkdir -p $DIR
fi
if [ ! -d "$logs" ]; then
        mkdir -p $logs
fi

echo $IN_VCF
echo $INDIR
echo $DIR
SECONDS=0

cd $DIR
convert2annovar.pl -format vcf4 $IN_VCF -includeinfo >${PREFIX}.avinput
table_annovar.pl  ${PREFIX}.avinput $ANNOVAR_DATA/hg38 \
	-buildver hg38 -out ${PREFIX} -remove \
	-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f \
	-nastring . -csvout -polish

duration=$SECONDS
echo "Annotation completed. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
```

The program <code>convert2annovar.pl</code> converts the VCF format into the ANNOVAR input format (.avinput file). For specification of the ANNOVAR input format, you may get more details at: [https://annovar.openbioinformatics.org/en/latest/user-guide/input/](https://annovar.openbioinformatics.org/en/latest/user-guide/input/).

```bash
convert2annovar.pl -format vcf4 $IN_VCF -includeinfo >${PREFIX}.avinput
table_annovar.pl  ${PREFIX}.avinput $ANNOVAR_DATA/hg38 \
	-buildver hg38 -out ${PREFIX} -remove \
	-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f \
	-nastring . -csvout -polish
```

The program <code>table_annovar.pl</code> is the core command for annotation. <code>-protocol</code> specifies the annotation sources to use:ExAC version 0.3 (referred as exac03), dbNFSP version 3.0a (referred to as dbnsfp30a), dbSNP version 147 with left-normalization (referred to as avsnp147). Each 'protocol' corresponds to one field in the <code>-option</code> argument: gene-based (referred to as g), gene-based with cross-reference annotation (referred to as gx), region-based (referred to as r) and filter-based (referred to as f). In the final output file, a tab-delimited file with many columns will contain the input columns and annotation columns, that correspond to the  combination of 'protocol' and 'option'. For more details, please refer to [https://annovar.openbioinformatics.org/en/latest/user-guide/startup/#table_annovarpl](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/#table_annovarpl).

**3\.** Now open the two output files for both software in Excel.

<code>PAP1_1708_01_T01_funcotator.maf<br>
PAP1_1708_01_T01.hg38_multianno.csv</code>

The Funcotator output is in the MAF format, the first several columns are the basic features of the mutation: genomic coordinates, reference and variant alleles, variant classifications and types. The variant classifications and types could be used to prioritize candidate driver mutations.

{% include image-modal.html link="practical_assets/Funcotator_1.PNG" %}

If we scroll to the right, there are annotations of genomic, cDNA and protein  changes.
**TODO: Mention the link to downstream analysis? Any other comments on the Funcotator results?**

{% include image-modal.html link="practical_assets/Funcotator_2.PNG" %}

The ANNOVAR output is in VCF format. The first several columns include basic information about the variants. The columns of 'Func.refGene' and 'ExonicFunc.refGene' could be used to prioritize candidate driver mutations. The column 'AAChange.refGene' provides similar information to the protein changes annotated by Funcotator.

{% include image-modal.html link="practical_assets/ANNOVAR_1.PNG" %}

For the next a few columns, the ExAC* columns are allele frequencies in all the samples and the subpopulations in the Exome Aggregation Consortium data sets. The column 'avsnp147' is the SNP ID in the dbSNP database. The other columns are prediction scores for the likelihood of non-synonymous variants to be deleterious using several popular tools, such as SIFT, PolyPhen2, HDIV,  LRT, and so on.

{% include image-modal.html link="practical_assets/ANNOVAR_2.PNG" %}

---  
