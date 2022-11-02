---
layout: page
permalink: sessions/session_1/practical
menubar_toc: true
---
<script src="{{ site.baseurl }}/assets/js/vanilla-back-to-top.min.js"></script>
<script>addBackToTop()</script>

<script src="{{ site.baseurl }}/assets/js/copyCodeSnippet.js" defer></script>
<script src="{{ site.baseurl }}/assets/js/copyCodeBlock.js" defer></script>

This practical section will focus first on connecting to a remote cluster
as well as using the linux command line. Then we will practice working with
some common bioinformatics file formats.

---

## Setup ssh connection

First, to use the computing clusters we need to establish an ssh connection
between your computer and the cluster.

**The steps to do this will differ
depending on whether you're using MacOS or Windows:**

#### MacOS

- Open the Terminal app, like so:
	- "Finder -> Applications -> Utilities -> Terminal"
- enter the following command: <code>ssh USERNAME@biowulf.nih.gov</code>,
replacing 'USERNAME' with your own username.

#### Windows

- Download <a href="https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html" target="_blank">PuTTY</a>. The version 'putty.exe' is a standalone package that does not require installation.

- Launch PuTTY. Under “Host Name (or IP address), type:
<code>USERNAME@biowulf.nih.gov</code>, replacing 'USERNAME' with your own username,
and click “Open”

{% include image-modal.html link="practical_assets/putty_openConnect.png" classes="center" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

- At the prompt, enter the account password. Please note that the password will not show on the screen, and 5 failed login attempts since the last successful login are allowed.

At this point you should be connected to the NIH Biowulf cluster, and your screen
should look something like this:

{% include image-modal.html link="practical_assets/biowulf_connected.png" classes="center" styles="display: block;margin-left: auto; margin-right: auto; max-width:85%" %}

### Request an interactive session

You're now logged into the login nodes that are shared by everyone who uses Biowulf. It is important that you do NOT run intensive computing jobs on the login nodes to make them overburdened. If you do, other users will be temporarily unable to use Biowulf. You should instead request an interactive session, then you will be assigned to the compute nodes.

**To do this, type this command:**

<code>sinteractive --mem=4g --cpus-per-task=2 --time=4:00:00</code>{% include code-snippet-copy.html %}

**and hit enter.**

You will see this screen when your request has been submitted and is waiting for computing resources:

{% include image-modal.html link="practical_assets/pending_interactive.png" classes="center" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

When you are allocated to a compute node you should see USERNAME@cn####, like this:

{% include image-modal.html link="practical_assets/granted_interactive.png" classes="center" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

When you want to close the interactive session use the command <code>exit</code>, and you will return to the login nodes.

### Locally mounting HPC System Directories
We will need to mount your personal <code>/data</code> drive to access files on Biowulf from your local computer at the end of this session. See the instructions here for mounting your data directory: <a target="_blank" href="https://hpc.nih.gov/docs/hpcdrive.html">(https://hpc.nih.gov/docs/hpcdrive.html).

---

## Basic Bash/Linux commands


Now let's practice some basic linux commands:

**1\.** First, let's check what directory we're working in on the command line. Simply type:

<code>pwd</code>{% include code-snippet-copy.html %}

and hit enter. This should output <code>/home/[your username]</code>.

**2\.** Since Biowulf home directories are relatively small we are going to be working in the /data/ directory, so let's make a folder there with the <code>mkdir</code> command:

<code>mkdir /data/[your username]/session1</code>{% include code-snippet-copy.html %}

To verify that your folder was created successfully, enter the following: <code>ls /data/[your username]</code>. This should list all contents of your /data/ directory, and among those contents should be the 'session1' folder.

**3\.** Let's move ourselves to that directory to begin working:

<code>cd  /data/[your username]/session1</code>{% include code-snippet-copy.html %}

Verify that you are in the right folder with the <code>pwd</code> command.

**4\.** Let's copy some files to practice with. We deposited two copies of our sample data in Biowulf and Github, respectively. You can follow either **4a** or **4b** to copy the data to your home directory.

**4a\.** Copy the data from the shared folder in Biowulf.
Before we do that, let's look at the contents of the course folder like so:

<code>ls -lh /data/classes/DCEG_Somatic_Workshop/Practical_session_1</code>{% include code-snippet-copy.html %}

Now let's copy that data over to our current directory:

<code>cp -r /data/classes/DCEG_Somatic_Workshop/Practical_session_1/* .</code>{% include code-snippet-copy.html %}

The option <code>-r</code> means copy directories recursively (i.e. including all the contents in sub-directories). If we omit this option, all the folders in this directory will be skipped.

The last two characters in the code above are special characters with specific meanings. The <code>*</code> in the code above is called a 'wildcard' and can be interpreted as 'anything'. In words, we're instructing the computer to copy all files in <code>/data/classes/DCEG_Somatic_Workshop/</code> matching the pattern 'anything', or in other words all files.

The <code>.</code> in the code means 'the current working directory', and is where we're copying the files to.

**4b\.** Alternatively, you can download the data from GitHub.

Let’s first look up the link from GitHub. Open the course resources in web browser: <a href="https://github.com/NCI-ITEB/tumor_epidemiology_approaches_materials" target="_blank">https://github.com/NCI-ITEB/tumor_epidemiology_approaches_materials</a>

Click the green button ‘code’. Click the button next to the https link to copy the link.

{% include image-modal.html classes="center" styles="display:block; max-width:50%; margin-left:auto; margin-right:auto" link="practical_assets/github_link.png" %}

With the link now copied to your clipboard, return to the Biowulf terminal and enter <code>git clone</code> followed by the link you've copied:

<code>git clone https://github.com/NCI-ITEB/tumor_epidemiology_approaches_materials.git</code>{% include code-snippet-copy.html %}

**5\.** Now we can switch to the directory of the sample input data.
<code>cd  sample_input_data</code>{% include code-snippet-copy.html %}

There is a BED file with transcript coordinates in chromosome 22 (gencode.hg38.chr22.bed). Let's check the first few lines of this file:

<code>head gencode.hg38.chr22.bed</code>{% include code-snippet-copy.html %}

By default, <code>head</code> will give you the first ten lines of a file to examine. If you want to see more, you can use the option <code>-n</code>, or examine the file manually with the <code>more</code> command.

{% include image-modal.html link="practical_assets/preview_bed.png" classes="center" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}


**6\.** Finally let's try sorting this BED file (gencode.hg38.chr22.bed) by genomic coordinates. Note that many commands/applications will require the input files (in BAM or BED formats) to be sorted by genomic cooridinates.

But before we do that let's first check the manual on the <code>sort</code> command:

<code>man sort</code>{% include code-snippet-copy.html %}

Use the arrow keys or your scroll wheel to read, and take note of the options <code>-k</code> and <code>-r</code>. To exit the manual enter <code>q</code>.

Then we can sort the file with this command:

<code>sort -k1,1 -k2,2n gencode.hg38.chr22.bed >gencode.hg38.chr22.sorted.bed</code>{% include code-snippet-copy.html %}

<code>-k1,1 -k2,2n</code>: sort first by the first field, then by the second field numerically. The sorting isn't saved automatically, so as before we use the <code>></code> to redirect the output to a new file.

---

## Working with FASTQ files
**7\.** Before we can work with our fastq files we need to decompress them using the decompression command <code>tar</code>:

<code>tar -xzvf Sample1.tar.gz</code>{% include code-snippet-copy.html %}

The <code>-xz</code> and <code>-f</code> options tell <code>tar</code> that we want to decompress our files and what file to extract from, respectively. The <code>-v</code> option is to print filenames as they're decompressed.

**8\.** Let's preview our files. We can do it in two ways:

<code>head Sample1_1.fq</code>{% include code-snippet-copy.html %}

We have used the <code>head</code> command before. This will return the first ten lines of the file by default.

Alternatively, we can use:

<code>more Sample1_1.fq</code>{% include code-snippet-copy.html %}

Use 'enter' or 'space' to scroll down, 'b' to scroll back, and '/' to search. To quit, either reach the end of the file or type <code>q</code>.

We have another compressed format for the same file (files with extension .fastq.gz). We can check the contents of these files directly without creating a plain text file, using this command:
<code>zless Sample1_1.fastq.gz</code>{% include code-snippet-copy.html %}
Same to the command <code>more</code>, use ‘space’ to scroll down, and type <code>q</code> to exit.

To check the number of paired-end reads, use the following commands:

<code>echo $(zcat Sample1_1.fastq.gz | wc -l)/4 | bc</code>{% include code-snippet-copy.html %}

<code>echo $(zcat Sample1_2.fastq.gz | wc -l)/4 | bc</code>{% include code-snippet-copy.html %}

You should see '300,000’ for both commands.

To briefly explain these commands:

- <code>zcat Sample1_1.fastq.gz | wc -l</code>: print the decompressed contents of the gzipped (‘.gz’) reads file using ‘zcat’, then count the number of lines by feeding (or piping, the ‘\|’ symbol) the output directly to ‘wc -l’
- <code>echo $(...)/4 | bc</code>: setup a division by 4 of the line count value from the <code>wc -l</code> command using <code>echo</code> (recall that each sequence in the fastq is composed of four lines). This mathematical expression is only evaluated after using the command <code>bc</code>.

**9\.** Let's subset our fastq file to a smaller collection of reads with the seqtk library. First load seqtk:

<code>module load seqtk</code>{% include code-snippet-copy.html %}

**10\.** then use the 'subseq' tool to extract the reads specified in the 'name.lst' file:

<code>seqtk subseq Sample1_1.fq name.lst >out.fq</code>{% include code-snippet-copy.html %}

‘name.lst’ is a list of identifiers of reads.

<code> head name.lst</code>{% include code-snippet-copy.html %}

Check the output file 'out.fq' with the <code>head</code> and/or <code>more</code> command.

---

## Working with BAM files

**11\.** We're going to be working with BAM files from 1000 genomes. The example data was downloaded from the [1000 Genomes Project](https://www.internationalgenome.org/home).

{% include image-modal.html link="practical_assets/1000_genomes_download.png" classes="center" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

We used the following command to download the example in this session:

<code>samtools view http://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/data/ERR3240193/HG00118.final.cram chr22:38952741-38992778</code>

And the CRAM file could be converted to BAM using this command:

<code>samtools view -h -b -@ 10 -T Homo_sapiens_assembly38.fasta -o HG00118__chr22@38952741@38992778.bam HG00118__chr22@38952741@38992778.cram</code>

**12\.** Now we rename the BAM file in our input data folder to reads.bam.

<code>mv HG00118__chr22@38952741@38992778.bam reads.bam</code>{% include code-snippet-copy.html %}

---

### Examine file format

**13\.** Let's look at the header of this BAM file. Try listing all the versions of samtools available on Biowulf:

<code>module spider samtools</code>{% include code-snippet-copy.html %}

Some modules will have different parameters in different versions. By default, usually it will load the latest version of the module. For our purposes any recent version will suffice, so enter <code>module load samtools</code> to load the default version of samtools. Then to view the BAM header, enter:

<code>samtools view -H reads.bam | more</code>{% include code-snippet-copy.html %}

The <code>view</code> mode of samtools is a tool to print sections of a BAM/SAM/CRAM file, and the <code>-H</code> option instructs samtools to print only the header. In this example we then feed the samtools output directly to the <code>more</code> command via the linux 'pipe' (the <code>|</code> symbol) so it's easier to read and browse.

Note that we've been using a BAM file which is in binary format, but the output is in readable text. Samtools has converted the output from BAM to SAM automatically.

If we want to extract specific fields in the header, we could use the <code>grep</code> command. For example, this command would extract the lines for read groups in the header.

<code>samtools view -H reads.bam|grep "^@RG"</code>{% include code-snippet-copy.html %}

The <code>^</code> symbol will match all lines starting with a specific pattern. In this example, it matches all lines starting with "@RG". Purely for your interest, these are the read group definition lines.

**14\.** Now let's view the first 20 lines of the aligned reads in the body section. The command will be very similar, but without the <code>-H</code> option samtools will ignore the header:

<code>samtools view reads.bam|head -20</code>{% include code-snippet-copy.html %}

As before, we use the pipe to feed the output into the <code>head -20</code> command so we can see just the first 20 lines.

---

### Extract, sort, index reads

**15\.** Let's isolate only unmapped reads, and reads with unmapped mates:

<code>samtools view -b -f4 reads.bam >unmapped.bam</code>{% include code-snippet-copy.html %}

<code>samtools view -b -f8 reads.bam >mate_unmapped.bam</code>{% include code-snippet-copy.html %}

The <code>-b</code> flag tells samtools to output in the compressed BAM format rather than SAM, and is very important when working with large alignment files. The <code>-f</code> flag requires all output reads to have the specified alignment flags, in this case flag '4' and flag '8'. These correspond to 'read unmapped' and 'mate unmapped', respectively.

A full list of SAM flags can be found at: <a href="https://broadinstitute.github.io/picard/explain-flags.html" target="_blank">https://broadinstitute.github.io/picard/explain-flags.html</a>

**16\.** Sort the original file by genomic coordinates and output into file reads_sorted.bam:

<code>samtools sort -o reads_sorted.bam reads.bam</code>{% include code-snippet-copy.html %}

Note that in the previous command we used <code>></code> to save a new file whereas in this case we've used <code>-o [file name]</code> to accomplish the same.

**17\.** Let's now index the sorted file. Indexing allows for more efficient lookup of reads and is required to run many bioinformatics algorithms:

<code>samtools index reads_sorted.bam</code>{% include code-snippet-copy.html %}

This should create a new index file <code>reads_sorted.bam.bai</code>.

**18\.** Take a look at some of the alignment statistics using:

<code>samtools flagstat reads_sorted.bam</code>{% include code-snippet-copy.html %}

**19\.** Find all reads mapping to chr22:38,700,000-39,300,000 and save to file chr22.bam:

<code>samtools view -h -b reads_sorted.bam chr22:38,700,000-39,300,000 > chr22.bam</code>{% include code-snippet-copy.html %}

The <code>-h</code> option will retain the original header in our output file.

<!--**20\.** We can visualize the alignment at XXXX using tview:

<code>samtools tview XXXX</code>{% include code-snippet-copy.html %}

Note that the program IGV is much more useful for this purpose with more features, but we will not cover it today.-->

---

## Working with BED files

**20.** Let's convert our chr2 BAM alignment file to BED format. First load bedtools and then use the 'bamtobed' mode:

<code>module load bedtools</code>{% include code-snippet-copy.html %}

<code>bedtools bamtobed -i chr22.bam > reads.bed</code>{% include code-snippet-copy.html %}

**21\.** For each gene that overlaps with alignments, report the base-pair overlap between the sequence alignment and genes. Here we can use the ‘reads.bed’ file to extract all the regions with alignments.

<code>bedtools intersect -a reads.bed -b gencode.hg38.chr22.sorted.bed >intersect_overlap.bed</code>{% include code-snippet-copy.html %}

<code>-a</code> and <code>-b</code> specify our two BED file inputs.

**22\.** For each gene that overlaps with alignments, report the full-length of the original gene. Report each gene with more than one hit only once.

<code>bedtools intersect -a gencode.hg38.chr22.sorted.bed -b reads.bed -wa -u >intersect_full_length_genes.bed</code>{% include code-snippet-copy.html %}

Compare to **21**, note that we have a <code>-wa</code> option in **22**. See the diagram below for the specifics on bedtools intersect.

{% include image-modal.html link="practical_assets/bedtools_intersect.png" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

With the option <code>-u</code>, each read with more than one hit will be reported only once.

**23\.** Report regions in genes that have no overlap with alignments (specified with <code>-v</code>):

<code>bedtools intersect -a gencode.hg38.chr22.sorted.bed -b reads.bed -v  >intersect_no_overlap.bed</code>{% include code-snippet-copy.html %}

**24\.** For a more advanced query, we can do the following: report all reads within 2000bp upstream or 1000bp downstream of genes.

<code>bedtools window -a reads.bed -b gencode.hg38.chr22.sorted.bed -l 2000 -r 1000 -u > intersect_reads_window.bed</code>{% include code-snippet-copy.html %}

The options <code>-l</code> and <code>-r</code> ('left' and 'right') define the base pair window upstream and downstream of the overlapping regions, respectively.

---

### Visualizing on UCSC

**25\.** We're going to visualize these reads on UCSC, and to do so we need to add some header lines to our BED file. Run the following series of commands (red text only):

- Configure browser:

{% include code-block-copy.html %}```
printf "browser position chr22:38,700,000-39,300,000\nbrowser hide all\n" > custom_UCSC_track.bed
```
<br>

- Add the track for overlapping regions:

{% include code-block-copy.html %}```
(printf "track name=\"overlap regions\" description=\"example for bedtools A intersect B\" visibility=1 color=0,0,255 useScore=1\n#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n" && cat intersect_overlap.bed)  >> custom_UCSC_track.bed
```
<br>

- Add the track for full length of genes:

{% include code-block-copy.html %}```
(printf "track name=\"original genes\" description=\"example for bedtools A intersect B -wa\" visibility=3 color=255,0,0 useScore=1\n#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n" && cat intersect_full_length_genes.bed)  >> custom_UCSC_track.bed
```
<br>

The <code>printf</code> and <code>cat</code> commands simply print some text which we then save to a file using <code>></code>. Take a look at this header we've just created using <code>head custom_UCSC_track.bed</code>{% include code-snippet-copy.html %}.

Note one very important detail in the previous commands: <code>>></code> will append text to the end of an existing file while <code>></code> will overwrite existing files. **When working with files of your own, be very careful of this difference or you could accidentally lose data!**

**26\.** Let's now visualize using the UCSC genome browser. Go to <a href="https://genome.ucsc.edu/"  target="_blank">https://genome.ucsc.edu/</a>. Under the "Genomes" tab, select "Human GRCh38/hg38" and then click the 'add custom tracks' button on the bottom of the genome browser.

{% include image-modal.html link="practical_assets/ucsc_select_genome.png" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

Next, upload the BED file via the "Choose File" button

{% include image-modal.html link="practical_assets/upload_BED.png" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

{% include image-modal.html link="practical_assets/upload_BED_2.png" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

and finally hit "go". An APOBEC3B homozygous deletion is highlighted below.

{% include image-modal.html link="practical_assets/custom_track_apobec.png" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}

<!--{% include image-modal.html link="practical_assets/custom_track_apobec2.png" styles="display: block;margin-left: auto; margin-right: auto; max-width:75%" %}-->
