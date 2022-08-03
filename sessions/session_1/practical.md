---
layout: page
permalink: sessions/session_1/practical
menubar_toc: true
---

This practical section will focus first on connecting to a remote cluster
as well as using the linux command line. Then we will practice working with
some common bioinformatics file formats.

## Setup ssh connection

First, to use the computing clusters we need to establish an ssh connection
between your computer and the cluster. The steps to do this will differ
depending on whether you're using MacOS or Windows:

#### MacOS

- Open the Terminal app, like so:
	- "Finder -> Applications -> Utilities -> Terminal"
- enter the following command: <code>ssh username@helix.nih.gov</code>,
replacing 'username' with your own username.

#### Windows

- Install [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)
- Launch PuTTY. Under “Host Name (or IP address), type:
<code>username@helix.nih.gov</code>, replacing 'username' with your own username,
and click “Open”
- At the prompt, enter the account password

At this point you should be connected to the NIH Helix cluster, and your screen
should look something like this:

<img src="practical_assets/login_screen.png">

---

## Basic Bash/Linux commands

Let's practice some basic linux commands. :

1\. First, let's check what directory we're working in on the command line. Simply type:

<code>pwd</code>

and hit enter. This should output <code>/home/[your username]</code>.

2\. Since Biowulf home directories are relatively small we are going to be working in the /data/ directory, so let's make a folder there with the <code>mkdir</code> command:

<code>mkdir /data/[your username]/session1</code>

. To verify that your folder was created successfully, enter the following: <code>ls /data/[your username]</code>. This should list all contents of your /data/ directory, and among those contents should be the 'session1' folder.

3\. Let's move ourselves to that directory to begin working:

<code>cd  /data/[your username]/session1</code>

. Verify that you are in the right folder with the <code>pwd</code> command.

4\. Let's copy some files to practice with. Before we do that, let's look at the contents of the course folder like so:

<code>ls -lh /data/course</code>

. Now let's copy that data over to our current directory:

<code>cp /data/course/* .</code>

. The last two characters in the code above are special characters with specific meanings. The <code>*</code> in the code above is called a 'wildcard' and can be interpreted as 'everything', in this case meaning all files. The <code>.</code> in the code means 'current directory', and is where we're copying the files to.

5\. For the sake of practice let's merge gencode.v19.og.bed and gencode.v19.tsg.bed into a single file gencode.v19.driver.bed:

<code>cat gencode.v19.og.bed gencode.v19.tsg.bed > gencode.v19.driver.bed
</code>

. The <code>cat</code> command will print the contents of one or more files to your screen, but using the <code>></code> character we can redirect the output to a new file. To verify that we've pasted the files together, try using the <code>wc -l</code> command on <code>gencode.v19.og.bed</code>, <code>gencode.v19.tsg.bed</code>, and <code>gencode.v19.driver.bed</code> individually to count how many lines each one has.

6\. We're going to try sorting our new file, so before we do that let's check the manual on the <code>sort</code> command:

<code>man sort</code>

. Use the arrow keys or your scroll wheel to read, and take note of the options <code>-k</code> and <code>-r</code>. To exit the manual enter <code>:q</code>.

7\. Finally let's sort the merged bed file by genomic coordinates:

<code>sort -k1,1 -k2,2n gencode.v19.driver.bed > genelist.bed</code>

 <code>-k1,1 -k2,2n</code>: sort first by the first field, then by the second field numerically. The sorting isn't saved automatically, so as before we use the <code>></code> to redirect the output to a new file.

---

## Working with Fastq files
8\. Before we can work with our fastq files we need to decompress them:

<code>tar -xvf sample1.fq.tar.gz</code>

. The <code>-x</code> and <code>-f</code> options tell tar to extract files and what file to extract from, respectively. The <code>-v</code> option is to print filenames as they're decompressed.

9\. Let's preview our files:

<code>head sample1.fq</code>

. By default, head will you give you the first ten lines of a file to examine. If you want to see more, you can use the option <code>-n</code>, or use the <code>more</code> command:

more sample1.fq

. Use 'enter' or 'space' to scroll down, 'b' to scroll back, and '/' to search. To quit, either reach the end of the file or type <code>:q</code>.

10\. Let's subset our fastq file to a smaller collection of reads with the seqtk library. First load seqtk:

<code>module load seqtk</code>

, then use the 'subseq' tool to extract the reads specified in the 'name.lst' file:

seqtk subseq sample1.fq name.lst > out.fq

. Check the output file 'out.fq' with the <code>head</code> and/or <code>more</code> command.

---

## Working with BAM files

---

## Working with BED files
