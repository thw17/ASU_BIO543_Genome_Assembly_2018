# Reproducible Genome Assembly and Variant Calling
A hands-on tutorial introducing users to reproducible reference-based genome assembly and variant calling.

This tutorial has been tested on a Mac and Linux operating systems and will assume you're working with one of them.  It assumes very basic command line knowledge.
If you're comfortable with the basics of changing directories (e.g., ```cd```), listing the contents of directories (e.g., ```ls```), and figuring out where you are
(e.g., ```pwd```), then you should be fine following along.  If you're not comfortable with these commands, the free command line tutorials at [Codecademy](https://www.codecademy.com/learn/learn-the-command-line) will get you up to speed.  [This two page cheat sheet](https://static1.squarespace.com/static/56a2965eab2810d951ae42d9/t/57e9b51b1b631b22b84d3bd2/1474934043330/UnixNotes.pdf) also contains all of the commands you'll need and many, many more as well.

## Notation
In this tutorial, all commands entered will follow a ```$ ```, while any associated output will follow this line in the
same box.  Do not include the ```$``` in your command, rather enter the command that follows.

## Setting up your environment
For today's tutorial, you'll need this repository and Anaconda.

#### Getting the repo
We'll use ```git``` to clone the repository (repo) for this tutorial onto your computer.  You can check to see if you have git installed by typing the command ```$ git```.  You should see some usage information.  If not, see [here](https://git-scm.com/) for information about installing.

Once git is installed, move to the directory on your computer or cluster where you'd like to work and type the command:

```
$ git clone https://github.com/thw17/ASU_BIO543_Genome_Assembly_2018
```
This should create a directory called ASU_BIO543_Genome_Assembly_2018 containing all of the contents of this repo.

Alternatively, if git isn't working for you, you can directly download the zipped directory via the green "Clone or download" button on [the repository's website](https://github.com/thw17/ASU_BIO543_Genome_Assembly_2018).

#### Setting up Anaconda
[Anaconda](https://docs.continuum.io/anaconda/) is an environment and package manager for the programming language Python and it makes installation, environment management, etc. simple without requiring root or administrator privileges.  Fortunately, its framework has been leveraged for a project called [Bioconda](https://bioconda.github.io/) that extends these capabilities to external programs as well.  All of the packages and programs we're using today can easily be installed with Anaconda/Bioconda with the following steps:

* First, install Python 3.6 version of Miniconda [available here](https://conda.io/miniconda.html), *OR* Anaconda [available here](https://www.continuum.io/downloads). Miniconda installs the conda framework without all of the Python packages installed with Anaconda (numpy, scipy, matplotlib, etc.). All of these packages can be easily installed later (via ```conda install <package name>```), so the decision is up to you. During installation, be sure to allow Miniconda/Anaconda to append to your .bashrc or .bash_profile (this will add it and all programs it installs to your PATH). *This means you'll have to pay attention to all prompts during installation!!*  If installation goes well, the command ``` which python ``` should result in something like ```/Users/<yourusername>/miniconda/bin/python ``` or ```/home/<yourusername>/miniconda/bin/python ```.

* Add bioconda channels to conda with the following commands:
  ```
  $ conda config --add channels defaults

  $ conda config --add channels conda-forge

  $ conda config --add channels bioconda
  ```
* Create the environment we'll be working in and install required packages with the command:

  ```
  $ conda create --name ASU_Genomics_Tutorial python=3.6 snakemake fastqc bwa samtools freebayes bcftools bioawk bedtools
  ```
* Load the new environment and add the samblaster package

  ```
  $ source activate ASU_Genomics_Tutorial

  $ conda install -c biobuilds samblaster
  ```

This will create a working environment called ASU_Genomics_Tutorial containing python 3.6 (python 3 is required for snakemake) and all of the tools listed in the command.  You can see the full list of programs available through bioconda [listed here](https://bioconda.github.io/) and the full list of python packages available through Anaconda [listed here](https://docs.continuum.io/anaconda/pkg-docs).

The reason for installing samblaster separately in the final command is that we have to get it from the biobuilds channel.  Bioconda does support samblaster, but only on Linux.  So to ensure that your environment works across both Mac and Linux operating systems (for the purpose of this tutorial) and to avoid confusion over the source of other programs that might be available from biobuilds and bioconda, we're installing the from biobuilds only in this single instance.

If you want to load our new environment, you can do so with the command:
```
$ source activate ASU_Genomics_Tutorial
```
and leave the environment with the command:
```
$ source deactivate
```

If you're in your environment, you can easily add additional programs (like we did for samblaster) and packages with the command:
```
$ conda install <program/package name>
```

For example, if we also want to take a look at Bowtie2, another read mapper (we'll use bwa today), we can easily add it by entering our environment ``` $ source activate ASU_Genomics_Tutorial ``` and entering the command ```$ conda install bowtie2 ```

## Reference-based genome assembly
There are, in general, two main flavors of genome assembly.  _De novo_ assembly involves taking raw sequencing reads and piecing them together into a genome assembly using only the information contained in the reads and their metadata (e.g., the sequences themselves, insert sizes, etc.).  While a number of _de novo_ assemblers exist and there's a great deal of work being done to improve algorithms, lengthen sequencing reads, develop methods to increase insert sizes, etc., _de novo_ assembly remains challenging, expensive (usually 100x or greater sequencing depth), and computationally demanding.  Fortunately, if we have a reference genome available to us, we can make do with much less sequencing (often 30x coverage or less; low coverage - 1-5x - are not uncommon for some purposes), and use tools that require far less memory and storage.

In this tutorial we'll walk through the basics of reference-based genome assembly.  While the dataset we're working with is tiny (we're using the human mitochondrial genome and a tiny subset of reads from the 1000 genomes project), you should be able to use this as a starting point for working with larger datasets down the road.

### What you you'll need
In the simplest cases, you need:
* **a reference genome assembly (in fasta format)**
* **sequencing reads (in fastq format - many processes are quicker if files are gzipped as well)**
* **a computing environment with adequate storage and memory, and required software installed**

#### Reference genome
This is the pre-prepared reference genome assembly that you're going to use to map the sequencing reads from your project/experiment.  In a perfect world, this assembly is of a high-quality, has a good set of annotations available (e.g., genes, functional elements, repeats, etc.), and is relatively closely related to the species that you're studying.  There are, for example, numerous reference genomes hosted at the [UCSC Genome Browser](http://hgdownload.soe.ucsc.edu/downloads.html) and [Ensembl](http://www.ensembl.org/info/data/ftp/index.html).  Accessing a reference genome probably won't be a problem if you're working with a model organism, but in other situations you'll have to consider whether a good assembly is available for a taxon evolutionarily close enough for your purposes (if not, you might need to think about assembling a reference for your project _de novo_).  For today, we're working with example sequencing reads from human samples and we have the human reference genome available, so we'll be fine.

Fasta format is a file format designed for DNA or protein sequences.  It looks something like (the first 10 lines of the ``` human_g1k_v37_MT.fasta``` file in the ```references``` directory:
  ```
  $ head reference/human_g1k_v37_MT.fasta

  >MT
  GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
  CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
  GCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATT
  ACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
  ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCA
  AACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAA
  ACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCAC
  TTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
  CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATA
  ```
Here, the name of the sequence, ```MT``` is given after ```>```.  The lines that follow contain the sequence itself.  Because most fasta files wrap lines every 50-80 characters (this isn't uniform across files, unfortunately), there will often be many lines composing each sequence.  Today's (```human_g1k_v37_MT.fasta```) file should only contain a single sequence, the 1000 genomes reference MT sequence, that's 16-17kb in length.  We can quickly check to make sure using (the very, very powerful) [bioawk](https://github.com/lh3/bioawk), which we installed earlier:
  ```
  $ bioawk -c fastx '{print ($name), length($seq)}' reference/human_g1k_v37_MT.fasta

  MT	16569
  ```
We see that we do indeed have a single sequence called "MT" that's 16,569 bases in length.  

#### Sequencing reads
There are a few types of sequencing technologies out there, but Illumina is still dominant, so that's what we'll focus on today.  We won't get into the nuts and bolts of sequencing itself (there are plenty of resources available online that describe it well like [this video](https://www.youtube.com/watch?v=fCd6B5HRaZ8) or [this review](http://www.nature.com/nrg/journal/v17/n6/abs/nrg.2016.49.html)).

If you've sent your samples to a core for sequencing, they'll likely return to you a series of FASTQ files.  If you used single-end sequencing, your files can be concatenated into a single FASTQ file per sample per lane (note that you can easily concatenate files, even if they are gzipped, using ```cat```). On the other hand, if you used paired-end sequencing, you'll end up with two FASTQ files per sample per lane - one for the forward read and one for the reverse read (see the video I linked to in the previous paragraph for more information about paired-end sequencing).  It's generally very important that paired reads are in the same order in both files.  If you're getting reads directly from a sequencing center, they're probably already organized this way.  However, you might have to sort and re-pair reads if you have, for example, stripped reads from a bam file containing an alignment (the README in the ```fastq``` directory explains how to do this, if needed).

For today's tutorial, we have paired-end reads from two individuals: ind1 and ind2.  Forward reads are in files ending in "1.fastq.gz" and reverse reads are in the two files ending in "2.fastq.gz".  We can take a look at the first two reads in the file ```ind1_1.fastq.gz``` using ```zcat``` on Linux or ```gzcat``` on Mac/Unix combined with ```head```:
  ```
$ zcat fastq/ind1_1.fastq.gz | head -n 8

@H7AGFADXX131213:1:2110:18963:43686
ATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGAC
+
>>=>>>???>??=?<?=?>>><@>=@<?>>7>><?><?<7?<@<7?>=?><8>@>@>?>=>>=A@>?>=B?????=B???@A@??>B@AB@@@?@A@?>@?>A@?;;@?@A@@@@>@>?=@@@?@@?<>?A@@A>@?A??@AA=@=>@A??@>>@A@>@BA>4@@AAAA??>@????@=B??@A?>??@?>>5>?:??AA??A?@>AB@?@A@@?A?@99C?@@@=@@9<B?@>?@?34@A?>@@@=199
@H7AGFADXX131213:2:1212:12308:62676
CATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTT
+
>?<>>>@?>?>?<>=?>>>>>>>?>>?A>;??????<?>?><??>?>?>>>=>??>A@>?@<?A???>@=>>??>>?@=@@=@=AA>@@=A>AA?A?@?A@>A@>A?@A@>?@??@=??@A>=?>@?@A?@??@>@A???@??@=@<@???@>@>@=A@A>>3?A?A???A??:A>@BA;>@9?;><4;@?@=<5:A=A?=@?@?A?@@AA@@A=@?@:@A?@>8A?8@@?=?AAAA?@7>A?B@?=770
```
In these files, each sequencing read is listed in a series of four lines:
* the sequence identifier line, beginning with @
* the sequence line (consisting of A, T, C, G, or N)
* a comment line (here, the comment lines only contain +)
* a quality score line (ASCII characters)

The sequence identifier can contain a lot of information (see [Illumina's description for more information](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm)), the combination of which will identify individual reads uniquely.  In this example, we don't have records for the instrument or run number, and instead have IDs organized as: flowcellID:lane:tile:x_position:y_position.

The sequence line here is 250 bases long and contains the nucleotide sequence corresponding to each read.

The comment line usually either is a lone + or a + and then the sequence identifier repeated.  You'll generally ignore this line, but it's good to be aware of what's on it in case you're counting the number of reads belonging to a certain tile, for example.  In this case if the sequence identifier is repeated and you're counting the number of times a tile number is present in a file using a command like ```grep```, your count will be double the actual number of reads coming from that tile.

The quality score line contains an ASCII character for every nucleotide in the sequence (i.e., it'll be 250 characters long for a 250 base read, 75 characters long for a 75 base read, etc.).  Each ASCII character can be converted into a integer PHRED quality score, ranging from 0 to 93, that indicates the probability that a particular base call is incorrect.  [See more details here - our Illumina scores use the 33 offset](http://www.drive5.com/usearch/manual/quality_score.html).

If we want to quickly take a look at a summary of the reads in our fastq, we can use ```fastqc```. To analyze all of the files in the ```fastq``` directory and output the results to the ```stats``` directory, enter the following command from our main directory:
  ```
  $ fastqc -o stats fastq/*.fastq.gz
  ```
and it will quickly generate a series of reports for each fastq file that are summarized in the .html files (one per fastq) that you can view using your browser (e.g., Firefox, Chrome, Safari, etc.).

We can also use ```bioawk``` to parse and analyze fastq files as well.  For example we can count all of the reads in each file:
  ```
  $ for i in fastq/*.fastq.gz; do echo $i; bioawk -c fastx 'END{print NR}' $i; done

  ind1_1.fastq.gz
  5000
  ind1_2.fastq.gz
  5000
  ind2_1.fastq.gz
  5000
  ind2_2.fastq.gz
  5000
  ```

This is an example of a "for loop" in BASH.  The ```echo``` command is simply printing each file name and then the ```bioawk``` command counts and prints the total number of records in each file.

We can also count the number of reads from tile 2111 in ```ind1_1.fastq.gz```:
  ```
  $ bioawk -c fastx '{print $name}' fastq/ind1_1.fastq.gz | grep ':2111:' | wc -l

  80
  ```
Or count the number of reads from each tile in ```ind1_1.fastq.gz```:
  ```
  $ bioawk -c fastx '{print $name}' fastq/ind1_1.fastq.gz | cut -d':' -f 3 | sort | uniq -c

  ```
#### Programs
Reference-based assembly requires a number of tools for steps including (but not limited to) read trimming, read mapping, bam processing (sorting, removing duplicates, adding read group information, etc.), variant calling, and variant filtering.  We installed everything we need for today's tutorial with the conda command above.  I'll discuss these steps in more detail below.

### Assembly: Step-by-step
For today's tutorial, the reference genome is in the ``` reference ``` directory of this repository, the sequencing reads are in two files for each sample in the ``` fastq ``` directory of this repository, and the "Setting Up Anaconda" section above should take care of the software requirements.  Because we're working with a small dataset today, we won't need too much in the way of memory/storage, but bigger projects will often require at minimum a high-memory computer, but more likely high-performance computing clusters, dedicated servers, or online services such as Amazon Web Services or [Galaxy](https://usegalaxy.org/).

So now we'll walk through the major steps of reference-based genome assembly and build a pipeline along the way.

#### Preparing your reference
Most reference genomes are quite large, so it's very inefficient to linearly search through, say, 3 billion characters spread across 25 (or MANY more) sequences.  So, many programs use various hashing/indexing strategies to more efficiently handle the data. We'll create two of the most commonly required index files - .dict and .fai - which will summarize the reference sequences.  We'll also create the required index for ```bwa```, our read mapper, which will allow ```bwa``` to quickly search for matches in the reference.  This is all we'll need for our purposes today, but check any additional tools you incorporate into your work down the line to see if they require additional indexing or processing of reference assemblies.  For example, each read mapper will likely require its own unique index, so your ```bwa``` indices won't work for, say, ```bowtie2```.

We'll first move into our ```reference``` directory:
  ```
  $ cd reference
  ```

and run the following three commands:
  ```
  $ samtools faidx human_g1k_v37_MT.fasta

  $ samtools dict -o human_g1k_v37_MT.dict human_g1k_v37_MT.fasta

  $ bwa index human_g1k_v37_MT.fasta
  ```
This will be quick and require very little memory on our small reference, but the bwa indexing in particular can take much longer on a full, human-sized reference and require ~4-6 GB of memory.

And that's it.  All of our reference files and indices are now contained in our reference directory.  Let's now move back to our main directory:

  ```
  $ cd ../
  ```

#### Fastq Quality Control
The first thing you should do when you get fastq (sequencing read) files is get a sense of their quality.  The quickest way to do this is to run ```fastqc``` and take a look at the reports.  We can run fastqc on every fastq file in ```fastq``` and output reports into ```stats``` with the command (from our main directory).  If you didn't earlier, run the following command:
  ```
  $ fastqc -o stats fastq/*.fastq.gz
  ```
This should take less than a minute to complete and output a .zip and .html file for each example fastq (it'll take longer for full-sized files).  You can open the .html locally on your computer in your browser.

If your sequences are of an unexpectedly low quality it might be worth contacting your sequencing center.  Otherwise, lower quality towards the ends of the reads, some PCR duplication, and sequence content that's slightly off are all pretty typical in sequencing experiments.  Also note that certain sequencing experiments, for example those targeting coding regions like exome sequencing and RNA seq, are likely to have some odd kmers.  If you're interested in trimming low-quality ends of reads or removing sequencing adapters, take a look at programs like [BBDuk - part of the BBTools suite](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) or [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).  Many read mappers like BWA handle adapters and low quality ends very well and downstream variant callers will build base quality into their genotype likelihoods, so it's often not necessary to do any preprocessing if your goal is variant calling.  For today's purposes, we'll take this approach and move onto our next step.  However, in your future experiments, you should carefully consider the costs/benefits of trimming.  Adapters, in particular, are often worth trimming and can be discovered with tools such as [BBMerge - also part of the BBTools suite](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/).

#### Read mapping
Our next step involves mapping our reads to our reference.  We'll use the [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) algorithm to do this, as it's among the most popular mappers and works very well mapping reads to a closely related reference.  Other popular mappers include [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Novoalign](http://www.novocraft.com/products/novoalign/), and [Stampy](http://www.well.ox.ac.uk/stampy) (Stampy, in particular, for mapping to a very evolutionary diverged reference genome).

The command line for bwa mem is quite straightforward.  From the main directory (/path/to/ASU_Genomics_Workshop_Summer2017), we can execute the following command:
  ```
  $ bwa mem -M reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz > bam/ind1.sam
  ```
This will map both sets of paired-end reads from ind1 to our reference genome and output the alignments in SAM format in our ```bam``` directory.  The ```-M``` flag ensures compatibility of the alignments with downstream tools.

If we take a look at the top of the file using ```$ head -n 8 bam/ind1.sam```, we see the following:
  ```
  @SQ     SN:MT   LN:16569
  @PG     ID:bwa  PN:bwa  VN:0.7.15-r1140 CL:bwa mem -M human_g1k_v37_MT.fasta ind1_1.fastq.gz ind1_2.fastq.gz
  H7AGFADXX131213:1:2110:18963:43686      99      MT      6969    60      250M    =       7187    468     ATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGAC      >>=>>>???>??=?<?=?>>><@>=@<?>>7>><?><?<7?<@<7?>=?><8>@>@>?>=>>=A@>?>=B?????=B???@A@??>B@AB@@@?@A@?>@?>A@?;;@?@A@@@@>@>?=@@@?@@?<>?A@@A>@?A??@AA=@=>@A??@>>@A@>@BA>4@@AAAA??>@????@=B??@A?>??@?>>5>?:??AA??A?@>AB@?@A@@?A?@99C?@@@=@@9<B?@>?@?34@A?>@@@=199      NM:i:0  MD:Z:250        AS:i:250        XS:i:0
  H7AGFADXX131213:1:2110:18963:43686      147     MT      7187    60      250M    =       6969    -468    ACACTTCCTCGGCCTAACCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCCCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGACCCCGTATACATA      #################?6=,?<:4<821+08=@???B9<>==77?=<6+<8'?=;9:9<==?@=<5<>==<>?>,*<=?=>>>?=<>@>?A?>?@@@>?=>><>@>>?<A@@?<<;7<=<,??=@<=@<?@@?<>>?>>7<?A8><>;87>B@6=>@@@=@===<+83(=?9598&0;91%<;=2<=79=<809;>78;3;:>><98>80;9:/45==<69<16<8;=>64;);<=*>?'<=0<:;<;9      NM:i:4  MD:Z:6T9T164T55A12      AS:i:230        XS:i:0
  H7AGFADXX131213:2:1212:12308:62676      99      MT      3728    60      250M    =       3978    500     CATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTT      >?<>>>@?>?>?<>=?>>>>>>>?>>?A>;??????<?>?><??>?>?>>>=>??>A@>?@<?A???>@=>>??>>?@=@@=@=AA>@@=A>AA?A?@?A@>A@>A?@A@>?@??@=??@A>=?>@?@A?@??@>@A???@??@=@<@???@>@>@=A@A>>3?A?A???A??:A>@BA;>@9?;><4;@?@=<5:A=A?=@?@?A?@@AA@@A=@?@:@A?@>8A?8@@?=?AAAA?@7>A?B@?=770      NM:i:0  MD:Z:250        AS:i:250        XS:i:0
  H7AGFADXX131213:2:1212:12308:62676      147     MT      3978    60      250M    =       3728    -500    CATAGCCGAATACACAAACATTATTATAATAAACCCCCTCACCACTACAATCTTCCTAGGAACAACATATGACGCCCTCTCCCCTGAACTCTACACAACATATTTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTATGAATTCGAACAGCATACCCCCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATA      @>=B?>95@>?@>@?@?@>??>;??>+@@@7@?8&?@@@?>?@@><=@A@B?@@><?A?>@@<:;7=826?=6//(;1;>>=>;>=?==>@>?@@?A@@>@@@???<:?=>?@@=0==>>>?@?>>@>@>A@@?@=@@?@@?=@A@AA6?A@@A?@?>;??@?5??A??9?><:9?@=>>>?@>@?>>??>>=?>>>??????>>??=>>>><7>>>8>>>>@=>?>=>>?>>>>>>>>>?<@>@><;<<      NM:i:2  MD:Z:34A40A174  AS:i:240        XS:i:0
  H7AGFADXX131213:2:1114:9857:73128       73      MT      15987   60      250M    =       15987   0       CACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTCCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATTCATCCAGTACATTACTGCCAGCCACCATGAATATTTTACGGTACCATAAATACTTGACCACCTGTAGTACATCAAAACCCACTCCACATCAAAACCCCCTCCCCATGCTTACTACCTAGTACATCAATCACTCCACAACTATCACAC      ;?;>?0??>?@<@>+>7<?>1/?@=??<A7>?>?>?:8/-B<2?1=-*((:+9??=A>@8=;&5;<A0:>?@9-3;7<=A>A@@*?>@?8@@>@A###########################################################################################################################################################      NM:i:16 MD:Z:38T63G0T2T0T0C29G36A8A30A1G1A6G6A0C2T12    AS:i:170        XS:i:0
  H7AGFADXX131213:2:1114:9857:73128       133     MT      15987   0       *       =       15987   0       ATGGCCCAAAAGAGAGGAGTGCAGGAAGATTTTGCGGGATAATGAATACCCGAAGAAGGGTGGACAAGGGATCCCTATCTCCGGTGGAACATACATAGGGCCGAGAAAGGACTTAACTGTAATGAGCTATGCATAATAGATAACTGTACAGTTCATCAAATTGTGAGGATGAGTATGATTATTTGTGTCCTCGAGGGAGAGGCGAGGACATGGATAAGCCGCTGAGTTGGGGACGTTGAAGGTTAATTGC      ##########################################################################################################################################################################################################################################################      AS:i:0  XS:i:0
  ```
As you can see, our first line (@SQ) contains information about the reference genome (there would be more lines if there were more sequences in our reference), our second line (@PG) contains information about our bwa commands, and the remaining lines contain information about each mapped read.  You can find more information about what information is contained in each record in the [SAM/BAM specificaions](https://samtools.github.io/hts-specs/SAMv1.pdf).

While it's exciting that we've successfully mapped our first sample, there are two potential problems with our command above.  First, while SAM format is convenient in that it's human-readable, alignment files are ENORMOUS, so we'll want to compress (BAM format) to minimize our data footprint.  Further, we want to ensure that all read-pairing, etc. didn't get lost in the mapping process.  We can easily add these steps to our pipeline by piping our output to ```samtools``` which handles this easily.  So our complete command now looks like this:
  ```
  $ bwa mem -M reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samtools fixmate -O bam - bam/ind1.bam
  ```
#### Processing our bam file: adding read groups, removing duplicates, sorting, and indexing.
The bam file we just created contains all of our alignments, but it's currently unordered, unlabeled, unfiltered, and unindexed.

##### Adding read groups
Read groups are very useful when you're working with multiple samples, sequencing lanes, flowcells, etc.  Importantly for us, it'll help our downstream variant caller label samples correctly and handle potential sequencing batch effects.  [Picard](https://broadinstitute.github.io/picard/command-line-overview.html) is very commonly used to add read groups to bam files, but ```bwa``` also has the ability to add read groups on the fly while mapping.  This latter option will save us time and space, so we'll add read groups to individual 1 with ```bwa``` by adding to our previous command:
  ```
  $ bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samtools fixmate -O bam - bam/ind1.bam
  ```
In a perfect world, we'd know more about our sample and could use these tags more appropriately.  ID is the name of a read group (containing a unique combination of the following tags).  SM is the sample name.  LB is the sequencing library. PU is the flowcell barcode.  PL is the sequencing technology.  A number of other options exist (see the [@RG section of the SAM/BAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information).

*Note that the fastq files we're using today actually contain reads from multiple lanes, etc., as I randomly grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.*

##### Removing Duplicates
During library preparation for sequencing, amplification steps can lead to PCR duplicates of reads.  The inclusion of duplicate reads can negatively affect our downstream analyses, so we need to remove them (this is the case with DNA sequencing - there's debate whether or not to do this in RNA-seq).  Most available tools, such as [Picard](https://broadinstitute.github.io/picard/command-line-overview.html), take a sorted bam file as input and output a sorted bam with duplicates either flagged or removed.  [Samblaster](https://github.com/GregoryFaust/samblaster), on the other hand, can take streaming output from bwa, which speeds up duplicate removal and allows us to produce one less bam file (saving us space). [Sambamba](http://lomereiter.github.io/sambamba/) is another nice option.  Note that at the time of writing this tutorial, [the developers of Samtools suggest against using it to remove duplicates](https://github.com/samtools/samtools/issues/614).  We'll opt for Samblaster here, and incorporate it into our bwa command like so:
  ```
  $ bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate -O bam - bam/ind1.rmdup.bam
  ```
This will flag, but not remove, duplicates in our bam.

A quick note here.  If the same sample library is sequenced across multiple lanes, you'll want to identify and mark/remove duplicates AFTER merging the bam files from all of the lanes.

*Note, again, that the fastq files we're using today actually contain reads from multiple lanes, etc., as I randomly grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.*

##### Sorting bam files
If we weren't piping data directly to Samblaster, and instead using a different tool for duplicate removal, we'd have to sort our bam (in genome coordinate order) first because duplicates are often identified as some combination of reads with identical sequences mapping to the exact same start and end points.  But because we were able to build duplicate removal and adding our read groups into our pipeline, sorting will come last.

Sorting doesn't require too much explanation.  Most genomic datasets are huge, so it's inefficient to move along unsorted bam files (we need access to all reads covering a given base for variant calling, for example).  ```samtools sort ``` is widely used, and that's what we'll employ here.  Like our previous tools, it handles streaming input, so we can simply add to our previous command to save space:
  ```
  $ bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind1.rmdup.sorted.bam -
  ```

So, now, after running our command, we can view the first 10 lines of ```ind1.rmdup.sorted.bam``` with the command:
  ```
  $ samtools view -h bam/ind1.rmdup.sorted.bam | head
  ```

We need ```samtools``` to view the bam file because of how it is compressed.  The ```-h``` flag ensures that the bam header is printed along with the records.

Anyway, that command will give the following output:
  ```
  @HD     VN:1.3  SO:coordinate
  @SQ     SN:MT   LN:16569
  @RG     ID:ind1 SM:ind1 LB:ind1 PU:ind1 PL:Illumina
  @PG     ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -M -R @RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz
  @PG     ID:SAMBLASTER   VN:0.1.23       CL:samblaster -i stdin -o stdout -M
  H7AGFADXX131213:1:1216:6248:58707       353     MT      1       60      216H34M =       1       206     GATCACAGGTCTATCACCCTATTAACCACTCACG      ;@8<B=A>>17@???A=@@B@@A@A>?A?A<95)      NM:i:0  MD:Z:34 AS:i:34 XS:i:0  RG:Z:ind1       SA:Z:MT,16354,+,216M34S,60,2;
  H7AGFADXX131213:1:1115:9742:18153       99      MT      1       60      122S128M        =       1       193     TCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATC      >>5?@>?7=>?>>?>>><?<?>>>>==>>>@?>??>>>?><A>=>>>=8@<?>>@=>?<=>@?=AA?@?>?@@@@AAA??@BAAA@???@>A>A>9??@@@@A??@A@AA??=@?@A>9B?@?A?@@>@?@??@@>@A=??BA?@>?>?@=@@@>:?@B@AA?A@@A@??A@@@>=???B@B@9??@@>>@@<&==>@A>9?:B?@?@=>A>@9B???7=@<<B?@A9>B?AA=A?A@??>?7B@?9968      NM:i:0  MD:Z:128        AS:i:128        XS:i:0  RG:Z:ind1       SA:Z:MT,16448,+,122M128S,60,1;  MQ:i:60
  H7AGFADXX131213:1:1203:19236:9314       99      MT      1       60      93S157M =       1       240     GCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGAATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCCATTAT      >>=?@@=@>??<@6:>>><6@<??=@:6=>==A><@>=><3>>?@?=>@>>@A??@>>?@=?>>@=9?A>==@@?@@@=>=<@=B??@<'@5=;@>?@=@==@=@=?=?;A?B??@?@==?>B?A588<A>A<>@<?@?==><?@=);<>@@B;&?;A>4;3/1&3<;5>;'27B<53:<=?/:4?&?=2=@=>@??A9@??=A???C@?=1=,:&9@,=?@>@>AAA@A@?==@>?AA@>?=>?1:77;      NM:i:2  MD:Z:71T79T5    AS:i:147        XS:i:0  RG:Z:ind1       SA:Z:MT,16477,+,93M157S,60,1;   MQ:i:60
  H7AGFADXX131213:1:2110:2821:77416       163     MT      1       60      89S161M =       1       225     AAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCCATTATTTAT      <=;=><8<A>@>>>=71<?>>>>>>?>>?><??>>>=>>>>>>?>>>>?>?>>>>@?>?<@=7??????A@?@?@?@?@=@??@=9B===A@@A>@>@@@A@@@A>?AB??@A@=>@=@?@=5@<B@AB?@??A?@@@@@@=?@@@B@A@:?@B@??8;;8<?@@A>:@:B?@@@A@@>@:B@A<7AA=@@=B?:@<?AA>A?B@@@@@7>>;???@A>@@AAA?A@@@?19=><@A@@?@A@??@?>==      NM:i:1  MD:Z:151T9      AS:i:156        XS:i:0  RG:Z:ind1       SA:Z:MT,16481,+,89M161S,60,1;   MQ:i:60
  H7AGFADXX131213:1:1211:16620:80940      353     MT      1       60      213H37M =       1       221     GATCACAGGTCTATCACCCTATTAACCACTCACGGGA   @A??@=@??@@A@>@A=??A@?A>@=?A>A@@=251/   NM:i:0  MD:Z:37 AS:i:37 XS:i:0  RG:Z:ind1       SA:Z:MT,16357,+,213M37S,60,1;
  ```
You can see that there are a few changes.  Most notably, there is a @HD line in the header indicating that the file is sorted in coordinate order ("SO:coordinate"), we have an additional @PG flag with information about the SAMBLASTER command, and we now have a @RG line in the header with our read group information and a RG tag on each record with its corresponding read group ID (here they're all ind1).

##### Indexing
Again, bam files can get pretty big and they're compressed, so we need to index them for other tools to use them.  Here's a simple command for indexing our bam (from our main directory):
  ```
  $ samtools index bam/ind1.rmdup.sorted.bam
  ```
This creates a file with a .bai extension that needs to remain in the same directory as its corresponding bam.

##### Summarizing the contents of your bam file
Now that we have our final bam file, we might be interested in knowing what it contains.  More specifically, we should get a sense of how many duplicates there were, how successful mapping was, and other measures along those lines.  Fortunately, [samtools offers tools to calculate summary statistics](http://www.htslib.org/doc/samtools.html).  We'll calculate stats using ```samtools stats``` because it provides a bit more detail, but you should have a look at ```samtools flagstat``` as well.  From our main directory, enter the command:
  ```
  $ samtools stats bam/ind1.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind1.rmdup.sorted.bam.stats
  ```
Because samtools stats offers a huge range of statistics including a number of very big tables in our output, we'll just grab the summary statistics.  This is what ```grep ^SN | cut -f 2-``` in our command does.

We can print the contents of each file to screen using the ```cat``` command.  Here's the command and its result for the samtools stats output for ```stats/ind1.rmdup.sorted.bam.stats```:
  ```
  $ cat stats/ind1.rmdup.sorted.bam.stats

  raw total sequences:	10000
  filtered sequences:	0
  sequences:	10000
  is sorted:	1
  1st fragments:	5000
  last fragments:	5000
  reads mapped:	9882
  reads mapped and paired:	9764	# paired-end technology bit set + both mates mapped
  reads unmapped:	118
  reads properly paired:	9616	# proper-pair bit set
  reads paired:	10000	# paired-end technology bit set
  reads duplicated:	12	# PCR or optical duplicate bit set
  reads MQ0:	0	# mapped and MQ=0
  reads QC failed:	0
  non-primary alignments:	117
  total length:	2490364	# ignores clipping
  bases mapped:	2460864	# ignores clipping
  bases mapped (cigar):	2411834	# more accurate
  bases trimmed:	0
  bases duplicated:	3000
  mismatches:	21935	# from NM fields
  error rate:	9.094738e-03	# mismatches / bases mapped (cigar)
  average length:	249
  maximum length:	250
  average quality:	25.9
  insert size average:	566.3
  insert size standard deviation:	902.9
  inward oriented pairs:	4669
  outward oriented pairs:	80
  pairs with other orientation:	2
  pairs on different chromosomes:	0
  ```
As you can see, it gives us a lot of information about number of reads, mapping, pairing, etc.  A quick glance shows us that our mapping was quite successful (9882 out of 10000 reads mapped).  We also had very few PCR duplicates (12 reads - note that this can be calculated because we flagged, but didn't remove duplicates using samblaster earlier), but this is probably because I randomly sampled read pairs to create our fastq files.

#### Variant calling
Now that our bam is processed and indexed, it's time to call variants!!  There are a few popular variant callers in use (e.g., [GATK's HaplotypeCaller](https://software.broadinstitute.org/gatk/), [Platypus](http://www.well.ox.ac.uk/platypus), and [VarScan](http://dkoboldt.github.io/varscan/)), but today we'll be using [Freebayes](https://github.com/ekg/freebayes) because it works well on both Mac and Linux and can be easily installed via conda.  Freebayes is extremely flexible (it allows a great deal of performance tuning, handles atypical ploidy, can take tumor/normal pairings, etc.), but it's simplest case is perfect for us today.  From our main directory, enter the command:
  ```
  $ freebayes -f reference/human_g1k_v37_MT.fasta bam/ind1.rmdup.sorted.bam > vcf/ind1.raw.vcf
  ```
The resulting list of variants is quite small.  Because the mitochondrial genome is haploid, we really only expect to see variant records for sites where individual 1 differs from the reference. However, because there are so many copies of mitochondria sequenced during normal genome sequencing, we might expect a bit of heteroplasmy (variation among different mitochondrial genomes within an individual) and the possibility of heterozygous calls.

We can print all heterozygous calls for which the site quality (a measure of our confidence that an allele different than the one found in the reference genome) divided by the total read depth is greater than 8 to the screen using ```bcftools``` with the following command (from the main directory):
  ```
  $ bcftools view -g het -i "(QUAL/DP) > 8" vcf/ind1.raw.vcf
  ```

You should see three heterozygous calls passing filters, along with the full VCF header.  While each caller, unfortunately, produces its own version of a VCF, the header usually contains a very detailed description of how to interpret the file.  So, for each VCF file you generate (or receive), you should spend some time reading the header to understand what measures are available for filtering.

That's all we're going to do in terms of filtering vcf files today. In case you're interested in doing more filtering, ```bcftools``` is a very flexible and powerful vcf parser.  You can find out more about what it can do [here](https://samtools.github.io/bcftools/bcftools.html).  ```SnpSift``` is another useful tool for working with vcf files ([check out the manual here](http://snpeff.sourceforge.net/SnpSift.html)). GATK's ```SelectVariants``` ([here]{https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php}) and ```VCFtools``` ([see here](https://vcftools.github.io/man_latest.html)) are other nice options.

### Putting it all together
To run our pipeline on our two samples, we can simply run the following commands:

```
$ bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind1.rmdup.sorted.bam -
```
```
$ samtools index bam/ind1.rmdup.sorted.bam
```
```
$ samtools stats bam/ind1.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind1.rmdup.sorted.bam.stats
```
```
$ bwa mem -M -R '@RG\tID:ind2\tSM:ind2\tLB:ind2\tPU:ind2\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind2_1.fastq.gz fastq/ind2_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind2.rmdup.sorted.bam -
```
```
$ samtools index bam/ind2.rmdup.sorted.bam
```
```
$ samtools stats bam/ind2.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind2.rmdup.sorted.bam.stats
```
And while we ran Freebayes on a single bam file before, it will just as easily take two files for joint calling:

```
$ freebayes -f reference/human_g1k_v37_MT.fasta bam/ind1.rmdup.sorted.bam bam/ind2.rmdup.sorted.bam > vcf/joint.raw.vcf
```

## Making your pipeline reproducible
So, now that we have our pipeline, it's time to make it as reproducible as possible.  There are a few reasons for this, including (but not limited to), allowing us to keep track of and control which versions of programs we're using, making it easy to share with collaborators, ensuring you know exactly what you did down the line, and making your research open so that reviewers and other researchers can reproduce your analyses and adapt them to their own research.

### Sharing your environment
Luckily, by working with Anaconda, we can easily share our environment (including exact versions of tools).  To do this, we can enter the command:
  ```
  $ conda env export > environment.yml
  ```
and we'll get a file that looks something like:
  ```
  name: ASU_Genomics_Tutorial
  channels: !!python/tuple
  - !!python/unicode
    'bioconda'
  - !!python/unicode
    'defaults'
  dependencies:
  - asn1crypto=0.22.0=py36_0
  - bioconda::bcftools=1.4.1=0
  - bioconda::bedtools=2.26.0=0
  - bioconda::bioawk=1.0=1
  - bioconda::bwa=0.7.15=1
  - bioconda::dropbox=5.2.1=py36_0
  - bioconda::fastqc=0.11.5=1
  - bioconda::filechunkio=1.6=py36_0
  - bioconda::freebayes=1.1.0=py36_1
  - bioconda::ftputil=3.2=py36_0
  - bioconda::java-jdk=8.0.92=1
  - bioconda::pysftp=0.2.9=py36_0
  - bioconda::samtools=1.4.1=0
  - bioconda::snakemake=3.13.0=py36_1
  - bioconda::urllib3=1.12=py36_0
  - bzip2=1.0.6=3
  - cffi=1.10.0=py36_0
  - cryptography=1.8.1=py36_0
  - curl=7.52.1=0
  - docutils=0.13.1=py36_0
  - idna=2.5=py36_0
  - mkl=2017.0.1=0
  - numpy=1.13.0=py36_0
  - openssl=1.0.2l=0
  - packaging=16.8=py36_0
  - pandas=0.20.2=np113py36_0
  - paramiko=2.1.2=py36_0
  - pip=9.0.1=py36_1
  - psutil=5.2.2=py36_0
  - pyasn1=0.2.3=py36_0
  - pycparser=2.17=py36_0
  - pyparsing=2.1.4=py36_0
  - python=3.6.1=2
  - python-dateutil=2.6.0=py36_0
  - pytz=2017.2=py36_0
  - pyyaml=3.12=py36_0
  - readline=6.2=2
  - requests=2.14.2=py36_0
  - setuptools=27.2.0=py36_0
  - six=1.10.0=py36_0
  - sqlite=3.13.0=0
  - tk=8.5.18=0
  - wheel=0.29.0=py36_0
  - wrapt=1.10.10=py36_0
  - xz=5.2.2=1
  - yaml=0.1.6=0
  - zlib=1.2.8=3
  prefix: /Users/thw/anaconda/envs/ASU_Genomics_Tutorial

  ```
This is perfect, except for the prefix line at the end (which will get in the way for future users that don't have the same prefix).  So, we'll need to delete that line.  You can do this in a text editor like ```nano``` or ```vi``` (both come standard in most Linux/UNIX environments).  If you haven't used ```vi``` or something similar before, you're probably going to have huge issues writing, saving, and exiting, so I would probably avoid it for now.  ```nano``` can be used by simply typing ```nano <filename>```, using your keyboard arrows to move down to the bottom of the file, manually deleting the last two lines, and following the instructions on the bottom of the screen to exit.

We can also delete the last line of the file easily with ```sed```.  The command for doing so is ```sed '$d' <filename>```. This will print the file, without its last line, to the screen.  If your file is like mine, with one blank line at the end, and the prefix statement on the second to last line, we can delete both and print to a new file with:
  ```
  $ sed '$d' environment.yml | sed '$d' > ASU_Genomics_Tutorial.yml
  ```
Now we can easily share our environment file, ```ASU_Genomics_Tutorial.yml```, with anyone.  They can then create an identical environment with the command:
  ```
  $ conda env create -f ASU_Genomics_Tutorial.yml
  ```
This will create an environment called ```ASU_Genomics_Tutorial``` on their computer (because of the "name" row in the .yml file).

### Sharing your pipeline

#### Bash script
Now that you're able to share your environment, what about your pipeline?  One option is that you can simply share your list of commands, and ask users to run them on their own.  You could make this easy for them by creating a short shell script.  I've included one in the main directory, ```example.sh```.  If we look inside (try ```$ cat example.sh```), we can see that it's very similar to our list of commands:
```
#!/usr/bin/env bash

# Prepare reference
samtools faidx reference/human_g1k_v37_MT.fasta
samtools dict -o human_g1k_v37_MT.dict human_g1k_v37_MT.fasta
bwa index reference/human_g1k_v37_MT.fasta

# Process sample 1 (ind1)
bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind1.rmdup.sorted.bam -
samtools index bam/ind1.rmdup.sorted.bam
samtools stats bam/ind1.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind1.rmdup.sorted.bam.stats

# Process sample 2 (ind2)
bwa mem -M -R '@RG\tID:ind2\tSM:ind2\tLB:ind2\tPU:ind2\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind2_1.fastq.gz fastq/ind2_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind2.rmdup.sorted.bam -
samtools index bam/ind2.rmdup.sorted.bam
samtools stats bam/ind2.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind2.rmdup.sorted.bam.stats

# Jointly call variants for both samples
freebayes -f reference/human_g1k_v37_MT.fasta bam/ind1.rmdup.sorted.bam bam/ind2.rmdup.sorted.bam > vcf/joint.raw.vcf
```
The one difference is the first line, ```#!/usr/bin/env bash```, which is called the shebang (the other lines beginning with ```#``` are just comments).  It allows us to run the script like a program from the command line.  For example, to run ```example.sh```, we could type (from the main directory):
```
$ ./example.sh
```
and you should see that our pipeline is runs exactly as it did earlier, only now using a single command.  Note that the ```./``` is required for your shell to recognize your script as a program.

One quick note: if we want others to be able to use our script (either to run it, or adapt it for their own purposes), we need to ensure that the script gives everyone permission to read, write, or execute it.  You can do this with the ```chmod``` command:
```
$ chmod 777 example.sh
```
This will allow _anyone_ to read, write, or run your script.  If you need to limit some aspects of reading/writing/executing, see [this page](http://ss64.com/bash/chmod.html) for more information on the different codes you can use.

#### Snakemake
So, for the general purpose of sharing our analyses, a bash script works well.  But what if we want to allow users to flexibly add more samples or change the reference genome?  What about if we need to run processes in parallel on a high performance cluster (this isn't a worry for our tiny example, but is often critical when working with real genomic datasets)?

One solution is [Snakemake](https://snakemake.readthedocs.io/en/stable/), a program that allows users to construct workflows in Python.  We won't cover Snakemake today, but I highly recommend checking out the [official tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html), [reading the documentation](https://snakemake.readthedocs.io/en/stable/), and browsing the [Google group](https://groups.google.com/forum/#!forum/snakemake). I've also included an example snakefile that will run all of our commands from today: ```example_snakefile```.
