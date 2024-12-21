# Bioinformatics workshop for beginners

## Medical Ain Shams Research Institute, Egypt, Jan 2025.

Note: Material used to prepare for the workshop was extracted from [HPRC Tutorial](https://github.com/human-pangenomics/hprc-tutorials/tree/GA-workshop)

Prepared by: **Mohammed Farahat**

## Table of Contents

* [Day1: Introduction to Unix/Linux OS](#day1-introduction-to-unix/linux-os)

  - [Introduction to Linux/Unix OS and command line](#introduction-to-linuxunix-os-and-command-line)
  - [Manipulating files, useful commands and tips](#manipulating-files-useful-commands-and-tips)
  - [Permissions, groups and control](#permissions-groups-and-control)
  - [Application of Linux](#application-of-linux)
  - [Environment variable](#environment-variable)
  - [Shell Scripting](#shell-scripting)
  - [Controlling tasks](#controlling-tasks)
  - [SSH into remote machine](#ssh-into-remote-machine)

* [Part 1: Familiarize Ourselves With The Data](#part-1-familiarize-ourselves-with-the-data)
    * [Subset Our Input Data](#subset-our-input-data)
    * [Data Comparison](#Now-lets-compare-the-data)
    * [Cleaning Data For Assembly](#cleaning-data-for-assembly)
        * [PacBio Adapter Trimming](#pacBio-adapter-trimming)
        * [ONT Read Length Filtering](#ont-read-length-filtering)
    * [Phasing Data](#phasing-data)
        * [Trio data: Meryl](#trio-data-meryl)
        * [Trio data: Yak](#trio-data-yak)
* [Part 2: Assembly](#part-2-assembly)
    * [Theoretical Walkthrough Of The Assembly Process](#theoretical-walkthrough-of-the-assembly-process)
    * [Exploring With Test Data](#Exploring-with-test-data)
        * [Run Hifiasm With Test Data](#run-hifiasm-with-test-data)
        * [Introduction To GFA Files](#introduction-to-gfa-files)
        * [View Hifiasm Test Assembly GFA in Bandage](#view-hifiasm-test-assembly-gfa-in-bandage)
    * [Run Verkko With Test Data](#run-verkko-with-test-data)
    * [Running on Full Data Sets](#running-on-full-data-sets)
    * [Building Graph Sense](#building-graph-sense)

# Day1: Introduction to Unix/Linux OS
## Introduction to Linux/Unix OS and command line
### What is Linux?
- UNIX is an Operating System (OS) initially developed in the 1960s.
- There are many different versions of UNIX, that share common similarities.
- The most popular varieties of UNIX are Solaris, Linux, and MacOS.
- UNIX systems have a graphical user interface (GUI) making them easier to use.
### Linux vs Unix
- Linux is a “clone” of the original Unix but doesn’t contain its code.
- Linux is free and open source, the original Unix is not (but some of its derivatives are).
- All command lines work the same on both.
### Why Linux?
- Linux is free and the most popular distributions are
Ubuntu, Fedora/Red Hat, Mandriva, etc.
- Low cost and very stable system.
- Most secure OS.
- Best multi-user and multi tasking OS.
- The world’s fastest super computers run Linux.
- Fast developing OS (many developers)
### Linux Distributions
Different Linux distributions are available. You can explore them here: [DistroWatch](http://distrowatch.com/).

- **Ubuntu**: A distribution that is easy and convenient to use for beginners.
- A simple guide to install Ubuntu on your machine can be found here: [Install Ubuntu Desktop](http://www.ubuntu.com/download/desktop/install-ubuntudesktop).

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/4a4c390b4c7222fd7a78966b371c3ef8ca202187/imgaes/Unix.png?raw=true" width="650"/>
</p>

### The Terminal
- A terminal refers to a wrapper program which runs a shell.
- There are many different Unix shells, the most popular shell for interactive use includes Bash: the default on most Linux installations.

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/imgaes/Terminal.png?raw=true" width="650"/>
</p>
Even though it is a command line interface, the mouse is still handy (scroll, copy, paste, etc.)

###  File-System Under Linux
<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/imgaes/FileSys.png?raw=true" width="650"/>
</p>

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/imgaes/FileSys2.png?raw=true" width="650"/>
</p>

### Home Directory and Working Directory
- When you first log in on a UNIX system, the working directory is your home directory.
- While working, you will be associated with one directory called the working directory or the current directory.
- An abbreviation of the working directory is displayed as part of the prompt on your terminal.
- The command `pwd` gives the absolute path of the working directory.

### What is a Path or a Pathname?

- A path locates a given file in the system hierarchy.
- An **absolute path** in the file system hierarchy for a given file or folder describes the parents all the way up to the root.
- A **relative path** describes the path to the file starting from the **current working directory**.

## Manipulating files, useful commands and tips
Content for this section.

## Permissions, groups and control
Content for this section.

## Application of Linux
Content for this section.

## Environment variable
Content for this section.

## Shell Scripting
Content for this section.

## Controlling tasks
Content for this section.

## SSH into remote machine
Content for this section.

# end of page


**Compute-node**
I recommend to reserve a compute-node to run the `bash` script and keep it alive for all the steps.

```
srun --time=4:00:00 --mem=60G --account=b188-cbio-037-ag --reservation=refgraph --cpus-per-task=32 --pty bash
```

**It will take a while to allocate the reourses**, so leave it and open another terminal window and continute.


---
### Conda environment
For the tools that are not preinstalled on ilifu, I have created a conda env that contains all the required tools and packages you will need. 

Download the yml file from here [refgraph.yml](https://github.com/Mo7ammedFarahat/Long-Read-Assembly/blob/main/refgraph.yml)  
Or use it directly from the workshop project dir `/cbio/projects/037/mohammed/condaEnv/refgraph.yml`

First, ask for interactive job **only to install conda env**, like:
```
sinteractive
```
as you see now your account changed from `username@slurm-login` to `username@compute-001`, that's the shared resources that you can run your light jobs, so let's create the conda on it while waiting the other resources to run the other tasks there.


To create the environment:
```bash
username@compute-001:~$ module load anaconda3
username@compute-001:~$ conda env create -f /cbio/projects/037/mohammed/condaEnv/refgraph.yml
username@compute-001:~$ conda activate refgraph
```


Some tools might be found on `\cbio\soft` and `\cbio\bin\`, **you need to add these paths to the compute-node**
```
export PATH=/cbio/bin:$PATH
```
```
export PATH=/cbio/soft:$PATH
```

Now let's check that the conda env and the paths are working fine:
```bash
username@compute-001:~$ mashmap --version
3.1.3
```
```
username@compute-001:~$ which verkko
username@compute-001:~$ /cbio/bin/verkko
```
As you see now you have the tools from PATHs and the tools installed in the conda env.
---

Now, go back to the previous **window where you asked for resources**, Once you got the resource, What we need now is:
* Activate the `refgraph` conda env in that compute-node (After you create it from the next step)
* * Export the two paths again as in the following step.
```
username@compute-xxx:~$ conda activate refgraph
username@compute-xxx:~$ export PATH=/cbio/bin:$PATH
username@compute-xxx:~$ export PATH=/cbio/soft:$PATH
```
You can also check again `which hapmap2` and `which verkko`.


---
## Dataset

  In this workshop we will be using data from HG002, which is a reference sample from the [Genome In A Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) consortium. The GIAB project releases benchmark data for genomic characterization, and you may have seen their benchmark variant calls and regions out in the wild. As part of their benchmarking material generation, they release datasets for their reference samples. We will be using those in this workshop.
    
    
    HG002 is actually part of a trio of reference samples. Below is the family listing, also known as the Ashkenazim trio:
    
    * HG002: Son
    * HG003: Father
    * HG004: Mother


You will find this dataset here:
`/cbio/projects/037/mohammed/T2T/HG002`

## Part 1: Familiarize Ourselves With The Data
**Create a directoty**
```bash
username@slurm-login:~$ mkdir ~/lra
username@slurm-login:~$ cd ~/lra
username@slurm-login:~$ mkdir day1_data
username@slurm-login:~$ cd day1_data
```
### Subset Our Input Data
In order to get a feel for the data, we only need a small portion of it. Pull the first few thousand reads of the **HiFi reads** and write them to new files.
```bash
zcat /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/deepconsensus/m64011_190830_220126.Q20.fastq.gz \
    | head -n 200000 \
    | pigz > hifi_50k_reads.fq.gz &
```
Next, downsample the **ONT UL reads**, too.

```bash
samtools fastq -@4 \
    /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/ont_ul/03_08_22_R941_HG002_1_Guppy_6.1.2_5mc_cg_prom_sup.bam \
    | head -n 20000 \
    | pigz > ont_ul_5k_reads.fq.gz &

```
### Now let's compare the data  
We are going to use a tool called NanoComp. This tool can take in multiple FASTQs (or BAMs) and will create summary statistics and nice plots that show things like read length and quality scores. NanoComp has nano in the name, and has some ONT-specific functionality, but it can be used with PacBio data just fine.

```bash
NanoComp --fastq \
    hifi_50k_reads.fq.gz \
    ont_ul_5k_reads.fq.gz \
    --names PacBio_HiFi ONT_UL \
    --outdir nanocomp_hifi_vs_ul
```
Once the run is complete (~2 minutes), navigate in your file browser to the folder that NanoComp just created and then click on the `NanoComp-report.html`

### Cleaning Data For Assembly
#### PacBio Adapter Trimming
PacBio's CCS software attempts to identify adapters and remove them. This process is getting better all the time, but some datasets (especially older ones) can have adapters remaining. If this is the case, adapters can find their way into the assemblies. 

Run CutAdapt to check for adapter sequences in the downsampled data that we are currently using. (The results will print to stdout on your terminal screen.)


```bash
module load cutadapt

cutadapt \
    -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" \
    -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" \
    --discard-trimmed \
    -o /dev/null \
    hifi_50k_reads.fq.gz \
    -j 0 \
    --revcomp \
    -e 0.05
```
Notice that we are writing output to `/dev/null`. We are working on a subset of these reads so the runtime is reasonable. There is no need to hold onto the reads that we are filtering on, it is just a subset of the data.

#### ONT Read Length Filtering
Hifiasm is often run with ONT data filtered to be over 50kb in length, so let's filter that data now to see how much of the data remains. 
```bash
module load seqkit/2.6.0
seqkit seq \
    -m 50000 \
    ont_ul_5k_reads.fq.gz \
    | pigz > ont_ul_5k_reads.50kb.fq.gz &
```
Now we can quickly check how many reads are retained.
```bash
NanoComp --fastq \
    ont_ul_5k_reads.fq.gz \
    ont_ul_5k_reads.50kb.fq.gz \
    --names ONT_UL ONT_UL_50kb \
    --outdir nanocomp_ul_vs_ul_50kb
```
### Phasing Data
Now that we've introduced the data that creates the graphs, it's time to talk about data types that can phase them in order to produce fully phased diploid assemblies.

At the moment the easiest and most effective way to phase human assemblies is with trio information. Meaning you sequence a sample, and then you also sequence its parents. You then look at which parts of the genome the sample inherited from one parent and not the other. This is done with kmer databases (DBs). In our case, we will use both Meryl (for Verkko) and yak (for hifiasm) so let's take a moment to learn about kmer DBs.

#### Trio data: Meryl
[Meryl](https://github.com/marbl/meryl) is a kmer counter that dates back to Celera. It creates kmer DBs, but it is also a toolset that you can use for finding kmers and manipulating kmer count sets. Meryl is to kmers what BedTools is to genomic regions.

Today we want to use Meryl in the context of creating databases from PCR-free Illumina readsets. These can be used both during the assembly process and during the post-assembly QC. 


#### Helpful Background

Meryl can be used to create hapmer DBs (haplotype + k-mer), which can be used as input for tools like Verkko and Merqury. Hapmer DBs are constructed from the *k*-mers that a child inherits from one parent and not the other. These *k*-mers are useful for phasing assemblies because if an assembler has two very similar sequences, it can look for maternal-specific *k*-mers and paternal-specific *k*-mers and use those to determine which haplotype to assign to each sequence.

<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/sequencing/meryl_venn.png?raw=true" width="350"/>
</p>

In the Venn diagram above, the maternal hapmer *k*-mers/DB are on the left-hand side (in the purple in red box). The paternal hapmer *k*-mers/DB are on the right-hand side (in the purple in blue box). 

**Wait, what is phasing?**

Phasing is the process of saying two things are on the same haplotype (<i>i.e.</i>, saying two blocks of sequence came from the maternal haplotype, or vice versa)

One way you will hear us talk about phasing in this workshop is in the context of ultra long reads. In this case, we may have two heterozygous regions separated by a homozygous region. When an assembler is walking this graph, if there is no external information about haplotype, then the assembler doesn't have a way of knowing that certain blocks of sequence came from the same sequence. For example, in the bottom image, the assembler might walk from the top left block, into the homozygous block, and then down to the <i>bottom right</i> block, switching between the two haplotypes.

However, if we can find a long read that maps to the top sequences in both, then we could say that these sequences come from the same haplotype. That is phasing.

<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/sequencing/ont_phasing.png?raw=true" width="550"/>
</p>

We can play a similar trick with parental data. If we find paternal markers in both sequences in the top, then we can say that they both come from the paternal haplotype. This is also phasing. 

<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/sequencing/trio_phasing2.png?raw=true" width="550"/>
</p>


#### Using Meryl

**Make sure you are in the right directory**

```bash
cd day1_data
```

**Now create a small file to work with**

```bash
zcat /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/ilmn/pat/HG003_HiSeq30x_subsampled_R1.fastq.gz \
    | head -n 20000000 \
    | pigz > HG003_HiSeq30x_5M_reads_R1.fastq.gz &
```    

**Create a *k*-mer DB from an Illumina read set**


```bash
module load merqury/1.3
module load meryl/1.4.1

meryl count \
    compress \
    k=30 \
    threads=4 \
    memory=8 \
    HG003_HiSeq30x_5M_reads_R1.fastq.gz \
    output paternal_5M_compress.k30.meryl
```
This should be pretty fast because we are just using a small amount of data to get a feel for the program. The output of Meryl is a folder that contains 64 index files and 64 data files. If you try and look at the data files you'll see that they aren't human readable. In order to look at the actual k-mers, you have to use meryl to print them.

**Look at the *k*-mers**


```bash
module load meryl/1.4.1

meryl print \
    greater-than 1 \
    paternal_5M_compress.k30.meryl \
    | head

```
The first column is the *k*-mer and the second column is the count of that *k*-mer in the dataset.

**Take a look at some statistics for the DB**


```bash
module load meryl/1.4.1

meryl statistics \
    paternal_5M_compress.k30.meryl \
    | head -n 20

```

We see a lot of *k*-mers missing and the histogram (frequency column) has a ton of counts at 1. This makes sense for a heavily downsampled dataset. Great. We just got a feel for how to use Meryl in general on subset data. Now let's actually take a look at how to create Meryl DBs for Verkko assemblies.
<details>
<summary>clipboard-question "How would we run Meryl for Verkko?"</summary>

    **Here is what the Slurm script would look like:**
    
    (Don't run this, it is slow! We have made these for you already.)

    
```bash
#!/bin/bash 

#SBATCH --job-name='meryl_run'
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=96G
#SBATCH --partition=Main
#SBATCH --output=testjob-%j-stdout.log
#SBATCH --error=testjob-%j-stderr.log


module load merqury/1.3
module load meryl/1.4.1
export MERQURY=$(dirname $(which merqury.sh))

## Create mat/pat/child DBs
meryl count compress k=30 \
    threads=32 memory=96 \
    maternal.*fastq.gz \
    output maternal_compress.k30.meryl

meryl count compress k=30 \
    threads=32 memory=96 \
    paternal.*fastq.gz \
    output paternal_compress.k30.meryl

meryl count compress k=30 \
    threads=32 memory=96    \
    child.*fastq.gz output    \
    child_compress.k30.meryl

## Create the hapmer DBs
$MERQURY/trio/hapmers.sh \
maternal_compress.k30.meryl \
paternal_compress.k30.meryl \
    child_compress.k30.meryl
```
</details>

#### Trio data: Yak

Yak (Yet-Another Kmer Analyzer) is the kmer counter that we need for Hifiasm assemblies and to QC assemblies made with either assembler so let's learn about how to make yak dbs. 

**In the Meryl section we subset R1, now subset R2 as well**


```bash
zcat /cbio/projects/037/mohammed/T2T/HG002/LRA/resources/ilmn/pat/HG003_HiSeq30x_subsampled_R2.fastq.gz \
    | head -n 20000000 \
    | pigz > HG003_HiSeq30x_5M_reads_R2.fastq.gz &
```  

**Look up yak's github and figure out how to make a count/kmer db for this data**

Yak won't work on our Jupyter instances, so create a slurm script that has 32 cores and 96GB of memory. That way it will work on our subset data and it will also work on full size data -- you'd just have to extend the time variable in slurm.

<details>
<summary>clipboard-question "Click here for the answer"</summary>

    Here is one way to call yak in a `yak.sbatch` script...

```bash
#!/bin/bash


#SBATCH --job-name='yak_run'
#SBATCH --cpus-per-task=32
#SBATCH --time=00:10:00
#SBATCH --mem=100G
#SBATCH --partition=Main
#SBATCH --output=testjob-%j-stdout.log
#SBATCH --error=testjob-%j-stderr.log

module load anaconda3
conda activate refgraph

#Export Paths
export PATH=/cbio/bin:$PATH

export PATH=/cbio/soft:$PATH





yak count \
    -t32 \
    -b37 \
    -o HG003_subset.yak \
    <(zcat HG003_HiSeq30x_5M_reads_R*.fastq.gz) \
    <(zcat HG003_HiSeq30x_5M_reads_R*.fastq.gz)

``` 

    Notice that for paired-end reads we have to stream both reads to yak twice!
</details>

When you are done you get out a non-human readable file. It doesn't need to be compressed or decompressed, and nothing else needs to be done in order to use it.

**Closing remarks on yak**

* If you have Illumina data for an entire trio (which we do), then you can use yak to make yak DBs for each parent separately to use hifiasm for trio assembly or yak trioeval for quality control (more on that later)
* You don't need to homopolymer compress yak dbs
* There is no need to create separate dbs for assembly and for QC
* yak can perform a variety of assembly QC tasks (as we will see) but it isn't really designed to play around with kmers like Meryl is

### Hi-C
Hi-C is a proximity ligation method. It takes intact chromatin and locks it in place, cuts up the DNA, ligates strands that are nearby and then makes libraries from them. It's easiest to just take a look at a cartoon of the process.

<p align="center"> <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/sequencing/hi-c-flow-2.png?raw=true" width="550"/> </p>

Given that Hi-C ligates molecules that are proximate (nearby) to each other, it can be used for spatial genomics applications. In assembly, we take advantage of the fact that most nearby molecules are on the same strand (or haplotype) of DNA. 

<details>
<summary>What are the advantages of trio phasing over Hi-C?</summary>

    Trio data is great for phasing because you can assign haplotypes to maternal and paternal bins. This has the added benefit of assigning all maternal contigs to the same assembly. Hi-C ensure that an entire chromosome is phased into one haplotype, but across chromosomes the assignment is random. 
</details>

<details>
<summary>So why wouldn't you always use trio data for phasing?</summary>

    It can be hard to get trio data. If a sample has already been collected it may be hard to go back and identify the parents and collect sample from them. In non-human samples, trios can also be difficult particularly with samples taken from the wild. 
</details>

<details>
<summary>Are there any difficulties in preparing Hi-C data?</summary>

    Yes! As you can see in the cartoon above Hi-C relies on having intact chromatin as an input, which means it needs whole, non-lysed cells. This means that cell lines are an excellent input source, but frozen blood is less good, for instance.
</details>

### Other (Phasing) Datatypes
We should also mention that there are other datatypes that can be used for phasing, though they are less common.

**Pore-C**  
Pore-C is a variant of Hi-C which retains the chromatin conformation capture aspect, but the sequencing is done on ONT. This allows long-read sequencing of concatemers. Where Hi-C typically has at most one "contact" per read, Pore-C can have many contacts per read. The libraries also do not need to be amplified, so Pore-C reads can carry base modification calls. 

**StrandSeq**  
StrandSeq is a technique that creates sparse Illumina datasets that are both cell- and strand-specific. Cell specificity is achieved by putting one cell per well into 384-well plates (often multiple). Strand specificity is achieved through selective fragmentation of nascent strands. (During DNA replication, BrdU is incorporated exclusively into nascent DNA strands. In the library preparation the BrdU strand is fragmented and only the other strand amplifies.) This strand specificity gives another way to identify haplotype-specific kmers and use them during assembly phasing.

*If you are interested in these phasing approaches, you can read more about them in the following articles:*

    - Lorig-Roach, Ryan, et al. "Phased nanopore assembly with Shasta and modular graph phasing with GFAse." *bioRxiv* (2023): 2023-02.
    
    - Porubsky, David, et al. "Fully phased human genome assembly without parental data using single-cell strand sequencing and long reads." *Nature biotechnology* 39.3 (2021): 302-308.

## Part 2: Assembly
 ***Here is a rundown of what we will do today:***

    * Learn about how Verkko creates assemblies
    * Run Hifiasm with test data (HiFi only)
    * Run Verkko with test data (HiFi + ONT)
    * See how to run both with full datasets
    * Compare the two assemblers
    * Talk about how much data we need 
    
At the end of the day you will hopefully have a feel for how to actually run each assembler, what data to give them, and when to choose one over the other.

## Theoretical Walkthrough Of The Assembly Process

### Verkko
In this section we will go over the rough outline of Verkko's approach to assembly. Hopefully this will help put the data types from yesterday in context. Knowing how each data type is used also helps you to make better decisions when planning your sequencing runs.

Both Verkko and Hifiasm can use a variety of data sources:



* PacBio HiFi: >10kbp, around 99.9% accuracy
* Oxford Nanopore Ultralong: >100kb, around 97% accuracy
* Phasing data from Hi-C or trio Illumina data

PacBio HiFi data is required for both assemblers. Other data types are optional&mdash;but they lead to much better assemblies. So let's jump ahead a bit and take a peek at how Verkko creates assemblies using figure 1 from the recent Verkko paper (Rautiainen, Mikko, et al.). It's ok if this is a bit confusing, you don't need to know the inner workings of Verkko in order to make great assemblies.

**PacBio HiFi is used to create the initial graph** 

Verkko's first theoretical task is to create an assembly graph from HiFi data, but it has to prepare the HiFi reads first. HiFi data is less accurate in homopolymer repeats and microsatellites, so before creating an assembly graph, the reads are "compressed" in these regions:
<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/verkko_process_hifi_repeat_compr.png?raw=true" width="550"/>
</p>

The reads are then error corrected. You don't have to worry about how this works (just know that it is very computationally intensive). Once that is done, a graph is created from the HiFi reads:

<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/verkko_process_hifi_graph.png?raw=true" width="550"/>
</p>

If you aren't familiar with what an assembly graph is, that is also ok! The annoying thing is that there are different ways to make assembly graphs, but they all have the common feature of linking together reads by their overlaps. In this graph the boxes (also called nodes) represent sequences and the lines (also called edges) represent the relationship between those overlaps.

Note that in the middle section the orange and blue lines represent the ONT reads aligned to the graph. 

**Oxford Nanopore Data (which is long) helps simplify the graph**

Using these alignments, nodes that are linked (or phased) by a read are combined. This "simplifies" the graph -- in other words, you get nice long nodes where you previously had shorter nodes. (Long nodes are good, they mean you have longer sequences that are assembled.)

<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/verkko_process_ont_resolved.png?raw=true" width="550"/>
</p>

**The graph can now be phased**

In the case of trio Illumina data, Verkko looks at the nodes and counts the number of maternal-specific or paternal-specific sequences of DNA (from meryl hapmer DBs). If it finds that a node has, for instance, a bunch of maternal specific sequences/*k*-mers and almost no paternal specific sequences/*k*-mers then the assembler will assign that node to be maternal. Nodes that are the same haplotype but are separated by a homozygous region are then merged.

<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/verkko_process_phasing.png?raw=true" width="550"/>
</p>

Finally the assembly graph can be converted into two contigs which represent maternal and paternal haplotypes.

<p align="center">
<img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/verkko_process_contigs.png?raw=true" width="550"/>
</p>

Maternal and paternal contigs for the entire assembly are then put into one diploid FASTA as well as two haploid FASTAs.

!!! info "How does Hi-C phasing work?"

Like using trio data, Hi-C phasing aims to find nodes that are near to each other and come from the same haplotype. To achieve this, Hi-C data is aligned to the graph and reads that are linked across nodes can be used to phase the graph as shown in this figure modified from (Garg, Shilpa):
<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/hic_phasing.png?raw=true" width="450"/>
</p>
One key difference between trio phasing and Hi-C is that Hi-C data cannot say that a set of nodes all come from the sample's mother or father, only that they come from the same haplotype. 



## Exploring With Test Data
Running assemblers is very computationally intensive and the output files can be big. Let's not jump straight into assembling human genomes. Instead we can use the test data that both assemblers provide as a way to both ensure that we know how to run the tool (which is easy) and we can start to get a feel for the process and outputs in real life. 

#### Run Hifiasm With Test Data

**Create a directory**


    
```bash
cd ~/lra
mkdir -p day2_assembly/hifiasm_test
cd day2_assembly/hifiasm_test
```

**Now download Hifiasm's test data**


```bash
wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz
```
This is HiFi data from about 2 million bases of chromosome 11. HiFi data is the only required data type for Hifiasm and Verkko. You can create assemblies from only HiFi data and you can add ONT and phasing later. Also notice that this data is in FASTA format! Presumably this is to make the file smaller since this is test data.


**And actually run the test data**


```bash
hifiasm \
    -o test \
    -t4 \
    -f0 \
    chr11-2M.fa.gz \
    2> test.log &
```
This should take around 3 minutes. Once the run is complete take a look at the top of the log:


```bash
head -n 60 test.log
```

Now check the [Hifiasm log interpretation](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#hifiasm-log-interpretation) section of the documentation to give that some context.

Now `ls` the directory to see what outputs are present. What do you see? No FASTA files, but there are a lot of files that  end in `gfa`. If you haven't seen these before, then we get to introduce you to another file format! 

#### Introduction To GFA Files
GFA stands for [Graphical Fragment Alignment](http://gfa-spec.github.io/GFA-spec/GFA1.html) and Hifiasm outputs assemblies in GFA files. GFAs aren't like bed or sam files which have one entry per line (or FASTA/Q that have 2/4 lines per entry). But this is bioinformatics, so you can rest assured that it is just a text file with a new file extension. It's easiest to just look at an example of a GFA file from the spec:

```bash
H    VN:Z:1.0
S   11  ACCTT
S   12  TCAAGG
S   13  CTTGATT
L   11  +   12  -   4M
L   12  -   13  +   5M
L   11  +   13  +   3M
P   14  11+,12-,13+ 4M,5M
```

> **Info:** "Here we see the following line types"

    * H (Header): File header. You get the idea.
        * The example here the header is just saying that the file follows GFA 1.0 spec. 
        * Notice that this line follows a TAG:TYPE:VALUE convention. Type in this case is Z which corresponds to printable string. 
    * S (Segment): A sequence of DNA
        * This is what we care about for the moment!
    * L (Link): Overlap between two segments
        * We can read the first Link line as saying that the end of Segment 11 (+) connects to the beginning of Segment 12 (-) and the overlap is 4 matching bases. In this case it would look like this:
        
    ```bash
        ACCTT    (Segment 11)
         ||||
         GGAACT  (Segment 12 -- reversed)
    ```

    * P (Path): Ordered list of segments (connected by links)
**So how do we get a FASTA from a GFA?**  
To get a FASTA we just pull the S lines from a GFA and print them to a file:

```bash
awk '/^S/{print ">"$2;print $3}' \
    test.bp.p_ctg.gfa \
    > test.p_ctg.fa 
```

> **Info:** "You can read this awk command as:"

    1. Give me all input lines that start with `S` 
    2. Then print the second column of those lines (which is the sequence ID) 
    3. Also print another line with the actual sequence

<details>
<summary>"Why does Hifiasm output GFAs and not FASTAs?"</summary>
  
    Hifiasm (and many other assemblers) use GFAs while they are actually assembling. The GFA represents/stores the assembly graph. Hifiasm probably doesn't output FASTAs just because everything in the FASTA is contained in the GFA, so why store it twice?
</details>

#### View Hifiasm Test Assembly GFA in Bandage
**This is OPTIONAL (in your own time)**  
  
We are going to take a look at the assembly GFA file in a browser called Bandage. Bandage provides a way to visualize something called unitig graphs.

!!! jupyter "Start Bandage" OR copy the files to your local machine using `scp` or mount it using `sshfs`
Example (from your local machine):
```
cd /path/to/local/dir
scp  username@transfer.ilifu.ac.za:/path/to/gfa .
```
OR
```
cd /path/to/local/dir
sshfs username@transfer.ilifu.ac.za:/path/to/your/files/on/ilifu .
```
    1. Open Jupyter Virtual Desktop [according to these instructions](../supplementary/supplementary_3.md)
    2. In the Virtual Desktop, click on the terminal emulator icon (in your toolbar at the bottom of your screen)
    3. Load the Bandage module with `module load Bandage`
    4. Type `Bandage &` to start Bandage

<p align="center">
    <img src="https://github.com/Mo7ammedFarahat/Long-Read-Assembly/blob/main/bandage_vdt.png" width="500"/>
</p>  
    **Load a unitig GFA**
    
    1. Click the *File* dropdown then *Load Graph*
    2. Navigate to our current folder (`day2_assembly/hifiasm_test`)
    3. Select the `test.bp.r_utg.noseq.gfa` file and press the **Open** icon
    4. Under **Graph Drawing** on the left-hand side click **Draw Graph**
    
    Ok, so what are we looking at? The thick lines are nodes&mdash;which in this case represent sequences. Since we loaded the unitig graph the sequences are unitigs. A unitig is a high confidence contig. It is a place where the assembler says "I know exactly what is going on here". The ends of unitigs are where it gets messy. At the ends, an assembler has choices to make about which unitig(s) to connect to next.
    
    **Now load a contig GFA**  
    Open the `test.bp.p_ctg.noseq.gfa` file to see how boring it is.
    
    In general, when using Bandage people look at the unitig GFAs (not contig GFAs). An assembly is a hypothesis, and the contigs output by the assembler are its best guess at the correct haplotype sequence. The contigs don't show much information about the decisions being made, however. They are the output. We view unitig GFAs so we can see the data structure at the point that the assembler was making tough decisions. 
    
    **Here are some things you can do with Bandage**
    
    1. Let's say you mapped a sample's ONT reads back onto that sample's *de novo* assembly and have identified a misjoin. You can open up bandage and find that  unitigs that went into the contig to see if it can be easily manually broken.
    2. If you have a phased diploid assembly with a large sequence that is missing, you can look at the unitig gfa, color the nodes by haplotype, and see which sequences are omitted. Those sequences can then be analyzed and manually added into the final assembly.
    3. You can label nodes with (HiFi) coverage and inspect regions with low quality too see if they have low coverage as well. If so, you might want to throw them out. (This does happen, in particular for small contigs that assemblers tend to output.)


### Run Verkko With Test Data
**Create a directory**

```bash
cd ~/lra
mkdir -p day2_assembly/verkko_test
cd day2_assembly/verkko_test
```

**Now download Verkko's test data**  


```bash
curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o hifi.fastq.gz
curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o ont.fastq.gz
```
You can see that this dataset is for *E. coli* and there is both HiFi and ONT data included.


**Create Slurm script for test Verkko run**

    
Start your favourite text editor
```bash
nano verkko_test.sbatch
```

And then paste in the following **(DO NOT FORGOT TO CHANGE CONDA ENV PATH TO YOURS!)**

```bash
#!/bin/bash


#SBATCH --job-name='test_verkko'
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --mem=24G
#SBATCH --partition=Main
#SBATCH --output=testjob-%j-stdout.log
#SBATCH --error=testjob-%j-stderr.log




module load anaconda3
conda activate refgraph

#Export Paths
export PATH=/cbio/bin:$PATH

export PATH=/cbio/soft:$PATH


## run verkko
verkko \
    -d assembly \
    --hifi ./hifi.fastq.gz \
    --nano ./ont.fastq.gz \

echo "Done."
date

```

**Run verkko test**
```bash
sbatch verkko_test.sbatch
```
This should only take a few minutes to complete.

You can keep track of the run with the `squeue` command.


```bash
squeue -u $USER
```

**How does Verkko run?**

It turns out that if you run Verkko more than once or twice you will have to know a bit about how it is constructed. Verkko is a program that reads in the parameters you gave it and figures out a few things about your verkko installation and then creates a configuration file (`verkko.yml`) and a shell script (`snakemake.sh`). The shell script is then automatically executed.

Take a look at the shell script that was created for your run


```bash
cat assembly/snakemake.sh
```

It is just a call to Snakemake!!! You can think of Verkko as a tool, but also as a pipeline because it is. This has some advantages. One is that if you know what Verkko is doing (which is somewhat achievable given that the Snakemake rules guide you through Verkko's logic), you can add to it, or even swap out how Verkko performs a given step for how you'd like to do it.  
It also means that you can restart a run at any given step (if you made a mistake or if the run failed). Lastly, and maybe most importantly, Snakemake supports Slurm as a backend. So if you have access to an HPC you could (and probably should) run Verkko and allow it to launch Slurm jobs for you. (This is in contrast to what we just did which was to run a Slurm job and just allow all jobs to run on the allocated resources that we requested for the entire run.)

**Now take a look at the jobs that were run**

You can view the stderr from the run in your Slurm logs, or in Snakemake's logs. Let's take a look at the top of the log:



```bash
head -n 35 assembly/.snakemake/log/*.log
```
This shows a list of Snakemake jobs that will get executed for this dataset. There are a few things to note. The first is that for larger datasets some jobs will get executed many times (hence the count column). This dataset is small, so most jobs have `count=1`.  
The second thing to note is that these jobs are sorted alphabetically, so we can get a feel for scale, but it's a bit hard to figure out what Verkko is really doing.

Open the logs and scroll through them


```bash
less assembly/.snakemake/log/*.log
```
You can see all of the Snakemake jobs, in order, that were run. Even for this tiny dataset there are many. Since there are a lot of jobs, there are a lot of outputs, and these are organized (roughly) by Snakemake rule. Take a look at the output folder in order to familiarize yourself with the layout.


```bash
ls -lh assembly
```

**Take a look at the initial HiFi graph**

Open the `assembly/1-buildGraph/hifi-resolved.gfa` file in Bandage. You will see that it is already pretty good. There are only three nodes.

**Now take a look at the ONT resolved graph**

Open the `assembly/5-untip/unitig-normal-connected-tip.gfa` file in Bandage. Now our three nodes have been resolved into one. 

<details>
<summary>question "Do we need to phase?"</summary>
     
We just worked through test data that was HiFi only (Hifiasm) then HiFi + ONT UL (Verkko) so you'd be forgiven for thinking that for these assemblers you don't need phasing data. In truth, if the HiFi and ONT data could produce contigs that are chromosome-level in scale, then trios/Hi-C wouldn't be needed, but assemblers aren't there yet. 

Let's look at some example data from Hifiasm which was produced with only HiFi for the Eastern Narrow Mouth Toad (*Gastrophryne carolinensis*). Below is a Merqury kmer count plot for the primary haplotype from this assembly:

**Hifiasm (primary) assembly without Purge Dups**

<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/merqury_kmer_primary_only.png?raw=true" width="300"/>
</p>
This is the primary assembly of a diploid, so we'd like to see the red section have two peaks -- one for heterozygous positions and one for homozygous positions. This is because the primary would ideally represent all of the homozygous regions as well as one set of alleles for heterozygous regions. We'd also like to see almost no blue section (since the blue section represents regions that are represented twice in the assembly), as these mean that a region is being represented twice. Unfortunately there are a lot of kmers seen twice in the assembly. When purge dups (a tool that looks for false duplications based on coverage information) is run, it gets better:

**Hifiasm assembly with Purge Dups**

<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/merqury_kmer_primary_only_purged.png?raw=true" width="300"/>
</p>
But you can see that even after purging duplications, there is a large blue area. Now look what happens when Hi-C phasing is used:

**Hifiasm assembly (one haplotype) with HiC phasing**
<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/merqury_kmer_hap1.png?raw=true" width="300"/>
</p>  

The blue section reduces in size drastically! So what is happening here?

Without phasing data such as trio or Hi-C, Hifiasm will create an assembly graph as you'd expect. Then it will try and figure out what path to walk in order to create a primary and an alternate assembly. It does this by just picking one side of the bubbles it encounters when walking the graph and assigning those to the primary assembly. The other sides gets put in the alternate assembly. In very heterozygous genomes, or very heterozygous regions of a genome, the primary assembly can still retain haplotigs from the alternate allele because the regions may look very different (and therefore don't create a nice bubble in the graph), even though they are still just representing alternate alleles at the same locus. This can result in the same genomic region incorrectly being represented twice in the primary assembly, which we call a false duplication.

Verkko for its part, is more conservative. Once Verkko gets to a bubble that it doesn't know how to phase it just stops and you get a break in your assembly.
</details>  

### Building Graph Sense

We&rsquo;re going to open Bandage and look around. This is a demonstration, so
you do not need to try to do this on your own machine. The goal is to develop
an intution for what the graph visualization is telling you.

* Overall feel for graphs and bandage
    * graphs: nodes, edges, and paths
        * simple examples: cycle & hairpin
    * bandage UI
    * loading CSVs and launching from CLI
* How does a typical human assembly look?
    * ~23 chunks
    * zippers: het-hom-het-etc.
    * non-zippers (rDNA tangles, sex chrs, breaks)
* Complex examples
    * hom with short hets w/o markers
    * long hom with long hets
    * &ldquo;scaffold&rdquo; from graph structure
    * complete fail

## Comparison of Computational Cost

### Hifiasm

Hifiasm is compiled into a single binary file, and, when executed, it manages all tasks and parallelism under one parent process. You can run it the same on a VM in the cloud or in an HPC. 

For a human sample with around 40X HiFi and 30X UL and either HiC or trio phasing Hifiasm can assemble with:

* 64 cores
* 240GB of memory (most samples will use less)
* Around 24 hours of total runtime

So Hifiasm takes about 1500 CPU hours to assemble this sample. On a cluster you can just execute the run command. If you are on a cloud and would like to take advantage of pre-emptible instances, you can break the run command into three parts (each take around 8 hours).


### Verkko

Verkko is written as a shell wrapper around a Snakemake pipeline. This has the advantages of easily restarting after failures and increased potential for parallelism in an HPC environment with multiple nodes available, but it is hard to profile all the individual tasks. If the cluster is not too busy a human assembly can finish in around a day. Most of the compute is done in the overlap and graph aligner jobs. So we can break the runtimes into steps that revolve around the big jobs. That looks something like this:

| **Step**            | **CPUs** | **Shards** | **Time/Shard (est)** | **Total CPU Hours** |
|---------------------|----------|------------|----------------------|---------------------|
| pre overlap         | 24       | 1          | 3                    | 72                  |
| overlap             | 8        | 600        | 1                    | 2400                |
| create graph        | 80       | 1          | 13                   | 1040                |
| graph aligner       | 12       | 100        | 2                    | 2400                |
| complete asm        | 64       | 1          | 12                   | 768                 |


This gives an estimate of around 9000 CPU hours for the same data as above. This is almost certainly an overestimate, but not by more than a factor of 2. 

Note that the runtime estimates for Hifiasm and Verkko don't consider the preparatory work of counting parental *k*-mers with yak or Meryl, which are necessary steps before running either in trio mode.

## Comparison of Outputs

Verkko and Hifiasm are both excellent assemblers. If you have a human sample with over 40X HiFi and over 15X ONT data over 100kb then the high level metrics that you will learn about tomorrow should be pretty comparable across the two assemblers. If you have a set of data that you spent a bunch of money on, and you are hoping to make a high-quality assembly, your best bet is to assemble with both and see which assembly is better. You could even stitch the good parts of each assembly together&mdash;though the people who have to do the actual work tend to flinch when they hear that.


### Verrko's Approach
As we saw previously, Verkko uses HiFi data to create a graph (in the case of Verkko it is a DeBruijn graph). ONT reads are aligned to the graph and the graph is then simplified. One thing that Verkko does is it outputs scaffolds&mdash;where Hifiasm only outputs contigs. With 40X+ HiFi, Verkko's scaffolding tends to add about 12 extra T2T chromosomes to a diploid human assembly. 

The way it is scaffolded comes from the graph (as shown below). One the left we see one haplotype and there is a tangle in the middle of the sequence. Verkko doesn't necessarily know how to walk through this tangle and it doesn't want to output incorrect sequence. So it just estimates the size of the nodes in the tangle and puts the corresponding number of N's into the final assembly. Similarly, on the right we have a gap in one haplotype. Verkko will infer the size of the missing sequence from the other haplotype and put that many N's into to top sequence.

<p align="center">
    <img src="https://github.com/human-pangenomics/hprc-tutorials/blob/GA-workshop/assembly/genomics_aotearoa/images/assembly/verkko_grapholds.png?raw=true" width="350"/>
</p>

This has led some people (well at least one person) to call this approach grapholding. 

### Hifiasm's Approach

Hifiasm creates string graphs from HiFi and ONT data separately (kind of) and then combines them. The argument here is that by creating a standalone ONT graph you don't risk losing information that may be missing in the HiFi-only graph. At the time moment (July 2023) Hifiasm does not include a scaffolding step. Though that will likely change in the coming months.

### How Should I Choose?

**It's not an easy choice, but here are some guidelines**

* If you can, use both
* If you have Hifi coverage under 40X: use Hifiasm
    * Verkko tends to perform less well at lower HiFi coverages
* If you have to pay for compute time: use Hifiasm (see the previous section)
    * Verkko is more expensive to run. If you are on an HPC that may be ok. If you are paying Amazon for your compute then Verkko assemblies can cost upwards of $300 (USD).
* If you want to assemble then fiddle with it to perfect the assembly: use both, then fix things with Verkko
    * Verkko allows you to see its inner workings. You can also make manual changes and then restart from that point in the assembly process. If you do things right, Verkko will take care of the rest. This was done, for instance, by the Verkko team on their version of the HG002 assembly: they manually resolved tangles in the graph.

## How Much Input Data Do I Need?

Let's do some back of the envelope math to see how much is an ideal amount of data that would go into an assembly...

**PacBio HiFi**

Computing overlaps isn't so different from calling variants. For each haplotype we probably want around 10X coverage in order to calculate good overlaps. So that would give about 20X total. But HiFi coverage is variable, and there are some well know regions (such as GA repeats) that drop out of HiFi data. In order to get as many regions as possible above the threshold for assembly we increase the value to 40X.

**ONT UL**

The answer to this depends on who you ask. In the case of Verkko you are often just looking for one or a handful of reads to span a tricky region in the graph. The Verkko team has shown that coverages over 15X of 100kb+ reads don't add much in terms of contiguity. (Though for advanced applications such as using ONT to fill-in missing parts of the graph the math may be different.)

**Trio and Hi-C**

In general the better your graph, the easier it is to phase. If you have only a few big bubbles in the graph, it is a lot easier to find Illumina data that will map to them in a haplotype-specific way. This hasn't been tested rigorously, but people tend to talk about 30X for these datasets.

    
    
