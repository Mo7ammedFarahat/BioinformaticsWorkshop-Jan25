* [Day3: Genomics Sequencing technologies and NGS Overview](Day3.md)
  
  - [Part 1: Introduction to DNA Sequencing](#part-1-introduction-to-dna-sequencing)
  - [Part 2: DNA Sequencing in the NGS Era](#part-2-dna-sequencing-in-the-ngs-era)
  - [Part 3: Overview of NGS Technologies](#part-3-overview-of-ngs-technologies)
  - [Part 4: DNA-Seq Protocol: Overview](#part-4-dna-seq-protocol-overview)
  - [Part 5: DNA-Seq Analysis Pipeline and File Formats](#part-5-dna-seq-analysis-pipeline-and-file-formats)

# Part 1: Introduction to DNA Sequencing

## What is DNA Sequencing?

DNA Sequencing is the process of reading the nucleotides present in DNA: determining the precise order of nucleotides within a DNA molecule.

DNA-Seq generally refers to any NGS method or technology that is used to determine the order of the four bases (A, T, C, G) in a strand of DNA.

## DNA Sequencing Technologies Overview

In fact, there are three main types of DNA sequencing technologies that are widely used today: **Sanger sequencing**, **Next-Generation Sequencing (NGS)**, and **Third-Generation Sequencing (3GS)**.

### Sanger Sequencing
Sanger sequencing, also known as **chain-termination sequencing**, is the traditional method of DNA sequencing. It is highly accurate for sequencing smaller DNA fragments, typically up to 900 base pairs. While it has been a gold standard for many years, its high cost and lower throughput make it less suitable for large-scale applications.

### Next-Generation Sequencing (NGS)
NGS encompasses a range of modern sequencing technologies that allow for massively parallel sequencing, enabling high-throughput, cost-effective, and rapid sequencing of entire genomes or specific regions of interest. NGS platforms like **Illumina sequencing**, **PacBio sequencing**, and **Oxford Nanopore sequencing** offer increased read lengths, higher accuracy, and the ability to sequence millions to billions of DNA fragments simultaneously.

NGS has enabled a wide range of applications, including whole-genome sequencing, RNA sequencing (RNA-Seq), and metagenomics. It has revolutionized genomics research, clinical diagnostics, and personalized medicine.

### Third-Generation Sequencing (3GS)
Third-generation sequencing refers to newer technologies that provide even longer reads and faster sequencing compared to second-generation NGS technologies. Key third-generation sequencing technologies include:

- **Pacific Biosciences (PacBio)**: Known for **Single Molecule, Real-Time (SMRT) sequencing**, PacBio provides very long reads (up to 100,000 base pairs) and high accuracy, especially for detecting structural variations, full-length transcripts, and de novo genome assembly.
- **Oxford Nanopore Technologies (ONT)**: ONT uses nanopore sequencing technology, allowing for ultra-long reads (up to millions of base pairs). It also enables real-time sequencing, providing rapid results that are ideal for field-based genomics or time-sensitive applications.

3GS technologies offer unique advantages in applications requiring long reads, such as assembling complex genomes, detecting structural variants, and characterizing repetitive regions that are difficult to sequence with short reads.

Third-generation sequencing is pushing the boundaries of genomics by enabling comprehensive and highly accurate genome assemblies, transcriptome analysis, and more, all with unprecedented speed and scalability.

## The Human Genome Project: Objectives

The **Human Genome Project (HGP)** was an international research effort that aimed to sequence the entire human genome. The project spanned 13 years (1990 – April 14, 2003) and was one of the most significant scientific endeavors in history. The goal was to map all 3 billion base pairs of human DNA.

### Key Facts:
- **Timeline**: 1990 to 2003
- **Cost**: Approximately $300 million
- **Leadership**: The project was primarily led by the U.S. Department of Energy (DoE) and the National Institutes of Health (NIH).
- **Collaborators**: The International Human Genome Sequencing Consortium (IHGSC), a group of publicly funded researchers, coordinated the global effort. 

At any given time, around **200 laboratories** in the United States were supporting the sequencing efforts, alongside contributions from more than **18 countries** across the globe.

The HGP provided the foundational reference for understanding human genetics, enabling major advancements in genomics research, medical science, and biotechnology.  

## The Human Genome Project: The Method (WGS)

The Human Genome Project (HGP) utilized two main sequencing approaches: **Hierarchical Genome Shotgun Sequencing** and **Whole-Genome Shotgun Sequencing (WGS)**.

### Hierarchical Genome Shotgun Sequencing:

1. **Shotgun Phase**:
   - The genome was fragmented into larger segments.
   - Clones of these fragments were created using vectors.
   - These clones were then sequenced.
   - Shotgun sequence data was assembled to reconstruct the genome.
   - This phase relied on a **physical map** of the human genome, which had been previously established.

2. **Finishing Phase**:
   - This phase involved **filling in gaps** in the sequence.
   - Ambiguous regions were resolved, ensuring the completion of the human genome sequence.

### Whole-Genome Shotgun Sequencing (Celera Genomics):

1. The genome was sheared randomly into **small fragments** of an appropriate size for sequencing.
2. These fragments were sequenced, and the data was used to **reassemble** the genome.

### Key Differences:
- **Hierarchical Genome Shotgun Sequencing** relied on a predefined physical map of the genome, while **Whole-Genome Shotgun Sequencing** took a more random approach by fragmenting the genome without the need for a physical map.

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/WGSSequencing.png" width="650"/>
</p>

## Sequencing Quality

Sequencing quality depends upon the average number of times each base in the genome is 'read' during the sequencing process.

For the **Human Genome Project (HGP)**:
- **'Draft sequence'**: Covering ~90% of the genome at ~99.9% accuracy.
- **'Finished sequence'**: Covering >95% of the genome at ~99.99% accuracy.

Producing truly high-quality 'finished' sequence by this definition is very expensive and labor-intensive.

Several releases of the human genome sequences...

<details>
  <summary>Finished Genome vs Draft Genome</summary>

### Variable Degrees of Completion of Published Genomes

#### Draft Sequencing
- High-throughput or shotgun phase (whole genome or clone-based approach).
- Assembly using specific algorithms (whole-genome or single-clone assembly).
  - Results in **lower accuracy** than finished sequence.
  - Some segments may be missing or in the wrong order or orientation.

#### Finishing
- Accuracy in bases identification + Quality Check + few if any gaps.
- Contiguous segments of sequence are ordered and linked to one another.
- No ambiguities or discrepancies about segment order and orientation.

#### Complete Genome
A genome represented by a single contiguous sequence with **no ambiguities**.
- The sequences available are finished to a certain high quality.

</details>

## The Human Genome Project: The Heritage

### Open Data Initiative of the Human Genome Project
- The HGP required that all human genome sequence information be **freely and publicly available**.  
- Existing DNA sequences have been stored in databases accessible for anyone to **exploit and analyze**.

### Key Databases
1. **Model Organism Databases**  
   - Sequences of known and hypothetical genes and proteins.  
   - Examples:  
     - **GenBank** (hosted by NCBI).  
     - **Ensembl**: [http://www.ensembl.org](http://www.ensembl.org)  
       - Offers additional data, annotations, visualization, and search tools.

2. **Non-Model Organism Databases**  
   - Community efforts for eukaryotic pathogens:  
     - **EuPathDB**: [http://eupathdb.org/eupathdb/](http://eupathdb.org/eupathdb/).

### Computational Tools
- Specialized computer programs have been developed to **analyze and interpret genomic data**, fostering research advancements.

### Genetic Differences: Single Nucleotide Polymorphisms (SNPs)
- SNPs are the **most common type of genetic variation** in a genome.
- These variations involve differences in individual bases.

### The Goal: Developing a Haplotype Map of the Human Genome
- The Haplotype Map aims to:  
  - **Identify and catalog** most of the millions of SNPs commonly occurring in the human genome.
  - Provide information about:
    - **Described variants** (type and location).  
    - Their **distribution within populations** and **across populations** globally.
  - **Link genetic variants** to the risk for specific diseases.

### The 1000 Genomes Project
- A milestone project designed to map human genetic diversity.  
- Progress has led to a **more complete and reliable dataset** with many **novel variants discovered**.

### 1000 Genomes Project
- [www.1000genomes.org](http://www.1000genomes.org/)  
- **Objectives**:
  - Identify **most genetic variants** with frequencies of at least **1%**.
  - Serve as a **freely accessible resource** for human genetic variation.  
- **Final Data Set**:
  - Covers data for **2,504 individuals** from **26 populations**.
  - Includes **low-coverage sequencing** and **exome sequence data** for all individuals.  
- **Ongoing Usability**:  
  - Managed by the **International Genome Sample Resource (IGSR)** to ensure accessibility of data generated by the project.

---

### UK10K Project
- [www.uk10k.org](http://www.uk10k.org/)  
- **Objectives**:
  - Focus on identifying **rare genetic variants** by studying the DNA of:
    - **4,000 individuals**.
    - **6,000 individuals** with documented diseases (focused on protein-coding regions).  
- **Outcomes**:
  - Establishing links between **genetic variants** and **rare diseases**.
### Part 2: DNA Sequencing in the NGS Era

## DNA Sequencing (DNA-Seq)

### Overview
- DNA-Seq is a highly effective sequencing strategy used today, thanks to **rapid DNA sequencing methods** that have significantly accelerated biological and medical research.
- Applications:
  - Determining the sequence of:
    - Individual genes.
    - Larger genetic regions.
    - Full chromosomes.
    - Entire genomes.
  - Addressing genome complexity with specialized technologies like **Methyl-Seq**, **ChIP-Seq**, **exome sequencing**, etc.

---

### Next-Generation Sequencing (NGS)
- **NGS** (also called **high-throughput sequencing**) describes various modern sequencing technologies offered by different platforms:
  - **Illumina (Solexa)** sequencing.
  - **Roche 454** sequencing.
  - **Ion Torrent**: Proton/PGM sequencing.
  - And more.
- **Advantages**:
  - Faster and cheaper compared to traditional **Sanger sequencing**.
  - Represents a **revolution** in genomics and molecular biology.



### Part 3: Overview of NGS Technologies

## NGS: Sequencing Technologies and Platforms

### Massively Parallel Sequencing

| **Technology** | **Company**            | **Support**            | **Chemistry**        |
|-----------------|------------------------|------------------------|----------------------|
| Solexa          | Illumina              | Bridge PCR on flowcell | Seq-By-Synthesis     |
| 454             | Roche Applied Science | emPCR on beads         | Pyrosequencing       |
| SOLiD           | AB / Life Technologies | emPCR on beads         | Seq-By-Ligation      |
| Ion Torrent     | Life Technologies     | emPCR on beads         | Proton detection     |

### Single Molecule Sequencing

| **Technology** | **Company**              | **Support**            | **Chemistry**        |
|-----------------|--------------------------|------------------------|----------------------|
| PacBio SMRT     | Pacific Biosciences     | Pol performance        | Real-time-Seq        |
| Nanopore        | Oxford Nanopore Tech / McNally | Translocation      | NA                  |

## Solexa Sequencing Workflow

### 1. Sample Preparation
- **Cleaving**: Input DNA is cleaved into short fragments.
- **Adaptor Ligation**: Fragments are ligated to adaptors and then annealed to a slide using the adaptors.
- **Denaturation**: DNA fragments are denatured into single strands for sequencing.

### 2. Sequencing Process
- **Modified Nucleotides**: Modified nucleotides are used, with each emitting a distinct colored light when excited by a laser.
- **Terminators**: Each nucleotide has a terminator, ensuring that only one base is added at a time.

### 3. Amplification
- **Bridge PCR**: The DNA fragments are amplified on the slide using bridge PCR in repeated cycles.

### 4. Imaging
- **Fluorescence Detection**: After each base addition, fluorescence signals are captured to determine the sequence of each DNA fragment.

### 5. Data Analysis
- **Image Processing**: The fluorescence signals are analyzed to reconstruct the DNA sequence.

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Solexa.png" width="650"/>
</p>


[Solexa / Illumina Sequencing Technology Comparison](https://emea.illumina.com/systems/sequencing-platforms/comparison-tool_msm_moved_msm_moved.html)

---
You can also read about the other technologies such as:
- 454 Technology
- SOLiD
- Ion Torrent
- PacBio
- Nanopore

### Four Main Advantages of NGS Over Classical Sanger Sequencing

#### 1. **Speed**
NGS is quicker than Sanger sequencing in two ways:
- The chemical reaction can be combined with the signal detection, whereas in Sanger sequencing, these are two separate processes.
- In Sanger sequencing, only 1 read can be taken at a time, whereas NGS is massively parallel.

#### 2. **Cost**
- The human genome sequence cost $300M with Sanger sequencing.
- Sequencing a human genome with Illumina allows approaching the $1,000 expected cost.

#### 3. **Sample Size**
NGS requires significantly less starting DNA/RNA compared to Sanger sequencing.

#### 4. **Accuracy**
- NGS covers more repeats than Sanger sequencing, resulting in greater coverage, higher accuracy, and sequence reliability. However, individual reads are less accurate for NGS.



### Part 4: DNA-Seq Protocol: Overview

Content for Part 4...

### Part 5: DNA-Seq Analysis Pipeline and File Formats

Content for Part 5...
