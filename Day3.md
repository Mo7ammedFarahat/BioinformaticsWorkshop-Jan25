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

![Hierarchical Genome Shotgun Sequencing](https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/WGSSequencing.png)

<p align="center">
<img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/unix.png?raw=true" width="650"/>
</p>
### Part 2: DNA Sequencing in the NGS Era

Content for Part 2...

### Part 3: Overview of NGS Technologies

Content for Part 3...

### Part 4: DNA-Seq Protocol: Overview

Content for Part 4...

### Part 5: DNA-Seq Analysis Pipeline and File Formats

Content for Part 5...
