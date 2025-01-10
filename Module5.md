# Next Generation Sequencing and its Applications (RNA-Seq)

## Tutorial Outline

### I. RNA-Seq Overview
- **RNA-Seq Background**
- **Experimental Design Considerations**

### II. RNA-Seq Analysis Workflow
- Mapping
- Expression Quantification
- Read Count Normalisation
- Differential Expression Analysis

### III. Downstream Analysis
- Functional Enrichment Analysis

---

## Learning Outcomes
By the end of this module, you should be able to:
- Understand the key factors to consider when designing an RNA-seq experiment.
- Describe and execute a bioinformatics analysis workflow for bulk RNA-seq, from sequenced reads to differentially expressed genes.

---

## Central Dogma

The **Central Dogma of Molecular Biology** describes the flow of genetic information within a cell, outlining how DNA directs the synthesis of proteins. It begins with **transcription**, where the DNA sequence is transcribed into messenger RNA (mRNA) by RNA polymerase in the nucleus. The mRNA then carries the genetic code to the cytoplasm, where it undergoes **translation**. During translation, ribosomes read the mRNA sequence in codons (three-nucleotide segments) and assemble the corresponding amino acids to form proteins, which perform various cellular functions.

This process is unidirectional, summarized as **DNA → RNA → Protein**, though exceptions exist, such as retroviruses using reverse transcription to convert RNA back into DNA. This framework is fundamental to understanding how genetic information is expressed and regulated in living organisms.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/centralDogma.png" width="650"/>
</p>

---

## What is the Transcriptome of a Cell?

The **transcriptome** of a cell refers to the complete set of RNA molecules transcribed from the genome at a specific time and under specific conditions. It includes all types of RNA, such as:
- Messenger RNA (mRNA)
- Transfer RNA (tRNA)
- Ribosomal RNA (rRNA)
- Non-coding RNAs (e.g., microRNAs and long non-coding RNAs)

The transcriptome reflects which genes are actively being expressed, providing a snapshot of the cell's functional state and how it responds to internal and external stimuli. Unlike the genome, which is relatively stable, the transcriptome is dynamic and varies between different cell types, developmental stages, and environmental conditions. Transcriptome analysis is commonly used to study gene expression and regulation, often through techniques like RNA sequencing (RNA-seq) or microarrays.

---

## How Can We Read the Whole Transcriptome?

### RNA Sequencing (RNA-Seq)
RNA sequencing (RNA-seq) is the application of next-generation sequencing technologies to cDNA molecules, which are obtained by reverse transcription from RNA, in order to get information about the RNA content of a sample. RNA-Seq has fast become the preferred method for measuring gene expression, providing an accurate proxy for absolute quantification of messenger RNA (mRNA) levels within a sample. It has rapidly matured in data handling, quality control (QC), and downstream statistical analysis methods.

### RNA-seq Applications
RNA-seq has numerous applications in biological and biomedical research. It is widely used for:
- **Gene Expression Profiling**: Quantifying gene activity under various conditions.
- **Differential Expression Analysis**: Identifying genes that are upregulated or downregulated between groups.
- **Alternative Splicing Events**: Discovering novel splicing events and isoforms.
- **Non-coding RNAs**: Studying miRNAs, lncRNAs, etc.
- **Single-cell Transcriptomics**: Revealing cellular heterogeneity.
- **Disease Mechanisms**: Exploring disease pathways, identifying RNA-based biomarkers.
- **Drug Development**: Evaluating treatment impacts on gene expression and identifying therapeutic targets.

---

## Best Experimental Design for an RNA-Seq Experiment

To ensure reliable results, a successful RNA-Seq study starts with a good experimental design. Here are key considerations:

### 1. Sample Selection and Experimental Criteria
- **Biological Replicates**:
  - **Biological Replicates vs Technical Replicates**:
    - **Biological Replicates**: Independent samples from biologically distinct individuals or units under the same conditions. They account for natural biological variability and ensure results reflect broader biological patterns.
    - **Technical Replicates**: Repeated measurements of the same biological sample to assess technical variability in processes like RNA extraction, library preparation, or sequencing.

#### Key Differences Between Biological and Technical Replicates

| Aspect                | Biological Replicates            | Technical Replicates            |
|-----------------------|-----------------------------------|---------------------------------|
| **Source**            | Different biological entities     | Same biological entity          |
| **Purpose**           | Measure biological variability    | Measure technical variability   |
| **Example**           | RNA from different individuals    | RNA sequenced twice from one sample |
| **Importance**        | Ensures generalizability          | Ensures protocol reliability    |

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/replica.png" width="650"/>
</p>

#### How Many Replicates Do I Really Need?
- **Increasing sequencing depth** enhances the ability to detect low-expression transcripts.
- **Returns diminish** beyond a certain sequencing depth.
- **Increasing biological replicates** improves the accuracy of expression levels, especially for low-expression transcripts, and reduces the coefficient of variation.

**General Rule of Thumb**:
- **Differential Gene Expression**: 10-30M reads (SE 50-75bp).
- **Alternative Splicing**: 50-100M reads (PE, 2x75bp).

---

- **Sample Diversity**: Ensure samples represent the population or condition of interest (e.g., age, sex, developmental stage).
- **Clear Hypothesis**: Define a research question to guide sample grouping (e.g., treatment vs. control, diseased vs. healthy).
- **Randomization**: Randomize sample collection to minimize bias.
- **Power Analysis**: Determine sample size based on effect size and variability.

#### Experimental Planning Considerations
- Were all RNA isolations performed on the same day?
- Were all library preparations performed on the same day?
- Did the same person perform the RNA isolation/library preparation?
- Did you use the same reagents and protocols for all samples?
- Were all samples processed in the same location?

If any answer is "No," you may have batch effects.

---

## RNA Extraction
- **Total RNA = mRNA + rRNA + tRNA + regulatory RNAs...**
- Ribosomal RNA (rRNA) can represent > 90% of total RNA.
- Be mindful of the protocol used (some will remove small RNAs).
- PCR amplification can reduce coverage of regions with high GC content; consider using amplification-free protocols.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/rnaSeq.png" width="650"/>
</p>

---

## RNA Fragmentation
- Fragment RNA into smaller pieces (~200 bp) to facilitate efficient library preparation and sequencing.

---

## Library Preparation
- Convert RNA to complementary DNA (cDNA) via reverse transcription.
- Add adaptors to cDNA fragments for sequencing library preparation.
- Amplify the library using PCR, minimizing amplification bias.

### What Type of Sequencing Library is Best for Your Experimental Question?

- **Stranded vs. Unstranded**: Stranded protocols are better for distinguishing antisense or overlapping transcripts.
- **Single vs. Paired-End**: Paired-end sequencing is better for de novo transcript discovery or isoform expression analysis. More than 55% of reads will span two or more exons.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/library.png" width="650"/>
</p>

---

## Sequencing
- Choose sequencing depth based on the goals of the experiment:
  - **Gene Expression Studies**: ~20 million reads/sample.
  - **Transcript Discovery/Alternative Splicing Analysis**: ~50 million reads/sample.
- **Paired-End Sequencing**: Provides better assembly and mapping accuracy.
- Optimize **multiplexing** to balance cost and coverage.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/overview.png" width="650"/>
</p>

---

## Data Analysis
### Upstream Analysis
Upstream analysis in any workflow refers to the initial steps of processing and analyzing raw sequencing data to prepare it for downstream, more advanced analyses. It focuses on cleaning, organizing, and aligning raw data to ensure it is of high quality and biologically meaningful for interpretation. The common steps in Upstream analysis are

#### Quality Control (QC):
In RNA sequencing (RNA-seq), quality control (QC) is a critical step to ensure the accuracy and reliability of downstream analyses. QC involves evaluating the quality of raw sequencing reads, identifying and removing low-quality reads, adapter contamination, and other artifacts. For short reads, tools like FastQC are widely used to assess quality metrics such as per-base quality scores, GC content, and sequence duplication levels. For long reads, tools like NanoPlot and pycoQC are utilized for visualizing read quality, read length distributions, and error rates specific to platforms like Oxford Nanopore or PacBio. 
The sequencing quality score of a given base, Q, is defined by the following equation:

 *Q = -10log10(e)*

where e is the estimated probability of the base call being wrong. Higher Q scores indicate a smaller probability of error. Lower Q scores can result in a significant portion of the reads being unusable. They may also lead to increased false-positive variant calls, resulting in inaccurate conclusions. A quality score of 20 (Q20) represents an error rate of 1 in 100 (meaning every 100 bp sequencing read may contain an error), with a corresponding call accuracy of 99%. When sequencing quality reaches Q30, virtually all of the reads will be perfect, with no errors or ambiguities. This is why Q30 is considered a benchmark for quality in next-generation sequencing (NGS).

    Q20: 99% accuracy (1 error in 100 bases)
    Q30: 99.9% accuracy (1 error in 1000 bases)

#### Trimming 
It is a critical preprocessing step in RNA-seq and other sequencing workflows. It involves removing unwanted or low-quality portions of sequencing reads to improve the quality and accuracy of downstream analyses, such as alignment or quantification. Trimming ensures that artifacts like adapter sequences, low-quality bases, or overly short reads do not introduce bias or errors in the analysis.

##### Why Trimming is Important:
1. **Adapter Removal**: During library preparation, adapter sequences may be partially or fully sequenced. These non-biological sequences can interfere with alignment.
2. **Low-Quality Bases**: Sequencing quality tends to decrease at the ends of reads. Removing these low-quality bases prevents alignment errors.
3. **Short Read Filtering**: Reads below a certain length may lack informative content and are often discarded to avoid misalignments.
4. **Poly-N Sequences**: Long stretches of "N" (unknown bases) can occur due to sequencing errors and are trimmed to improve data quality.

##### Tools for Trimming:
1. **Cutadapt**: Specifically designed for adapter trimming, offering flexibility for custom adapter sequences.
2. **Trimmomatic**: A versatile tool for quality and adapter trimming with customizable parameters.

#### Indexing
**Indexing** in bioinformatics is the process of creating a structured and searchable representation of a reference genome or transcriptome to facilitate rapid and efficient alignment of sequencing reads. This step speeds up the mapping process by organizing the reference into smaller, manageable units, such as k-mers, or advanced data structures like suffix arrays or Burrows-Wheeler Transform (BWT), to reduces thr computational overhead by allowing alignment tools to quickly find matches for input reads without searching the entire genome. For RNA-seq, splice-aware indexing is crucial, as it incorporates exon-exon junctions (derived from annotations like GTF files), enabling accurate alignment of spliced reads. Tools like **HISAT2** abd **STAR** are commonly used for indexing in RNAseq, with each optimized for specific applications. Overall, indexing is a critical preprocessing step that ensures the scalability, speed, and accuracy of downstream alignment in sequencing workflows.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/indexing.png" width="650"/>
</p>

#### Mapping
Mapping in RNA-seq refers to the process of aligning sequencing reads to a reference genome or transcriptome to identify their origin. It serves multiple purposes, including locating where in the genome the reads originated, assessing the overall quality of RNA-seq data, and exploring the structure and expression of genes of interest. Mapping provides critical insights into the organization of transcripts and allows researchers to quantify gene expression accurately.

**Mapping Options**  
Mapping can be performed against a reference genome or transcriptome, depending on the research goals. Alternatively, transcriptomes can be reconstructed directly from RNA-seq data. This reconstruction can be **genome-guided**, leveraging a reference genome, or **de novo**, without using a reference, which is particularly useful for non-model organisms lacking well-annotated genomes.

**Challenges in Mapping**  
Mapping RNA-seq data presents several challenges. Genomic variations, sequencing errors, and non-unique sequences can complicate accurate read alignment. Additionally, the presence of introns in eukaryotic genes—absent in mature mRNA—requires special handling, as conventional mapping tools may fail to account for spliced regions.

**The Need for Splice-Aware Algorithms**  
Given the complexity of eukaryotic gene structures, **splice-aware mapping algorithms** are essential. These tools, such as HISAT2 and STAR, are designed to accurately map reads spanning exon-exon junctions, ensuring reliable representation of spliced transcripts. These algorithms are indispensable for analyzing RNA-seq data from eukaryotic organisms, where correct handling of introns and splicing patterns is critical for accurate transcriptomic analysis.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/mapping-1.png" width="650"/>
</p>

#### Feature counting
It is the final step in the upstream analysis, following mapping step in RNA-seq, where the number of reads aligned to specific genomic features (e.g., genes, exons) is quantified to measure expression levels. This process typically uses the SAM or BAM file, which contains information about the mapped reads, as input. Feature counting tools like **HTSeq-count** or **FeatureCount** determine how many reads align to each feature, accounting for overlaps and complexities. When reads align to or overlap multiple features, different strategies are used to handle ambiguity. The **union** mode assigns reads to a feature if they overlap with any part of it, **intersect-strict** counts only those reads that fall entirely within a single feature, and **intersect-nonempty** counts reads overlapping with any feature but requires partial overlap.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/featurecount.png" width="650"/>
</p>


---
### Downstream analysis
Downstream analysis refers to the suite of computational and statistical processes performed after the initial steps of quality control, read alignment, and feature counting which are aimed at interpreting the data to gain biological insights into gene expression patterns, regulatory mechanisms, or functional roles. In this analysis we transform raw expression data into meaningful results, guiding hypotheses and conclusions.

#### Normalization

**Normalization** in RNA-seq can adjusts raw gene expression data to account for technical biases and differences between samples to allow for accurate comparisons across conditions. So, without normalization, raw read counts would be influenced by factors like sequencing depth, gene length, and technical variation, leading to misleading interpretations. Normalization methods adjust for these biases by transforming data so that the expression levels can be reliably compared across samples or experimental conditions.

### Ok Why Normalization is Needed?:
1. **Sequencing Depth Bias**: Samples with higher sequencing depth (more total reads) will naturally have more reads mapped to each gene, even if the actual gene expression levels are the same. This introduces a bias, making samples with higher read depth seem like they have higher gene expression.

- A sample contains a large amount of RNA
- Each gene contributes a portion of the total amount of reads
- We assume the reads are randomly distributed in the sample
- The sample is an approximation of the true proportions
- Genes from higher depth samples appear to have higher expression
- We need to compensate for different sequencing depths
- Account for more uncertainty in samples with lower depths

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/normalization-seq.png" width="650"/>
</p>

3. **Gene Length Bias**: Longer genes tend to have more reads mapped to them simply due to their larger size. This bias affects the comparison of gene expression levels between genes of varying lengths.


<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/normalization-gene.png" width="650"/>
</p>

**Normalization methods**
Normalization methods correct for these biases by scaling the data by adjusting the counts based on the total read count (sequencing depth) or the length of the genes. The common normalization methods include:
- **TPM (Transcripts per Million)**
TPM (Transcripts per Million) is a normalization method used in RNA-seq to quantify gene expression levels, taking into account both sequencing depth and gene length.

TPM = read counts / the length of each gene in kilobases
- **RPKM (Reads per Kilobase of transcript per Million mapped reads)**
It adjusts two primary biases: sequencing depth and gene length.

<p align="center">
  <img src="https://github.com/Mo7ammedFarahat/MASRI-Jan25/blob/main/images/Module5/tpm-rkpm.png" width="650"/>
</p>

- **FPKM (Fragments per Kilobase of transcript per Million mapped reads)**
 takes into account that two reads can map to one fragment (and so it doesn’t count this  fragment twice)

### Differential expression analysis
It aims to identify genes that are expressed at significantly different levels between different conditions, experimental groups, or treatments to understand the biological processes underlying these conditions. Statistical models are used to identify genes whose expression levels significantly differ between experimental groups. Commonly used tools include:
     - **DESeq2**: Uses a negative binomial distribution and a generalized linear model to estimate the distribution of counts and test for differential expression.
     - **edgeR**: Similar to DESeq2 but uses a method based on empirical Bayes.
     - **limma-voom**: A method for RNA-seq data based on linear models, which applies a variance-stabilizing transformation (voom).
These tools typically generate a **p-value** for each gene, which reflects the likelihood that the observed differences in gene expression are due to random chance.
   - **Multiple Testing Correction**: Since many genes are tested simultaneously, a correction for multiple comparisons is applied (e.g., **Benjamini-Hochberg** false discovery rate (FDR)) to control for Type I errors.

**Result Interpretation**:
   - **Volcano Plots**: These plots show the relationship between the log fold change (x-axis) and the -log10 p-value (y-axis) for each gene. They help identify genes that are both significantly differentially expressed and biologically relevant.
     
   - **Heatmaps**: Heatmaps can visualize the expression patterns of differentially expressed genes across samples. Genes are usually clustered based on similar expression patterns, highlighting potential biological similarities between groups.
   - **MA Plots**: MA plots show the relationship between the average gene expression level and the fold change. Genes with high fold changes are typically more interesting, but both low and high expression levels are analyzed.
   - **Gene Ontology (GO) and Pathway Analysis**: To interpret the biological significance of differentially expressed genes, **GO enrichment analysis** or **pathway analysis** (e.g., using tools like **DAVID**, **Reactome**, or **KEGG**) is performed. This helps to identify biological processes, molecular functions, or cellular components that are overrepresented among the differentially expressed genes.
   - **PCA (Principal Component Analysis)**: PCA is often used to visualize the overall variance in gene expression data. It can help identify outliers or batch effects that might need to be corrected for.
   - **List of Differentially Expressed Genes**: This includes genes that show significant differences in expression between experimental conditions. The list typically includes:
   -  **Gene IDs**
   - **Log fold change (LFC)**: Indicates how much more or less a gene is expressed in one condition compared to another.
   - **P-value**: Indicates the probability that the observed difference is due to chance.
   - **Adjusted p-value (FDR)**: The p-value corrected for multiple testing.
   - Genes with a **low p-value** (typically < 0.05 or adjusted < 0.1) are considered statistically significant.
   - Genes with **high fold change** (either upregulated or downregulated) are often the most biologically relevant, particularly if they are involved in important pathways or processes.

3. **Enrichment Analysis**:
   - **Enriched Pathways**: Pathway analysis identifies biological pathways or processes that are overrepresented among the differentially expressed genes, providing insights into the underlying biology.
   - **Gene Ontology (GO)** terms can help elucidate the molecular functions, biological processes, and cellular components pathways associated with the differentially expressed genes.

---
## Practical Session
### **1. Install Miniconda:**
To begin, we'll install Miniconda, which is a minimal installer for Conda.

- Open a terminal (`Ctrl+T`) and run the following command:
  ```bash
  curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  ```
- Navigate to your home directory and execute the installation script:
  ```bash
  sh Miniconda3-latest-Linux-x86_64.sh
  ```
  - Approve all prompts by typing `Y`.

---

### **2. Set Up the Conda Environment:**
- Create and activate an environment for RNA-Seq analysis:
  ```bash
  conda create -n rnaseq -y
  conda activate rnaseq
  ```
- Configure Conda channels:
  ```bash
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda config --set channel_priority strict
  ```

---

### **3. Install Required Tools:**
Install the essential bioinformatics tools for QC, trimming, alignment, and counting reads.

```bash
conda install fastqc multiqc --channel conda-forge --channel bioconda --strict-channel-priority
conda install -y bioconda::trimmomatic
conda install star samtools --channel conda-forge --channel bioconda --strict-channel-priority
conda install subread --channel conda-forge --channel bioconda --strict-channel-priority
```

---

### **4. Download the Data:**
Download the necessary data files for RNA-Seq analysis:

This data presentsan experiment profiling Drosophila cells after the depletion of a regulatory gene. In the study by Brooks et al. (2011), the authors used RNA-Seq to identify genes and pathways regulated by the Pasilla gene (the Drosophila equivalent of mammalian splicing regulators Nova-1 and Nova-2). They depleted the Pasilla (PS) gene in Drosophila melanogaster through RNA interference (RNAi). Total RNA was extracted and used to prepare both single-end and paired-end RNA-Seq libraries for both treated (PS depleted) and untreated samples. The libraries were sequenced to generate RNA-Seq reads for each sample. By comparing the RNA-Seq data from the treated and untreated samples, the effects of Pasilla gene depletion on gene expression were assessed.

These datasets are subsampled from the original dataset

```bash
curl -O https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger
curl -O https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger
curl -O https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger
curl -O https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf.gz

```
OUr datasets are:
- GSM461177 (untreated): GSM461177_1 and GSM461177_2
- GSM461180 (treated): GSM461180_1 and GSM461180_2

Please rename the extention to be (_1 or _2 .fastq.gz)
Make sure that your samples are placeed in a Directory called (datqasets), and the reference genome and annotation files are in (ref) Directory.

---

## **Step 1: Quality Control (QC)**

Quality control ensures that your sequencing data is free from contaminants, has high base quality, and is suitable for downstream analysis. FastQC generates QC reports, and MultiQC aggregates multiple reports into a summary.

**Create the QC script:**

- Run the following command to create the `qc.sh` script:

```bash
nano qc.sh
```

- Copy and paste the following script into your terminal:

```bash
#!/bin/bash
# Step 1: Quality Control (QC)

# Directories
DATASET_DIR="datasets"
QC_OUTPUT_DIR="qc_reports"

mkdir -p $QC_OUTPUT_DIR

# Run FastQC for all FASTQ files
for file in $DATASET_DIR/*.fastq.gz; do
    fastqc -o $QC_OUTPUT_DIR $file
done

# Aggregate QC reports using MultiQC
multiqc -o $QC_OUTPUT_DIR $QC_OUTPUT_DIR
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command in your terminal:

```bash
sh qc.sh
```

---

## **Step 2: Trimming**

Trimming removes low-quality bases and sequencing adapters, ensuring accurate read alignment and downstream analysis. Trimmomatic is used for trimming paired-end reads.

**Create the trimming script:**

- Run the following command to create the `trimming.sh` script:

```bash
nano trimming.sh
```

- Copy and paste the following script:

```bash
#!/bin/bash
# Step 2: Trimming

# Directories
DATASET_DIR="datasets"
TRIMMED_DIR="trimmed_reads"

mkdir -p $TRIMMED_DIR

# Trimming paired-end reads
for sample in $(ls $DATASET_DIR/*_1.fastq.gz | sed 's/_1.fastq.gz//' | xargs -n 1 basename); do
    TrimmomaticPE -phred33 \
        $DATASET_DIR/${sample}_1.fastq.gz $DATASET_DIR/${sample}_2.fastq.gz \
        $TRIMMED_DIR/${sample}_1.paired.fastq.gz $TRIMMED_DIR/${sample}_1.unpaired.fastq.gz \
        $TRIMMED_DIR/${sample}_2.paired.fastq.gz $TRIMMED_DIR/${sample}_2.unpaired.fastq.gz \
        ILLUMINACLIP:primers_and_adabters/primers.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:3:25 MINLEN:36
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh trimming.sh
```

---

## **Step 3: Alignment**

Reads are aligned to the reference genome to determine their genomic origin. STAR is a fast aligner that supports RNA-Seq data.

**Create the mapping script:**

- Run the following command to create the `mapping.sh` script:

```bash
nano mapping.sh
```

- Copy and paste the following script:

```bash
#!/bin/bash
# Step 3: Alignment

# Directories
REFERENCE_DIR="reference"
TRIMMED_DIR="trimmed_reads"
ALIGNMENT_DIR="alignments"

mkdir -p $ALIGNMENT_DIR

# Generate STAR index for the reference genome
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $REFERENCE_DIR --genomeFastaFiles $REFERENCE_DIR/genome_reference.fasta --sjdbGTFfile $REFERENCE_DIR/genome_annotation.gff --sjdbOverhang 100

# Align paired-end reads to the reference genome
for sample in $(ls $TRIMMED_DIR/*_1.paired.fastq.gz | sed 's/_1.paired.fastq.gz//' | xargs -n 1 basename); do
    STAR --runThreadN 8 --genomeDir $REFERENCE_DIR --readFilesIn $TRIMMED_DIR/${sample}_1.paired.fastq.gz $TRIMMED_DIR/${sample}_2.paired.fastq.gz --outFileNamePrefix $ALIGNMENT_DIR/${sample}_ --outSAMtype BAM SortedByCoordinate
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh mapping.sh
```

---

## **Step 4: Post-Alignment Processing**

Sorting and indexing the BAM files makes them suitable for downstream read counting. SAMtools is used for sorting and indexing BAM files.

**Create the BAM indexing script:**

- Run the following command to create the `bamindex.sh` script:

```bash
nano bamindex.sh
```

- Copy and paste the following script:

```bash
#!/bin/bash
# Step 4: Post-Alignment Processing

# Directories
ALIGNMENT_DIR="alignments"

# Sort and index BAM files
for bamfile in $ALIGNMENT_DIR/*.bam; do
    samtools sort -o ${bamfile%.bam}.sorted.bam $bamfile
    samtools index ${bamfile%.bam}.sorted.bam
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh bamindex.sh
```

---

## **Step 5: Read Counting**

FeatureCounts is used for counting reads that align to genomic features, such as exons or genes.

**Create the read counting script:**

- Run the following command to create the `counting.sh` script:

```bash
nano counting.sh
```

- Copy and paste the following script:

```bash
#!/bin/bash
# Step 5: Read Counting

# Directories
ALIGNMENT_DIR="alignments"
ANNOTATION_FILE="reference/genome_annotation.gff"
COUNT_DIR="read_counts"

mkdir -p $COUNT_DIR

# Count reads for each BAM file
for bamfile in $ALIGNMENT_DIR/*.sorted.bam; do
    sample=$(basename $bamfile .sorted.bam)
    featureCounts -T 8 -a $ANNOTATION_FILE -o $COUNT_DIR/${sample}_counts.txt $bamfile
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh counting.sh
```

---

### **Step 6: Differential Expression Analysis**

Once you have your count data, you can perform differential gene expression (DGE) analysis to identify genes that are differentially expressed between experimental conditions.

### **In R install DESeq2**

First, make sure you have **R**, **DESeq2**, and the required packages installed. If not, you can install them from the CRAN repository and Bioconductor.

- **Install Packages:**
  Open R or RStudio and run the following command to install the `DESeq2` package and the other packages from Bioconductor and CRAN:

  ```r
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
  install.packages(ggplot2)
  install.packages(pheatmap)
  install.packages(dplyr)
  ```

---

### **Prepare Count Data for DESeq2**

To perform DGE analysis using DESeq2, you'll need to format your count data in a matrix where rows are genes, and columns are samples in a tsv file. Make sure you have the metadata also.

### **Start Your Analysis**

```r
# Step 1: Install and load required libraries
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(dplyr)

# Step 2: Load count matrix and metadata
# Assume you have a CSV file for counts and a CSV file for metadata
counts <- read.csv("/path-to/PRJNA143369.tsv", row.names = 1, sep = "\t")  # Counts matrix (genes in rows, samples in columns)
metadata <- read.csv("/path-to/Meta.tsv", row.names = 1, sep = "\t")    # Metadata (samples in rows, experimental conditions in columns)

# Step 3: Data filtering
# Filter out genes with low counts across all samples
# Keep genes with at least 10 counts in at least 3 samples (adjust threshold as needed)
filtered_counts <- counts[rowSums(counts >= 10) >= 3, ]

# Step 4: Create DESeqDataSet
# Create DESeqDataSet object from the filtered counts matrix and metadata
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = metadata,
                              design = ~ CASE)  # Adjust 'condition' to your experimental condition column name

# Step 5: Normalization (using DESeq2's internal normalization)
dds <- DESeq(dds)

# Step 6: Differential Expression Analysis
res <- results(dds)

# Step 7: Filter significant results
# Adjust for multiple testing (padj < 0.05)
res_sig <- res[which(res$padj < 0.05), ]

# Print the significant results
print(head(res_sig))

# Step 8: Create Volcano Plot
volcano_data <- data.frame(log2FoldChange = res$log2FoldChange,
                           pvalue = res$pvalue,
                           padj = res$padj)

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(padj < 0.05, 
                                ifelse(log2FoldChange > 0, "red", "blue"), "gray")), 
             alpha = 0.7) +
  scale_color_manual(values = c("gray" = "gray", "red" = "red", "blue" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme(legend.position = "none")

# Step 9: Create Heatmap for Top 100 Significant Genes
top_genes <- rownames(res_sig[order(res_sig$padj),])[1:100]
heatmap_data <- assay(dds)[top_genes, ]

# Generate heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row", 
         main = "Heatmap of Top 100 Significant Genes", 
         annotation_col = metadata)

# Step 10: Create PCA Plot
# Perform PCA
rld <- rlog(dds, blind = FALSE)  # Log-transformation of the data
pca_data <- plotPCA(rld, intgroup = "CASE", returnData = TRUE)  # Adjust 'condition' if necessary

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = CASE)) + 
  geom_point(size = 3) + 
  labs(title = "PCA Plot") + 
  theme_minimal()

# Step 11: Create MA Plot
# Plot the MA plot (log2 fold change vs mean normalized counts)
# Prepare the data for ggplot2
res_df <- as.data.frame(res)
res_df$color <- ifelse(res_df$padj < 0.05, 
                       ifelse(res_df$log2FoldChange > 0, "red", "blue"), 
                       "gray")

# Create the MA plot using ggplot2
ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = color)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_identity() +  # Use the color assigned in the data frame
  scale_x_log10() +  # Log scale for baseMean
  labs(title = "MA Plot", x = "Mean of Normalized Counts (log scale)", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(plot.title = element_text(hjust = 0.5))


# Subset results for significant genes (padj < 0.05)
sig_genes <- subset(res, padj < 0.05)

# Export the significant genes to a CSV file
write.csv(sig_genes, file = "significant_genes.csv", row.names = TRUE, quote = F)

# You can perform functional enrichment analysis for your genes using ShinyGO
```
---

### **That's it! You've now completed an RNA-Seq data analysis pipeline.**

