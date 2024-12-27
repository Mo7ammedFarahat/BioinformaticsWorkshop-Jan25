- [Variant Calling Tutorial](#variant-calling-tutorial)
  - [1. Install Miniconda](#1-install-miniconda)
  - [2. Set Up the Conda Environment](#2-set-up-the-conda-environment)
  - [3. Install Required Tools](#3-install-required-tools)
  - [4. Download the Data](#4-download-the-data)
  - [Step 1: Quality Control (QC)](#step-1-quality-control-qc)
  - [Step 2: Trimming](#step-2-trimming)
  - [Step 3: Alignment](#step-3-alignment)
  - [Step 4: Post-Alignment Processing](#step-4-post-alignment-processing)
  - [Step 5: Variant Calling](#step-5-variant-calling)
  - [Step 6: Annotating Variants](#step-6-annotating-variants)
  - [BAM File Structure](#structure-of-the-bam-file)
  - [VCF File Structure](#structure-of-the-vcf-file)
  - [Annotation Details](#annotation)


# **Variant Calling Tutorial**

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
- Create and activate an environment for genome assembly:
  ```bash
  conda create -n variantcall -y
  conda activate variantcall
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
Install the essential bioinformatics tools for QC, trimming, alignment, and variant calling:
```bash
conda install fastqc multiqc --channel conda-forge --channel bioconda --strict-channel-priority
conda install -y bioconda::trimmomatic
conda install bwa samtools --channel conda-forge --channel bioconda --strict-channel-priority
conda install igv --channel conda-forge --channel bioconda --strict-channel-priority
conda install snpeff --channel conda-forge --channel bioconda --strict-channel-priority
```

---

### **4. Download the Data:**
Download the necessary data files for *Mycobacterium tuberculosis*:

```bash
curl -O https://zenodo.org/record/3960260/files/004-2_1.fastq.gz
curl -O https://zenodo.org/record/3960260/files/004-2_2.fastq.gz
curl -O https://zenodo.org/record/3960260/files/MTB_ancestor_reference.fasta
curl -O https://zenodo.org/record/3960260/files/Mycobacterium_tuberculosis_h37rv.ASM19595v2.45.chromosome.Chromosome.gff3
```

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

Trimming removes low-quality bases and sequencing adapters, ensuring accurate read alignment and variant calling. Trimmomatic is used for trimming paired-end reads.

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

Reads are aligned to the reference genome to determine their genomic origin. BWA is a fast and accurate aligner suitable for this task.

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

# Index the reference genome
bwa index $REFERENCE_DIR/MTB_ancestor_reference.fasta

# Align paired-end reads to the reference genome
for sample in $(ls $TRIMMED_DIR/*_1.paired.fastq.gz | sed 's/_1.paired.fastq.gz//' | xargs -n 1 basename); do
    bwa mem $REFERENCE_DIR/MTB_ancestor_reference.fasta \
        $TRIMMED_DIR/${sample}_1.paired.fastq.gz $TRIMMED_DIR/${sample}_2.paired.fastq.gz \
        | samtools view -Sb - > $ALIGNMENT_DIR/${sample}.bam
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh mapping.sh
```

---

### **Structure of the BAM File**

BAM (Binary Alignment/Map) files store the results of sequence alignments. A BAM file contains binary representations of the alignment results and includes:

1. **Header**: Information about the reference genome, sequence read groups, and alignment tools.
2. **Alignments**: The actual alignments, including sequence name, position, CIGAR string (indicating how the read maps), mapping quality, and more.

Example BAM file structure:

```
@SQ SN:chr1 LN:248956422
@RG ID:group1 SM:sample1
r001 99 chr1 7 60 8M 4= 15 8 GCTTAGAC 60
r002 0 chr1 16 60 8M 4= 20 8 CTGGAATA 60
```

---

## **Step 4: Post-Alignment Processing**

Sorting and indexing the BAM files makes them suitable for downstream variant calling. SAMtools is used for sorting and indexing BAM files.

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

## **Step 5: Variant Calling**

Variant calling identifies genetic variations relative to the reference genome. SAMtools `mpileup` and `bcftools call` are commonly used for this task.

**Create the variant calling script:**

- Run the following command to create the `variant_calling.sh` script:

```bash
nano variant_calling.sh
```

- Copy and paste the following script:

```bash
#!/bin/bash
# Step 5: Variant Calling

# Directories
REFERENCE_DIR="reference"
ALIGNMENT_DIR="alignments"
VCF_DIR="vcf_files"

mkdir -p $VCF_DIR

# Call variants for each sample
for sorted_bam in $ALIGNMENT_DIR/*.sorted.bam; do
    sample=$(basename $sorted_bam .sorted.bam)
    bcftools mpileup -f $REFERENCE_DIR/MTB_ancestor_reference.fasta $sorted_bam \
        | bcftools call -mv -Ov -o $VCF_DIR/${sample}.vcf
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh variant_calling.sh
```

---

### **Structure of the VCF File**

The VCF (Variant Call Format)

 file contains information about genetic variants such as SNPs, insertions, deletions, and structural variations. It includes:

1. **Header**: Describes the VCF file format version, the reference genome, and annotations.
    #CHROM: The chromosome name where the variant is located.
    POS: The 1-based position of the variant on the chromosome.
    ID: A unique identifier for the variant, often taken from a database (e.g., rsID). A . means no ID is available.
    REF: The reference base(s) at the variant position on the reference genome.
    ALT: The alternative base(s) observed at the variant position.
    QUAL: A quality score for the variant call, representing the confidence of the call (higher is better).
    FILTER: Information about any filters applied to the variant (e.g., PASS if it passed all filters, or specific filter names if it did not). A . indicates no filtering information.
    INFO: Additional information about the variant in key-value pairs, separated by semicolons.
    FORMAT: Specifies the data fields included for each sample in the subsequent columns.
    alignments/004-2.sorted.bam: The sample-specific information. The column name matches the BAM file name or sample ID.
3. **Body**: Contains variant information for each chromosome position, including the reference and alternate alleles, genotype, and quality metrics.

Example VCF file:

```
##fileformat=VCFv4.2
##reference=MTB_ancestor_reference.fasta
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
chr1    1000    .       A       G       60      PASS    .       GT      0/1
chr1    1005    .       T       C       60      PASS    .       GT      1/1
```
For example:
```
Chromosome	2532	.	C	T	225.417	.	DP=171;VDB=0.982727;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,67,69;MQ=60	GT:PL	1/1:255,255,0
```
1. **`Chromosome`**: The chromosome where the variant is located (e.g., `Chromosome`).
2. **`2532`**: The position of the variant on the chromosome.
3. **`.`**: No unique ID is available for this variant.
4. **`C`**: The reference base at this position is Cytosine.
5. **`T`**: The alternative base at this position is Thymine.
6. **`225.417`**: The quality score of the variant call.
7. **`.`**: No filter applied to the variant.
8. **`INFO`**: Provides additional information:
   - `DP=171`: Depth of reads covering the position (171 reads).
   - `VDB=0.982727`: Variant Distance Bias quality metric.
   - `SGB=-0.693147`: Segregation-based quality metric.
   - `MQSBZ=0`: Mapping quality strand bias z-score.
   - `MQ0F=0`: Fraction of reads with a mapping quality of zero.
   - `AC=2`: Allele count for the alternative allele.
   - `AN=2`: Total number of alleles (2 for diploid organisms).
   - `DP4=0,0,67,69`: Breakdown of base counts:
     - `0,0`: Reference allele forward and reverse read counts.
     - `67,69`: Alternate allele forward and reverse read counts.
   - `MQ=60`: Mean mapping quality of the reads covering the position.
9. **`FORMAT`**: Specifies the format of the sample-specific fields:
   - `GT`: Genotype.
   - `PL`: Phred-scaled likelihoods for genotypes.
10. **`1/1:255,255,0`**: Sample-specific data:
    - `1/1`: The genotype is homozygous for the alternate allele (`T` in this case).
    - `255,255,0`: Genotype likelihoods for `0/0`, `0/1`, and `1/1`.

---

### **Step 6: Annotating Variants**

The final step is annotating variants with biological information, such as gene effects. SNPeff is a tool that annotates VCF files.

**Create the annotation script:**

- Run the following command to create the `annotation.sh` script:

```bash
nano annotation.sh
```

- Copy and paste the following script:

```bash
#!/bin/bash
# Step 6: Annotation

# Directories
VCF_DIR="vcf_files"
ANNOTATED_VCF_DIR="annotated_vcf"

mkdir -p $ANNOTATED_VCF_DIR

# Annotate variants using SNPeff
for vcf in $VCF_DIR/*.vcf; do
    sample=$(basename $vcf .vcf)
    snpEff ann -v MTB_ancestor_reference $vcf > $ANNOTATED_VCF_DIR/${sample}.annotated.vcf
done
```

- Save the script by pressing `Ctrl+X`, then `Y`, and hit `Enter`.
- To run the script, use the following command:

```bash
sh annotation.sh
```

### INFO Field:

The `INFO` field contains a variety of annotations

- **DP=171**: The total read depth at this position is 171, meaning 171 reads support this position.
- **VDB=0.982727**: Variant detection bias. A higher value suggests the variant is detected confidently with less bias.
- **SGB=-0.693147**: Strand bias score, which reflects whether there is any strand-specific bias in the variant calls.
- **MQSBZ=0**: Mapping quality bias score for the variant (value of 0 indicates no bias).
- **MQ0F=0**: The fraction of reads with mapping quality 0.
- **AC=2**: Allele count; there are 2 occurrences of this alternate allele in the sample.
- **AN=2**: Allele number; there are 2 total alleles considered in this sample.
- **DP4=0,0,67,69**: Depth of bases in four categories: reference strand (+) and (-), alternate strand (+) and (-).
- **MQ=60**: Mapping quality score, which is the Phred-scaled likelihood of the alignment being correct.
  
---

### Annotation (`ANN`):

The annotation (`ANN`) provides detailed information about the functional impact of the variant and its relation to nearby genes. Here’s a breakdown:

1. **T|synonymous_variant|LOW|dnaN|Rv0002|transcript|CCP42724|protein_coding|1/1|c.481C>T|p.Leu161Leu|481/1209|481/1209|161/402||WARNING_REF_DOES_NOT_MATCH_GENOME**
   - **T**: Indicates this annotation is for the alternate allele "T".
   - **synonymous_variant**: The variant is a synonymous mutation, meaning it does not change the amino acid of the protein (Leu161Leu).
   - **LOW**: Indicates the impact is considered low (synonymous variants often do not affect protein function).
   - **dnaN**: The gene name where the variant is located (likely a gene related to DNA replication or repair).
   - **Rv0002**: The gene ID in the database.
   - **transcript**: Refers to the transcript annotated.
   - **CCP42724**: The transcript ID.
   - **protein_coding**: The type of gene (protein-coding).
   - **c.481C>T**: The exact nucleotide change at the DNA level (c.481, C to T).
   - **p.Leu161Leu**: The protein-level change (Leucine to Leucine, indicating no actual amino acid change).
   - **WARNING_REF_DOES_NOT_MATCH_GENOME**: A warning suggesting that the reference genome used may not match the actual reference sequence for this variant.

2. **T|upstream_gene_variant|MODIFIER|recF|Rv0003|transcript|CCP42725|protein_coding**: 
   - **upstream_gene_variant**: The variant is located upstream (before) the gene "recF", which may affect its expression.
   - **MODIFIER**: This type of variant may have little to no functional impact.
   - **recF**: The gene affected by this variant.
   - **Rv0003**: Gene ID.
   - **transcript**: The transcript ID.
   - **CCP42725**: The transcript ID for this gene.
   - **protein_coding**: Gene type.
   - No coding change here, but the variant may affect gene regulation.

3. **Other annotations for genes like gyrB, gyrA, dnaA**:
   - These annotations describe variants in genes related to DNA replication and repair, all with **upstream_gene_variant** or **downstream_gene_variant** effects, suggesting that these variants could affect gene expression or regulation without directly altering the protein coding sequence.

---

### Summary of the Variant:

- **Location**: Chromosome 2532, with a C>T variant.
- **Genetic Impact**: 
  - The variant is **synonymous** in the **dnaN** gene, meaning it does not alter the protein sequence (Leu161Leu).
  - There are several **upstream** and **downstream** gene variants in other genes (recF, gyrB, gyrA, dnaA), which might affect gene regulation.
- **Genotype**: The sample is homozygous for the alternate allele (T), with a high-quality score (225.417) and high read depth (171).
- **Warnings**: There are warnings indicating that the reference genome used does not match the current genome sequence, and there are some missing start codons in certain transcripts.
