# Bulk RNA-Seq Data Analysis – Hypoxia Response in LNCaP and PC3 Cell Lines

## 🎯 Project Overview

This project demonstrates step-by-step analysis of Bulk RNA-Seq data, starting from raw SRA files to obtaining differentially expressed genes (DEGs).<br>
The dataset help us to investigate how hypoxia (low oxygen conditions) alters gene expression in two prostate cancer cell lines:

**1.LNCaP (Lymph Node Carcinoma of the Prostate)**<br>
Origin: Derived from a lymph node metastasis of a human prostate adenocarcinoma (prostate cancer spread to the lymph node).

**2.PC3 (Prostate Cancer 3)**<br>
Origin: Derived from a bone metastasis of a grade IV prostate adenocarcinoma.

Here we comparing Normoxia vs Hypoxia in each cell line.

## 📌 Objectives

- To learn the basic workflow of Bulk RNA-Seq analysis.
- To perform quality control, alignment, quantification, and differential expression analysis.
- To visualize results using plots such as PCA, heatmaps, and volcano plots and pathway enrichment analysis of DEGs.
- To document the process in a way that is clear for beginners to reproduce.

## 📂 Dataset Information

Dataset Source: GSE106305 (NCBI GEO)
 (Guo et al., Nature Communications, 2019)

## How the raw data is organized in the database

- Each biological sample (e.g., LNCaP Normoxia replicate 1) is not stored as one file.
- Instead, it is split into multiple technical runs (SRR IDs) on the SRA database.

**For example:**<br>
LNCaP Normoxia Replicate 1 (GSM3145509) is split across 4 SRR files:
**SRR7179504, SRR7179505, SRR7179506, SRR7179507.**

This is common in sequencing experiments, where one sample is sequenced in multiple lanes/runs.

**What we downloaded**

- In total, we downloaded 20 SRR files.
```
SRR7179504, SRR7179505, SRR7179506, SRR7179507
SRR7179508, SRR7179509, SRR7179510, SRR7179511
SRR7179520, SRR7179521, SRR7179522, SRR7179523
SRR7179524, SRR7179525, SRR7179526, SRR7179527
SRR7179536, SRR7179537, SRR7179540, SRR7179541
```
- These represent 8 biological samples (2 replicates for each condition).


## ⚙️ Steps to be followed 

### 1️⃣ Data Download & Conversion
- **Install** SRA Toolkit  
- **Download** raw `.sra` files from NCBI using accession numbers  
- **Convert** `.sra` → compressed `.fastq.gz` files  
- *(Optional)* Automate multiple downloads with a Python script  

**Output:** FASTQ files (raw sequencing reads)
<img width="1237" height="144" alt="image" src="https://github.com/user-attachments/assets/64fd647b-8579-4470-836d-3da17cddb781" />

---

### 2️⃣ Quality Control (QC)
- Run **FastQC** on all FASTQ files to check:  
  - Per-base quality scores  
  - GC content  
  - Adapter contamination  
  - Overrepresented sequences  
- Use **MultiQC** to summarize results into one report  

**Output:** HTML reports showing sequencing quality
<img width="1737" height="250" alt="Screenshot 2025-08-24 151055" src="https://github.com/user-attachments/assets/b4980b8f-09ee-4c0f-8a2f-38d3a3479d96" />

---

### 3️⃣ Read Trimming (if needed)
- Use **Trimmomatic** or **fastp** to remove:  
  - Low-quality bases  
  - Adapter contamination  
- Re-run FastQC to confirm improvements  

**Output:** Cleaned FASTQ files

<img width="368" height="45" alt="image" src="https://github.com/user-attachments/assets/5f93cfc9-a75c-43b5-8e71-edb665c26ef3" />

---

### 4️⃣ Sample Preparation
- **Merge** multiple technical replicates (several SRR files per sample) into one file  
- **Rename** files clearly (e.g., `LNCAP_Normoxia_S1.fastq.gz`)  

**Output:** 8 final FASTQ files (one per replicate)
<img width="983" height="280" alt="22222222222222" src="https://github.com/user-attachments/assets/a72b7de8-f5cf-4231-ad8e-b6bb95de383e" />

---

### 5️⃣ Alignment to Reference Genome
- **Download** HISAT2 prebuilt genome index for GRCh38 (Human)  
- **Install** HISAT2 + Samtools  
- **Align** each FASTQ to the reference genome → SAM/BAM files  
- **Sort & Index** BAM files for downstream use  

**Output:** Aligned, sorted, and indexed BAM files
<img width="1831" height="629" alt="Screenshot 2025-08-24 152117" src="https://github.com/user-attachments/assets/b4099a07-32ec-4786-8b6b-eb2b8f635a7b" />

---

### 6️⃣ Quantification of Gene Expression
- Use **featureCounts** (Subread package) to assign reads to genes  
- **Input:** BAM files + GRCh38 GTF annotation  
- **Output:** Count matrix (genes × samples)  

**Result:**  raw counts of all the samples which is then used to generate final count matrix using python script.
<img width="1614" height="194" alt="Screenshot 2025-08-24 152222" src="https://github.com/user-attachments/assets/d9839eae-b8c7-4c5d-b097-0ffb1db688bd" />

---

### 7️⃣ Post-alignment Quality Check
- Use **Qualimap RNA-seq** to evaluate:  
  - Mapping statistics  
  - Coverage  
  - Strand specificity  
- Ensure all samples pass QC before moving forward  

**Output:** QC reports for aligned reads
<img width="1257" height="311" alt="Screenshot 2025-08-23 202801" src="https://github.com/user-attachments/assets/b78ccb89-5668-41a8-a4f4-a895bc959c74" />

---

### 8️⃣ Differential Expression Analysis (in R)
- Import counts matrix into **R (DESeq2 workflow)**  
- **Normalise** read counts (to remove library size bias)  
- Perform **differential expression testing**:  
  - Normoxia vs Hypoxia (for each cell line separately)  
- Visualisation and **pathway enrichment** 

**Output:**  
- List of significantly up/downregulated genes  
- plots (PCA plots,Volcano plots,Heatmaps) 
- Pathway enrichment analysis  

------

## 🗂️ Recommended Folder Layout
```
Bulk_RNA_Seq_Analysis/
├─ fastq/ # FASTQs (raw + merged/renamed)
├─ fastqc_results/ # FastQC HTMLs
├─ multiqc_report/ # MultiQC summary
├─ alignedreads/ # BAM + BAI
├─ quants/ # featureCounts outputs
├─ rnaseq_qc_results/ # Qualimap outputs
```

## 🧭 End-to-End Workflow — Commands + Explanations


### Step 0 : Creating environment and downloading SRA Toolkit

It is recommended to work in a **clean environment** so that all tools are properly installed and version conflicts are avoided.  
You can use **Conda (preferred)** or your system package manager (APT) to install the required tools.

#### 📦 Using Conda (recommended)
```bash
# (optional) create and activate a dedicated environment
conda create -n rnaseq python=3.9 -y
conda activate rnaseq

# go to your project directory
mkdir Bulk_RNA_Seq_Analysis
cd Bulk_RNA_Seq_Analysis

# install SRA Toolkit
conda install -c conda-forge -c bioconda sra-tools=3.2.1

# verify installation
which prefetch
which fastq-dump
````

#### 🖥️ Using APT (if Conda is not available)

```bash
sudo apt update
sudo apt install sra-toolkit
```

#### ℹ️ Why this step?

* `prefetch` → downloads `.sra` files from NCBI
* `fastq-dump` → converts `.sra` → `.fastq.gz` for downstream analysis



### Step 1: Download and Convert SRA Files to FASTQ

In this step, we will download sequencing data from the NCBI SRA database and convert them to FASTQ format for downstream analysis.

### 1.1 Manual Download Example (Single SRR)
To download and convert a single SRA run (e.g., `SRR7179504`), use the following commands:

```bash
# Download the SRA file
prefetch SRR7179504

# Convert SRA to FASTQ format
fastq-dump --outdir fastq --gzip --skip-technical --readids \
--read-filter pass --dumpbase --split-3 --clip \
SRR7179504/SRR7179504.sra
````

> **Note:** This command generates compressed FASTQ files in the `fastq/` directory and ensures only high-quality reads are retained.

### 1.2 Automate for Multiple SRR IDs

Processing multiple runs manually is time-consuming. Use a Python script to automate the download and conversion for all SRR IDs.

1. **Check Python installation**:

```bash
python3 --version || (sudo apt update && sudo apt install python3 -y)
```

2. **Run the automation script**:

```bash
python3 fastq_download.py
```

> **Why automation matters:** Looping through all SRR IDs saves time and ensures that all sequencing runs (e.g., 20 runs) are processed efficiently without manual intervention.



## Step 2: Raw Read Quality Control (FastQC → MultiQC)

Before alignment, it’s critical to check the quality of raw reads to identify issues such as adapter contamination, low-quality bases, or GC bias.

### 2.1 Install Tools
I recommend installing via **Conda** (preferred), but APT or pip can also be used.

**Option 1: Conda (Recommended)**
```bash
conda install -c bioconda fastqc
conda install -c bioconda multiqc
````

**Option 2: APT (System-wide)**

```bash
sudo apt update && sudo apt install fastqc multiqc
```

**Option 3: pip (Latest MultiQC)**

```bash
pip install multiqc
```


### 2.2 Run Quality Control

1. **Run FastQC on all FASTQ files**

```bash
mkdir -p fastqc_results
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```

2. **Aggregate results with MultiQC**

```bash
multiqc fastqc_results/ -o multiqc_report/
```

### Why this step?

- ✅ Confirms sequencing quality
- ✅ Detects adapter contamination
- ✅ Identifies GC bias and other systematic issues
- ✅ Provides both per-sample (FastQC) and summarized (MultiQC) reports



## Step 3: (Optional) Read Trimming with Trimmomatic

Trimming can remove low-quality bases or adapter contamination. In this project, trimming was tested but since the aligner used here can better handle this things,untrimmed reads were used for downstream analysis.  


### 3.1 Install Trimmomatic
Trimmomatic requires Java. Install Java first, then Trimmomatic:

```bash
# Check Java or install if missing
java -version || (sudo apt update && sudo apt install default-jre -y)

# Install Trimmomatic via Conda
conda install -c bioconda trimmomatic
````

### 3.2 Example: Trim a Single File

Here I demonstrate trimming on one FASTQ file (`SRR7179504_pass.fastq.gz`):

```bash
trimmomatic SE -threads 4 -phred33 \
  fastq/SRR7179504_pass.fastq.gz \
  fastq/SRR7179504_trimmed.fastq.gz \
  TRAILING:10
```

### 3.3 Re-check QC

After trimming, re-run FastQC on the trimmed file:

```bash
fastqc fastq/SRR7179504_trimmed.fastq.gz
```

### Why this step?

- 🔹 Removes low-quality bases from read ends
- 🔹 Can reduce adapter contamination
- 🔹 Ensures cleaner input for alignment

> **Note:** In our dataset, trimming did not materially improve QC



## Step 4: Merge Technical Runs & Rename Samples

Concatenate LNCaP technical runs → replicates and rename PC3 runs:
```
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
# After verifying the new files exist and sizes look correct:
rm -rf SRR*
```



### 📑 Merging 20 SRR Runs → 8 Final FASTQ Samples gives below result

| Cell line | Condition | Replicate | SRR runs merged            | Final FASTQ                     |
|-----------|-----------|-----------|----------------------------|---------------------------------|
| LNCaP     | Normoxia  | Rep1      | 504, 505, 506, 507         | LNCAP_Normoxia_S1.fastq.gz      |
| LNCaP     | Normoxia  | Rep2      | 508, 509, 510, 511         | LNCAP_Normoxia_S2.fastq.gz      |
| LNCaP     | Hypoxia   | Rep1      | 520, 521, 522, 523         | LNCAP_Hypoxia_S1.fastq.gz       |
| LNCaP     | Hypoxia   | Rep2      | 524, 525, 526, 527         | LNCAP_Hypoxia_S2.fastq.gz       |
| PC3       | Normoxia  | Rep1      | 536                        | PC3_Normoxia_S1.fastq.gz        |
| PC3       | Normoxia  | Rep2      | 537                        | PC3_Normoxia_S2.fastq.gz        |
| PC3       | Hypoxia   | Rep1      | 540                        | PC3_Hypoxia_S1.fastq.gz         |
| PC3       | Hypoxia   | Rep2      | 541                        | PC3_Hypoxia_S2.fastq.gz         |


### Why merging?
- ✅ Ensures that all reads from the same biological replicate are analyzed together  
- ✅ Prevents splitting of replicate data across multiple files  
- ✅ Simplifies downstream processing (alignment, quantification, etc.)  

After merging, we also **renamed the files** for clarity and consistency.




## Step 5: Reference Genome & Annotation

To perform alignment and gene quantification, we need a **reference genome index** (for HISAT2) and a **gene annotation file** (GTF).


### 5.1 Download Prebuilt HISAT2 Index

HISAT2 requires a prebuilt index of the reference genome. Download and extract:

```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
````


### 5.2 Install Alignment & BAM Tools

Install HISAT2 (aligner) and Samtools (for handling BAM files):

```bash
sudo apt install hisat2
sudo apt install samtools
```


### 5.3 Download Gene Annotation (Ensembl GTF)

Download and extract the **Ensembl GTF annotation file**:

```bash
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
gunzip Homo_sapiens.GRCh38.114.gtf.gz
```

---

### Why this step?

* ✅ **HISAT2** requires a reference genome index to align reads
* ✅ **featureCounts** and other quantification tools need a gene annotation (GTF) to assign reads to genes



## Step 6: Alignment (HISAT2 → Samtools)

In this step, we align the quality-checked FASTQ reads to the reference genome using **HISAT2**, and then process the results into sorted and indexed BAM files with **Samtools**.


### 6.1 Single-Sample Example
The following command aligns one FASTQ file (`LNCAP_Hypoxia_S1.fastq.gz`) and outputs a sorted and indexed BAM file:

```bash
hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1.fastq.gz | \
  samtools sort -o alignedreads/LNCAP_Hypoxia_S1.bam

samtools index alignedreads/LNCAP_Hypoxia_S1.bam
````


### 6.2 Automating the Alignment (Multiple Samples)

For multiple samples, use the provided helper script:

```bash
chmod +x hisat2_alignment1.sh
./hisat2_alignment1.sh
```

This script automates alignment for 7 samples; one additional sample was processed manually.


### Notes

⚠️ If `samtools index` fails when piped, run it separately after sorting completes.
✅ Result: Sorted and indexed BAM files stored in the `alignedreads/` directory, ready for downstream analysis.


## Step 7: Quantification (featureCounts)

After alignment, we need to count how many reads map to each gene.  
We use **featureCounts** (part of the Subread package) to generate a **gene-by-sample count matrix**, which is essential for differential expression (DE) analysis.


### 7.1 Install Subread
Install Subread and check that featureCounts is available:

```bash
sudo apt-get update
sudo apt-get install subread
featureCounts -v
````


### 7.2 Prepare Output Folder

Create a folder for storing quantification results:

```bash
mkdir -p quants
```

### 7.3 Count Reads (Single BAM Example)

Run featureCounts on one BAM file (`tmp.bam`) with strand-specific mode (`-S 2`):

```bash
featureCounts -S 2 -a Homo_sapiens.GRCh38.114.gtf \
  -o quants/featurecounts.txt tmp.bam
```

### 7.4 Automate for All BAM Files

To automate the process of quantification, use the provided script:

```bash
chmod +x featurecounts.sh
./featurecounts.sh
```

### Why this step?

- ✅ Converts aligned reads into **per-gene read counts**
- ✅ Produces a **gene × sample count matrix**
- ✅ Provides the featurecounts which can be used to create count matrix ,which is the input required for downstream **differential expression analysis (DEA)**



## Step 8: QC of Aligned Reads (Qualimap RNA-seq)

After alignment and quantification, it is important to assess the quality of aligned reads.  
We use **Qualimap** to evaluate mapping quality, coverage, and biases.


### 8.1 Install Qualimap
Download manually or install via Conda:

**Option 1: Manual Download**
```bash
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.2.zip
unzip qualimap_v2.2.2.zip
cd qualimap_v2.2.2
chmod +x qualimap
````

**Option 2: Conda Install**

```bash
conda install bioconda::qualimap
```

### 8.2 Per-Sample QC Example

Run QC on one BAM file (`LNCAP_Hypoxia_S1.bam`):

```bash
./qualimap_v2.3/qualimap rnaseq \
  -bam alignedreads/LNCAP_Hypoxia_S1.bam \
  -gtf Homo_sapiens.GRCh38.114.gtf \
  -outdir rnaseq_qc_results \
  --java-mem-size=8G
```

### 8.3 Automate QC for All Samples

For batch QC, use the helper script:

```bash
chmod +x run_qualimap.sh
./run_qualimap.sh
```

### Why this step?

- ✅ Validates **mapping percentage**
- ✅ Checks **coverage distribution**
- ✅ Detects **gene-body bias**
- ✅ Confirms **library strandedness**



## Step 9: Build the Counts Matrix

After quantification, each sample has its own featureCounts output.  
In this step, we merge all per-sample counts into a **single counts matrix** (genes × samples).


### 9.1 Merge Counts
Use the provided Jupyter Notebook to combine featureCounts outputs:

```python
countsmatrix.ipynb
```

This notebook consolidates all per-sample counts into one matrix where:

**Rows = genes**
**Columns = samples**

---

### 9.2 DESeq2_analysis using R/Rstudio

Once the counts matrix is built, proceed to **R / RStudio**, and folllow up the .Rmd file provided for downstream analysis with **DESeq2**:

- ✅ Normalization
- ✅ Differentially expressed genes (DEGs)
- ✅ PCA plots, Volcano plots, Heatmap
- ✅ Pathway ebnrichment analysis


## 🛠️ Tools & Software Used

- SRA Toolkit – Download sequencing data
- FastQC + MultiQC – Quality control
- Trimmomatic / fastp – Read trimming (optional)
- HISAT2 – Alignment to reference genome
- Samtools – BAM sorting and indexing
- featureCounts – Read quantification
- Qualimap – QC for aligned reads
- R (DESeq2, ggplot2, pheatmap) – Statistical analysis ,and visualization

## 📚 Learning Outcomes

- Learned how to process raw RNA-Seq data step by step
- Understood the importance of QC at every stage
- Experienced handling large datasets and merging technical replicates
- Performed differential expression analysis using DESeq2
- Gained insights into how hypoxia impacts prostate cancer cell lines

### 👩‍🏫 Mentor<br>
This work was carried out under the guidance of "Smriti Arora"

### For questions or collaborations:<br>

**📬 Contact**<br>
**Author:** Gayatri Sunil Samant<br>
**Email:** gayatrisamant05@gmail.com






