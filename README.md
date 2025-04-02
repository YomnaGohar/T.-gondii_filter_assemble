# T.-gondii_filter_assemble: Snakemake Pipeline for Filtering and Assembly of Toxoplasma gondii Nanopore Data

This repository contains a **Snakemake pipeline** designed to process Oxford Nanopore sequencing data for *Toxoplasma gondii*. It performs:

### 🧬 Pipeline Overview

This pipeline processes Oxford Nanopore sequencing data for *Toxoplasma gondii* through the following stages:

1. **Host contamination filtering**  
   Each FASTQ file is filtered individually to remove reads originating from the host (human or mouse).  
   ⚠️ **Important:** Filtered FASTQ files are saved in the same folder as the original FASTQ file.  
   ➤ Ensure you have write permissions to these directories.

2. **Sample-level merging**  
   After filtering, all FASTQ files for each sample are concatenated into a single combined FASTQ file.

3. **Mitochondrial read filtering**  
   The concatenated FASTQ files are then filtered to retain only reads mapping to the *T. gondii* mitochondrial genome.

4. **Genome assembly**  
   The resulting mitochondrial reads are assembled into contigs to reconstruct the mitochondrial genome.

---

## 📁 Repository Structure

```bash
T.-gondii_filter_assemble/
├── workflow/             
│   ├── Snakefile    
├── rules/            
├── config/               
├── scripts/              
├── resources/            
├── .gitignore
├── README.md

### ⚙️ Dependencies

- Python ≥ 3.6
- snakemake 
- `minimap2`
- `samtools`
- `seqkit`
- `flye`
---

### 📝 Configuration

Before running the pipeline, you must edit the configuration file:

- A sample configuration template is available at: `config/config.yaml.sample`
- To use it, copy and rename it to `config/config.yaml`

### 🚀 Running the Pipeline

```bash
git clone https://github.com/YomnaGohar/T.-gondii_filter_assemble.git
cd T.-gondii_filter_assemble

To run the pipeline, you must navigate to the `workflow/` directory:

```bash
snakemake --cores 4 filter
```
This step generates a folder named <filtered_with_pipeline> inside the directory containing each input FASTQ file. The folder includes FASTQ files corresponding to human or mouse contamination, as well as the remaining Toxoplasma gondii reads. 

```bash
snakemake --cores 4 filter_mit
```
This step generates a file named <toxo.fastq> inside the <combined_with_pipeline> folder in the specified output directory. This file contains the combined, filtered FASTQ reads for each sample. Additionally, a FASTQ file named <non_mit_reads.fastq> will be generated in the same directory, containing T. gondii reads with mitochondrial sequences removed.

```bash
snakemake --cores 4 assemble
```
This step generates a folder named <assembly_without_mit_removal> that contains an assembly created from the <toxo.fastq> file located in the <combined_with_pipeline> folder within the specified output directory.

```bash
snakemake --cores 4 assemble_after_mit_removal
```
This step generates a folder named <assembly_with_mit_removal> that contains an assembly created from the <non_mit_reads.fastq> file located in the <combined_with_pipeline> folder within the specified output directory.






