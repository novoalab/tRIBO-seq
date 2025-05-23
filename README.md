# Ribo-embedded Nano-tRNAseq Analysis Pipeline

This repository contains scripts and workflows for analyzing ribo-embedded Nano-tRNAseq data. 
We apply both filtered and non-filtered analysis methods to assess reproducibility, differential expression, and differential modification.

## Table of Contents

- [Experimental conditions](#Experimental-conditions)
- [Analysis Workflow](#AnalysisWorkflow)
	- [1. Reproducibility between libraries and replicates](#1.-Reproducibility-between-libraries-and-replicates)
	- [2. Principal Component Analysis](#2.-Principal-Component-Analysis)
	- [3. Differential Expression Analysis](#3.-Differential-Expression-Analysis)
 - 	[4. Differential Modification Analysis](#3.-Differential-Modification-Analysis)
- [Expected Output](#Expected-output)
- [Pre-filtering of the BAM files](#Pre-filtering-of-the-BAM-files)
- [Dependencies and Versions](#Dependencies-and-Versions)
- [Citation](#Citation) 
- [Contact](#Contact)

### Structure of the directory  (this I would remove later on) 

@Hasan, please edit this as needed!)

📂 tRNA/
 ├── 📂 scripts/               # Analysis scripts
 │    ├── filter_script.py     # Script for filtering BAM files
 │    ├── analysis_pipeline.R  # Main analysis script
 │    ├── scatterplot.R        # Library reproducibility plot
 │    ├── corrplot.R           # Correlation heatmap script
 │    ├── differential_abundance_and_pca.R # Differential expression and Principal component analysis in codon level
 │    ├── differential_abundance_and_pca_aminoacid.R # Differential expression and Principal component analysis in amino acid level
 │    ├── differential_modification.R # Differential modification analysis
 ├── 📂 data/                  # Processed data and demo files
 │    ├── demo_data.bam        # A small subset BAM file for testing
 │    ├── metadata.tsv         # Metadata for samples
 ├── 📂 results/               # Output figures and processed results
 │    ├── 📂 filtered/         # Results using filtering
 │    ├── 📂 non_filtered/     # Results without filtering
 ├── README.md                 # This README file
 ├── .gitignore                # Ignores large BAM files and unnecessary files

Notes: BAM files should NOT be committed to avoid bloating the repository. Instead, use .tsv output files for processed data. Consider including a small subset BAM file (demo_data.bam) as an example dataset for testing the filtering script.


## Experimental Conditions

Note: Each condition has a designated prefix (Again, @Hasan, please edit this as needed!)

* **Arsenite** treatment: 'Arsenite_'
* **Methionine** starvation: 'MetStarve6h_' or 'MetStarve16h_'
* **Arginine** starvation: 'ArgStarve_'
* **Leucine** starvation: 'LeuStarve_'

In Differential Abundance Analysis four suffix assigned
* **CONTROL** : Ribo-Embedded tRNAs vs Total tRNAs in Control Samples
* **(TREATMENT)** : This suffix changes with respect to non-control condition and it reflects Ribo-Embedded tRNAs vs Total tRNAs in non-Control Samples
* **RE** : Non-Control vs Control samples in Ribo-Embedded tRNAs
* **TOTAL** : Non-Control vs Control samples in Total tRNAs
  
For each condition, we perform analysis using both filtered and non-filtered data. Please ensure that the filtering script is included in scripts/filter_script.py.

## Analysis Workflow

For both filtered and non-filtered datasets, we conduct the following analyses:

```
Rscript script_name.R --help (Helper function for all scripts below)
```

### 1. Reproducibility Between libraries and replicates

#### 1.1. Scatterplots comparing read counts between replicates

```
Rscript scatterplot.R --counts /path/to/counts.tsv --prefix YourExperiment_ --type RE or Total --output /path/to/out
```

#### 1.2. Correlation plots across all libraries in an experiment

```
Rscript corrplot.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --output /path/to/out
```

### 3. Differential Expression and Principal Component Analysis

#### 3.1 Codon Level Analysis
* Using the models defined (~SeqType, ~Trt_Total, ~Trt_Ribo, ~SeqType*Trt)
* Output: Volcano plots of differentially expressed tRNAs
* 
```
Rscript differential_abundance_and_pca.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --output /path/to/out
```
#### 3.2 Amino Acid Level Analysis
* Using the models defined (~SeqType, ~Trt_Total, ~Trt_Ribo, ~SeqType*Trt)
* Output: Volcano plots of differentially expressed tRNAs
* 
```
Rscript differential_abundance_and_pca_aminoacid.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --output /path/to/out
```

### 4. Differential Modification Analysis

#### tRNA heatmap
```
Rscript differential_modification.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --fasta /path/to/tRNAs/alignment/database.aln.fa --output /path/to/out
```

## Expected Output

Results of the pipeline using the demo data can be found in the [results](https://github.com/novoalab/Ribo-embedded/results) folder. 

[Comments - remove this later 
* Figures: Saved as .pdf for easy viewing and publication use
* Result files: Stored as .tsv for further inspection and replotting ]

## Pre-filtering of the BAM files

```
python3 filter_script.py -i /path/to/input_list.tsv -o /path/to/out -c
```

## Dependencies and versions 

Software | Version
--- | ---
R | 4.4.2
python | 3.12.4
pandas | 2.2.3
pysam | 0.22.1

## Citation

If you find this work useful, please cite: 

Hasan Yilmaz*, Mie Monti*, Alessia Del Piano, Michele Arnoldi, Isabelle Bonomo, Laia Llovera, Eva Maria Novoa# and Massimiliano Clamer#. Revealing Stress-Specific Differences Between Riboembedded and Total tRNA Populations (in preparation).

## Contact

If you have any issues running this code, please go first over previous [issues](https://github.com/novoalab/Ribo-embedded/issues). If you still can't figure it out based on the prior responses/issues raised, please open a new issue. Thanks!   


 
