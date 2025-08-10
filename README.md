# tRIBO-seq

This repository contains scripts and workflows for analyzing tRIBO-seq datasets. tRBIBO-seq is a method to selectively profile ribosome-embedded native tRNAs using nanopore sequencing technologies. 

- **What's included**

This repository includes scripts to perform the following steps: 
* quality metrics and reproducibility
* read filtering
* differential expression
* differential modification
* fragmentation of tRNAs


- **If you find this work useful**

Please cite: Hasan Yilmaz*, Mie Monti*, Alessia Del Piano, Michele Arnoldi, Isabelle Bonomo, Laia Llovera, Massimiliano Clamer# and Eva Maria Novoa#. **Selective profiling of translationally active tRNAs and their dynamics under stress** _(manuscript in preparation)_

---

## Table of Contents

- [Introduction](#introduction)
- [Analysis Workflow](#analysis-workflow)  
  - [1. Processing alignment files](#1-processing-alignment-files)  
  - [2. Reproducibility between libraries and replicates](#2-reproducibility-between-libraries-and-replicates)  
  - [3. Differential expression analysis](#3-differential-expression-analysis)  
  - [4. Differential modification analysis](#4-differential-modification-analysis)  
  - [5. Differential fragmentation analysis](#5-differential-fragmentation-analysis)  
- [Expected Output](#expected-output)
- [Dependencies and Versions](#dependencies-and-versions)
- [Citation](#citation)
- [Contact](#contact)

---

## Introduction

This repository contains all the necessary scripts to reproduce the analyses presented in the *tRIBO-seq* manuscript. A dummy dataset is provided in the `Example_data` directory to illustrate the workflow and enable testing of the pipeline.

### File Structure

- Files prefixed with `ribo_nanotrna_` correspond to ribosome-associated tRNA (ribo-tRNA) counts.
- Files prefixed with `nanotrna_` correspond to total tRNA counts.
- BAM file inputs are required for this analysis pipeline.

### Prefix Configuration

For each experimental treatment, define a **prefix** that will be used throughout the pipeline.

Example: For Arsenite treatment, set the prefix to `'Arsenite_'`.

This prefix will guide the naming of output files and the interpretation of comparisons.

### Differential analysis

The `Differential_Abundance_Analysis` module includes comparisons across four categories:

- **CONTROL**: Ribo-tRNAs vs total tRNAs in *control* samples  
- **(TREATMENT)**: Ribo-tRNAs vs total tRNAs in *treated* samples (based on the defined prefix)  
- **RE**: Ribo-tRNAs in *treated* vs *control* conditions  
- **TOTAL**: Total tRNAs in *treated* vs *control* conditions  

Each analysis is conducted using both **filtered** and **non-filtered** datasets.

### Filtering Parameters (Default)

The filtering criteria applied to reads are:

- Mapping Quality (MapQ) ≥ 10  
- 3′ Splint Adapter must be present (see `three_prime_adapter_check` in `nanoCount.py`)  
- Full-Length Reads only (`-c` flag in `nanoCount.py`)


---

## Analysis Workflow

```
Rscript script_name.R --help (Helper function for all scripts below)
python3 nanoCount.python --help
```

### 1. Processing alignment files

```
python3 nanoCount.py -i /path/to/input_list.txt -o /path/to/output/folder -f /path/to/reference.fasta \\
    -a /path/to/reference.aln.fa [optional] -c [optional] -nc
```
Arguments -c and -nc are switch options to get only **full_length** or **fragments**. If none of them used pipeline
will result **all** reads.

**Example Data**

Nano-tRNAseq All Reads (where All = full length *and* fragmented reads)
```
python3 nanoCount.py -i Example/Data/nano_trnaseq_input_list.txt -o Example_Data/Result/nanotrnaseq/All -f Example_Data/reference/hg38.fa \\
    -a Example_Data/reference/hg38.aln.fa
```
Nano-tRNAseq only Full Length Reads
```
python3 nanoCount.py -i Example/Data/nano_trnaseq_input_list.txt -o Example_Data/Result/nanotrnaseq/Full_Length -f Example_Data/reference/hg38.fa \\
    -a Example_Data/reference/hg38.aln.fa -c
```
Nano-tRNAseq only Fragmented Reads
```
python3 nanoCount.py -i Example/Data/nano_trnaseq_input_list.txt -o Example_Data/Result/nanotrnaseq/Fragment -f Example_Data/reference/hg38.fa \\
    -a Example_Data/reference/hg38.aln.fa -nc
```
Ribo Nano-tRNAseq All Reads
```
python3 nanoCount.py -i Example/Data/ribo_nano_trnaseq_input_list.txt -o Example_Data/Result/ribo_nanotrnaseq/All -f Example_Data/reference/hg38.fa \\
    -a Example_Data/reference/hg38.aln.fa
```
Ribo Nano-tRNAseq only Full Length Reads
```
python3 nanoCount.py -i Example/Data/ribo_nano_trnaseq_input_list.txt -o Example_Data/Result/ribo_nanotrnaseq/Full_Length -f Example_Data/reference/hg38.fa \\
    -a Example_Data/reference/hg38.aln.fa -c
```
Ribo Nano-tRNAseq only Fragmented Reads;
```
python3 nanoCount.py -i Example/Data/ribo_nano_trnaseq_input_list.txt -o Example_Data/Result/ribo_nanotrnaseq/Fragment -f Example_Data/reference/hg38.fa \\
    -a Example_Data/reference/hg38.aln.fa -nc
```

For both filtered and non-filtered datasets, we conduct the following analyses:

### 2. Reproducibility between libraries and replicates

#### 2.1. Proportion of Aligments between treatments

Generates density plots comparing read alignment proportion (read alignment area length / reference length) between control and treatment samples.

```
Rscript proportion_plot.R --input /path/to/coverage_proportion.tsv --prefix YourExperiment_ --suffix RE or Total --output path/to/out
```

**Example Data**
After this point outputs of nanoCount.py will be used. Please follow read.me

```
Rscript proportion_plot.R --input Example_Data/Result/nanotrnaseq/All/coverage_proportion.tsv --prefix Example_ --suffix Total --output Example_Data/Result/QC
Rscript proportion_plot.R --input Example_Data/Result/ribo_nanotrnaseq/All/coverage_proportion.tsv --prefix Example_ --suffix RE --output Example_Data/Result/QC
```

#### 2.2. Scatterplots between reps

Generates scatterplots comparing read counts between biological replicates.

```
Rscript scatterplot.R --ribo path/to/counts.tsv --total path/to/counts.tsv --prefix YourExperiment_ --samplefile /path/to/samplefile.tsv --control YourControlCondition --treatment YourTreatmentCondition --axes AxisLimits --output /path/to/out
```

**Example Data**
```
Rscript scatterplot.R --ribo Example_Data/Result/ribo_nanotrnaseq/All/counts.tsv --total Example_Data/nanotrnaseq/All/counts.tsv --prefix Example_ --samplefile Example_Data/samplefile.tsv --control Control --treatment Treatment --axes 140000 --output Example_Data/Result/QC
```

#### 2.3. Correlation across all libraries in experiment

Generates a correlation matrix across all samples in an experiment.

```
Rscript corrplot.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --output /path/to/out
```

**Example Data**
```
Rscript corrplot.R --total Example_Data/Result/nanotrnaseq/All/counts.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/counts.tsv --prefix Example_ --output Example_Data/Result/QC
```

#### 2.4. Barplots of tRNAs/mt-tRNAs and rRNAs in experiments

Generates a barplot for experiment.

```
Rscript biotype_plot.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --output /path/to/out
```

**Example Data**
```
Rscript biotype_plot.R --total Example_Data/Result/nanotrnaseq/All/counts.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/counts.tsv --prefix Example_ --output Example_Data/Result/QC
```

#### 2.5. MapQ density plots across samples

Generates a density plot of mapq of each reads.

```
Rscript mapq_plot.R --total /path/to/totaltRNAs/quality.tsv --re /path/to/riboembeddedtRNAs/quality.tsv --prefix YourExperiment_ --output /path/to/out
```

**Example Data**
```
Rscript mapq_plot.R --total Example_Data/Result/nanotrnaseq/All/quality.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/quality.tsv --prefix Example_ --output Example_Data/Result/QC
```

### 3. Differential expression analysis

Differential expression is performed using the DESeq2 framework and outputs include PCA plots, volcano plots, and significance tables.
Analysis is conducted at both codon and amino acid levels.

#### 3.1 Codon level analysis

Model terms include: ~SeqType, ~Trt_Total, ~Trt_Ribo, ~SeqType*Trt

```
Rscript differential_abundance_and_pca.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/samplefile.tsv --output /path/to/out
```

**Example Data**
```
Rscript differential_abundance_and_pca.R --total Example_Data/Result/nanotrnaseq/All/counts.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/counts.tsv --prefix Example_ --sample_list Example_Data/samplefile.tsv --output Example_Data/Result/Differential_Abundance/codon_level
```

#### 3.2 Amino Acid Level Analysis

Same model structure as codon-level analysis. This script aggregates counts at the amino acid level.

```
Rscript differential_abundance_and_pca_aminoacid.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/samplefile.tsv --output /path/to/out
```

**Example Data**
```
Rscript differential_abundance_and_pca_aminoacid.R --total Example_Data/Result/nanotrnaseq/All/counts.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/counts.tsv --prefix Example_ --sample_list Example_Data/samplefile.tsv --output Example_Data/Result/Differential_Abundance/amino_level
```

### 4. Differential Modification Analysis

#### 4.1 tRNA heatmap in Condition Level

Visualizes modifications across tRNAs in Condition level (Treatment vs. Control in Ribo-nanotRNAseq and nano-tRNAseq)

```
Rscript differential_modification.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --fasta /path/to/tRNAs/alignment/database.aln.fa --output /path/to/out
```

**Example Data**
```
Rscript differential_modification.R --total Example_Data/Result/nanotrnaseq/All/counts_single_nucleotide.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/counts_single_nucleotide.tsv --prefix Example_ --sample_list Example_Data/samplefile.tsv --fasta Example_Data/reference/hg38.aln.fa --output Example_Data/Result/Differential_Modification/
```

#### 4.2 tRNA heatmap in Type Level

Visualizes modifications across tRNAs in Type level (Ribo-nanotRNAseq vs. nano-tRNAseq in Control and Treatment)

```
Rscript differential_modification_exp.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --fasta /path/to/tRNAs/alignment/database.aln.fa --control YourControlCondition --treatment YourTreatmentCondition --output /path/to/out
```

**Example Data**
```
Rscript differential_modification_exp.R --total Example_Data/Result/nanotrnaseq/All/counts_single_nucleotide.tsv --re Example_Data/Result/ribo_nanotrnaseq/All/counts_single_nucleotide.tsv --prefix Example_ --sample_list Example_Data/samplefile.tsv --fasta Example_Data/reference/hg38.aln.fa --control Control --treatment Treatment --output Example_Data/Result/Differential_Modification/
```

### 5. Differential Fragmentation Analysis

#### 5.1 Bulk fragmentation

```
Rscript fragmentation.R --re_fl /path/to/ribo_trnas/Full_Length/counts.tsv --re_fr /path/to/ribo_trnas/Fragment/counts.tsv --total_fl /path/to/total_trnas/Full_Length/counts.tsv --total_fr /path/to/total_trnas/Fragment/counts.tsv --sample_list /path/to/sample_list.tsv --prefix YourExperiment --condition YourControlCondition<Control,Treatmen> --type YourControlType<Total,RE> --output /path/to/out
```

**Example Data**
```
Rscript fragmentation.R --re_fl Example_Data/Result/ribo_nanotrnaseq/Full_Length/counts.tsv --re_fr Example_Data/Result/ribo_nanotrnaseq/Fragment/counts.tsv --total_fl Example_Data/Result/nanotrnaseq/Full_Length/counts.tsv --total_fr Example_Data/Result/nanotrnaseq/Fragment/counts.tsv --sample_list Example_Data/samplefile_fragmentation.tsv --prefix Example_ --condition CONTROL <Control,Treatmen> --type TOTAL <Total,RE> --output /Example_Data/Result/Fragmentation
```


#### 5.2 Site-specific fragmentation
To generate site-specific fragmentation profiles from BAM files, run:
```
python3 bam2ends.py -i /path/to/input.bam -o /path/to/output/bam2ends.tsv
```

**Example Data**
```
python3 bam2ends.py -i Example/Data/*.bam -o Example_Data/per_site_fragmentation_bam2ends.tsv 
```

The resulting .tsv file contains per-site fragmentation data and can be used for downstream statistical testing. This output was tested using the **edgeR** framework ([edgeR User’s Guide – page 104](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)).

The following linear model was used:
`model.matrix(~ type + condition + type:condition, data = sample_info)`

Where:

- `type`: Indicates whether the read has a fragment end (`end`) or not (`no_end`)
- `condition`: Experimental condition (e.g., `control` or `treatment`)

## Dependencies and versions 

Software | Version
--- | ---
python | 3.12.4
pandas | 2.2.3
pysam | 0.22.1
numpy | 1.26.4
R | 4.4.2
pacman | 0.5.1
argsparse | 2.2.5
optparse | 1.7.5
BiocManager | 1.30.25
ComplexHeatmap | 2.22.0
circlize | 0.4.16
corrplot | 0.95
DESeq2 | 1.46.0
data.table | 1.16.4
dplyr | 1.1.4
EnhancedVolcano | 1.24.0
glue | 1.8.0
gsubfn | 0.7
ggplot2 | 3.5.1
GGally | 2.2.1
ggrastr | 1.0.2
ggridges | 0.5.6
ggstatsplot | 0.13.0
grDevices | 4.4.2
ggpubr | 0.6.0
patchwork | 1.3.0
rtracklayer | 1.66.0
rlang | 1.1.5
seqinr | 4.2-36
tibble | 3.2.1
tidyr | 1.3.1
wesanderson | 0.3.7


## Citation

Hasan Yilmaz*, Mie Monti*, Alessia Del Piano, Michele Arnoldi, Isabelle Bonomo, Laia Llovera, Massimiliano Clamer# and Eva Maria Novoa#. **Selective profiling of translationally active tRNAs and their dynamics under stress** _(manuscript in preparation)_[DOI to be shared once available]


## Contact

If you have any issues running this code, please go first over previous [issues](https://github.com/novoalab/tRIBO-seq/issues/). If you still can't figure it out based on the prior responses/issues raised, please open a new issue. Thanks!   
