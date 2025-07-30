# RiboNano-tRNAseq

This repository contains scripts and workflows for analyzing RiboNano-tRNAseq data.  
We apply both filtered and non-filtered analysis methods to assess reproducibility, differential expression, differential modification, and fragmentation of tRNAs.

If you find this work useful, please cite:

Hasan Yilmaz*, Mie Monti*, Alessia Del Piano, Michele Arnoldi, Isabelle Bonomo, Laia Llovera, Eva Maria Novoa#, and Massimiliano Clamer#.
Selective profiling of translationally active tRNAs and their dynamics under stress _(manuscript in preparation)_
[DOI to be shared once available]

---

## Table of Contents

- [Parameter Setting](#parameter-setting)
- [Analysis Workflow](#analysis-workflow)  
  - [1. Reproducibility between libraries and replicates](#1-reproducibility-between-libraries-and-replicates)  
  - [2. Differential expression analysis](#2-differential-expression-analysis)  
  - [3. Differential modification analysis](#3-differential-modification-analysis)  
  - [4. Differential fragmentation analysis](#4-differential-fragmentation-analysis)  
- [Expected Output](#expected-output)
- [Dependencies and Versions](#dependencies-and-versions)
- [Citation](#citation)
- [Contact](#contact)

---

## Parameter setting

For your treatment of interest, please define a prefix.  
**Example**: For Arsenite treatment use the prefix `'Arsenite_'`.

In the `Differential Abundance Analysis` module, four suffixes are used:

- **CONTROL**: Ribo-tRNAs vs total tRNAs in control samples  
- **(TREATMENT)**: Changes based on the prefix (e.g., Ribo-tRNAs vs total tRNAs in treated samples)  
- **RE**: Treated vs control (Ribo-tRNAs)  
- **TOTAL**: Treated vs control (Total tRNAs)  

For each condition, the analysis is performed using both filtered and non-filtered data.

The default filtering parameters are:  
@Hasan to complete here

---

## Analysis Workflow

```
Rscript script_name.R --help (Helper function for all scripts below)
```

The command completed for filtering is:

```
python3 filter_script.py -i /path/to/input_list.tsv -o /path/to/out -c
```

For both filtered and non-filtered datasets, we conduct the following analyses:

### 1. Reproducibility between libraries and replicates

#### 1.1. Scatterplots between reps

Generates scatterplots comparing read counts between biological replicates.

```
Rscript scatterplot.R --counts /path/to/counts.tsv --prefix YourExperiment_ --type RE or Total --output /path/to/out
```

#### 1.2. Correlation across all libraries in experiment

Generates a correlation matrix across all samples in an experiment.

```
Rscript corrplot.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --output /path/to/out
```

### 2. Differential expression analysis

Differential expression is performed using the DESeq2 framework and outputs include PCA plots, volcano plots, and significance tables.
Analysis is conducted at both codon and amino acid levels.

#### 2.1 Codon level analysis

Model terms include: ~SeqType, ~Trt_Total, ~Trt_Ribo, ~SeqType*Trt

```
Rscript differential_abundance_and_pca.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --output /path/to/out
```

#### 2.2 Amino Acid Level Analysis

Same model structure as codon-level analysis. This script aggregates counts at the amino acid level.

```
Rscript differential_abundance_and_pca_aminoacid.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --output /path/to/out
```

### 3. Differential Modification Analysis

#### 3.1 tRNA heatmap

Visualizes modifications across tRNAs

```
Rscript differential_modification.R --total /path/to/totaltRNAs/counts.tsv --re /path/to/riboembeddedtRNAs/counts.tsv --prefix YourExperiment_ --sample_list /path/to/sample_list.tsv --fasta /path/to/tRNAs/alignment/database.aln.fa --output /path/to/out
```

### 4. Differential Fragmentation Analysis

#### 4.1 Bulk fragmentation

```
Rscript TBC
```

#### 4.2 Site-specific fragmentatio
```
python3 bam2ends.py -i /path/to/input.bam -o /path/to/output/bam2ends.tsv
```
This returns a table where each row is XYZ @Mie to complete. 

## Expected Output

Results of the pipeline using the demo data can be found in the [results] folder. 

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

@Hasan there are a lot of libraries missing here no?

## Citation

If you find this work useful, please cite: 

Hasan Yilmaz*, Mie Monti*, Alessia Del Piano, Michele Arnoldi, Isabelle Bonomo, Laia Llovera, Eva Maria Novoa# and Massimiliano Clamer#. Revealing Stress-Specific Differences Between Riboembedded and Total tRNA Populations (in preparation).

## Contact

If you have any issues running this code, please go first over previous [issues]. If you still can't figure it out based on the prior responses/issues raised, please open a new issue. Thanks!   
