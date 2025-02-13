# tRNA Analysis Pipeline

This repository contains scripts and workflows for harmonizing tRNA analysis for the paper. The data is organized by experimental conditions, and we apply both filtered and non-filtered analysis methods to assess reproducibility, differential expression, and differential modification.

# Structure of the directory 

@Hasan, please edit this as needed!)

📂 tRNA/
 ├── 📂 scripts/               # Analysis scripts
 │    ├── filter_script.py     # Script for filtering BAM files
 │    ├── analysis_pipeline.R  # Main analysis script
 │    ├── scatterplot.R        # Library reproducibility plot
 │    ├── corrplot.R           # Correlation heatmap script
 │    ├── differential_abundance_and_pca.R # Differential expression and PC analysis
 │    ├── differential_modification.R # Differential modification analysis
 ├── 📂 data/                  # Processed data and demo files
 │    ├── demo_data.bam        # A small subset BAM file for testing
 │    ├── metadata.tsv         # Metadata for samples
 ├── 📂 results/               # Output figures and processed results
 │    ├── 📂 filtered/         # Results using filtering
 │    ├── 📂 non_filtered/     # Results without filtering
 ├── README.md                 # This README file
 ├── .gitignore                # Ignores large BAM files and unnecessary files

#Experimental Conditions

Each condition has a designated prefix:
Again, @Hasan, please edit this as needed!

Arsenite treatment: 'Arsenite_'
Methionine starvation: 'MetStarve6h_' or 'MetStarve16h_'
Arginine starvation: 'ArgStarve_'
Leucine starvation: 'LeuStarve_'
For each condition, we perform analysis using both filtered and non-filtered data. Please ensure that the filtering script is included in scripts/filter_script.py.

#Analysis Workflow

For both filtered and non-filtered datasets, we conduct the following analyses:

1. Reproducibility Between Libraries
* Scatterplots comparing read counts between libraries
* Correlation plots across all libraries in an experiment
2. Differential Expression and Principal Component Analysis (PCA)
* Using the models defined (~SeqType, ~Trt_Total, ~Trt_Ribo, ~SeqType*Trt)
* Output: Volcano plots of differentially expressed tRNAs
* Output: PCA visualization of sample clustering
3. Differential Modification Analysis
* Analyzing tRNA modification levels
* Output: Heatmaps of differential modifications
4. Output Format
* Figures: Saved as .pdf for easy viewing and publication use
* Result files: Stored as .tsv for further inspection and replotting

#Notes

BAM files should NOT be committed to avoid bloating the repository. Instead, use .tsv output files for processed data.
Consider including a small subset BAM file (demo_data.bam) as an example dataset for testing the filtering script.



