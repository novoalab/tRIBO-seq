# -----------------------------------------------------------------------------
# Title:         Bio-type Plot
# Author:        Hasan YILMAZ
# Date:          25/03/25
# Location:      Trento, Italy
# Description:   The script is to create barplot plot of biotypes on 
#nano-tRNASeq data.
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, ggplot2, ggpubr, tibble, tidyr, dplyr, glue
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Quality Control Check")
parser$add_argument("--total", type = "character", help = "Counts of Total tRNASeq file path")
parser$add_argument("--re", type = "character", help = "Counts of Ribo-Embedded tRNASeq file path")
parser$add_argument("--prefix", type = "character", help = "Analysis Prefix")
parser$add_argument("--output", type = "character", help = "output directory")
args <- parser$parse_args()


# -------------------------- Helper Functions ------------------------------ #
# 1. Collapse iso-acceptors to single tRNA
collapsetRNA <- function(data){
  
  data <- data[!grepl('oligo|unique|antisense', data$reference),]
  
  for (i in 1:dim(data)[1]){
    splitted <- strsplit(data$reference[i], '-', fixed = T)[[1]]
    if (length(splitted)>1){
      data$reference[i] <- paste0(splitted[1], '-', splitted[2])
    }
    else{
      data$reference[i] <- splitted[1]
    }
  }
  data <- aggregate(. ~ reference, data, sum)
  
  return(data)
}

# 2. Convert wide data to long data
convert_df <- function(data){
  # Classify RNA types
  data <- data %>%
    mutate(Type = case_when(
      grepl("^mt", reference) ~ "mt-tRNA",
      grepl("^RDN", reference) ~ "rRNA",
      TRUE ~ glue("tRNA"))
    )
  # Reshape to long format
  data <- data %>%
    pivot_longer(cols = -c(reference, Type),
                 names_to = "Sample",
                 values_to = "Count")
  
  summary_data <- data %>%
    group_by(Sample, Type) %>%
    summarise(TotalCount = sum(Count), .groups = "drop")
  
  sample_levels <- summary_data %>%
    distinct(Sample) %>%
    arrange(!grepl("^Control", Sample)) %>%
    pull(Sample)
  
  summary_data <- summary_data %>%
    mutate(Sample = factor(Sample, levels = sample_levels))
  
  return(summary_data)
}

# 3. Read and merge the data
process <- function(arguments){
  message('Processing Started...')
  data_re <- read.table(arguments$re, sep='\t', header=T)
  data_re <- collapsetRNA(data_re)
  data_re <- convert_df(data_re)
  
  plotting(data_re, arguments, 'RE')
  
  data_total <- read.table(arguments$total, sep='\t', header=T)
  data_total <- collapsetRNA(data_total)
  data_total <- convert_df(data_total)
  
  plotting(data_total, arguments, 'Total')
}


plotting <- function(data, arguments, exp){
  plot <- ggplot(data, aes(x = Sample, y = TotalCount, fill = Type)) +
    geom_bar(stat = "identity", position='fill') +
    labs(x = "Sample",
         y = "Count",
         fill = "RNA Type") +
    theme_pubr(border = T, legend = 'right')
  ggsave(filename = glue('{arguments$prefix}{exp}_RNA_types.pdf'),
         device = 'pdf', plot = plot, path = arguments$output,
         height = 8, width = 10)
}

process(arguments = args)
