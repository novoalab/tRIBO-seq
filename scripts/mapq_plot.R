# -----------------------------------------------------------------------------
# Title:         MAP-Q Plot
# Author:        Hasan YILMAZ
# Date:          25/03/25
# Location:      Trento, Italy
# Description:   The script is to create density plot of mapq values on 
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


plotter <- function(data, arguments, exp){
  plot <- data %>% ggplot(aes(x=mapq, fill=sample, color=sample, group = sample)) +
    geom_density(alpha=0.4) +
    labs(x='Mapping Quality (Map-Q)', y='Density') +
    theme_pubr(border = T, legend = 'right') +
    theme(legend.title=element_blank())
  
  ggsave(filename = glue('{arguments$prefix}{exp}_mapq.pdf'),
         device = 'pdf', plot = plot, path = arguments$output,
         height = 8, width = 10)
}


main <- function(arguments){
  data_re <- read.table(arguments$re, sep='\t', header=T)
  sample_levels <- data_re %>%
    distinct(sample) %>%
    arrange(!grepl("^Control", sample)) %>%
    pull(sample)
  data_re <- data_re %>%
    mutate(sample = factor(sample, levels = sample_levels))
  
  plotter(data_re, arguments, 'RE')
  
  data_total <- read.table(arguments$total, sep='\t', header=T)
  sample_levels <- data_total %>%
    distinct(sample) %>%
    arrange(!grepl("^Control", sample)) %>%
    pull(sample)
  data_total <- data_total %>%
    mutate(sample = factor(sample, levels = sample_levels))
  
  plotter(data_total, arguments, 'Total')
}

main(arguments = args)