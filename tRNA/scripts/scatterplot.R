# -----------------------------------------------------------------------------
# Title:         Scatter Plot
# Author:        Hasan YILMAZ
# Date:          11/02/25
# Location:      Trento, Italy
# Description:   The script is to create scatter plots to compare samples on 
#nano-tRNASeq data.
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, ggplot2, ggpmisc, ggpubr, ggrepel, tibble, dplyr, tidyr,
  rlang, glue
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Library Reproducibility Check")
parser$add_argument("--counts", type = "character", help = "Counts file path")
parser$add_argument("--prefix", type = "character", help = "Analysis Prefix")
parser$add_argument("--type", type = "character", help = "Sample Type Name (RE|Total)")
parser$add_argument("--output", type = "character", help = "output directory")
args <- parser$parse_args()


# -------------------------- Helper Functions ------------------------------ #
# 1. Collapse iso-acceptors to single tRNA
collapsetRNA <- function(arguments){
  
  data <- read.table(arguments$counts, sep='\t', header = T)
  
  data <- data[!grepl('mt|oligo|RDN|unique|antisense', data$reference),]
  
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

# ----------------------------- Main Pipeline ------------------------------ #
scatter_samples <- function(arguments, color_map){
  
  counts_data <- collapsetRNA(arguments)
  
  counts_data <- tibble::column_to_rownames(counts_data ,'reference')
  
  samples <- colnames(counts_data)
  
  for (i in 1:dim(counts_data)[2]){
    if (i == dim(counts_data)[2]){break}
    for (j in (i+1):dim(counts_data)[2]){
      temp_data = counts_data[,c(i,j)]
      temp_data <- sweep(temp_data, 2, colSums(temp_data), "/")
      temp_data <- tibble::rownames_to_column(temp_data, 'reference')
      temp_data <- temp_data %>% 
        separate(reference, into = c("AminoAcid", NA), sep = "-", remove = FALSE)
      plot <- ggplot(temp_data, aes(x=!!sym(samples[i]), y=!!sym(samples[j])))+
        geom_point(size = 3, aes(colour = AminoAcid)) +
        scale_color_manual(values = color_map) +
        geom_text_repel(aes(label = reference), size = 3, max.overlaps = 15) +
        stat_cor(method = "spearman", cor.coef.name = 'rho') +
        stat_smooth(method='lm', formula = 'y ~ x', se = F, color='black') +
        theme_pubr(border=TRUE, legend = 'right')
      ggsave(filename = glue('{arguments$prefix}{arguments$type}_{samples[i]}.vs.{samples[j]}.pdf'),
             device = 'pdf', plot = plot, path = arguments$output, height = 8, width = 10)
    }
  }
  
}

color_map <- c('Ala' = '#e6194B', 'Arg' = '#3cb44b', 'Asn' = '#ffe119',
               'Asp' = '#4363d8', 'Cys' = '#f58231', 'Gln' = '#911eb4',
               'Glu' = '#42d4f4', 'Gly' = '#f032e6', 'His' = '#bfef45',
               'Ile' = '#fabed4', 'Leu' = '#469990', 'Lys' = '#dcbeff',
               'Met' = '#9A6324', 'Phe' = '#fffac8', 'Pro' = '#800000',
               'SeC' = '#aaffc3', 'Ser' = '#808000', 'Thr' = '#ffd8b1',
               'Trp' = '#000075', 'Tyr' = '#a9a9a9', 'Val' = '#000000',
               'iMet' = '#00441b')

# ---------------------------- Execute Script ----------------------------- #
scatter_samples(arguments = args, color_map = color_map)
