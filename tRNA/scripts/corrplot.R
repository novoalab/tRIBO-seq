# -----------------------------------------------------------------------------
# Title:         Correlation Plot
# Author:        Hasan YILMAZ
# Date:          11/02/25
# Location:      Trento, Italy
# Description:   The script is to create correlation plot to compare samples on 
#nano-tRNASeq data.
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, ggplot2, GGally, corrplot, ggpubr, grDevices,
  tibble, tidyr, dplyr, glue
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Correlation Check")
parser$add_argument("--total", type = "character", help = "Counts of Total tRNASeq file path")
parser$add_argument("--re", type = "character", help = "Counts of Ribo-Embedded tRNASeq file path")
parser$add_argument("--prefix", type = "character", help = "Analysis Prefix")
parser$add_argument("--output", type = "character", help = "output directory")
args <- parser$parse_args()


# -------------------------- Helper Functions ------------------------------ #
# 1. Collapse iso-acceptors to single tRNA
collapsetRNA <- function(data){
  
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

# 2. Read and merge the data
preprocess <- function(arguments){
  message('Preprocessing Started...')
  data_re <- read.table(arguments$re, sep='\t', header=T)
  data_re <- collapsetRNA(data_re)
  data_re <- data_re %>% rename_with(~ paste0(.x, "_RE"), .cols = 2:ncol(data_re))
  data_total <- read.table(arguments$total, sep='\t', header=T)
  data_total <- collapsetRNA(data_total)
  data_total <- data_total %>% rename_with(~ paste0(.x, "_Total"), .cols = 2:ncol(data_total))
  
  data_return <- data_total %>% full_join(data_re, by='reference')
  data_return <- tibble::column_to_rownames(data_return,'reference')
  data_return <- sweep(data_return, 2, colSums(data_return), "/")
  
  return(data_return)
}

# 3. Create correlation Plot
correlation_plot <- function(counts_data, arguments){
  
  message('Correlation plots are preparing...')
  
  plot <- ggpairs(counts_data, showStrips = F, proportions = 'auto') + 
    theme_pubr(border=T) +
    font('xy.text', size=6)
  
  ggsave(filename = glue('{arguments$prefix}correlation.pdf'),
         device = 'pdf', plot = plot, path = arguments$output,
         height = 12, width = 12)
}

# 4.Create correlatoin matrix plot
correlation_matrix <- function(counts_data, arguments){

  message('Correlation heatmaps are preparing...')
  
  corr_matrix <- cor(counts_data)
  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  
  grDevices::cairo_pdf(file = glue('{arguments$output}/{arguments$prefix}correlation_matrix.pdf'),
                       height = 12, width = 12)
  
  corrplot(corr_matrix, method = "color", col=col(200), 
           addCoef.col = "black", tl.col="black", tl.srt=45, is.corr = F, col.lim = c(0,1))
  
  dev.off()
}


# ---------------------------- Execute Script ----------------------------- #
data <- preprocess(arguments = args)
correlation_plot(counts_data = data, arguments = args)
correlation_matrix(counts_data = data, arguments = args)
