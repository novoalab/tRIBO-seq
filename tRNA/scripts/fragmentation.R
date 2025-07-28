# -----------------------------------------------------------------------------
# Title:         Fragmentation Analysis
# Author:        Hasan YILMAZ
# Date:          11/02/25
# Location:      Trento, Italy
# Description:   The script is to analyse fragment data with full-length data on 
#nano-tRNASeq data.
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, BiocManager, DESeq2, EnhancedVolcano, tibble, tidyr, dplyr, rlang,
  ggplot2, ggpubr, glue, ggrastr
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Fragmentation Analysis")
parser$add_argument("--full", type = "character", help = "Counts of Full-Length tRNAs")
parser$add_argument("--fragment", type = "character", help = "Counts of Fragmented tRNAs")
parser$add_argument("--samples", type = "character", help = "sample list")
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

# 2.Two factor differential abundance analysis
two_factor_analysis <- function(arguments){
  
  message('Data preprocessing for two factor analysis...')
  data_full <- read.table(arguments$full, sep='\t', header=T)
  data_full <- collapsetRNA(data_full)
  data_full <- data_full %>% rename_with(~ paste0(.x, "_Full"), .cols = 2:ncol(data_full))
  data_fragment <- read.table(arguments$fragment, sep='\t', header=T)
  data_fragment <- collapsetRNA(data_fragment)
  data_fragment <- data_fragment %>% rename_with(~ paste0(.x, "_Fragment"), .cols = 2:ncol(data_fragment))
  
  data <- data_full %>% full_join(data_fragment, by='reference')
  data <- tibble::column_to_rownames(data,'reference')
  
  sample_list <- read.table(arguments$samples, sep='\t', header = T)
  sample_list$Sample <- paste0(sample_list$Sample, '_', sample_list$Type)
  
  data[data == 0] = 0L
  
  message('Running two factor analysis...')
  dds = DESeqDataSetFromMatrix(
    countData = data,
    colData = sample_list,
    design =~ Type*Condition)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- dds[rowMedians(counts(dds)) >= 1]
  dds$Condition <- relevel(dds$Condition, ref = "Control")
  dds$Type <- relevel(dds$Type, ref = "Full")
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[4], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x = result, file = glue('{arguments$output}/{arguments$prefix}Fragmentation.tsv'),
              sep = "\t", quote = FALSE)
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  message('Plots are preparing...')
  pca <- plotPCA(vsd, intgroup = "Type") +
    geom_label_repel(aes(label=name)) +
    theme_pubr(border = T) +
    theme(legend.position = "none", aspect.ratio = 1)
  ggsave(glue('{arguments$prefix}Fragmentation_pca.pdf'), plot = pca,
         device = 'pdf', path = arguments$output, width = 8, height = 8)
  
  result <- result[!is.na(result$padj),]
  result[result$padj == 0, "padj"] = 2.225074e-307
  
  absMax <- ceiling(max(abs(result$log2FoldChange)))
  
  volcano <- EnhancedVolcano(result,
                             lab = rownames(result),
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = '',
                             subtitle = '',
                             xlim = c(-absMax, absMax),
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             pointSize = 5.0,
                             labSize = 6,
                             legendLabSize = 12,
                             drawConnectors = TRUE,
                             ylim = c(0,max(-log10(result$padj))+0.3),
                             col=c('#848484', '#13BCCF', '#A6D113', '#EB9200'),
                             colAlpha = 1, raster = T) + theme_pubr(border = T)
  
  ggsave(glue('{arguments$prefix}Fragmentation_volcano.pdf'), plot = volcano,
         device = 'pdf', path = arguments$output, width = 8, height = 12)
  
  
}

# ---------------------------- Execute Script ----------------------------- #

two_factor_analysis(arguments = args)


