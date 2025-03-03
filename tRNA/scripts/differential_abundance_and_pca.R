# -----------------------------------------------------------------------------
# Title:         Differential Expression and PCA Plot
# Author:        Hasan YILMAZ
# Date:          13/02/25
# Location:      Trento, Italy
# Description:   The script is to create differential expression and PCA plot
#to compare samples on nano-tRNASeq data.
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, BiocManager, DESeq2, EnhancedVolcano, tibble, tidyr, dplyr, rlang,
  ggplot2, ggpubr, glue, ggrastr
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Differential Expression and PC Analysis")
parser$add_argument("--total", type = "character", help = "Counts of Total tRNASeq file path")
parser$add_argument("--re", type = "character", help = "Counts of Ribo-Embedded tRNASeq file path")
parser$add_argument("--prefix", type = "character", help = "Analysis Prefix")
parser$add_argument("--sample_list", type = 'character', help = "Sample list file path")
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

# 2.Single factor differential abundance analysis
single_factor_analysis <- function(arguments, analysis, seq_type){
  
  message('Data preprocessing for single factor analysis...')
  data_re <- read.table(arguments$re, sep='\t', header=T)
  data_re <- collapsetRNA(data_re)
  data_re <- data_re %>% rename_with(~ paste0(.x, "_RE"), .cols = 2:ncol(data_re))
  data_total <- read.table(arguments$total, sep='\t', header=T)
  data_total <- collapsetRNA(data_total)
  data_total <- data_total %>% rename_with(~ paste0(.x, "_TOTAL"), .cols = 2:ncol(data_total))
  
  data <- data_total %>% full_join(data_re, by='reference')
  data <- tibble::column_to_rownames(data,'reference')
  
  sample_list <- read.table(arguments$sample_list, sep='\t', header = T)
  sample_list$Sample <- paste0(sample_list$Sample, '_', sample_list$Type)
  
  if (analysis == 'Condition'){
    sample_list <- sample_list[sample_list$Type == toupper(seq_type),]
    data <- data[,colnames(data) %in% sample_list$Sample]
    dsgn <- formula(~ Condition)
  }
  else if (analysis == 'Type'){
    sample_list <- sample_list[sample_list$Condition == seq_type,]
    data <- data[,colnames(data) %in% sample_list$Sample]
    dsgn <- formula(~ Type)
  }
  

  data[data == 0] = 0L
  
  message('Running single factor analysis...')
  dds = DESeqDataSetFromMatrix(
    countData = data,
    colData = sample_list,
    design = dsgn)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- dds[rowMedians(counts(dds)) >= 1]
  if (analysis == 'Condition'){
    dds$Condition <- relevel(dds$Condition, ref = "Control")
  }
  else if (analysis == 'Type'){
    dds$Type <- relevel(dds$Type, ref = "TOTAL")
  }
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[2], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x = result, file = glue('{arguments$output}/{arguments$prefix}{toupper(seq_type)}_diffAbundances.tsv'),
              sep = "\t", quote = FALSE)
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  message('Plots are preparing...')
  
  pca <- plotPCA(vsd, intgroup = analysis) +
    geom_label_repel(aes(label=name)) +
    theme_pubr(border = T) +
    theme(legend.position = "none", aspect.ratio = 1)
  ggsave(glue('{arguments$prefix}{toupper(seq_type)}_pca.pdf'), plot = pca,
         device = 'pdf', path = arguments$output, width = 8, height = 8)
  
  result <- result[!is.na(result$padj),]
  
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
  
  ggsave(glue('{arguments$prefix}{toupper(seq_type)}_volcano.pdf'), plot = volcano,
         device = 'pdf', path = arguments$output, width = 8, height = 12)
  
}

# 3.Two factor differential abundance analysis
two_factor_analysis <- function(arguments){
  
  message('Data preprocessing for two factor analysis...')
  data_re <- read.table(arguments$re, sep='\t', header=T)
  data_re <- collapsetRNA(data_re)
  data_re <- data_re %>% rename_with(~ paste0(.x, "_RE"), .cols = 2:ncol(data_re))
  data_total <- read.table(arguments$total, sep='\t', header=T)
  data_total <- collapsetRNA(data_total)
  data_total <- data_total %>% rename_with(~ paste0(.x, "_TOTAL"), .cols = 2:ncol(data_total))
  
  data <- data_total %>% full_join(data_re, by='reference')
  data <- tibble::column_to_rownames(data,'reference')
  
  sample_list <- read.table(arguments$sample_list, sep='\t', header = T)
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
  dds$Type <- relevel(dds$Type, ref = "TOTAL")
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[4], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x = result, file = glue('{arguments$output}/{arguments$prefix}Type*Condition_diffAbundances.tsv'),
              sep = "\t", quote = FALSE)
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  message('Plots are preparing...')
  pca <- plotPCA(vsd, intgroup = "Type") +
    geom_label_repel(aes(label=name)) +
    theme_pubr(border = T) +
    theme(legend.position = "none", aspect.ratio = 1)
  ggsave(glue('{arguments$prefix}Type*Condition_pca.pdf'), plot = pca,
         device = 'pdf', path = arguments$output, width = 8, height = 8)
  
  result <- result[!is.na(result$padj),]
  
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
  
  ggsave(glue('{arguments$prefix}Type*Condition_volcano.pdf'), plot = volcano,
         device = 'pdf', path = arguments$output, width = 8, height = 12)
  
  
}


# ---------------------------- Execute Script ----------------------------- #
conditions <- read.table(args$sample_list, sep ='\t', header = T)
conditions <- unique(conditions$Condition)

single_factor_analysis(arguments = args, analysis = 'Condition', seq_type = 're')
single_factor_analysis(arguments = args, analysis = 'Condition', seq_type = 'total')
single_factor_analysis(arguments = args, analysis = 'Type', seq_type = conditions[1])
single_factor_analysis(arguments = args, analysis = 'Type', seq_type = conditions[2])
two_factor_analysis(arguments = args)

