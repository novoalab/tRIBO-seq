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
parser$add_argument("--all", type = "character", help = "Counts of All tRNAs")
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

diff_exp <- function(arguments, data, sample_list, t){
  data <- tibble::column_to_rownames(data,'reference')
  
  data[data == 0] = 0L
  
  message('Running one factor analysis...')
  dds = DESeqDataSetFromMatrix(
    countData = data,
    colData = sample_list,
    design =~ Type)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- dds[rowMedians(counts(dds)) >= 1]
  dds$Type <- relevel(dds$Type, ref = "All")
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[2], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x = result, file = glue('{arguments$output}/{arguments$prefix}{t}_Fragmentation.tsv'),
              sep = "\t", quote = FALSE)
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  message('Plots are preparing...')
  pca <- plotPCA(vsd, intgroup = "Type") +
    geom_label_repel(aes(label=name)) +
    theme_pubr(border = T) +
    theme(legend.position = "none", aspect.ratio = 1)
  ggsave(glue('{arguments$prefix}{t}_Fragmentation_pca.pdf'), plot = pca,
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
  
  ggsave(glue('{arguments$prefix}{t}_Fragmentation_volcano.pdf'), plot = volcano,
         device = 'pdf', path = arguments$output, width = 8, height = 12)
}

diff_exp_two_fact <- function(arguments, data_ct, data_tr,
                              sample_list_ct, sample_list_tr, t){
  
  data <- data_ct %>% left_join(data_tr, by ='reference')
  
  sample_list <- rbind(sample_list_ct, sample_list_tr)
  sample_list$Condition <- factor(sample_list$Condition, ordered = F)
  sample_list$Type <- factor(sample_list$Type, ordered = F)
  
  data <- tibble::column_to_rownames(data,'reference')
  
  data[data == 0] = 0L
  
  message('Running two factor analysis...')
  dds = DESeqDataSetFromMatrix(
    countData = data,
    colData = sample_list,
    design =~ Type)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- dds[rowMedians(counts(dds)) >= 1]
  dds$Condition <- relevel(dds$Condition, ref = "Control")
  dds$Type <- relevel(dds$Type, ref = "All")
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[2], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x = result, file = glue('{arguments$output}/{arguments$prefix}{t}_Fragmentation.tsv'),
              sep = "\t", quote = FALSE)
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  message('Plots are preparing...')
  pca <- plotPCA(vsd, intgroup = "Type") +
    geom_label_repel(aes(label=name)) +
    theme_pubr(border = T) +
    theme(legend.position = "none", aspect.ratio = 1)
  ggsave(glue('{arguments$prefix}{t}_Fragmentation_pca.pdf'), plot = pca,
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
  
  ggsave(glue('{arguments$prefix}{t}_Fragmentation_volcano.pdf'), plot = volcano,
         device = 'pdf', path = arguments$output, width = 8, height = 12)
}

# 2.Two factor differential abundance analysis
process_data <- function(arguments){
  
  message('Data preprocessing for differential expression analysis...')
  data_full <- read.table(arguments$full, sep='\t', header=T)
  data_full <- collapsetRNA(data_full)
  data_fragment <- read.table(arguments$all, sep='\t', header=T)
  data_fragment <- collapsetRNA(data_fragment)
  
  sample_list <- read.table(arguments$samples, sep='\t', header = T)
  sample_list_ct <- sample_list[sample_list$Condition == 'Control',]
  sample_list_tr <- sample_list[sample_list$Condition != 'Control',]
  
  if (dim(sample_list_ct)[1] == 0){
    data_ct = data.frame()
    sample_list_full_tr <- sample_list_tr[sample_list_tr$Type == 'Full',]
    data_full_tr <- data_full[,c('reference', sample_list_full_tr$Sample)]
    sample_list_fragment_tr <- sample_list_tr[sample_list_tr$Type == 'All',]
    data_fragment_tr <- data_fragment[,c('reference', sample_list_fragment_tr$Sample)]
    data_tr <- data_full_tr %>% left_join(data_fragment_tr, by ='reference',
                                          suffix = c("_Full", "_All"))
    diff_exp(arguments = arguments, data = data_tr, sample_list = sample_list_tr, t = 'Treatment')
  }
  if (dim(sample_list_tr)[1] == 0){
    data_tr = data.frame()
    sample_list_full_ct <- sample_list_ct[sample_list_ct$Type == 'Full',]
    
    data_full_ct <- data_full[,c('reference', sample_list_full_ct$Sample)]
    
    sample_list_fragment_ct <- sample_list_ct[sample_list_ct$Type == 'All',]
    
    data_fragment_ct <- data_fragment[,c('reference', sample_list_fragment_ct$Sample)]
    
    
    data_ct <- data_full_ct %>% left_join(data_fragment_ct, by ='reference',
                                          suffix = c("_Full", "_All"))
    diff_exp(arguments = arguments, data = data_ct, sample_list = sample_list_ct, t = 'Control')
  }
  
  if (dim(sample_list_tr)[1] != 0 && dim(sample_list_ct)[1] != 0){
    sample_list_full_tr <- sample_list_tr[sample_list_tr$Type == 'Full',]
    data_full_tr <- data_full[,c('reference', sample_list_full_tr$Sample)]
    sample_list_fragment_tr <- sample_list_tr[sample_list_tr$Type == 'All',]
    data_fragment_tr <- data_fragment[,c('reference', sample_list_fragment_tr$Sample)]
    data_tr <- data_full_tr %>% left_join(data_fragment_tr, by ='reference',
                                          suffix = c("_Full", "_All"))
    diff_exp(arguments = arguments, data = data_tr, sample_list = sample_list_tr, t = 'Treatment')
    sample_list_full_ct <- sample_list_ct[sample_list_ct$Type == 'Full',]
    
    data_full_ct <- data_full[,c('reference', sample_list_full_ct$Sample)]
    
    
    sample_list_fragment_ct <- sample_list_ct[sample_list_ct$Type == 'All',]
    
    data_fragment_ct <- data_fragment[,c('reference', sample_list_fragment_ct$Sample)]
    
    
    data_ct <- data_full_ct %>% left_join(data_fragment_ct, by ='reference',
                                          suffix = c("_Full", "_All"))
    diff_exp(arguments = arguments, data = data_ct, sample_list = sample_list_ct, t = 'Control')
    
    
    diff_exp_two_fact(arguments = arguments, data_ct = data_ct, data_tr = data_tr,
                      sample_list_ct = sample_list_ct, sample_list_tr = sample_list_tr,
                      t = 'Treatment.vs.Control')
  }
}

# ---------------------------- Execute Script ----------------------------- #

process_data(arguments = args)

