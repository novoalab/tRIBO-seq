#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Title:         Find Fragmented tRNAs
# Author:        Hasan YILMAZ
# Date:          29/04/25
# Location:      Trento, Italy
# Description:   This script is to find fragmented tRNAs from count tables
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(argparse, BiocManager, DESeq2, EnhancedVolcano, tibble, tidyr,
               dplyr, rlang, ggplot2, ggpubr, glue, ggrastr, data.table)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Fragmented tRNA Detection")
parser$add_argument(
  "--re_fl",
  type = "character",
  help = "Path to RE Full Length counts.tsv")
parser$add_argument(
  "--re_fr",
  type = "character",
  help = "Path to RE Fragments counts.tsv")
parser$add_argument(
  "--total_fl",
  type = "character",
  help = "Path to Total Full Length counts.tsv")
parser$add_argument(
  "--total_fr",
  type = "character",
  help = "Path to Total Fragments counts.tsv")
parser$add_argument(
  "--sample_list",
  type = "character",
  help = "Path to sample list.tsv file"
)
parser$add_argument(
  "--output",
  type = "character",
  help = "Path to output directory")
parser$add_argument(
  "--prefix",
  type = "character",
  help = "Analysis Prefix")
parser$add_argument(
  "--condition",
  type = "character",
  help = "Base condition of condition (eg. CONTROL or TREATMENT)")
parser$add_argument(
  "--type",
  type = "character",
  help = "Base condition of type (eg. TOTAL or RE)")
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

# 2.Preprocessing
preprocess <- function(data_re_fl, data_re_fr,
                       data_total_fl, data_total_fr,
                       sample_list,
                       process = "standard"){
  message('Data preprocessing for analysis...')
  
  if (process == "standard"){
    
    sample_dict <- setNames(sample_list$Condition, sample_list$Sample)
    
    return_list = list()
    for (data_type in c("RE", "TOTAL")){
      if (data_type == "RE"){
        data_fl <- collapsetRNA(data_re_fl)
        data_fl <- data_fl %>% rename_with(~ paste0(.x, "_", sample_dict[.x] ,"_FL"), .cols = 2:ncol(data_fl))
        data_fr <- collapsetRNA(data_re_fr)
        data_fr <- data_fr %>% rename_with(~ paste0(.x, "_", sample_dict[.x], "_FR"), .cols = 2:ncol(data_fr))
        
        temp_data <- data_fl %>% full_join(data_fr, by='reference')
        temp_data <- tibble::column_to_rownames(temp_data,'reference')
        return_list[[data_type]] = temp_data
      }
      else {
        data_fl <- collapsetRNA(data_total_fl)
        data_fl <- data_fl %>% rename_with(~ paste0(.x, "_", sample_dict[.x], "_FL"), .cols = 2:ncol(data_fl))
        data_fr <- collapsetRNA(data_total_fr)
        data_fr <- data_fr %>% rename_with(~ paste0(.x, "_", sample_dict[.x], "_FR"), .cols = 2:ncol(data_fr))
        
        temp_data <- data_fl %>% full_join(data_fr, by='reference')
        temp_data <- tibble::column_to_rownames(temp_data,'reference')
        return_list[[data_type]] = temp_data
      }
    }
  }
  else{
    #Get sample names to split data
    control_samples <- unique(sample_list[sample_list$Condition == "CONTROL", "Sample"])
    treatment_samples <- unique(sample_list[sample_list$Condition != "CONTROL", "Sample"])
    
    data_re_fl <- collapsetRNA(data_re_fl)
    data_re_fr <- collapsetRNA(data_re_fr)
    data_total_fl <- collapsetRNA(data_total_fl)
    data_total_fr <- collapsetRNA(data_total_fr)
    
    return_list = list()
    for (data_type in c("CONTROL", "TREATMENT")){
      
      if (data_type == "CONTROL"){
        temp_re_fl <- data_re_fl[, colnames(data_re_fl) %in%
                                   c("reference", control_samples)]
        
        temp_re_fl <- temp_re_fl %>% 
          rename_with(~ paste0(.x, "_RE", "_FL"), .cols = 2:ncol(temp_re_fl))
        
        temp_re_fr <- data_re_fr[, colnames(data_re_fr) %in%
                                   c("reference", control_samples)]
        
        temp_re_fr <- temp_re_fr %>%
          rename_with(~ paste0(.x, "_RE", "_FR"), .cols = 2:ncol(temp_re_fl))
        
        temp_total_fl <- data_total_fl[, colnames(data_total_fl) %in%
                                         c("reference", control_samples)]
        
        temp_total_fl <- temp_total_fl %>%
          rename_with(~ paste0(.x, "_TOTAL", "_FL"), .cols = 2:ncol(temp_total_fl))
        
        temp_total_fr <- data_total_fr[, colnames(data_total_fr) %in%
                                         c("reference", control_samples)]
        temp_total_fr <- temp_total_fr %>%
          rename_with(~ paste0(.x, "_TOTAL", "_FR"), .cols = 2:ncol(temp_total_fr))
        
        
        temp_data <- temp_total_fl %>%
          full_join(temp_total_fr, by = "reference") %>%
          full_join(temp_re_fl, by = "reference") %>%
          full_join(temp_re_fr, by = "reference")
        
        temp_data <- tibble::column_to_rownames(temp_data,'reference')
        
        return_list[[data_type]] = temp_data
      }
      else {
        temp_re_fl <- data_re_fl[, colnames(data_re_fl) %in%
                                   c("reference", treatment_samples)]
        
        temp_re_fl <- temp_re_fl %>% 
          rename_with(~ paste0(.x, "_RE", "_FL"), .cols = 2:ncol(temp_re_fl))
        
        temp_re_fr <- data_re_fr[, colnames(data_re_fr) %in%
                                   c("reference", treatment_samples)]
        
        temp_re_fr <- temp_re_fr %>%
          rename_with(~ paste0(.x, "_RE", "_FR"), .cols = 2:ncol(temp_re_fl))
        
        temp_total_fl <- data_total_fl[, colnames(data_total_fl) %in%
                                         c("reference", treatment_samples)]
        
        temp_total_fl <- temp_total_fl %>%
          rename_with(~ paste0(.x, "_TOTAL", "_FL"), .cols = 2:ncol(temp_total_fl))
        
        temp_total_fr <- data_total_fr[, colnames(data_total_fr) %in%
                                         c("reference", treatment_samples)]
        temp_total_fr <- temp_total_fr %>%
          rename_with(~ paste0(.x, "_TOTAL", "_FR"), .cols = 2:ncol(temp_total_fr))
        
        temp_data <- temp_total_fl %>%
          full_join(temp_total_fr, by = "reference") %>%
          full_join(temp_re_fl, by = "reference") %>%
          full_join(temp_re_fr, by = "reference")
        
        temp_data <- tibble::column_to_rownames(temp_data,'reference')
        
        return_list[[data_type]] = temp_data
      }
    }
  }
  message(glue("Preprocessing Done for {process} analysis"))
  return(return_list)
}

# 3.Two factor differential abundance analysis
two_factor_analysis <- function(data, condition, sub_out, arguments){
  
  dir.create(glue("{arguments$output}/{sub_out}"), recursive = TRUE)
  
  data <- copy(data)
  
  design_matrix <- data.frame(SampleID=colnames(data))
  design_matrix <- design_matrix %>%
    separate(SampleID,
             into = c("Sample", "ID", "Condition", "Type"),
             sep = "_",
             remove = FALSE) %>%
    select(-Sample, -ID)
  
  data[data == 0] = 0L
  
  message('Running two factor analysis...')
  dds = DESeqDataSetFromMatrix(
    countData = data,
    colData = design_matrix,
    design =~ Type*Condition)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- dds[rowMedians(counts(dds)) >= 1]
  dds$Condition <- relevel(dds$Condition, ref = condition)
  dds$Type <- relevel(dds$Type, ref = "FL")
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[4], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x = result, file = glue('{arguments$output}/{sub_out}/{arguments$prefix}Type*Condition_diffAbundances.tsv'),
              sep = "\t", quote = FALSE)
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  message('Plots are preparing...')
  pca <- plotPCA(vsd, intgroup = "Type") +
    geom_label_repel(aes(label=name)) +
    theme_pubr(border = T) +
    theme(legend.position = "none", aspect.ratio = 1)
  ggsave(glue('{arguments$prefix}Type*Condition_pca.pdf'), plot = pca,
         device = 'pdf', path = glue("{arguments$output}/{sub_out}"),
         width = 8, height = 8)
  
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
  
  ggsave(glue('{arguments$prefix}Type*Condition_volcano.pdf'), plot = volcano,
         device = 'pdf', path = glue("{arguments$output}/{sub_out}"),
         width = 8, height = 12)
  
  
}

# ---------------------------- Execute Script ----------------------------- #

sample_list <- read.table(args$sample_list, sep='\t', header = T)
sample_list <- sample_list %>%
  mutate(across(c(Type, Condition), ~ toupper(.)))

data_re_full_length <- read.table(args$re_fl, sep = '\t', header = T)
data_re_fragment <- read.table(args$re_fr, sep = '\t', header = T)
data_total_full_length <- read.table(args$total_fl, sep = '\t', header = T)
data_total_fragment <- read.table(args$total_fr, sep = '\t', header = T)

condition_data <- preprocess(data_re_fl = data_re_full_length,
                             data_re_fr = data_re_fragment,
                             data_total_fl = data_total_full_length,
                             data_total_fr = data_total_fragment,
                             sample_list = sample_list,
                             process = "standard")

type_data <- preprocess(data_re_fl = data_re_full_length,
                        data_re_fr = data_re_fragment,
                        data_total_fl = data_total_full_length,
                        data_total_fr = data_total_fragment,
                        sample_list = sample_list,
                        process = "type")

two_factor_analysis(data = condition_data[["RE"]], condition = args$condition,
                    sub_out = "RE", arguments = args)

two_factor_analysis(data = condition_data[["TOTAL"]], condition = args$condition,
                    sub_out = "TOTAL", arguments = args)

two_factor_analysis(data = type_data[["CONTROL"]], condition = args$type,
                    sub_out = "CONTROL", arguments = args)

two_factor_analysis(data = type_data[["TREATMENT"]], condition = args$type,
                    sub_out = "TREATMENT", arguments = args)
