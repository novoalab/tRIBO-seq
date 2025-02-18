# -----------------------------------------------------------------------------
# Title:         Differential Modification Plot
# Author:        Hasan YILMAZ
# Date:          13/02/25
# Location:      Trento, Italy
# Description:   The script is to create differential modification heatmaps
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, BiocManager, DESeq2, EnhancedVolcano, ComplexHeatmap, circlize,
  tibble, tidyr, dplyr, rlang, gsubfn, ggplot2, ggpubr, glue, ggrastr,
  rtracklayer, seqinr, data.table
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Optimized nano-tRNASeq Analysis")
parser$add_argument("--total", type = "character", help = "Total tRNASeq file path")
parser$add_argument("--re", type = "character", help = "Ribo-Embedded tRNASeq file path")
parser$add_argument("--prefix", type = "character", help = "Analysis Prefix")
parser$add_argument("--sample_list", type = 'character', help = "Sample list file path")
parser$add_argument("--fasta", type = 'character', help = "Reference Fasta file path")
parser$add_argument("--output", type = "character", help = "Output directory")
args <- parser$parse_args()

# -------------------------- Helper Functions ------------------------------ #
# 1. Build Reference Data
ref_data_builder <- function(arguments){
  fasta <- read.fasta(arguments$fasta, as.string = T)
  
  tRNAs <- setdiff(getName(fasta), c("reference_annotation", "secondary_structure"))
  tRNAs <- tRNAs[!grepl('mt', tRNAs)]
  
  reference <- getAnnot(fasta$reference_annotation)
  reference <- strsplit(strsplit(reference, " ")[[1]][2], ",")[[1]]
  reference <- data.table(sprinzl_number = reference)
  reference[(.N-1), sprinzl_number := 'C2']
  reference[(.N-2), sprinzl_number := 'C1']
  
  reference_vector <- reference$sprinzl_number
  
  final_reference <- expand.grid(reference = unlist(tRNAs), sprinzl_number = reference$sprinzl_number)
  final_reference <- final_reference[order(final_reference$reference), ]
  rownames(final_reference) <- NULL
  
  single_reference <- rbindlist(lapply(fasta, function(item) {
    if (grepl('mt|reference_annotation|secondary_structure', getName(item))) return(NULL)
    temp_ref <- copy(reference)
    temp_ref[, bases := getSequence(item)]
    temp_ref <- temp_ref[bases != '-', .(reference = getName(item), sprinzl_number)]
    temp_ref
  }))
  
  return(list(
    'single_reference' = single_reference,
    'reference' = final_reference,
    'reference_vector' = reference_vector
    ))
  }

# 2. Calculate Error Counts
error_count <- function(counts, sample_list){
  
  return_data <- counts[,c('reference', 'position', 'sprinzl_number', 'base')]
  
  for (sample in sample_list$Sample){
    for (i in 1:dim(counts)[1]){
      bases <- c('A', 'G', 'C', 'T', 'deletion', 'insertion')
      base <- counts$base[i]
      colBase <- paste0(sample, '.', base)
      columns <- bases[-grep(base, bases)]
      columns <-  paste0(sample, '.', columns)
      if (base == ''){
        return_data[i, paste0(sample, '-notmodified')] <- NA
        return_data[i, paste0(sample, '-modified')] <- NA
      }
      else{
        return_data[i, paste0(sample, '-notmodified')] <- counts[i, colBase]
        return_data[i, paste0(sample, '-modified')] <- sum(counts[i, columns])
      }
    }
  }
  return(return_data)
}

deseq2_single <- function(count, design, arguments, seq_type){
  
  count[count == 0] = 0L
  
  dds = DESeqDataSetFromMatrix(
    countData = count,
    colData = design,
    design = ~Condition*Mod)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- dds[rowMedians(counts(dds)) >= 1]
  
  dds$Condition <- relevel(dds$Condition, ref = "Control")
  dds$Mod <- relevel(dds$Mod, ref = "notmodified")
  
  dds = DESeq(dds)
  
  result = results(dds, name = resultsNames(dds)[4], alpha = 0.05, pAdjustMethod = 'BH')
  result = as.data.frame(result)
  
  write.table(x=result, file = glue('{arguments$output}/{arguments$prefix}{toupper(seq_type)}_mod.tsv'),
              sep='\t', quote = F, col.names = T, row.names = T)
  
  volcano <- EnhancedVolcano(result,
                             lab = rownames(result),
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = '',
                             subtitle = '',
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             pointSize = 5.0,
                             labSize = 6,
                             legendLabSize = 12,
                             drawConnectors = TRUE,
                             ylim = c(0,max(-log10(result$padj))+0.3),
                             col=c('#848484', '#13BCCF', '#A6D113', '#EB9200'),
                             colAlpha = 1, raster = T) + theme_pubr(border = T)
  
  ggsave(glue('{arguments$prefix}{toupper(seq_type)}_mod.pdf'), plot = volcano,
         device = 'pdf', path = arguments$output, width = 8, height = 12)
  
  
  return(result)
}

# 4. Generate Heatmaps
heatmap <- function(data, reference, reference_vector, arguments, seq_type){
  
  data <- reference %>% left_join(data, by = c('reference', 'sprinzl_number'))
  
  data$sprinzl_number <- factor(data$sprinzl_number, levels = reference_vector, ordered = TRUE)
  matrix <- spread(data[,c('reference', 'sprinzl_number', 'log2FoldChange')],
                        key = sprinzl_number, value = log2FoldChange)
  matrix <- column_to_rownames(matrix, 'reference')
  
  maxValue <- max(ceiling(max(matrix, na.rm = T)), abs(floor(min(matrix, na.rm = T))))
  
  col_fun = colorRamp2(
    c(-maxValue, -0.5, 0, 0.5, maxValue),
    c("#2b8cbe","#FFFFF0", "#FFFFF0","#FFFFF0", "#e34a33")
    )
  
  to_change <- data$padj > 0.05 & !is.na(data$log2FoldChange)
  filtered_matrix <- data
  filtered_matrix$log2FoldChange[to_change] <- 0
  filtered_matrix <- spread(filtered_matrix[,c('reference', 'sprinzl_number', 'log2FoldChange')],
                            key = sprinzl_number, value = log2FoldChange)
  filtered_matrix <- column_to_rownames(filtered_matrix, 'reference')
  
  
  pdf(glue('{arguments$output}/{arguments$prefix}{toupper(seq_type)}_mod_heatmap.pdf'),
      width = 16, height = 10)
  full_map <- Heatmap(
    as.matrix(matrix),
    name = "log2FC",
    col = col_fun,
    show_row_names = TRUE,
    show_column_names = TRUE,
    na_col = 'grey',
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  draw(full_map)
  dev.off()
  
  pdf(glue('{arguments$output}/{arguments$prefix}{toupper(seq_type)}_mod_filtered_heatmap.pdf'),
      width = 16, height = 10)
  filtered_map <- Heatmap(
    as.matrix(filtered_matrix),
    name = "log2FC",
    col = col_fun,
    show_row_names = TRUE,
    show_column_names = TRUE,
    na_col = 'grey',
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  draw(filtered_map)
  dev.off()
}

# ----------------------------- Main Pipeline ------------------------------ #
process_data_single <- function(arguments){
  message("Building Reference Data...")
  list[single_ref, reference, ref_vector] <- ref_data_builder(arguments = arguments)
  
  sample_list <- read.table(arguments$sample_list, sep='\t', header=T)
  sample_list_re <- sample_list[sample_list$Type == toupper('re'),]
  sample_list_total <- sample_list[sample_list$Type == toupper('total'),]
  
  re <- read.table(arguments$re, sep='\t', header=T)
  re <- re[re$sprinzl_number != '',]
  re <- re[!grepl('mt', re$reference),]
  total <- read.table(arguments$total, sep='\t', header=T)
  total <- total[total$sprinzl_number != '',]
  total <- total[!grepl('mt', total$reference),]
  
  message("Calculating Errors for RE...")
  error_re <- error_count(counts = re, sample_list = sample_list_re)
  error_re$reference <- paste0(error_re$reference, '_', error_re$sprinzl_number)
  error_re <- subset(error_re, select=-c(position, base, sprinzl_number))
  rownames(error_re) <- NULL
  error_re <- column_to_rownames(error_re, 'reference')
  
  message("Calculating Errors for Total...")
  error_total <- error_count(counts = total, sample_list = sample_list_total)
  error_total$reference <- paste0(error_total$reference, '_', error_total$sprinzl_number)
  error_total <- subset(error_total, select=-c(position, base, sprinzl_number))
  rownames(error_total) <- NULL
  error_total <- column_to_rownames(error_total, 'reference')
  
  design_re <- sample_list_re %>%
    mutate(Mod = "notmodified") %>%
    bind_rows(sample_list_re %>% mutate(Mod = "modified")) %>%
    mutate(
      Sample = paste0(Sample, "-", Mod)
    ) %>%
    select(Sample, Condition, Mod) %>%
    arrange(factor(Condition, levels = unique(sample_list_re$Condition)),
            gsub("-.*", "", Sample), 
            factor(Mod, levels = c("notmodified", "modified")))
  
  design_total <- sample_list_total %>%
    mutate(Mod = "notmodified") %>%
    bind_rows(sample_list_total %>% mutate(Mod = "modified")) %>%
    mutate(
      Sample = paste0(Sample, "-", Mod)
    ) %>%
    select(Sample, Condition, Mod) %>%
    arrange(factor(Condition, levels = unique(sample_list_total$Condition)),
            gsub("-.*", "", Sample), 
            factor(Mod, levels = c("notmodified", "modified")))
  
  message("Running DESeq2 for RE...")
  results_re <- deseq2_single(count = error_re, design = design_re, arguments = arguments, seq_type = 're')
  results_re <- rownames_to_column(results_re, 'reference')
  results_re <- results_re %>% separate(reference, sep = '_', into = c('reference', 'sprinzl_number'))
  results_re <- single_ref %>% left_join(results_re, by = c('reference', 'sprinzl_number'))
  results_re <- results_re %>%
    mutate(
      baseMean = ifelse(is.na(baseMean), 0, baseMean),
      log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
      lfcSE = ifelse(is.na(lfcSE), 0, lfcSE),
      stat = ifelse(is.na(stat), 0, stat),
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      padj = ifelse(is.na(padj), 1, padj)
    )
  
  message("Running DESeq2 for Total...")
  results_total <- deseq2_single(count = error_total, design = design_total, arguments = arguments, seq_type = 'total')
  results_total <- rownames_to_column(results_total, 'reference')
  results_total <- results_total %>% separate(reference, sep = '_', into = c('reference', 'sprinzl_number'))
  results_total <- single_ref %>% left_join(results_total, by = c('reference', 'sprinzl_number'))
  results_total <- results_total %>%
    mutate(
      baseMean = ifelse(is.na(baseMean), 0, baseMean),
      log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
      lfcSE = ifelse(is.na(lfcSE), 0, lfcSE),
      stat = ifelse(is.na(stat), 0, stat),
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      padj = ifelse(is.na(padj), 1, padj)
    )
  
  message("Creating Heatmaps...")
  heatmap(data = results_re, reference = reference,
          reference_vector = ref_vector, arguments = arguments, seq_type = 're')
  heatmap(data = results_total, reference = reference,
          reference_vector = ref_vector, arguments = arguments, seq_type = 'total')
  
  message("Analysis Complete.")
  
}

# ---------------------------- Execute Script ----------------------------- #
process_data_single(arguments = args)
