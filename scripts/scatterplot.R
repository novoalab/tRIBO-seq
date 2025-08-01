# -----------------------------------------------------------------------------
# Title:         Scatter Plot 
# Author:        Hasan YILMAZ
# Date:          11/02/25
# Location:      Trento, Italy
# Description:   The script is to create scatter plots to compare experiments on 
#nano-tRNASeq data.
# -----------------------------------------------------------------------------

# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  optparse, ggplot2, ggpmisc, ggpubr, ggrepel, tibble, dplyr, tidyr,
  rlang, glue
)

# ----------------------------- Parse Arguments ----------------------------- #
option_list <- list(
  make_option(c("-r", "--ribo"), type="character",
              help="Count File of Ribo-nanotRNASeq", metavar="file"),
  make_option(c("-t", "--total"), type="character",
              help="Count File of nano-tRNASeq", metavar="file"),
  make_option(c("-p", "--prefix"), type="character",
              help="Prefix for experiments", metavar="character"),
  make_option(c("-s", "--samplefile"), type="character",
              help="Sample file", metavar="file"),
  make_option(c("-c", "--control"), type="character", 
              help="Name for control condition", metavar="character"),
  make_option(c("-y", "--treatment"), type="character",
              help="Name for treatment condition", metavar="character"),
  make_option(c("-x", "--axes"), type="integer",
              help="Axes limit", metavar="integer"),
  make_option(c("-o", "--output"), type="character",
              help="Output Directory", metavar="directory")
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

# -------------------------- Helper Functions ------------------------------ #

# 0. Check Arguments
check_args <- function(arguments){
  if (is.null(arguments$ribo) | is.null(arguments$total) |
      is.null(arguments$axes) | is.null(arguments$output)){
    print_help(opt_parser)
    stop("
    One or more mandatory arguments were not provided!
    Please provide, --ribo, --total, --axes, --output
         ")
  }
}

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

# 2. Data Preprocessing
rpm_normalization <- function(counts) {
  # Compute library size per sample (column)
  library_sizes <- colSums(counts)
  
  # Scale counts by library size and multiply by 1 million
  rpm <- sweep(counts, 2, library_sizes, FUN = "/") * 1e6
  
  return(rpm)
}

merge_data <- function(data_re, data_total, rowname){
  
  data_re$mean_re <- rowMeans(data_re)
  data_total$mean_total <- rowMeans(data_total)
  
  data_re <- rownames_to_column(data_re, rowname)
  data_total <- rownames_to_column(data_total, rowname)
  
  data_re <- data_re[,c('reference', 'mean_re')]
  data_total <- data_total[,c('reference', 'mean_total')]
  
  data_join <- full_join(data_re, data_total, by = rowname)
  
  return(data_join)
}

merge_data_all <- function(data_re, data_total){
  data_re <- tibble::rownames_to_column(data_re, 'reference')
  data_total <- tibble::rownames_to_column(data_total, 'reference')
  
  data_join <- full_join(data_re, data_total, by='reference')
  data_join <- tibble::column_to_rownames(data_join, 'reference')
  return(data_join)
}

process_data <- function(arguments, list_of_samples){
  
  data_re <- read.table(arguments$ribo, sep='\t', header = T)
  data_total <- read.table(arguments$total, sep='\t', header = T)
  
  data_re <- collapsetRNA(data_re)
  data_total <- collapsetRNA(data_total)
  
  data_re <- tibble::column_to_rownames(data_re ,'reference')
  data_total <- tibble::column_to_rownames(data_total ,'reference')
  
  norm_re <- rpm_normalization(data_re)
  norm_total <- rpm_normalization(data_total)
  
  if (arguments$axes < max(norm_re, norm_total)){
    message(glue("Max count in data is higher than {arguments$axes}"))
    message(glue("Please increase axes limits higher than {max(norm_re, norm_total)}"))
    stop("Axes limit smaller than data limits!")
  }
  
  if (list_of_samples){
    sample_list <- read.table(arguments$samplefile, sep = '\t', header = T)
    control_names <- sample_list[sample_list$Condition == arguments$control, "Sample"]
    treatment_names <- sample_list[sample_list$Condition == arguments$treatment, "Sample"]
    
    control_re <- norm_re[colnames(norm_re) %in% control_names]
    control_total <- norm_total[colnames(norm_total) %in% control_names]
    
    treatment_re <-  norm_re[colnames(norm_re) %in% treatment_names]
    treatment_total <- norm_total[colnames(norm_total) %in% treatment_names]
    
    
    control_data <- merge_data(control_re, control_total, 'reference')
    treatment_data <- merge_data(treatment_re, treatment_total, 'reference')
    
    return(list("control" = control_data,
                "treatment" = treatment_data))
  }else{
    return(merge_data(norm_re, norm_total, 'reference'))
  }
}


scatter_samples <- function(arguments, color_map, list_of_samples){
  
  if (list_of_samples == TRUE){
    message("Data Processing..")
    data_list <- process_data(arguments = arguments, list_of_samples = list_of_samples)
    
    for (name in c('control', 'treatment')){
      data <- as.data.frame(data_list[name][[1]])
      message(glue("Plots for {name} preparing.."))
      temp_data <- data %>% 
        separate(reference, into = c("AminoAcid", NA), sep = "-", remove = FALSE)
      plot <- ggplot(temp_data, aes(x=mean_total, y=mean_re))+
        geom_point(size = 3, aes(colour = AminoAcid)) +
        scale_color_manual(values = color_map) +
        lims(x = c(0, arguments$axes), y = c(0, arguments$axes)) +
        xlab('Normalized Counts\nNano-tRNAseq') +
        ylab('Normalized Counts\nRibo-nano-tRNAseq') +
        geom_text_repel(aes(label = reference), size = 3, max.overlaps = 15) +
        stat_cor(method = "spearman", cor.coef.name = 'rho') +
        stat_smooth(method='lm', formula = 'y ~ x', se = F, color='black') +
        theme_pubr(border=TRUE, legend = 'right')
      ggsave(filename = glue('{arguments$prefix}{name}_Total.vs.Ribo.pdf'),
             device = 'pdf', plot = plot, path = arguments$output, height = 8, width = 10)
    }
  }else{
    
    message("Data Processing..")
    data <- process_data(arguments = arguments, list_of_samples = list_of_samples)
    message(glue("Preparing the plot.."))
    temp_data <- data %>% 
      separate(reference, into = c("AminoAcid", NA), sep = "-", remove = FALSE)
    plot <- ggplot(temp_data, aes(x=mean_total, y=mean_re))+
      geom_point(size = 3, aes(colour = AminoAcid)) +
      scale_color_manual(values = color_map) +
      lims(x = c(0, arguments$axes), y = c(0, arguments$axes)) +
      xlab('Normalized Counts\nNano-tRNAseq') +
      ylab('Normalized Counts\nRibo-nano-tRNAseq') +
      geom_text_repel(aes(label = reference), size = 3, max.overlaps = 15) +
      stat_cor(method = "spearman", cor.coef.name = 'rho') +
      stat_smooth(method='lm', formula = 'y ~ x', se = F, color='black') +
      theme_pubr(border=TRUE, legend = 'right')
    ggsave(filename = glue('{arguments$prefix}Total.vs.Ribo.pdf'),
           device = 'pdf', plot = plot, path = arguments$output, height = 8, width = 10)
  }
}


scatter_by_sample <- function(arguments, color_map){
  
  counts_re <- read.table(arguments$ribo, sep='\t', header=T)
  counts_total <- read.table(arguments$total, sep='\t', header= T)
  
  counts_re <- collapsetRNA(counts_re)
  counts_total <- collapsetRNA(counts_total)
  
  counts_re <- tibble::column_to_rownames(counts_re ,'reference')
  counts_total <- tibble::column_to_rownames(counts_total ,'reference')
  
  norm_re <- rpm_normalization(counts_re)
  norm_total <- rpm_normalization(counts_total)
  
  if (arguments$axes < max(norm_re, norm_total)){
    message(glue("Max count in data is higher than {arguments$axes}"))
    message(glue("Please increase axes limits higher than {max(norm_re, norm_total)}"))
    stop("Axes limit smaller than data limits!")
  }
  
  counts_data <- merge_data_all(data_re = norm_re, data_total = norm_total)

  samples <- colnames(counts_data)
  
  for (i in 1:dim(counts_data)[2]){
    if (i == dim(counts_data)[2]){break}
    for (j in (i+1):dim(counts_data)[2]){
      temp_data = counts_data[,c(i,j)]
      temp_data <- tibble::rownames_to_column(temp_data, 'reference')
      temp_data <- temp_data %>% 
        separate(reference, into = c("AminoAcid", NA), sep = "-", remove = FALSE)
      
      plot <- ggplot(temp_data, aes(x=!!sym(samples[i]), y=!!sym(samples[j])))+
        geom_point(size = 3, aes(colour = AminoAcid)) +
        lims(x = c(0,arguments$axes), y = c(0,arguments$axes)) +
        scale_color_manual(values = color_map) +
        labs(x = glue("Normalized Counts\n {samples[i]}"),
             y = glue("Normalized Counts\n {samples[j]}")
        ) +
        geom_text_repel(aes(label = reference), size = 3, max.overlaps = 15) +
        stat_cor(method = "spearman", cor.coef.name = 'rho') +
        stat_smooth(method='lm', formula = 'y ~ x', se = F, color='black') +
        theme_pubr(border=TRUE, legend = 'right')
      ggsave(filename = glue('{arguments$prefix}{samples[i]}.vs.{samples[j]}.pdf'),
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
check_args(arguments = args)
if (is.null(args$samplefile)){
  message("samplefile not provided! Experiment-wise analysis initiated!")
  scatter_samples(arguments = args, color_map = color_map, list_of_samples = FALSE)
}else{
  message("samplefile provided! Condition-wise analysis initiated!")
  if (is.null(args$control) | is.null(args$treatment)){
    stop("One or more mandatory arguments were not provided!\
         Please provide, --control, --treatment")
  }
  scatter_samples(arguments = args, color_map = color_map, list_of_samples = TRUE)
}

scatter_by_sample(arguments = args, color_map = color_map)



