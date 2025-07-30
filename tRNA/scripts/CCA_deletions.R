# -----------------------------------------------------------------------------
# Title:         Deletion Bar Plot 
# Author:        Hasan YILMAZ
# Date:          14/07/25
# Location:      Trento, Italy
# Description:   The script is to create barplot for deletions in CCA for nanotRNAseq
# -----------------------------------------------------------------------------


# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(optparse, data.table, dplyr, tibble, ggplot2, ggpubr, glue, gsubfn
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
              help="Sample File for experiment", metavar="file"),
  make_option(c("-o", "--output"), type="character",
              help="Output Directory", metavar="directory")
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

# ----------------------------- Helper Functions ---------------------------- #

create_list_columns <- function(sample_list){
  base_list <- c('A', 'G', 'C', 'T', 'insertion')
  combinations <- expand.grid(sample_list, base_list)
  include_list <- paste(combinations$Var1, combinations$Var2, sep =' ')
  deletion_list <- paste(sample_list, 'deletion', sep=' ')
  return(list('include_list' = include_list,
              'deletion_list' = deletion_list))
}


read_n_process <- function(arguments){
  
  sample_list <- fread(arguments$samplefile)
  
  ct_re_list <- sample_list[Condition == 'Control' & Type == 'RE', Sample]
  tr_re_list <- sample_list[Condition != 'Control' & Type == 'RE', Sample]
  ct_tt_list <- sample_list[Condition == 'Control' & Type == 'TOTAL', Sample]
  tr_tt_list <- sample_list[Condition != 'Control' & Type == 'TOTAL', Sample]
  
  list[ct_re_inc, ct_re_del] = create_list_columns(ct_re_list)
  list[tr_re_inc, tr_re_del] = create_list_columns(tr_re_list)
  list[ct_tt_inc, ct_tt_del] = create_list_columns(ct_tt_list)
  list[tr_tt_inc, tr_tt_del] = create_list_columns(tr_tt_list)
  
  data_re <- fread(arguments$ribo)
  data_re <- data_re[sprinzl_number %in% c('C1', 'C2', 'A'),]
  data_re <- data_re[!grepl('^mt', data_re$reference),]
  data_total <- fread(arguments$total)
  data_total <- data_total[sprinzl_number %in% c('C1', 'C2', 'A'),]
  data_total <- data_total[!grepl('^mt', data_total$reference),]

  
  
  data_re_ct <- data_re[, row_sum := rowSums(.SD), .SDcols = ct_re_inc]
  data_re_ct[, deletion_sum := rowSums(.SD), .SDcols = ct_re_del]
  data_re_ct[, deletion_percent := deletion_sum*100/row_sum]
  data_re_ct$type <- 'Ribo-nanotRNAseq'
  
  data_re_tr <- data_re[, row_sum := rowSums(.SD), .SDcols = tr_re_inc]
  data_re_tr[, deletion_sum := rowSums(.SD), .SDcols = tr_re_del]
  data_re_tr[, deletion_percent := deletion_sum*100/row_sum]
  data_re_tr$type <- 'Ribo-nanotRNAseq'
  
  
  data_tt_ct <- data_total[, row_sum := rowSums(.SD), .SDcols = ct_tt_inc]
  data_tt_ct[, deletion_sum := rowSums(.SD), .SDcols = ct_tt_del]
  data_tt_ct[, deletion_percent := deletion_sum*100/row_sum]
  data_tt_ct$type <- 'nano-tRNAseq'
  
  data_tt_tr <- data_total[, row_sum := rowSums(.SD), .SDcols = tr_tt_inc]
  data_tt_tr[, deletion_sum := rowSums(.SD), .SDcols = tr_tt_del]
  data_tt_tr[, deletion_percent := deletion_sum*100/row_sum]
  data_tt_tr$type <- 'nano-tRNAseq'
  
  
  data_plot_ct <- rbind(data_re_ct[,c('reference', 'sprinzl_number',
                                'deletion_percent', 'type')],
                        data_tt_ct[,c('reference', 'sprinzl_number',
                                   'deletion_percent', 'type')]
  )
  data_plot_ct$sprinzl_number <- factor(
    data_plot_ct$sprinzl_number, levels = c('C1', 'C2', 'A')
    )
  data_plot_ct$type <- factor(data_plot_ct$type,
                              levels = c('Ribo-nanotRNAseq', 'nano-tRNAseq'))
  
  
  data_plot_tr <- rbind(data_re_tr[,c('reference', 'sprinzl_number',
                                      'deletion_percent', 'type')],
                        data_tt_tr[,c('reference', 'sprinzl_number',
                                      'deletion_percent', 'type')]
  )
  data_plot_tr$sprinzl_number <- factor(
    data_plot_tr$sprinzl_number, levels = c('C1', 'C2', 'A')
  )
  data_plot_tr$type <- factor(data_plot_tr$type,
                              levels = c('Ribo-nanotRNAseq', 'nano-tRNAseq'))
  
  return(list('control_plot' = data_plot_ct,
              'treatment_plot' = data_plot_tr)
  )
  
}

plotting <- function(data_plot, arguments, t){
  
  y1 = ggplot(data_plot, aes(x=reference, y=deletion_percent, fill=sprinzl_number)) +
    ylim(c(0,100)) +
    theme_classic() +
    labs(y = 'Percent Deletions', x='', fill = 'Positions') +
    geom_bar(position=position_dodge(), stat="identity", colour='black')
  
  y1 <- y1 +
    facet_wrap(~type, nrow=2, strip.position = "bottom") +
    theme_pubr(border = TRUE) +
    theme(strip.placement = "outside", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  data_plot_2 <-data_plot[,
                          .(sum_deletion_percents = mean(deletion_percent)),
                          by=c('sprinzl_number', 'type')]
  data_plot_2$sprinzl_number <- factor(data_plot_2$sprinzl_number, levels = c('C1', 'C2', 'A'))
  
  y2 = ggplot(data_plot_2, aes(x=sprinzl_number, y=sum_deletion_percents, fill=sprinzl_number)) +
    ylim(c(0,100)) +
    theme_classic() +
    labs(y = 'Percent Deletions', x='', fill = 'Positions') +
    geom_bar(position=position_dodge(), stat="identity", colour='black')
  
  y2 <- y2 +
    facet_wrap(~type, nrow=2, strip.position = "bottom") +
    theme_pubr(border = TRUE) +
    theme(strip.placement = "outside")
  
  ggsave(filename = glue("{arguments$prefix}{t}_CCA_deletions_aminoacid.pdf"),
         plot = y1,
         device = 'pdf',
         path = arguments$output,
         height = 10, width = 20)
  
  ggsave(filename = glue("{arguments$prefix}{t}_CCA_deletions_sum.pdf"),
         plot = y2,
         device = 'pdf',
         path = arguments$output,
         height = 10, width = 8)
}

list[control_plot, treatment_plot] <- read_n_process(args)
plotting(data_plot = control_plot, arguments = args, t='Control')
#plotting(data_plot = treatment_plot, arguments = args, t='Treatment')

