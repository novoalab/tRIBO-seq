
# ----------------------------- Load Libraries ----------------------------- #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  argparse, tidyr, dplyr, tibble, plyr, ggplot2, patchwork, ggridges, ggpubr, ggstatsplot, glue, wesanderson
)

# ----------------------------- Parse Arguments ----------------------------- #
parser <- ArgumentParser(description = "Fragmentation Check")
parser$add_argument("--input", type = "character", help = "Proportion data")
parser$add_argument("--prefix", type = "character", help = "Analysis Prefix")
parser$add_argument("--suffix", type = "character", help = "Experiment Type")
parser$add_argument("--output", type = "character", help = "output directory")
args <- parser$parse_args()


data <- read.table(args$input, sep='\t', header=T)
output = args$output

data <- data[!grepl('mt', data$reference),]
data$condition <- as.factor(data$condition)

message("Density plot creating...")
control_data <- data[data$condition == 'Control',]
mean_control <- mean(control_data$proportion)
max_density_control <- round(max(density(control_data$proportion)$y))
not_control_data <- data[data$condition != 'Control',]
mean_not_control <- mean(not_control_data$proportion)
max_density_not_control <- round(max(density(not_control_data$proportion)$y))

#ymax <- max(max_density_control, max_density_not_control)
ymax <- 10
control_color <- wes_palette("GrandBudapest1", 2, type = "discrete")[1]
not_control_color <- wes_palette("GrandBudapest1", 2, type = "discrete")[2]

#Plotting
plot_control <- control_data %>% ggplot(aes(x=proportion)) +
  geom_density(fill=control_color, color=control_color, alpha = 0.8) +
  geom_vline(xintercept = mean(control_data$proportion), color = control_color, linetype = 'dashed') +
  ylim(0, ymax) +
  annotate("text",
           x = mean_control - 0.05,
           y = ymax - (ymax * 0.09),  # Top of the plot
           label = paste0("Mean = ", round(mean_control, 3)),
           vjust = -0.5,
           hjust = 0.5,
           size = 4,
           color = control_color) +
  labs(title = unique(control_data$condition),
       subtitle = 'Proportion = Length of Read / Length of tRNA',
       x = 'Proportion', y = 'Density') +
  theme_pubr(border = T)

plot_not_control <- not_control_data %>% ggplot(aes(x=proportion)) +
  geom_density(fill=not_control_color, color=not_control_color, alpha = 0.8) +
  geom_vline(xintercept = mean(not_control_data$proportion), color = not_control_color, linetype = 'dashed') +
  ylim(0, ymax) +
  annotate("text",
           x = mean_not_control - 0.05,
           y =  ymax - (ymax * 0.09),  # Top of the plot
           label = paste0("Mean = ", round(mean_not_control, 3)),
           vjust = -0.5,
           hjust = 0.5,
           size = 4,
           color = not_control_color) +
  labs(title = unique(not_control_data$condition),
       x = 'Proportion', y = 'Density') +
  theme_pubr(border = T)

plot <- plot_control / plot_not_control

ggsave(plot = plot, filename = glue('{args$prefix}{args$suffix}_proportion_density.pdf'), device = 'pdf',
       height = 8, width = 12, path = output)