# Load the required packages
library(tidyverse)

# Read the SNP density data file
snpden <- read.table("guavagwas_hetsites.snpden", header = TRUE)

# Define the order of the chromosomes to be used in the visualization
target <- c("C-Scaffold-1", "C-Scaffold-2", "C-Scaffold-3", "C-Scaffold-4", "C-Scaffold-5",
            "C-Scaffold-6", "C-Scaffold-7", "C-Scaffold-8", "C-Scaffold-9", "C-Scaffold-10", "C-Scaffold-11")

# Reorder the chromosome column of the data frame according to the target order
snpden$CHROM <- factor(snpden$CHROM, levels = target)

# Create a ggplot object for all individuals combined
snpden_plot <- snpden %>%
  ggplot(aes(x = BIN_START, y = CHROM, fill = groups)) + 
  geom_tile() +
  labs(x = 'Chromosome Length', 
       y = 'Chromosome', 
       title = "Combined Heterozygous SNP Densities for All Individuals", 
       subtitle = "Subtitle text here") + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.y = unit(0.15, "lines"),
        plot.title = element_text(hjust = .5, size = 15),
        plot.subtitle = element_text(hjust = .5, size = 13, color = "gray")) +
  scale_fill_manual(values = c("#000081", "#0000f3", "#004dff", "#00b3ff", "#29ffce", 
                               "#7bff7b", "#ceff29", "#ffc600", "#ff6800", "#f30900", "brown", "#800000", "black", "pink", "cyan", "purple", "magenta"),
                    name = "Variants/kb",
                    labels = c("<0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25",
                               "0.25-0.50","0.50-0.75","0.75-1.0","1.0-1.25","1.25-1.5",
                               "1.5-1.75","1.75-2.0","2.0-2.25","2.25-2.5")) +
  scale_x_continuous(name = 'Chromosome length', 
                     labels = c('0Mb', '50Mb', '100Mb', '150Mb', '200Mb', '250Mb'),
                     breaks = c(0, 50000000, 100000000, 150000000, 200000000, 250000000), 
                     expand = c(0, 0))

# Save the plot
ggsave(filename = "Combined_SNP_density_heatmap.png", 
       plot = snpden_plot, 
       device = 'png',
       dpi = 600, 
       units = 'cm', 
       width = 28, 
       height = 18, 
       bg = "white")
