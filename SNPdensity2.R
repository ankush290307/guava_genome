# Load the required packages
library(ggplot2)
library(dplyr)

# Set working directory (adjust path as needed)
setwd("H:/GuavaGWAS/new")

# Read the SNP density data file
snpden <- read.table("guavagwas_hetsites.snpden", header = TRUE)

# Define the order of the chromosomes to be used in the visualization
target <- c("C-Scaffold-1", "C-Scaffold-2", "C-Scaffold-3", "C-Scaffold-4", "C-Scaffold-5",
            "C-Scaffold-6", "C-Scaffold-7", "C-Scaffold-8", "C-Scaffold-9", "C-Scaffold-10", "C-Scaffold-11")

# Reorder the chromosome column of the data frame according to the target order
snpden$CHROM <- factor(snpden$CHROM, levels = target)

# Remove rows with NA values in VARIANTS.KB
snpden <- snpden[!is.na(snpden$VARIANTS.KB), ]

# Create a ggplot object for overall SNP density heatmap with spacing between chromosomes
snpden_plot <- ggplot(snpden, aes(x = BIN_START, y = CHROM, fill = VARIANTS.KB)) + 
  geom_tile() +
  labs(x = 'Chromosome Length', 
       y = 'Chromosome', 
       title = "Overall Heterozygous SNP Densities", 
       subtitle = "Across All Individuals") + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = .5, size = 15),
        plot.subtitle = element_text(hjust = .5, size = 13, color = "gray"),
        strip.text.y = element_text(angle = 0, size = 8),
        panel.spacing = unit(0.5, "lines")) +  # Increase spacing between chromosomes
  scale_fill_gradient(low = "pink", high = "darkgreen", name = "Variants/kb") +
  scale_x_continuous(name = 'Chromosome length', 
                     labels = c('0Mb', '50Mb', '100Mb', '150Mb', '200Mb', '250Mb'),
                     breaks = c(0, 50000000, 100000000, 150000000, 200000000, 250000000), 
                     expand = c(0, 0))

# Save the plot
ggsave(filename = "Overall_SNP_density_heatmap.png", 
       plot = snpden_plot, 
       device = 'png',
       dpi = 600, 
       units = 'cm', 
       width = 28, 
       height = 18, 
       bg = "white")

