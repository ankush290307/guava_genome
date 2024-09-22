#PCA Plotting Script (plot_pca.R)
# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2", repos='http://cran.us.r-project.org')
if (!require("dplyr")) install.packages("dplyr", repos='http://cran.us.r-project.org')
if (!require("cowplot")) install.packages("cowplot", repos='http://cran.us.r-project.org')

# Load libraries
library(ggplot2)
library(dplyr)
library(cowplot)

# Read PCA data
pca_data <- read.table("results/plink/guava_pca.eigenvec", header=FALSE)
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:20))

# Merge with sample information
samples <- read.table("samples.txt", header=TRUE)
pca_data <- merge(pca_data, samples, by.x="IID", by.y="ID")

# Plot PCA
pca_plot <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Group, shape=Group)) +
  geom_point(size=3, alpha=0.8) +
  theme_cowplot() +
  labs(title="PCA of Guava Varieties", x="Principal Component 1", y="Principal Component 2") +
  scale_color_brewer(palette="Set3") +
  theme(text = element_text(size=12), legend.position="right")

# Save the plot
ggsave("results/plink/PCA_plot.png", plot=pca_plot, width=10, height=8, dpi=300)

#ADMIXTURE Cross-Validation Error Plotting Script (plot_admixture_cv.R)
# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2", repos='http://cran.us.r-project.org')

# Load libraries
library(ggplot2)

# Read CV error data
cv_data <- read.table("results/admixture/cv_error.txt", header=FALSE)
colnames(cv_data) <- c("K", "CV_Error")

# Plot CV error
cv_plot <- ggplot(cv_data, aes(x=K, y=CV_Error)) +
  geom_line(size=1) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="ADMIXTURE Cross-Validation Error", x="Number of Clusters (K)", y="CV Error") +
  theme(text = element_text(size=12))

# Save the plot
ggsave("results/admixture/CV_Error_plot.png", plot=cv_plot, width=8, height=6, dpi=300)
#ADMIXTURE Ancestry Proportions Plotting Script (plot_admixture_ancestry.R)

# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2", repos='http://cran.us.r-project.org')
if (!require("reshape2")) install.packages("reshape2", repos='http://cran.us.r-project.org')

# Load libraries
library(ggplot2)
library(reshape2)

# Read Q file for the best K (replace K with the chosen number, e.g., 3)
ancestry <- read.table("results/admixture/guava.3.Q")

# Add sample information
samples <- read.table("samples.txt", header=TRUE)
ancestry$IID <- samples$ID
ancestry <- merge(ancestry, samples, by="IID")

# Melt data for plotting
ancestry_melted <- melt(ancestry, id.vars=c("IID", "Group"))

# Plot ancestry proportions
ancestry_plot <- ggplot(ancestry_melted, aes(x=IID, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  facet_wrap(~Group, scales="free_x") +
  theme_minimal() +
  labs(title="Ancestry Proportions", x="Individuals", y="Ancestry Proportion") +
  scale_fill_brewer(palette="Set3") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=90, hjust=1))

# Save the plot
ggsave("results/admixture/Ancestry_plot.png", plot=ancestry_plot, width=12, height=6, dpi=300)


library(dplyr)
library(ggplot2)

# Load data
mydata <- read.table("nucleotide_diversity/nucleotide_diversity.windowed.pi", header = TRUE)
mydata1 <- na.omit(mydata)

# Filter for chromosomes C-Scaffold-1 to C-Scaffold-11
chromosomes <- paste0("C-Scaffold-", 1:11)
mydata1_filtered <- mydata1 %>% filter(CHROM %in% chromosomes)

# Define chromosome lengths
chrom_lengths <- data.frame(
  CHROM = chromosomes,
  length = c(62699497, 58864030, 58031391, 56684224, 55492731, 50572956, 49542903, 44952543, 44696232, 44064185, 39768966)
)

# Calculate cumulative positions
cumulative_lengths <- cumsum(chrom_lengths$length)
names(cumulative_lengths) <- chrom_lengths$CHROM

mydata1_filtered <- mydata1_filtered %>%
  mutate(
    cumulative_start = lag(cumulative_lengths[CHROM], default = 0),
    adjusted_position = BIN_START + cumulative_start
  )

# Calculate the 95th percentile of PI as the threshold
threshold_pi <- quantile(mydata1_filtered$PI, 0.95)

# Label positions for chromosomes
label_positions <- mydata1_filtered %>%
  group_by(CHROM) %>%
  summarize(label_pos = (max(adjusted_position) + min(adjusted_position)) / 2, .groups = 'drop')

# Plot
ggplot(mydata1_filtered, aes(x = adjusted_position, y = PI, color = CHROM)) +
  geom_line(size = 0.5) +
  geom_hline(yintercept = threshold_pi, color = "red", linetype = "dashed") +
  scale_color_manual(values = rainbow(length(chromosomes))) +
  scale_x_continuous(
    breaks = label_positions$label_pos,
    labels = chromosomes
  ) +
  labs(x = "Chromosome", y = "Nucleotide Diversity (π)", title = "Nucleotide Diversity Across Chromosomes") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# Save the plot
ggsave("nucleotide_diversity_plot.png", width = 12, height = 6)


#Fst Plotting Script (plot_fst.R)

# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2", repos='http://cran.us.r-project.org')
if (!require("dplyr")) install.packages("dplyr", repos='http://cran.us.r-project.org')

# Load libraries
library(ggplot2)
library(dplyr)

# Read Fst data
fst1 <- read.table("results/fst/group1_group2_fst.weir.fst", header=TRUE)
fst2 <- read.table("results/fst/group1_group3_fst.weir.fst", header=TRUE)
fst3 <- read.table("results/fst/group2_group3_fst.weir.fst", header=TRUE)

# Combine data
fst_data <- rbind(
  mutate(fst1, Comparison="Group1 vs Group2"),
  mutate(fst2, Comparison="Group1 vs Group3"),
  mutate(fst3, Comparison="Group2 vs Group3")
)

# Plot Fst
fst_plot <- ggplot(fst_data, aes(x=BIN_START, y=WEIR_AND_COCKERHAM_FST, color=Comparison)) +
  geom_line() +
  theme_minimal() +
  labs(title="Fst Values across Genomic Windows", x="Genomic Position", y="Fst") +
  scale_color_manual(values=c("blue", "green", "red")) +
  theme(text = element_text(size=12))

# Save the plot
ggsave("results/fst/Fst_plot.png", plot=fst_plot, width=12, height=6, dpi=300)

#Selective Sweep Detection Plotting Script (plot_sweed.R)

# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2", repos='http://cran.us.r-project.org')

# Load libraries
library(ggplot2)

# Read SweeD results
sweed_data <- read.table("results/sweed/guava_chr1.SweeD_Report", header=TRUE)  # Example for chromosome 1

# Plot SweeD results
sweed_plot <- ggplot(sweed_data, aes(x=Position, y=Likelihood)) +
  geom_line() +
  theme_minimal() +
  labs(title="Selective Sweep Detection", x="Genomic Position", y="Likelihood") +
  theme(text = element_text(size=12))

# Save the plot
ggsave("results/sweed/SweeD_plot.png", plot=sweed_plot, width=12, height=6, dpi=300)


#Tajima's D Plotting Script (plot_tajimas_d.R)

# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2", repos='http://cran.us.r-project.org')

# Load libraries
library(ggplot2)

# Read Tajima's D data
tajimas_d_data <- read.table("tajimas_d/tajimas_d.Tajima.D", header=TRUE)

# Plot Tajima's D
tajimas_d_plot <- ggplot(tajimas_d_data, aes(x=BIN_START, y=TajimaD, color=CHROM)) +
  geom_line() +
  theme_minimal() +
  labs(title="Tajima's D across Genomic Windows", x="Genomic Position", y="Tajima's D") +
  theme(text = element_text(size=12), legend.position="none")

# Save the plot
ggsave("tajimas_d/Tajimas_D_plot.png", plot=tajimas_d_plot, width=12, height=6, dpi=300)

tajimas_d_plot
# Load necessary libraries
library(dplyr)
library(circlize)
library(ggplot2)

# Read the nucleotide diversity data
mydata <- read.table("nucleotide_diversity/nucleotide_diversity.windowed.pi", header = TRUE)
mydata1 <- na.omit(mydata)

# Filter for chromosomes C-Scaffold-1 to C-Scaffold-11
chromosomes <- paste0("C-Scaffold-", 1:11)
mydata1_filtered <- mydata1 %>% filter(CHROM %in% chromosomes)

# Read Tajima's D data
tajimas_d_data <- read.table("tajimas_d/tajimas_d.Tajima.D", header = TRUE)
tajimas_d_filtered <- tajimas_d_data %>% filter(CHROM %in% chromosomes)

# Define chromosome lengths
chrom_lengths <- data.frame(
  CHROM = chromosomes,
  length = c(62699497, 58864030, 58031391, 56684224, 55492731, 50572956, 49542903, 44952543, 44696232, 44064185, 39768966)
)

# Prepare data for Circos plot
mydata1_filtered <- mydata1_filtered %>%
  left_join(chrom_lengths, by = "CHROM") %>%
  group_by(CHROM) %>%
  mutate(cumulative_start = lag(cumsum(c(0, head(length, -1))), default = 0)) %>%
  ungroup() %>%
  mutate(adjusted_position = BIN_START + cumulative_start)

tajimas_d_filtered <- tajimas_d_filtered %>%
  left_join(chrom_lengths, by = "CHROM") %>%
  group_by(CHROM) %>%
  mutate(cumulative_start = lag(cumsum(c(0, head(length, -1))), default = 0)) %>%
  ungroup() %>%
  mutate(adjusted_position = BIN_START + cumulative_start)

# Calculate the 95th percentile of PI as the threshold
threshold_pi <- quantile(mydata1_filtered$PI, 0.95)

# Initialize Circos plot
circos.clear()
circos.par("track.height" = 0.1, gap.degree = 1)
circos.initialize(factors = mydata1_filtered$CHROM, x = mydata1_filtered$adjusted_position)

# Add tracks for nucleotide diversity
circos.track(ylim = range(mydata1_filtered$PI), panel.fun = function(x, y) {
  circos.axis(labels.cex = 0.6, major.tick.percentage = 0.2)
})

circos.trackLines(factors = mydata1_filtered$CHROM, x = mydata1_filtered$adjusted_position, y = mydata1_filtered$PI, col = "blue")

# Add a cutoff line for the 95th percentile of PI
circos.trackLines(factors = mydata1_filtered$CHROM, x = mydata1_filtered$adjusted_position, y = rep(threshold_pi, nrow(mydata1_filtered)), col = "red", lty = 2)

# Add tracks for Tajima's D
circos.track(ylim = range(tajimas_d_filtered$TajimaD), panel.fun = function(x, y) {
  circos.axis(labels.cex = 0.6, major.tick.percentage = 0.2)
})

circos.trackLines(factors = tajimas_d_filtered$CHROM, x = tajimas_d_filtered$adjusted_position, y = tajimas_d_filtered$TajimaD, col = "green")

# Add chromosome labels
chrom_labels <- chrom_lengths %>%
  mutate(mid_point = cumsum(length) - length / 2)

for (i in 1:nrow(chrom_labels)) {
  circos.text(
    x = chrom_labels$mid_point[i],
    y = max(c(range(mydata1_filtered$PI), range(tajimas_d_filtered$TajimaD))),
    labels = chrom_labels$CHROM[i],
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.5),
    cex = 0.8
  )
}

# Save the plot
circos.clear()

# Load the required packages
library(ggplot2)
library(dplyr)
library(data.table)
library(gridExtra)

# Set working directory (adjust path as needed)
setwd("H:/GuavaGWAS/new")

# Define the order of the chromosomes to be used in the visualization
target <- c("C-Scaffold-1", "C-Scaffold-2", "C-Scaffold-3", "C-Scaffold-4", "C-Scaffold-5",
            "C-Scaffold-6", "C-Scaffold-7", "C-Scaffold-8", "C-Scaffold-9", "C-Scaffold-10", "C-Scaffold-11")

# Read the nucleotide diversity data
nucleotide_diversity_data <- fread("nucleotide_diversity/nucleotide_diversity.windowed.pi", header = TRUE) %>%
  na.omit() %>%
  filter(CHROM %in% target)

# Read Tajima's D data
tajimas_d_data <- fread("tajimas_d/tajimas_d.Tajima.D", header = TRUE) %>%
  na.omit() %>%
  filter(CHROM %in% target)

# Reorder the chromosome column of the data frames according to the target order
nucleotide_diversity_data$CHROM <- factor(nucleotide_diversity_data$CHROM, levels = target)
tajimas_d_data$CHROM <- factor(tajimas_d_data$CHROM, levels = target)

# Define threshold lines (example values, adjust as needed)
tajimas_d_threshold <- 2
nucleotide_diversity_threshold <- 0.02

# Function to create individual plots for each chromosome
create_plots <- function(data, y_var, y_label, threshold, title, color) {
  plots <- list()
  for (chr in target) {
    p <- ggplot(data %>% filter(CHROM == chr), aes(x = BIN_START, y = !!sym(y_var))) + 
      geom_line(color = color) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      labs(x = 'Chromosome Length', y = y_label, title = chr) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "grey80")
      )
    plots[[chr]] <- p
  }
  return(plots)
}

# Create plots for Nucleotide Diversity
nucleotide_diversity_plots <- create_plots(
  data = nucleotide_diversity_data,
  y_var = "PI",
  y_label = "Nucleotide Diversity (π)",
  threshold = nucleotide_diversity_threshold,
  title = "Nucleotide Diversity Across Chromosomes",
  color = "blue"
)

# Create plots for Tajima's D
tajimas_d_plots <- create_plots(
  data = tajimas_d_data,
  y_var = "TajimaD",
  y_label = "Tajima's D",
  threshold = tajimas_d_threshold,
  title = "Tajima's D Across Chromosomes",
  color = "green"
)

# Combine the plots for Nucleotide Diversity into a single figure
combined_nucleotide_diversity_plot <- marrangeGrob(grobs = nucleotide_diversity_plots, nrow = 4, ncol = 3)
ggsave(filename = "Nucleotide_Diversity_Combined.png", combined_nucleotide_diversity_plot, width = 20, height = 15)

# Combine the plots for Tajima's D into a single figure
combined_tajimas_d_plot <- marrangeGrob(grobs = tajimas_d_plots, nrow = 4, ncol = 3)
ggsave(filename = "Tajimas_D_Combined.png", combined_tajimas_d_plot, width = 20, height = 15)


# Load the required packages
library(ggplot2)
library(dplyr)
library(data.table)

# Set working directory (adjust path as needed)
setwd("H:/GuavaGWAS/new")

# Define the order of the chromosomes to be used in the visualization
target <- c("C-Scaffold-1", "C-Scaffold-2", "C-Scaffold-3", "C-Scaffold-4", "C-Scaffold-5",
            "C-Scaffold-6", "C-Scaffold-7", "C-Scaffold-8", "C-Scaffold-9", "C-Scaffold-10", "C-Scaffold-11")

# Read the nucleotide diversity data
nucleotide_diversity_data <- fread("nucleotide_diversity/nucleotide_diversity.windowed.pi", header = TRUE) %>%
  na.omit() %>%
  filter(CHROM %in% target)

# Reorder the chromosome column of the data frame according to the target order
nucleotide_diversity_data$CHROM <- factor(nucleotide_diversity_data$CHROM, levels = target)

# Calculate the threshold for nucleotide diversity (e.g., 95th percentile)
nucleotide_diversity_threshold <- quantile(nucleotide_diversity_data$PI, 0.95, na.rm = TRUE)
print(paste("Nucleotide Diversity Threshold (95th percentile):", nucleotide_diversity_threshold))

# Create a ggplot object for combined nucleotide diversity
nucleotide_diversity_plot <- ggplot(nucleotide_diversity_data, aes(x = CHROM, y = PI, fill = CHROM)) + 
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = nucleotide_diversity_threshold, linetype = "dashed", color = "red") +
  labs(x = 'Chromosome', y = 'Nucleotide Diversity (π)', title = "Nucleotide Diversity Across Chromosomes") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80"),
    legend.position = "right"
  ) +
  scale_fill_manual(values = rainbow(length(target)))

# Save the combined nucleotide diversity plot
ggsave(filename = "Combined_Nucleotide_Diversity_Plot.png", plot = nucleotide_diversity_plot, width = 15, height = 10)

