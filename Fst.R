# Load the required packages
library(ggplot2)
library(dplyr)
library(data.table)
library(gridExtra)



# Read Fst data for each comparison
fst_group1_vs_group2 <- fread("guavagwas_group1_vs_group2.weir.fst", header = TRUE)
fst_group1_vs_group3 <- fread("guavagwas_group1_vs_group3.weir.fst", header = TRUE)
fst_group2_vs_group3 <- fread("guavagwas_group2_vs_group3.weir.fst", header = TRUE)

# Inspect the structure of each file
str(fst_group1_vs_group2)
str(fst_group1_vs_group3)
str(fst_group2_vs_group3)

# Check the first few rows of each file
head(fst_group1_vs_group2)
head(fst_group1_vs_group3)
head(fst_group2_vs_group3)
# Define the order of the chromosomes to be used in the visualization
target <- c("C-SCAFFOLD-1", "C-SCAFFOLD-2", "C-SCAFFOLD-3", "C-SCAFFOLD-4", "C-SCAFFOLD-5",
            "C-SCAFFOLD-6", "C-SCAFFOLD-7", "C-SCAFFOLD-8", "C-SCAFFOLD-9", "C-SCAFFOLD-10", "C-SCAFFOLD-11")
# Calculate Fst threshold (e.g., 95th percentile)
calculate_fst_threshold <- function(fst_data) {
  quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
}

fst_threshold1 <- calculate_fst_threshold(fst_group1_vs_group2)
fst_threshold2 <- calculate_fst_threshold(fst_group1_vs_group3)
fst_threshold3 <- calculate_fst_threshold(fst_group2_vs_group3)

# Print calculated thresholds
print(paste("Fst threshold for Group 1 vs Group 2:", fst_threshold1))
print(paste("Fst threshold for Group 1 vs Group 3:", fst_threshold2))
print(paste("Fst threshold for Group 2 vs Group 3:", fst_threshold3))

# Function to create individual plots for each chromosome
create_fst_plots <- function(data, y_label, threshold, color) {
  plots <- list()
  for (chr in target) {
    p <- ggplot(data %>% filter(CHROM == chr), aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) + 
      geom_line(color = color) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      labs(x = 'Chromosome Length', y = y_label, title = chr) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "grey80"),
        panel.spacing = unit(0.5, "lines")  # Increase spacing between panels
      )
    plots[[chr]] <- p
  }
  return(plots)
}

# Create plots for Fst comparisons with calculated thresholds
fst_group1_vs_group2_plots <- create_fst_plots(
  data = fst_group1_vs_group2,
  y_label = "Fst (Group 1 vs Group 2)",
  threshold = fst_threshold1,
  color = "blue"
)

fst_group1_vs_group3_plots <- create_fst_plots(
  data = fst_group1_vs_group3,
  y_label = "Fst (Group 1 vs Group 3)",
  threshold = fst_threshold2,
  color = "green"
)

fst_group2_vs_group3_plots <- create_fst_plots(
  data = fst_group2_vs_group3,
  y_label = "Fst (Group 2 vs Group 3)",
  threshold = fst_threshold3,
  color = "purple"
)

# Combine the plots for each comparison into single figures
combined_fst_group1_vs_group2_plot <- marrangeGrob(grobs = fst_group1_vs_group2_plots, nrow = 4, ncol = 3, top = "Fst (Group 1 vs Group 2) Across Chromosomes")
ggsave(filename = "Fst_Group1_vs_Group2_Combined.png", combined_fst_group1_vs_group2_plot, width = 20, height = 15)

combined_fst_group1_vs_group3_plot <- marrangeGrob(grobs = fst_group1_vs_group3_plots, nrow = 4, ncol = 3, top = "Fst (Group 1 vs Group 3) Across Chromosomes")
ggsave(filename = "Fst_Group1_vs_Group3_Combined.png", combined_fst_group1_vs_group3_plot, width = 20, height = 15)

combined_fst_group2_vs_group3_plot <- marrangeGrob(grobs = fst_group2_vs_group3_plots, nrow = 4, ncol = 3, top = "Fst (Group 2 vs Group 3) Across Chromosomes")
ggsave(filename = "Fst_Group2_vs_Group3_Combined.png", combined_fst_group2_vs_group3_plot, width = 20, height = 15)

# Plot all Fst comparisons in single combined plots

# Function to create combined plot for Fst comparisons
create_combined_fst_plot <- function(data1, data2, data3, y_label1, y_label2, y_label3, color1, color2, color3) {
  combined_data <- bind_rows(
    data1 %>% mutate(comparison = y_label1),
    data2 %>% mutate(comparison = y_label2),
    data3 %>% mutate(comparison = y_label3)
  )
  
  ggplot(combined_data, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = comparison)) + 
    geom_line() +
    geom_hline(yintercept = max(fst_threshold1, fst_threshold2, fst_threshold3), linetype = "dashed", color = "red") +
    labs(x = 'Chromosome Length', y = "Fst", title = "Fst Comparisons Across Chromosomes") + 
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey80"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c(color1, color2, color3)) +
    facet_wrap(~CHROM, ncol = 1)
}

# Create combined plot for all comparisons
combined_fst_plot <- create_combined_fst_plot(
  data1 = fst_group1_vs_group2,
  data2 = fst_group1_vs_group3,
  data3 = fst_group2_vs_group3,
  y_label1 = "Group 1 vs Group 2",
  y_label2 = "Group 1 vs Group 3",
  y_label3 = "Group 2 vs Group 3",
  color1 = "blue",
  color2 = "green",
  color3 = "purple"
)

# Save the combined Fst plot
ggsave(filename = "Combined_Fst_Plot.png", plot = combined_fst_plot, width = 15, height = 20)

