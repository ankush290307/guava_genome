# Install and load necessary libraries
install.packages("factoextra")
library(factoextra)
library(tidyverse)



# Read PCA data
pca <- read_table2("guavapca.eigenvec", col_names = FALSE)
eigenval <- scan("guavapca.eigenval")

# Set names for individuals and PCs
names(pca)[1] <- "FID"
names(pca)[2] <- "IID"
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

# Organize individuals into groups
group1 <- c("AC10-5", "AC10-7", "AC10-8", "AC11-9", "AC1-4", "AC14-4", "AC1-5", "AC18-4", "AC2-1", "AC2-3", "AC3-4", "AC5-2", "AC5-5", "AC5-6", "AC5-7", "AC6-2", "AC6-4", "AC8-2", "AC8-4", "AC8-5", "AC8-6", "AC-9", "CISH-G1", "CISH-G5", "CISH-G6", "Strawberry_guava")
group2 <- c("50-50", "AICRP10-4", "Arka_kiran", "HB17-16", "HB8-8", "HB9-2", "Hisar_surkha", "Lalit", "Pear_shaped", "Portugal", "Punjab_kiran", "Punjab_pink", "Taiwan_pink")
group3 <- c("2y2", "AICRP10-5", "Allahabad_safeda", "Arka_amulya", "AS11-3", "AS12-2", "AS13-2", "AS15-3", "AS1-6", "AS1-7", "AS2-8", "AS5-8", "AS7-3", "AS7-7", "AS7-8", "Barf_khana", "BDS", "Behat_coconut", "BP-1", "BP-2", "G-Bilas", "GVP", "HB2-6", "HB4-8", "HB5-5", "HB5-8", "HB6-3", "HB9-7", "Hisar_safeda", "KMG", "L-49", "Pant_prabhat", "Punjab_safeda", "R2P4", "Sabir", "Sabir-2", "Sabir-3", "Seedless", "Shweta", "Thai_guava", "TS-1", "VNR", "Lemon_guava")
group4 <- c("Purple_local")

# Add group information to the PCA dataframe
pca <- pca %>%
  mutate(group = case_when(
    IID %in% group1 ~ "Group 1",
    IID %in% group2 ~ "Group 2",
    IID %in% group3 ~ "Group 3",
    IID %in% group4 ~ "Group 4",
    TRUE ~ "Unknown"
  ))

# Calculate variance explained by each PC
pve <- eigenval / sum(eigenval) * 100

# Function to create PCA plot with ellipses
create_pca_plot <- function(pca_data, x, y, pve, color_palette) {
  ggplot(pca_data, aes_string(x = x, y = y, color = "group")) + 
    geom_point(size = 3) +
    stat_ellipse(type = "norm") +
    scale_color_manual(values = color_palette) +
    coord_fixed(ratio = 1) +
    theme_classic() +
    xlab(paste0(x, " (", signif(pve[as.numeric(substr(x, 3, 3))], 3), "%)")) + 
    ylab(paste0(y, " (", signif(pve[as.numeric(substr(y, 3, 3))], 3), "%)")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

# Define color palette
color_palette <- c("Group 1" = "red", "Group 2" = "lightpink", "Group 3" = "green", "Group 4" = "purple")

# Create PCA plots
pca_plot_pc1_pc2 <- create_pca_plot(pca, "PC1", "PC2", pve, color_palette)
pca_plot_pc1_pc3 <- create_pca_plot(pca, "PC1", "PC3", pve, color_palette)
pca_plot_pc2_pc3 <- create_pca_plot(pca, "PC2", "PC3", pve, color_palette)

# Save the PCA plots
ggsave("PCA_plot_PC1_PC2.png", pca_plot_pc1_pc2, width = 10, height = 8)
ggsave("PCA_plot_PC1_PC3.png", pca_plot_pc1_pc3, width = 10, height = 8)
ggsave("PCA_plot_PC2_PC3.png", pca_plot_pc2_pc3, width = 10, height = 8)

# Combine PCA plots into one figure for publication
combined_pca_plot <- ggarrange(pca_plot_pc1_pc2, pca_plot_pc1_pc3, pca_plot_pc2_pc3, 
                               ncol = 1, nrow = 3, common.legend = TRUE, legend = "right")

# Save the combined PCA plot
ggsave("Combined_PCA_Plot.png", combined_pca_plot, width = 10, height = 24)

