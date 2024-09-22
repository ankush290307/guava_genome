# Load required libraries
library(ggplot2)
library(reshape2)

# Define groups based on your provided list
group1 <- c("AC10-5", "AC10-7", "AC10-8", "AC11-9", "AC1-4", "AC14-4", "AC1-5", "AC18-4", "AC2-1", "AC2-3", "AC3-4",
            "AC5-2", "AC5-5", "AC5-6", "AC5-7", "AC6-2", "AC6-4", "AC8-2", "AC8-4", "AC8-5", "AC8-6", "AC-9",
            "CISH-G1", "CISH-G5", "CISH-G6", "Strawberry_guava")
group2 <- c("50-50", "AICRP10-4", "Arka_kiran", "HB17-16", "HB8-8", "HB9-2", "Hisar_surkha", "Lalit", "Pear_shaped",
            "Portugal", "Punjab_kiran", "Punjab_pink", "Taiwan_pink")
group3 <- c("2y2", "AICRP10-5", "Allahabad_safeda", "Arka_amulya", "AS11-3", "AS12-2", "AS13-2", "AS15-3", "AS1-6",
            "AS1-7", "AS2-8", "AS5-8", "AS7-3", "AS7-7", "AS7-8", "Barf_khana", "BDS", "Behat_coconut", "BP-1",
            "BP-2", "G-Bilas", "GVP", "HB2-6", "HB4-8", "HB5-5", "HB5-8", "HB6-3", "HB9-7", "Hisar_safeda", "KMG",
            "L-49", "Pant_prabhat", "Punjab_safeda", "R2P4", "Sabir", "Sabir-2", "Sabir-3", "Seedless", "Shweta",
            "Thai_guava", "TS-1", "VNR", "Lemon_guava")
group4 <- c("Purple_guava", "Purple_guava_rep2")

# Define colors for groups


# Assign group labels to individuals
assign_group_labels <- function(individuals) {
  group_labels <- rep(NA, length(individuals))
  group_labels[individuals %in% group1] <- "Group 1"
  group_labels[individuals %in% group2] <- "Group 2"
  group_labels[individuals %in% group3] <- "Group 3"
  group_labels[individuals %in% group4] <- "Group 4"
  return(group_labels)
}

# Actual individual names
individual_names <- c("AC10-5", "AC10-7", "AC10-8", "AC11-9", "AC1-4", "AC14-4", "AC1-5", "AC18-4", "AC2-1", "AC2-3", 
                      "AC3-4", "AC5-2", "AC5-5", "AC5-6", "AC5-7", "AC6-2", "AC6-4", "AC8-2", "AC8-4", "AC8-5", 
                      "AC8-6", "AC-9", "CISH-G1", "CISH-G5", "CISH-G6", "Strawberry_guava", "50-50", "AICRP10-4", 
                      "Arka_kiran", "HB17-16", "HB8-8", "HB9-2", "Hisar_surkha", "Lalit", "Pear_shaped", "Portugal", 
                      "Punjab_kiran", "Punjab_pink", "Taiwan_pink", "2y2", "AICRP10-5", "Allahabad_safeda", 
                      "Arka_amulya", "AS11-3", "AS12-2", "AS13-2", "AS15-3", "AS1-6", "AS1-7", "AS2-8", "AS5-8", 
                      "AS7-3", "AS7-7", "AS7-8", "Barf_khana", "BDS", "Behat_coconut", "BP-1", "BP-2", "G-Bilas", 
                      "GVP", "HB2-6", "HB4-8", "HB5-5", "HB5-8", "HB6-3", "HB9-7", "Hisar_safeda", "KMG", "L-49", 
                      "Pant_prabhat", "Punjab_safeda", "R2P4", "Sabir", "Sabir-2", "Sabir-3", "Seedless", "Shweta", 
                      "Thai_guava", "TS-1", "VNR", "Lemon_guava", "Purple_guava", "Purple_guava_rep2")

# Read the Q file for K=5
Q_file <- "guavaplinkbed.5.Q"
if (!file.exists(Q_file)) {
  stop(paste("File", Q_file, "does not exist"))
}
admixture_data <- read.table(Q_file, header = FALSE)
colnames(admixture_data) <- paste0("Cluster", 1:ncol(admixture_data))

# Add individual names and group labels to the data
admixture_data$Individual <- individual_names
admixture_data$Group <- assign_group_labels(individual_names)

# Melt the data frame for ggplot2
admixture_data_melt <- melt(admixture_data, id.vars = c("Individual", "Group"))

# Calculate the number of individuals in each group for proper scaling
group_counts <- table(admixture_data_melt$Group)
group_scale_factors <- c(Group_1 = group_counts["Group 1"] / max(group_counts), 
                         Group_2 = group_counts["Group 2"] / max(group_counts), 
                         Group_3 = group_counts["Group 3"] / max(group_counts), 
                         Group_4 = 1)

# Create the plot
p <- ggplot(admixture_data_melt, aes(x = Individual, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Group, scales = "free_x", nrow = 1, 
             labeller = labeller(Group = c("Group 1" = "Group 1", "Group 2" = "Group 2", "Group 3" = "Group 3", "Group 4" = "Group 4"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, margin = margin(t = 20, r = 5), size = 4),  # Adjust angle, spacing, and font size
        strip.background = element_rect(fill = c("red", "lightpink", "green", "purple")[match(unique(admixture_data_melt$Group), names(cols))], color = NA),
        strip.text = element_text(color = "black", face = "bold"),
        axis.ticks.x = element_blank()) +  # Remove ticks
  labs(x = "Individual", y = "Admixture Proportion", fill = "Ancestry Component") +
  ggtitle("Admixture Proportions for K = 5")

ggsave("admixture_plot_high_res.png", p, width = 10, height = 6, dpi = 300)

