# Load required packages
library(ggplot2)
library(dplyr)

# Set working directory (if necessary)
setwd("H:/GuavaGWAS/new/LDdecay.stat")

# Read the PopLDdecay output file without headers to manually assign correct column names
ld_data <- read.table("LDdecay.stat", header = FALSE, sep = "\t")

# Assign correct column names
colnames(ld_data) <- c("Dist", "Mean_r^2", "Mean_D'", "Sum_r^2", "Sum_D'", "NumberPairs")

# Display the structure of the data
str(ld_data)

# Display the first few rows to verify the column names
head(ld_data)

# Plotting the LD decay using Dist and Mean_r^2 columns
ld_plot <- ggplot(ld_data, aes(x = Dist, y = `Mean_r^2`)) +
  geom_point(alpha = 0.6, color = "pink", size = 1) +  # Points with some transparency and color
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 1) +  # Smooth line to show trend
  labs(title = "Linkage Disequilibrium Decay",
       subtitle = "Decay of linkage disequilibrium (LD) with physical distance",
       x = "Distance (bp)",
       y = expression(LD~(r^2)),
       caption = "Data source: PopLDdecay output") +  # Adding a title, subtitle, and axis labels
  theme_minimal(base_size = 15) +  # A clean theme with larger base font size for better readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Save the plot as a PNG file
ggsave("LDDecayPlot.png", ld_plot, width = 10, height = 6)

# Print the plot
print(ld_plot)
