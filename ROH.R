# Load libraries
library(tidyverse)
library(ggrepel)

#Read data with read_delim() for better control over input file parsing
clouded_roh <- read_delim("guava_all_individuals.roh.edited.txt", delim = "\t", skip = 1, col_names = c("Sample", "Chromosome", "RoH_length"))

#Compute NROH and Sum of ROH (SROH)
clouded_nroh <- clouded_roh %>% 
  group_by(Sample) %>% 
  summarize(NROH = n())

clouded_sroh <- clouded_roh %>% 
  group_by(Sample) %>% 
  summarize(SROH = sum(RoH_length))

#Compute FROH
clouded_froh <- inner_join(clouded_nroh, clouded_sroh, by = "Sample") %>% 
  mutate(FROH = NROH / SROH)

#Create a table with NROH, SROH, and FROH for each sample
summary_table <- clouded_froh

#Display the table
print(summary_table)

#Save the table to a CSV file
write_csv(summary_table, "summary_table.csv")

#Plot NROH vs. SROH and save the plot to a file
froh_plot <- inner_join(clouded_nroh, clouded_sroh, by = "Sample") %>% 
  ggplot(aes(x = SROH, y = NROH)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "NROH vs. SROH")

ggsave("froh_plot.png", froh_plot, width = 8, height = 6, dpi = 300)

# Load libraries
library(tidyverse)
library(ggrepel)

#Read data with read_delim() for better control over input file parsing
clouded_roh <- read_delim("guava_all_individuals.roh.edited.txt", delim = "\t", skip = 1, col_names = c("Sample", "Chromosome", "RoH_length"))

#Compute NROH and Sum of ROH (SROH)
clouded_nroh <- clouded_roh %>% 
  group_by(Sample) %>% 
  summarize(NROH = n())

clouded_sroh <- clouded_roh %>% 
  group_by(Sample) %>% 
  summarize(SROH = sum(RoH_length))

#Compute FROH
clouded_froh <- inner_join(clouded_nroh, clouded_sroh, by = "Sample") %>% 
  mutate(FROH = NROH / SROH)

#Create a table with NROH, SROH, and FROH for each sample
summary_table <- clouded_froh

#Display the table
print(summary_table)

#Save the table to a CSV file
write_csv(summary_table, "summary_table.csv")

#Plot NROH vs. SROH and save the plot to a file
froh_plot <- inner_join(clouded_nroh, clouded_sroh, by = "Sample") %>% 
  ggplot(aes(x = SROH, y = NROH)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "NROH vs. SROH")

ggsave("froh_plot.png", froh_plot, width = 8, height = 6, dpi = 300)

# Create a boxplot of RoH lengths for each sample and save the plot to a file
roh_boxplot <- clouded_roh %>% 
  ggplot(aes(x = Sample, y = RoH_length)) +
  geom_boxplot() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  # Rotate x-axis labels by 90 degrees
  labs(title = "RoH Lengths per Sample")

# Save plot as PNG
ggsave("roh_boxplot.png", roh_boxplot, width = 8, height = 6, dpi = 300)




#Save plot as PNG
ggsave("roh_boxplot.png", roh_boxplot, width = 8, height = 6, dpi = 300)

