# Load necessary libraries
library(dplyr)
library(readr)

# Load the data
fst_g1_vs_g2 <- read_csv("positive_selection_fst_g1_vs_g2.csv")
annotation <- read_csv("annotation.csv")

# Standardize chromosome identifiers
fst_g1_vs_g2 <- fst_g1_vs_g2 %>% mutate(CHROM = toupper(CHROM))
annotation <- annotation %>% mutate(Pseudochromosome = toupper(Pseudochromosome))

# Define the upstream region in base pairs
upstream_region <- 1000

# Function to find meaningful loci
find_meaningful_loci <- function(fst_df, annotation_df, upstream_region) {
  meaningful_loci <- list()
  for (i in 1:nrow(fst_df)) {
    chrom <- fst_df$CHROM[i]
    pos <- fst_df$POS[i]
    gene_rows <- annotation_df %>%
      filter(Pseudochromosome == chrom & Start - upstream_region <= pos & pos <= End)
    if (nrow(gene_rows) > 0) {
      for (j in 1:nrow(gene_rows)) {
        meaningful_loci <- append(meaningful_loci, list(c(
          chrom,
          pos,
          fst_df$WEIR_AND_COCKERHAM_FST[i],
          gene_rows$`Gene ID`[j],
          gene_rows$Description[j]
        )))
      }
    }
  }
  meaningful_loci_df <- do.call(rbind, meaningful_loci) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("CHROM", "POS", "WEIR_AND_COCKERHAM_FST", "Gene ID", "Description"))
  return(meaningful_loci_df)
}

# Find meaningful loci for g1 vs g2 comparison
meaningful_loci_g1_vs_g2 <- find_meaningful_loci(fst_g1_vs_g2, annotation, upstream_region)

# Save the result to a CSV file
write_csv(meaningful_loci_g1_vs_g2, "meaningful_loci_g1_vs_g2.csv")

print(head(meaningful_loci_g1_vs_g2))

# Repeat the process for g1_vs_g3 and g2_vs_g3 comparisons
fst_g1_vs_g3 <- read_csv("positive_selection_fst_g1_vs_g3.csv")
fst_g2_vs_g3 <- read_csv("positive_selection_fst_g2_vs_g3.csv")

meaningful_loci_g1_vs_g3 <- find_meaningful_loci(fst_g1_vs_g3, annotation, upstream_region)
meaningful_loci_g2_vs_g3 <- find_meaningful_loci(fst_g2_vs_g3, annotation, upstream_region)

write_csv(meaningful_loci_g1_vs_g3, "meaningful_loci_g1_vs_g3.csv")
write_csv(meaningful_loci_g2_vs_g3, "meaningful_loci_g2_vs_g3.csv")

print(head(meaningful_loci_g1_vs_g3))
print(head(meaningful_loci_g2_vs_g3))





# Load necessary libraries
library(dplyr)
library(tidyr)

# Read the input files
annotation <- read.csv("annotation.csv", stringsAsFactors = FALSE)
positive_tajima <- read.csv("positive_selection_tajima_d.csv", stringsAsFactors = FALSE)
negative_tajima <- read.csv("negative_selection_tajima_d.csv", stringsAsFactors = FALSE)
purifying_diversity <- read.csv("purifying_selection_nucleotide_diversity.csv", stringsAsFactors = FALSE)

# Define a function to annotate regions
annotate_regions <- function(data, annotation, upstream_threshold = 1000) {
  annotated <- data %>%
    rowwise() %>%
    mutate(
      annotation_info = list(
        annotation %>%
          filter(
            CHROM == data$CHROM &
              ((POS >= Start & POS <= End) |
                 (POS < Start & POS >= Start - upstream_threshold))
          )
      )
    ) %>%
    unnest(annotation_info) %>%
    select(CHROM, POS, everything())
  
  return(annotated)
}

# Annotate the regions for each selection type
annotated_positive_tajima <- annotate_regions(positive_tajima, annotation)
annotated_negative_tajima <- annotate_regions(negative_tajima, annotation)
annotated_purifying_diversity <- annotate_regions(purifying_diversity, annotation)
annotated_purifying_diversity
# Save the annotated data
write.csv(annotated_positive_tajima, "annotated_positive_tajima_d.csv", row.names = FALSE)
write.csv(annotated_negative_tajima, "annotated_negative_tajima_d.csv", row.names = FALSE)
write.csv(annotated_purifying_diversity, "annotated_purifying_diversity.csv", row.names = FALSE)
