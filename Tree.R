# Load necessary libraries
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(igraph)
library(tidyr)
library(dplyr)

# Read VCF file and population data
vcf <- read.vcfR("guavagwas.vcf")
pop.data <- read.table("pop", sep = "\t", header = TRUE)

# Print the list of sample names from the VCF file
sample_names <- colnames(vcf@gt)[-1]
print(sample_names)

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

# Ensure the sample IDs in the VCF match those in the population data
all(sample_names == pop.data$AccessID)

# Convert VCF to genlight object
gl <- vcfR2genlight(vcf)
ploidy(gl) <- 2

# Create a group vector based on the sample names

group_vector[sample_names %in% group1] <- "Group 1"
group_vector[sample_names %in% group2] <- "Group 2"
group_vector[sample_names %in% group3] <- "Group 3"
group_vector[sample_names %in% group4] <- "Group 4"

# Add group information to the genlight object
pop(gl) <- group_vector

# Calculate distance and create a tree
x.dist <- poppr::bitwise.dist(gl)
tree <- aboot(gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

# Perform PCA
pca <- glPca(gl, nf = 3)
barplot(100 * pca$eig / sum(pca$eig), col = heat.colors(0), main = "PCA Eigenvalues")
title(ylab = "Percent of variance\nexplained", line = 2)
title(xlab = "Eigenvalues", line = 1)

# Prepare PCA scores for plotting
pca.scores <- as.data.frame(pca$scores)
pca.scores$pop <- pop(gl)

# Define colors for groups
cols <- c("Group 1" = "red", "Group 2" = "lightpink", "Group 3" = "green", "Group 4" = "purple")

# Set the output file to save the plot as a PNG
png("phylogenetic_tree.png", width = 1200, height = 800, res = 150)

# Plot phylogenetic tree with improved spacing and appearance
par(mar = c(5, 4, 4, 15) + 0.1)  # Adjust margins to give more space for labels
plot.phylo(tree, cex = 0.2, font = 1, adj = 0, tip.color = cols[group_vector], label.offset = 0.02, srt = 0)
nodelabels(tree$node.label, adj = c(1.6, -0.5), frame = "n", cex = 0.3, font = 3, xpd = TRUE)
legend('topleft', legend = names(cols), fill = cols, border = FALSE, bty = "n", cex = 1.2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

# Close the device to finalize the plot
dev.off()

# Create Minimum Spanning Network (MSN)
dist <- bitwise.dist(gl)
msn <- poppr.msn(gl, dist, showplot = FALSE, include.ties = TRUE)
node.size <- rep(2, times = nInd(gl))
names(node.size) <- indNames(gl)
vertex.attributes(msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(gl, msn, palette = cols, gadj = 70)

# Plot PCA results using ggplot2
p1 <- ggplot(pca.scores, aes(x = PC1, y = PC2, colour = pop)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(title = "PCA - PC1 vs PC2")

p2 <- ggplot(pca.scores, aes(x = PC1, y = PC3, colour = pop)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(title = "PCA - PC1 vs PC3")

p3 <- ggplot(pca.scores, aes(x = PC2, y = PC3, colour = pop)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(title = "PCA - PC2 vs PC3")

# Display PCA plots
print(p1)
print(p2)
print(p3)

# Perform Discriminant Analysis of Principal Components (DAPC)
pnw.dapc <- dapc(gl, n.pca = 3, n.da = 2)
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

# Create compoplot for DAPC results
compoplot(pnw.dapc, col = cols, posi = 'top')

# Prepare DAPC results for plotting
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)

