if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
install.packages(c("ape", "dplyr", "tidyr", "ggplot2"))
install.packages("tibble")
install.packages("ggnewscale")
install.packages("RColorBrewer")
# Load required libraries
library(ggtree)
library(ape)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)

# Load data
filtered_genes <- read.csv("/Users/hannah/Tech/Read_Lab/filtered_genes_starting_with_7.csv")
isl_data <- read.csv("/Users/hannah/Tech/Read_Lab/S0015_ISL_25_1.csv", sep = " ")
tree <- read.tree("/Users/hannah/Tech/Read_Lab/pangenome_NJ.nwk") 

# Process ISL data
isl_data_split <- isl_data %>%
  separate(label, into = c("date_collected", "location", "unknown1", "unknown2"), sep = "-|\\s")

# Clean tree tip labels
tree$tip.label <- gsub("'", "", tree$tip.label)

# Prepare heatmap data (convert to binary presence/absence matrix)
heatmap_data <- filtered_genes %>%
  select(-Gene, -`Non.unique.Gene.name`, -Annotation) %>%  
  mutate_all(~ ifelse(!is.na(.) & . != "", 1, 0))         

# Convert heatmap data to factors for discrete color mapping
heatmap_data <- as.data.frame(lapply(heatmap_data, function(x) {
  factor(x, levels = c(0, 1)) # Ensure levels are explicitly defined
}))

# Set row names to gene identifiers
rownames(heatmap_data) <- filtered_genes$Gene

# Verify alignment between tree tips and heatmap columns
missing_labels <- setdiff(tree$tip.label, colnames(heatmap_data))
if (length(missing_labels) > 0) {
  stop(paste("Missing labels in heatmap data:", paste(missing_labels, collapse = ", ")))
}
heatmap_data_transposed <- t(heatmap_data)

# Merge heatmap data with ISL data
heatmap_data_with_location <- heatmap_data_transposed %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(isl_data_split, by = "SampleID") %>%  
  column_to_rownames("SampleID")

heatmap_data_with_location$location <- factor(heatmap_data_with_location$location)
# Ensure date_collected is a factor
heatmap_data_with_location$date_collected <- factor(
  heatmap_data_with_location$date_collected,
  levels = unique(heatmap_data_with_location$date_collected)
)
# Create pastel color vectors for location
location_colors <- setNames(
  c("#a1c9f4", "#ffb482", "#8de5a1", "#ff9f9b", "#d0bbff", 
    "#debb9b", "#fab0e4", "#cfcfcf", "#fffea3"),  # Brighter pastel colors
  levels(heatmap_data_with_location$location)
)

date_collected_colors <- setNames(
  c("#023eff", "#ff7c00", "#1ac938", "#e8000b", "#8b2be2", 
    "#9f4800", "#f14cc1", "#a3a3a3", "#ffc400"),  # Bright colors
  levels(heatmap_data_with_location$date_collected)
)
# For the gene presence/absence heatmap, use softer colors
presence_absence_colors <- c("#E0E0E0", "#A1CAF1")  # Soft gray and pastel blue
# Initialize the tree plot with adjusted tip labels
p <- ggtree(tree, branch.length = "none",color = "black", size = 1.5) + 
  geom_tiplab(size = 10, align = TRUE, offset = 1, fontface = "bold") 

# Add location heatmap (first layer)
p <- gheatmap(
  p, heatmap_data_with_location["location"], 
  offset = 4, width = 0.1,
  colnames_position = "top", color = FALSE,
  colnames_angle = 90,
  colnames_offset_y = 2.3, font.size = 12
) +
  scale_fill_manual(
    values = location_colors,
    name = "Location",
    guide = guide_legend(order = 1, title.position = "top")
  )

# Add new scale fill before the date_collected heatmap
p <- p + ggnewscale::new_scale_fill()

# Add date_collected heatmap (second layer)
p <- gheatmap(
  p, heatmap_data_with_location["date_collected"], 
  offset = 6, width = 0.1,  # Adjusted spacing
  colnames_position = "top", color = FALSE,
  colnames_angle = 90,
  colnames_offset_y = 3.3, font.size = 12
) +
  scale_fill_manual(
    values = date_collected_colors,
    name = "Date Collected",
    guide = guide_legend(order = 2, title.position = "top")
  )

# Add new scale fill before the presence/absence heatmap
p <- p + ggnewscale::new_scale_fill()

# Add presence/absence heatmap (third layer)
p <- gheatmap(
  p, heatmap_data_transposed, 
  offset = 8, width = 2,  # Increased width significantly to match image
  colnames_angle = 90, colnames_position = "top", color = FALSE,
  colnames_offset_y = 2.3, font.size = 12
) +
  scale_fill_manual(
    values = presence_absence_colors,
    name = "Gene Presence/Absence",
    labels = c("Absent", "Present"),
    guide = guide_legend(order = 3, title.position = "top")
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 20, b = 20),
    plot.title = element_text(size = 35, face = "bold", hjust = .5), # Larger title
    legend.text = element_text(size = 30, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 30, face = "bold"), # Larger legend titles
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = margin(12, 10, 12, 12)
  ) +
  ggtitle("Gene Presence/Absence Heatmap") +
  geom_segment(
    data = p$data,
    aes(xend = max(p$data$x), yend = y),
    linetype = "solid", color = "black"
  )
p <- p + 
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), limits = c(0, max(tree$Nnode) * 1.25))

# Save plot to file with larger dimensions to accommodate larger text
ggsave("/Users/hannah/Tech/Read_Lab/gene_heatmap_tree.png", 
       p, width = 49, height = 30, dpi = 400)

# Display plot
print(p)