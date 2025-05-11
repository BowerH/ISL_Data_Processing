#install.packages(c( "ggtree", "ape","dplyr", "tibble", "tidyr", "ggplot2", "ggnewscale", "RColorBrewer"))
#BiocManager::install("ggtree")
# Load required libraries
install.packages("cowplot")
library(ggtree)
library(ape)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(cowplot)

# List all *_ISL_25_1.csv files
isl_csv_files <- list.files(
  path = base_path,
  pattern = "_ISL_25_1\\.csv$",
  full.names = TRUE
)

for (csv_file in isl_csv_files) {
  # Extract ISL name from filename
  isl_name <- gsub("_ISL_25_1\\.csv$", "", basename(csv_file))
  isl_folder <- file.path(base_path, paste0(isl_name, "_ISL_25_1"))
  
  message("Processing ", isl_name)
  
  # Load data
  isl_data <- read.csv(csv_file)
  filtered_genes_file <- file.path(isl_folder, paste0("filtered_", isl_name, "_ISL_25_1_assemblies_genes.tsv"))
  tree_file <- file.path(base_path, "GGCaller", paste0("ggcaller_out_", isl_name, "_ISL_25_1_assemblies"), "pangenome_NJ.nwk")
  
  filtered_genes <- read.csv(filtered_genes_file, sep = "\t")
  tree <- read.tree(tree_file)
  
  # Process ISL data
  isl_data_split <- isl_data %>%
    separate(label, into = c("date_collected", "location", "unknown1", "unknown2"), 
             sep = "-|\\s", extra = "drop", fill = "right")
  
  # Clean tree tip labels
  tree$tip.label <- gsub("'", "", tree$tip.label)
  tree$tip.label <- gsub("_[^_]+$", "", tree$tip.label)
  
  # Prepare heatmap data (convert to binary presence/absence matrix)
  heatmap_data <- filtered_genes %>%
    select(-Gene)
  
  heatmap_data <- as.data.frame(lapply(heatmap_data, function(x) {
    factor(x, levels = c(0, 1))
  }))
  
  # Strip suffixes from column names so they match tree tip labels
  colnames(heatmap_data) <- gsub("_[^_]+$", "", colnames(heatmap_data))
  
  # Set gene names as row names
  rownames(heatmap_data) <- filtered_genes$Gene
  
  # Reorder heatmap to match tree tip order
  heatmap_data <- heatmap_data[, tree$tip.label]
  heatmap_data_transposed <- t(heatmap_data)
  
  # Merge heatmap data with ISL data
  heatmap_data_with_location <- heatmap_data_transposed %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    left_join(isl_data_split, by = "SampleID") %>%
    column_to_rownames("SampleID")
  
  heatmap_data_with_location$location <- factor(heatmap_data_with_location$location)
  heatmap_data_with_location$date_collected <- factor(
    heatmap_data_with_location$date_collected,
    levels = unique(heatmap_data_with_location$date_collected)
  )
  
  # Color vectors
  location_colors <- setNames(
    c("#FFB3BA", # Vibrant pastel pink
      "#FFDFBA", # Vibrant pastel orange
      "#FFBAED", # Vibrant pastel magenta
      "#BAFFC9", # Vibrant pastel green
      "#BAE1FF", # Vibrant pastel blue
      "#D5BAFF", # Vibrant pastel purple
      "#FFFFBA", # Vibrant pastel yellow
      "#Bafffa"  # Vibrant pastel cyan
    ),
    levels(heatmap_data_with_location$location)
  )
  date_collected_colors <- setNames(
    c("#FFD6E0", # Vibrant blush pink
      "#FFB7B2", # Vibrant coral
      "#FFDAC1", # Vibrant apricot
      "#E2F0CB", # Vibrant mint
      "#B5ead7", # Vibrant seafoam
      "#C7CEEA", # Vibrant periwinkle
      "#B5B2FF", # Vibrant lavender
      "#A0E7E5", # Vibrant aqua
      "#FFB5E8"  # Vibrant baby pink
    ),
    levels(heatmap_data_with_location$date_collected)
  )
  
  num_genes <- ncol(heatmap_data_transposed)
  
  # Sizing logic based on number of genes
  if (num_genes <= 50) {
    presence_width <- 1.2
    plot_width <- 20
    plot_height <- 25
    tip_label_size <- 8
    font_size <- 6
  } else if (num_genes <= 150) {
    presence_width <- 1.8
    plot_width <- 28
    plot_height <- 30
    tip_label_size <- 9
    font_size <- 7
  } else {
    presence_width <- 2.5
    plot_width <- 36
    plot_height <- 40
    tip_label_size <- 10
    font_size <- 8
  }
  
  presence_absence_colors <- c("#E0E0E0", "#A1CAF1")
  
  # ---------------------------------------------------------------
  # Step 1: Create individual plot components
  # ---------------------------------------------------------------
  
  # A. Create the tree plot
  tree_plot <- ggtree(tree, color = "black", size = 1.5) +
    geom_tiplab(size = tip_label_size, align = TRUE, offset = 0.1, fontface = "bold") +
    theme_tree2() +
    xlab("Branch length") +
    theme(
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size = 22)
    ) +
    theme(plot.margin = margin(0, 10, 0, 100)) +
    coord_cartesian(clip = "off")
  # B. Create location annotation heatmap
  location_df <- heatmap_data_with_location["location"] %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = tree$tip.label)) %>%
    pivot_longer(cols = -SampleID)
  
  location_heatmap <- ggplot(location_df, aes(x = name, y = SampleID, fill = value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = location_colors,
                      name = "Location",         # Custom legend title
                      na.value = "grey" ) +
    labs(x = "Location") +
    theme_void() +
    theme(
      axis.text.x = element_blank(),  # <--- remove the second label
      axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18, face = "bold"),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14)
    )
  
  # C. Create date collected annotation heatmap
  date_df <- heatmap_data_with_location["date_collected"] %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = tree$tip.label)) %>%
    pivot_longer(cols = -SampleID)
  
  date_heatmap <- ggplot(date_df, aes(x = name, y = SampleID, fill = value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = date_collected_colors,
                      name = "Date Collected",   # Custom legend title
                      na.value = "grey") +
    labs(x = "Date Collected") +
    theme_void() +
    theme(
      axis.text.x = element_blank(),  # <--- remove the second label
      axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18, face = "bold"),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14)
    )
  
  # D. Create gene presence/absence heatmap
  gene_df <- heatmap_data_transposed %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = tree$tip.label)) %>%
    pivot_longer(cols = -SampleID)
  
  gene_heatmap <- ggplot(gene_df, aes(x = name, y = SampleID, fill = value)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_manual(values = presence_absence_colors,
                      name = "Gene Presence",    # Custom legend title
                      labels = c("Absent", "Present")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 18, face = "bold"),  # <-- This makes gene labels bold
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      legend.key.size = unit(1.2, "cm")
    )
  
  # 1. Create plots with legends ON for extraction
location_heatmap_legend <- location_heatmap + theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.key.size = unit(2, "cm"),         # Make legend keys larger
  legend.key.width = unit(2, "cm"),        # Optional: make wider
  legend.key.height = unit(1, "cm"),       # Optional: make taller
  legend.title = element_text(size = 24),  # Bigger legend title
  legend.text  = element_text(size = 20),  # Bigger legend text
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.background = element_rect(fill = "white", color = NA)
)
date_heatmap_legend     <- date_heatmap     + theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.key.size = unit(2, "cm"),         # Make legend keys larger
  legend.key.width = unit(2, "cm"),        # Optional: make wider
  legend.key.height = unit(1, "cm"),       # Optional: make taller
  legend.title = element_text(size = 24),  # Bigger legend title
  legend.text  = element_text(size = 20),  # Bigger legend text
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.background = element_rect(fill = "white", color = NA)
)
gene_heatmap_legend     <- gene_heatmap     + theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.key.size = unit(2, "cm"),         # Make legend keys larger
  legend.key.width = unit(2, "cm"),        # Optional: make wider
  legend.key.height = unit(1, "cm"),       # Optional: make taller
  legend.title = element_text(size = 24),  # Bigger legend title
  legend.text  = element_text(size = 20),  # Bigger legend text
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.background = element_rect(fill = "white", color = NA)
)

# 2. Extract legends
get_legend_35 <- function(plot) {
  legends <- cowplot:::get_plot_component(plot, "guide-box", return_all = TRUE)
  nonzero <- vapply(legends, function(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  if (length(idx) > 0) {
    return(legends[[idx]])
  } else {
    return(legends[[1]])
  }
}
legend_location <- get_legend_35(location_heatmap_legend)
legend_date     <- get_legend_35(date_heatmap_legend)
legend_gene     <- get_legend_35(gene_heatmap_legend)


all_legends <- plot_grid(
  legend_location,
  legend_date,
  legend_gene,
  nrow = 5,
  align = "h"
)

# 3. Use legend-less plots in the main grid
location_heatmap_nolegend <- location_heatmap + theme(legend.position = "none")
date_heatmap_nolegend     <- date_heatmap     + theme(legend.position = "none")
gene_heatmap_nolegend     <- gene_heatmap     + theme(legend.position = "none")

# 4. Combine plots and legends
combined_plot <- plot_grid(
  tree_plot, spacer, location_heatmap_nolegend, date_heatmap_nolegend, gene_heatmap_nolegend,
  nrow = 1, align = "h", rel_widths = c(3, 0.5, 0.2, 0.2, 4), axis = "tb"
)
all_legends <- plot_grid(legend_location, legend_date, legend_gene, nrow = 1, align = "h")

  
  # ---------------------------------------------------------------
  # Step 3: Add legends and finalize
  # ---------------------------------------------------------------

final_plot <- plot_grid(
  ggdraw() + draw_label("Gene Presence/Absence Heatmap", fontface = "bold", size = 45, hjust = 0.5),
  combined_plot,
  all_legends,
  ncol = 1,
  rel_heights = c(0.08, 1, 0.12)
)
  # ---------------------------------------------------------------
  # Step 4: Save the plot
  # ---------------------------------------------------------------
  output_path <- file.path(isl_folder, "gene_heatmap_tree.png")
  ggsave(output_path, final_plot,
         width = plot_width * 1.5,
         height = plot_height,
         dpi = 300,
         limitsize = FALSE,
         bg = "white") 
}