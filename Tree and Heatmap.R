#!/usr/bin/env Rscript

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(argparse)
  library(ggtree)
  library(ape)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggnewscale)
  library(RColorBrewer)
  library(cowplot)
})

# ---- Parse arguments ----
parser <- ArgumentParser()
parser$add_argument("--gene_data", required=TRUE, help="Path to filtered genes TSV")
parser$add_argument("--tree", required=TRUE, help="Path to Newick tree file")
parser$add_argument("--metadata", required=TRUE, help="Path to metadata CSV")
parser$add_argument("--output", required=TRUE, help="Path to output PNG")
args <- parser$parse_args()

# ---- Load data ----
filtered_genes <- read.csv(args$gene_data, sep="\t")
tree <- read.tree(args$tree)
metadata <- read.csv(args$metadata)

# ---- Parse metadata ----
metadata <- metadata %>%
  separate(
    label, 
    into = c("date_collected", "location"), 
    sep = "-", 
    extra = "merge"  # Handle remaining fields
  )
if(!file.exists(args$metadata)) {
  stop(paste("Metadata file not found:", args$metadata))
}
  
  # ---- Clean tree tip labels ----
  tree$tip.label <- gsub("'", "", tree$tip.label)
  tree$tip.label <- gsub("_[^_]+$", "", tree$tip.label)
  
  # ---- Prepare heatmap data ----
  heatmap_data <- filtered_genes %>% select(-Gene)
  heatmap_data <- as.data.frame(lapply(heatmap_data, function(x) factor(x, levels = c(0, 1))))
  colnames(heatmap_data) <- gsub("_[^_]+$", "", colnames(heatmap_data))
  rownames(heatmap_data) <- filtered_genes$Gene
  
  # ---- Reorder heatmap to match tree tip order ----
  heatmap_data <- heatmap_data[, tree$tip.label]
  heatmap_data_transposed <- t(heatmap_data)
  
  # ---- Merge heatmap data with metadata ----
  heatmap_data_with_meta <- heatmap_data_transposed %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    left_join(metadata, by = "SampleID") %>%
    column_to_rownames("SampleID")
  
  heatmap_data_with_meta$location <- factor(heatmap_data_with_meta$location)
  heatmap_data_with_meta$date_collected <- factor(
    heatmap_data_with_meta$date_collected,
    levels = unique(heatmap_data_with_meta$date_collected)
  )
  
  # ---- Color vectors ----
  location_colors <- setNames(
    c("#FFB3BA", "#FFDFBA", "#FFBAED", "#BAFFC9", "#BAE1FF", "#D5BAFF", "#FFFFBA", "#Bafffa"),
    levels(heatmap_data_with_meta$location)
  )
  date_collected_colors <- setNames(
    c("#FFD6E0", "#FFB7B2", "#FFDAC1", "#E2F0CB", "#B5ead7", "#C7CEEA", "#B5B2FF", "#A0E7E5", "#FFB5E8"),
    levels(heatmap_data_with_meta$date_collected)
  )
  presence_absence_colors <- c("#E0E0E0", "#A1CAF1")
  
  # ---- Plot sizing ----
  num_genes <- ncol(heatmap_data_transposed)
  if (num_genes <= 50) {
    presence_width <- 1.2; plot_width <- 20; plot_height <- 25; tip_label_size <- 8; font_size <- 6
  } else if (num_genes <= 150) {
    presence_width <- 1.8; plot_width <- 28; plot_height <- 30; tip_label_size <- 9; font_size <- 7
  } else {
    presence_width <- 2.5; plot_width <- 36; plot_height <- 40; tip_label_size <- 10; font_size <- 8
  }
  
  # ---- Tree plot ----
  tree_plot <- ggtree(tree, color = "black", size = 1.5) +
    geom_tiplab(size = tip_label_size, align = TRUE, offset = 0.1, fontface = "bold") +
    theme_tree2() +
    xlab("Branch length") +
    theme(axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 22),
          plot.margin = margin(0, 10, 0, 100)) +
    coord_cartesian(clip = "off")
  
  # ---- Location annotation heatmap ----
  location_df <- heatmap_data_with_meta["location"] %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = tree$tip.label)) %>%
    pivot_longer(cols = -SampleID)
  
  location_heatmap <- ggplot(location_df, aes(x = name, y = SampleID, fill = value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = location_colors, name = "Location", na.value = "grey") +
    labs(x = "Location") +
    theme_void() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18, face = "bold"),
          plot.margin = margin(0, 0, 0, 0),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))
  
  # ---- Date collected annotation heatmap ----
  date_df <- heatmap_data_with_meta["date_collected"] %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = tree$tip.label)) %>%
    pivot_longer(cols = -SampleID)
  
  date_heatmap <- ggplot(date_df, aes(x = name, y = SampleID, fill = value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = date_collected_colors, name = "Date Collected", na.value = "grey") +
    labs(x = "Date Collected") +
    theme_void() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18, face = "bold"),
          plot.margin = margin(0, 0, 0, 0),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))
  
  # ---- Gene presence/absence heatmap ----
  gene_df <- heatmap_data_transposed %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = factor(SampleID, levels = tree$tip.label)) %>%
    pivot_longer(cols = -SampleID)
  
  gene_heatmap <- ggplot(gene_df, aes(x = name, y = SampleID, fill = value)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_manual(values = presence_absence_colors, name = "Gene Presence", labels = c("Absent", "Present")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18, face = "bold"),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(1.2, "cm"))
  
  # ---- Legends ----
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
  
  legend_location <- get_legend_35(location_heatmap + theme(legend.position = "bottom"))
  legend_date <- get_legend_35(date_heatmap + theme(legend.position = "bottom"))
  legend_gene <- get_legend_35(gene_heatmap + theme(legend.position = "bottom"))
  
  all_legends <- plot_grid(legend_location, legend_date, legend_gene, nrow = 1, align = "h")
  
  # ---- Remove legends for main plots ----
  location_heatmap_nolegend <- location_heatmap + theme(legend.position = "none")
  date_heatmap_nolegend <- date_heatmap + theme(legend.position = "none")
  gene_heatmap_nolegend <- gene_heatmap + theme(legend.position = "none")
  
  # ---- Combine plots ----
  spacer <- ggplot() + theme_void()
  combined_plot <- plot_grid(
    tree_plot, spacer, location_heatmap_nolegend, date_heatmap_nolegend, gene_heatmap_nolegend,
    nrow = 1, align = "h", rel_widths = c(3, 0.5, 0.2, 0.2, 4), axis = "tb"
  )
  
  final_plot <- plot_grid(
    ggdraw() + draw_label("Gene Presence/Absence Heatmap", fontface = "bold", size = 45, hjust = 0.5),
    combined_plot,
    all_legends,
    ncol = 1,
    rel_heights = c(0.08, 1, 0.12)
  )
  
  # ---- Save plot ----
  ggsave(args$output, final_plot,
         width = plot_width * 1.5,
         height = plot_height,
         dpi = 300,
         limitsize = FALSE,
         bg = "white")
}