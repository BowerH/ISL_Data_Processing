#!/usr/bin/env 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import os
import glob
from Bio import Phylo
import argparse


# ===== Helper Functions =====

def additional_gene_filtering(df, save_dir, isl_label):
    gene_data = df.iloc[:, 1:]
    sample_counts = gene_data.shape[1]

    initial_gene_count = df.shape[0]
    filtered = []

    for idx, row in gene_data.iterrows():
        present_count = row.sum()
        presence_ratio = present_count / sample_counts

        if sample_counts > 50:
            filtered.append(0.2 <= presence_ratio <= 0.8)
        elif sample_counts > 30:
            filtered.append(presence_ratio < 0.9)
        else:
            filtered.append(True)

    filtered_df = df[filtered]
    final_gene_count = filtered_df.shape[0]
    removed_count = initial_gene_count - final_gene_count

    # Save summary
    summary_path = os.path.join(save_dir, f"filtered_{isl_label}_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("[Filtering Summary]\n")
        f.write(f"Initial number of genes:\t{initial_gene_count}\n")
        f.write(f"Number of genes after filtering:\t{final_gene_count}\n")
        f.write(f"Number of genes removed:\t{removed_count}\n")
    print(f"Saved summary: {summary_path}")

    return filtered_df

def get_tree_order(tree_path):
    try:
        tree = Phylo.read(tree_path, "newick")
        return [term.name for term in tree.get_terminals()]
    except Exception as e:
        print(f"❌ Could not parse tree at {tree_path}: {e}")
        return []

def plot_and_save(heatmap_df, filename_prefix, new_column_names, sorted_columns, save_dir, cmap):
    heatmap_df = heatmap_df.set_index('Gene').rename(columns=new_column_names)
    valid_columns = [col for col in sorted_columns if col in heatmap_df.columns]
    heatmap_df = heatmap_df[valid_columns]

    # Save TSV
    tsv_path = os.path.join(save_dir, f"{filename_prefix}_genes.tsv")
    heatmap_df.to_csv(tsv_path, sep='\t')
    print(f"Saved TSV: {tsv_path}")

    # Set up figure and axis
    num_genes = heatmap_df.shape[0]
    fig_height = max(8, min(num_genes * 0.2, 60))
    fig, ax = plt.subplots(figsize=(20, fig_height))

    # Plot heatmap
    sns.heatmap(
        heatmap_df,
        cmap=cmap,
        cbar=False,
        yticklabels=True,
        linewidths=0.1,
        linecolor='lightgray',
        ax=ax
    )

    if num_genes > 50:
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)

    ax.set_xlabel("Isolates", fontsize=14)
    ax.set_ylabel("Genes", fontsize=14)
    ax.set_title("Binary Gene Presence/Absence Heatmap (GGCaller)", fontsize=14)
    ax.tick_params(axis='x', rotation=90)

    # Force layout calculation
    # Draw layout before positioning
    fig.tight_layout()
    fig.canvas.draw()

    # Get axis position and compress downward slightly to make room above
    box = ax.get_position()
    new_box = [box.x0, box.y0, box.width, box.height * 0.88]
    ax.set_position(new_box)

    # Add legend just below the title (above the heatmap)
    peach_patch = mpatches.Patch(color='#f4a582', label='Present')
    gray_patch = mpatches.Patch(color='#d9d9d9', label='Absent')

    ax.legend(
        handles=[peach_patch, gray_patch],
        title='Gene Status',
        loc='upper right',
        bbox_to_anchor=(.94, 1.07),  # top-right of the axes
        ncol=2,
        frameon=True,
        edgecolor='black',
        fancybox=False,
        framealpha=1,
        borderpad=0.5
    )

    # Save and close
    png_path = os.path.join(save_dir, f"{filename_prefix}_heatmap.png")
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f"Saved Heatmap: {png_path}")


# ===== Iterate through GGCaller outputs =====
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Path to gene_presence_absence.csv')
parser.add_argument('--output', required=True, help='Output directory path')
parser.add_argument('--tree', required=True, help='Path to pangenome_NJ.nwk')
args = parser.parse_args()

# Remove directory loop and process single input
def main():
    gene_csv = args.input
    save_dir = args.output
    tree_file = args.tree
    
    os.makedirs(save_dir, exist_ok=True)
    
    try:
        df = pd.read_csv(gene_csv)
        # Remove metadata file dependency - handle in Nextflow
        def process_data(df):
            isolate_columns = df.columns[3:]
            presence_absence = df[isolate_columns].notna().astype(int)
            presence_absence.insert(0, "Gene", df["Gene"])
            return presence_absence
        presence_absence = process_data(df)
        plot_and_save(presence_absence, save_dir, tree_file)
        
    except Exception as e:
        print(f"❌ Processing failed: {e}")

if __name__ == "__main__":
    main()