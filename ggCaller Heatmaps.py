import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
# Load the gene presence/absence file
df = pd.read_csv("/Users/hannah/Tech/Read_Lab/gene_presence_absence.csv")

# Load metadata
ISLdata = pd.read_csv("/Users/hannah/Tech/Read_Lab/S0015_ISL_25_1.csv", sep=" ")
ISLdata[['Date_Collected', 'OS', 'Cn1', 'Cn2', 'CL']] = ISLdata['label'].str.split('-', expand=True)


# Identify isolate columns
isolate_columns = df.columns[3:]

# Create binary presence/absence matrix
presence_absence = df[isolate_columns].notna().astype(int)
presence_absence.insert(0, "Gene", df["Gene"])

# Create custom colormap
custom_cmap = ListedColormap(["#e69597", "#a0ced9"])

# Select data for heatmap
data_only = presence_absence.iloc[:, 1:]

# Filter to keep only genes with at least one absence (0)
genes_with_absence = presence_absence[data_only.apply(lambda x: (x == 0).any(), axis=1)]

# Check if we found any genes with absence
if genes_with_absence.empty:
    print("No genes with absence found. Using all genes instead.")
    genes_to_plot = presence_absence
else:
    print(f"Found {len(genes_with_absence)} genes with absence patterns.")
    genes_to_plot = genes_with_absence

# Create isolate time mapping using the original format (with periods)
isolate_time_mapping = {}
for index, row in ISLdata.iterrows():
    sample_id = row['SampleID']
    # Keep the original format with periods to match column headers
    isolate_time_mapping[sample_id] = row['Date_Collected']

# Get the isolate columns from your presence_absence matrix
isolate_columns = genes_to_plot.columns[1:]  # Skip the Gene column

# Create new column names with time point suffix
new_column_names = {}
for col in isolate_columns:
    # Extract the sample ID part from the column name
    time_point = isolate_time_mapping.get(col, 'Unknown')
    new_column_names[col] = f"{col}_{time_point}"

# Rename the columns in your dataframe
genes_with_absence_renamed = genes_to_plot.iloc[:, 1:].rename(columns=new_column_names)

# Define time order
time_order = {'III': 0, 'En': 1, '3M': 2, '6M': 3, '9M': 4, 'Unknown': 5}

# Sort columns by time point
sorted_columns = sorted(new_column_names.values(),
                       key=lambda x: time_order.get(x.split('_')[-1], 6))

# Use only columns that exist in the dataframe
valid_columns = [col for col in sorted_columns if col in genes_with_absence_renamed.columns]

# Create heatmap with renamed isolates
plt.figure(figsize=(20, 20))
sns.heatmap(genes_with_absence_renamed[valid_columns], cmap=custom_cmap, cbar=True,
            yticklabels=genes_to_plot['Gene'])

# Formatting
plt.xlabel("Isolates")
plt.ylabel("Genes")
plt.title("Binary Gene Presence/Absence Heatmap-GG Caller)")
plt.xticks(rotation=90)  # Make sure labels are readable
plt.tight_layout()
plt.savefig("/Users/hannah/Tech/Read_Lab/ggcaller_heatmap.png", dpi=300)
plt.show()
