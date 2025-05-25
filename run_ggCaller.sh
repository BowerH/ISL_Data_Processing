conda create -n ggcaller_env python=3.9
conda activate ggcaller_env
conda install ggcaller

#!/bin/bash

# Get a list of all *_assemblies folders
ALL_DIRS=(*_assemblies)

# Loop starting from index 1 to skip the first directory
for ((i=1; i<${#ALL_DIRS[@]}; i++)); do
    ISL_DIR="${ALL_DIRS[$i]}"
    OUTDIR="ggcaller_out_${ISL_DIR}"

    # Skip if already completed
    if [ -d "$OUTDIR" ] && [ -f "$OUTDIR/gene_annotations.gff" ]; then
        echo "⏩ Skipping $ISL_DIR — already processed."
        continue
    fi

    echo "📦 Processing $ISL_DIR..."

    # Step 1: Gunzip all .fna.gz files
    echo "🧬 Unzipping .fna.gz files..."
    gunzip -f "$ISL_DIR"/*.fna.gz

    # Step 2: Create input.txt for this ISL folder
    INPUT_FILE="${ISL_DIR}_input.txt"
    echo "📄 Creating $INPUT_FILE..."
    > "$INPUT_FILE"
    for fna in "$ISL_DIR"/*.fna; do
        realpath "$fna" >> "$INPUT_FILE"
    done

    # Step 3: Run GgCaller
    echo "🚀 Running GgCaller on $ISL_DIR..."
    ggcaller --refs "$INPUT_FILE" --out "$OUTDIR" --save --threads 12

    # Step 4: Clean up gz files (already unzipped)
    echo "🧼 Cleaning up .fna.gz files..."
    rm -f "$ISL_DIR"/*.fna.gz

    echo "✅ Done with $ISL_DIR"
    echo
done
