nextflow.enable.dsl=2

params.base_dir = "/Users/hannah/Tech/Read_Lab/ISL_Data_Processing"
params.ggcaller_dir = ""
params.prokka_input = "${params.base_dir}/prokkafile/input.txt"
params.metadata_file = "${params.base_dir}/metadata.csv"
params.out_dir = "${params.base_dir}/results"
params.scripts_dir = "${projectDir}" 
params.run_ggcaller = true



process run_ggcaller {
    conda "environment.yml"
    publishDir "${params.out_dir}/ggcaller", mode: 'copy'
    
    input:
    path prokka_input
    
    output:
    path "ggcaller_out", emit: results
    
    script:
    """
    mkdir -p ggcaller_out
    ggcaller --refs ${prokka_input} --out ggcaller_out/ --save
    """
}

process generate_heatmaps {
    conda "environment.yml"
    publishDir "${params.out_dir}/heatmaps", mode: 'copy'
    
    input:
    path ggcaller_dir
    
    output:
    path "heatmap_results", emit: data
    
        script:
    """
    python ${params.scripts_dir}/ggCaller_Heatmaps.py \
        --input ${ggcaller_dir}/gene_presence_absence.csv \
        --tree ${ggcaller_dir}/pangenome_NJ.nwk \
        --output heatmap_results/
    """
}

process visualize_tree {
    conda "environment.yml"
    publishDir "${params.out_dir}/visualizations", mode: 'copy'
    
    input:
    path heatmap_data
    path ggcaller_dir
    path metadata_file
    
    output:
    path "*.png", emit: plots
    
    script:
    """
    Rscript ${params.scripts_dir}/Tree-and-Heatmap.R \
        --gene_data ${heatmap_data}/filtered_genes.tsv \
        --tree ${ggcaller_dir}/pangenome_NJ.nwk \
        --metadata ${metadata_file} \
        --output final_plot.png
    """
}

workflow {
    main_input = params.run_ggcaller 
        ? run_ggcaller(params.prokka_input) 
        : Channel.fromPath(params.ggcaller_dir, checkIfExists: true)
    
    heatmap_ch = main_input | generate_heatmaps
    
    visualize_tree(heatmap_ch, main_input, params.metadata_file)
}