# ISL_Data_Processing


Installs for nextflow:

Added environment.yml script:

```
conda env create -f environment.yml

conda activate PA_analysis_env

#verify

conda list

```

## Currently NOT tested running ggcaller in script.

to run ggcaller script, It should be in a folder with all the assembly files (.fna) including the reference genome. 

``` ./run_ggcaller.sh```

## Running nextflow script

Installs for nextflow:

Added environment.yml script:

```
conda env create -f environment.yml

conda activate PA_analysis_env

#verify

conda list

#Now actually to run it:

nextflow run working_PA_pipeline.nf  --base_dir /Users/hannah/Tech/Read_Lab/ISL_Data_Processing/  --run_ggcaller false  --ggcaller_dir (Directory that contains ggcaller output.  --metadata_file (the file with all of the data collected data, should be given.)


```

