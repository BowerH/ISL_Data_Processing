conda create -n ggcaller_env python=3.9
conda activate ggcaller_env
conda install ggcaller

touch ggCall_out/
##Input file is .fna files from PROKKA. Cat them all to a .txt to run.

ggcaller --refs prokkafile/input.txt --out ggCall_out/ --save
