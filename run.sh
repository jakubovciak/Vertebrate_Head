#!/bin/bash

set -e

echo "
Started annotation workflow ...
"

Rscript R_demo.R

echo "
Finished annotation workflow, check figures in Results directory
"

cd Transitions

echo "
Started celltype transitions workflow ...
"

Rscript -e "rmarkdown::render(input = '01_Integration.Rmd', output_dir = '../Results', clean = TRUE)"

Rscript 02_Convert_to_anndata.R

conda activate ctk_pseudotime_env

python 03_Cytotrace_pseudotime.py

Rscript -e "rmarkdown::render(input = '04_Urd_transition_matrix.Rmd', output_dir = '../Results', clean = TRUE)"

Rscript -e "rmarkdown::render(input = '05_Transition_graphs.Rmd', output_dir = '../Results', clean = TRUE)"

echo "
Finished celltype transitions workflow, check rendered reports in Results directory
"

cd ..

echo "
Started SAMap workflow ...
"

conda deactivate

cd SAMap

cd input_data/Amp_h5ad
Rscript export_individual_h5ad.R

cd ../../analysis

conda activate SAMap

sh ./run_1to1_comparisons.sh

chown -R rstudio:rstudio *

echo "
Finished SAMap workflow, check rendered reports in individual folders
"

cd ..
