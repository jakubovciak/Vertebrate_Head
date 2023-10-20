#!/bin/sh

#### Markos et al 2023
#### Script used for demonstration of annotation workflow

# This script only runs timepoint annotation workflow, which outputs raw pdf files used for main figure #1.
# Note that orientation of UMAP projections can change slightly due to randomness of the embedding process.

Rscript AmphiG4stage.R
Rscript AmphiN0stage.R
Rscript AmphiN2stage.R
Rscript AmphiN5stage.R
