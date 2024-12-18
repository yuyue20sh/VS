#!/bin/bash

source activate r433
Rscript 00_seu2mtx.R
source activate sc
python 01_create_adata.py
python 02_qc.py
python 03_concatenate.py
python 04a_h5ad2mtx.py
source activate r433
Rscript 04b_integrate.R
source activate sc

echo "Done!"
