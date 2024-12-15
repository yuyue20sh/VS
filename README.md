### Env

```shell
conda create -n sc python=3.12
conda activate sc
pip install 'scanpy[leiden]'
pip install harmonypy

conda create -n r433 python=3.10
conda activate r433
conda install conda-forge::r-base=4.3.3
conda install conda-forge::r-seurat
```

### Run

```shell
conda activate r433
Rscript 00_seu2mtx.R
conda activate sc
python 01_create_adata.py
```
