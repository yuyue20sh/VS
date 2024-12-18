### Env

```shell
conda create -n sc python=3.12
conda activate sc
pip install 'scanpy[leiden]'

conda create -n r433 python=3.10
conda activate r433
conda install conda-forge::r-base=4.3.3
conda install conda-forge::r-seurat
```

### Run

You can run [run.sh](./run.sh) to reproduce the results.

Or you can also run the following commands one by one.

```shell
conda activate r433
Rscript 00_seu2mtx.R
conda activate sc
python 01_create_adata.py
python 02_qc.py
python 03_concatenate.py
python 04a_h5ad2mtx.py
conda activate r433
Rscript 04b_integrate.R
conda activate sc

```

### Commit

```
commit bf1b9e9b1c1763a7e968368961134bfa9975152d (HEAD -> main, origin/main)
Author: yuyue <yuyue20sh@163.com>
Date:   Sun Dec 15 20:32:56 2024 +0800

    add create adata
```

next commit message: add integrate
