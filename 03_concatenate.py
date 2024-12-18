"""Concatenate adatas of each sample.
"""


import scanpy as sc
from tqdm import tqdm

from pathlib import Path
import warnings


if __name__ == '__main__':
    
    data_dir = './outs/adata/qc/'
    out_dir = './outs/adata/'
    
    # load adatas
    adatas = {}
    for f in tqdm(Path(data_dir).glob('./*.h5ad'), total=len(list(Path(data_dir).glob('./*.h5ad')))):
        id = f.stem
        adatas[id] = sc.read_h5ad(f)

    # concatenate
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.concat(adatas, label='sample', merge='same')
    adata.obs_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=10)

    # save
    adata.write('%s/concatenated.h5ad' % out_dir)

    print('Done!')
