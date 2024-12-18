"""Read scRNA data and create adata for each sample.
"""


import scanpy as sc
import pandas as pd

from pathlib import Path
import warnings

import utils.io as io


if __name__ == '__main__':
    
    data_dir = './data/'
    mtx_dir = './outs/mtx/ajp/'
    save_dir = './outs/adata/raw/'


    # create directory if not exists
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    # load gse datasets
    for f in Path(data_dir).glob('./*/*'):
        basename = f.name
        ext = f.suffix

        # only read scRNA data
        if 'matrix' not in basename:
            continue
        
        print(f)

        # get sample information
        gsm, sample = basename.split('-')[0].split('_')[:2]
        gse = f.parent.name
        
        if gse == 'GSE216784':
            dataset = 'NC'
        elif gse == 'GSE230375':
            dataset = 'NO'
        elif gse == 'GSE250061':
            dataset = 'BJC'
        else:
            raise ValueError('Unknown dataset: %s for %s' % (gse, f))

        # read data
        if ext == '.gz':
            adata = io.read_mtx(str(f))
        elif ext == '.h5':
            adata = io.read_h5(str(f))
        else:
            raise ValueError('Unknown file extension: %s for %s' % (ext, f))
        
        # set sample information
        adata.obs['gsm'] = gsm
        adata.obs['sample'] = sample
        adata.obs['gse'] = gse
        adata.obs['dataset'] = dataset

        # save
        adata.write('%s/%s.h5ad' % (save_dir, sample))

    # load ajp datasets
    for d in Path(mtx_dir).glob('./*'):
        # read data
        adata = io.read_mtx('%s/matrix.mtx' % d)
        adata.obs = adata.obs.loc[:, ['gsm', 'sample', 'gse', 'dataset']]
        
        # save
        adata.write('%s/%s.h5ad' % (save_dir, adata.obs['sample'].values[0]))

    print('Done!')
