"""Read scRNA data and create adata for each sample
"""


import scanpy as sc
import pandas as pd

from pathlib import Path
import warnings


def read_mtx(mtx_file):
    """
    Read mtx format data

    Args:
        mtx_file: str, file path

    Returns:
        adata: AnnData, scRNA data

    """
    adata = sc.read_mtx(mtx_file).T
    
    # get file names
    filename = Path(mtx_file).name
    parent_dir = Path(mtx_file).parent
    barcodes_file = '%s/%s' % (parent_dir, filename.replace('matrix', 'barcodes').replace('mtx', 'tsv'))
    features_file = '%s/%s' % (parent_dir, filename.replace('matrix', 'features').replace('mtx', 'tsv'))
    metadata_file = '%s/%s' % (parent_dir, filename.replace('matrix', 'metadata').replace('mtx', 'tsv'))

    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
    features = pd.read_csv(features_file, sep='\t', header=None)

    # add metadata
    adata.obs_names = barcodes[0]
    if features.shape[1] > 1:
        adata.var_names = features[1]
        adata.var['ens_id'] = features[0].to_list()
        adata.var['symbol'] = features[1].to_list()
    else:
        adata.var_names = features[0]
        adata.var['ens_id'] = '-'
        adata.var['symbol'] = features[0].to_list()
    
    if Path(metadata_file).exists():
        metadata = pd.read_csv(metadata_file, sep='\t', index_col=0, header=0)
        adata.obs = metadata.loc[adata.obs_names]

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    return adata


def read_h5(h5_file):
    """
    Read h5 format data

    Args:
        h5_file: str, file path

    Returns:
        adata: AnnData, scRNA data

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.read_10x_h5(h5_file)

    # select columns
    adata.var = adata.var.iloc[:, 0].to_frame()
    adata.var.columns = ['ens_id']
    adata.var['symbol'] = adata.var_names

    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    
    return adata


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
            adata = read_mtx(str(f))
        elif ext == '.h5':
            adata = read_h5(str(f))
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
        adata = read_mtx('%s/matrix.mtx' % d)
        adata.obs = adata.obs.loc[:, ['gsm', 'sample', 'gse', 'dataset']]
        
        # save
        adata.write('%s/%s.h5ad' % (save_dir, adata.obs['sample'].values[0]))

    print('Done!')
