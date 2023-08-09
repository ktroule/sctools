from joblib import Parallel, delayed
import multiprocessing
#import scanpy as sc
#import pandas as pd
import glob
import os

def soupy_ratio(data_dir, num_cores = -1):
    
    '''
    Calculates the ratio of gene counts associated to empty droplets to all droplets.
    
    This aims to identify soupy genes, those whose expression ratio is high in the empty droplets.
    A higher score (i.e close to 1), the more likely the gene is soupy.
    
    Parameters:
        - data_dir (int): Folder containing all the libraries associated to a given experiment.
        - num_cores (int): Number of cores to use, default all unless specified.
            
    Parameters:
        - a dataframe with the ratio of soupiness associated to every gene.
    
    '''
    
    # -- List all subfolders containing raw and filtered matrices
    dirs_dropl = glob.glob(os.path.join(data_dir, '*', 'raw_feature_bc_matrix'))
    dirs_cells = glob.glob(os.path.join(data_dir, '*', 'filtered_feature_bc_matrix'))

    # -- Function to load adata and modify cell barcodes
    def load_adata(obj_dir):
        adata = scanpy.read_10x_mtx(obj_dir)
        adata.obs.index = os.path.basename(obj_dir.split('/')[-2]) + '_' + adata.obs.index
        return adata


    # -- Set number of cores for parallelization
    if num_cores <= 0: 
        num_cores = multiprocessing.cpu_count()

    # -- Read all matrices and concateante into a single object for raw and filtered
    if __name__ == "__main__":
        print('Reading raw_feature_bc_matrix objects')
        dropl_holder = Parallel(n_jobs = num_cores)(delayed(load_adata)(obj_dir = i)
                                                    for i in dirs_dropl)

        print('Reading filtered_feature_bc_matrix objects')
        cells_holder = Parallel(n_jobs = num_cores)(delayed(load_adata)(obj_dir = i)
                                                    for i in dirs_cells)

    # -- Concatenate objets
    adata_dropl = dropl_holder[0].concatenate(dropl_holder[1:], join = 'outer')
    del dropl_holder
    
    adata_cells = cells_holder[0].concatenate(cells_holder[1:], join = 'outer')
    del cells_holder
    
    # -- From droplets remove barcodes associtate to cells
    adata_dropl = adata_dropl[~adata_dropl.obs.index.isin(list(adata_cells.obs.index))]

    # -- Calculate number of gene counts 
    counts_dropl = adata_dropl.X.sum(axis = 0).A1
    counts_cells = adata_cells.X.sum(axis = 0).A1

    # -- Join counts in a dataframe and discard genes whose counts is 0 in both objects
    soupy_ratio = pd.DataFrame({'droplets' : counts_dropl, 'cells' : counts_cells},
                               index = adata_dropl.var.index)
    
    # -- Remove genes 0 counts for droplets and cells
    soupy_ratio = soupy_ratio[soupy_ratio[['droplets', 'cells']].sum(axis = 1) != 0]
    
    # -- Calculate ratio of soupiness
    soupy_ratio['soupy_ratio'] = soupy_ratio['droplets']/soupy_ratio[['droplets', 'cells']].sum(axis = 1)

    return(soupy_ratio)
