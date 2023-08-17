import scanpy as sc

def calculate_pca(adata, max_scaling_value = 10, n_comps = 100):
    '''
    Calculates the PCA for downstream analyses. This function is meant to avoid scaling
    gene expression on the original adata.
    
    adata = calculate_pca(adata, max_scaling_value = 10, n_comps = 100)
    
    Parameters:
    Output:
    '''
    # -- subset adata to highly variable genes
    bdata = adata[:, adata.var.highly_variable]
    
    # -- scale bdata and calculate PCA
    sc.pp.scale(bdata,
                max_value = max_scaling_value)

    sc.tl.pca(bdata,
            n_comps = n_comps,
            svd_solver = 'arpack')

    # -- fill NaNs with False so that subsetting to HVGs is possible
    adata.var['highly_variable'].fillna(value = False,
                                        inplace = True)

    # -- transfer PCA outputs from bdata to adata
    adata.obsm['X_pca'] = bdata.obsm['X_pca'].copy()
    adata.uns['pca'] = bdata.uns['pca'].copy()

    adata.varm['PCs'] = np.zeros(shape = (adata.n_vars,
                                          n_comps))
    adata.varm['PCs'][adata.var['highly_variable']] = bdata.varm['PCs']

    # -- print results
    sc.pl.pca_variance_ratio(adata,
                            n_pcs = n_comps,
                            log = True)
    
    return(adata)
    