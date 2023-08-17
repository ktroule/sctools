from operator import itemgetter

def pseudobulk_matrix(adata, bulk_variable, n_cells_replica):

    '''
    Creates a dataframe with the the aggregsted counts within each bulk_variable. 
    This is aimed for differential expression methodologies requiring pseudobulk input.

    Parameters:
        - adata (obj): Scanpy object, in obs should contain the bulk_variable annotation.
        - bulk_variable (str): Column with strings indicating which cell belongs to which group.
        - n_cells_replica (int): Number of cells to sumple from each bulk_variable. If the number of cells for a given
            bilk_variable is lower than n_cells_replica, all cells used.
    
    Output:
        - A dataframe with the aggregated count expression. Dataframe format is in gene x bulk_variable.

    '''
    
    # -- Select n_cells_replica per bulk_variable, if the number of cells is lower than
    # -- n_cells_replica then all cells are used
    pseudo_bulk_metadata = adata.obs.groupby([bulk_variable]) \
    .apply(lambda x: x.sample(n = n_cells_replica if len(x)>n_cells_replica else len(x)))

    pseudobulk_index = list(map(itemgetter(1), pseudo_bulk_metadata.index))

    adata = adata[pseudobulk_index].copy()

    
    # -- Iterate over bulk_varaible and 
    counts = {}
    for deg_sample in set(adata.obs[bulk_variable]):
        
        group_samples = adata.obs[adata.obs[bulk_variable] == deg_sample].index
        counts[deg_sample] = adata[group_samples].X.sum(axis = 0).tolist()[0]

    ## -- Create final pesudobulk matrix
    mtx = pd.DataFrame(counts,
               index = list(adata.var.index))
    
    return(mtx)
