import pandas as pd
import scanpy as sc

def colorize_labels(adata, colorize, ncols):

    '''
    Plots UMAPs depiciting the label of interest.

    Parameters:
        - adata (obj): a scanpy object.
        - colorize (str): a column of `adata.obs`.
        - ncols (int): number of UMAPs to display per row.
    Output:
        - a set of umaps colores individually per label of interst
    '''
    # -- Copy adata of interest
    x = adata.copy()

    # -- Codify column of interest as dummy variables
    annot_dummy = pd.get_dummies(x.obs[colorize])

    n_cells = [ '(' + str(i) + ')' for i in list(annot_dummy.sum(axis = 0)) ]
    col_names = [c + '\n' + n for c, n in zip(annot_dummy.columns, n_cells)]
    annot_dummy.columns = col_names

    # -- Concat dummy-coded variables toa data.obs
    x.obs = pd.concat([x.obs, annot_dummy],
                    axis = 1)
    
    # -- Replace original dummy coding for a and b
    for col in col_names:
        x.obs[col] = x.obs[col].map({ 0 : 'b',
                                      1 : 'a'})

   # -- Plot columns of interest
    sc.pl.umap(
        x,
        color = col_names,
        frameon = False,
        ncols = ncols,
        colorbar_loc = None,
        legend_loc = None,
        alpha = 0.5,
        palette = {'a' : 'red',
                   'b' : '0.9'})