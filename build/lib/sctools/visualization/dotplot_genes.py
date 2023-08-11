import scanpy as sc
import panas as pd

def dotplot_genes(adata, colorize, groupby, group = 'All', dot_min = 0, dendrogram = False):
    
    '''
    This function takes a dictionary or a dataframe and plot a dotplot, also checks if one of the selected genes
    is not presented in the adata, discards it from the plot and warns about it.

    Parameters:
        - adata (obj): scanpy object.
        - colorizer (dict/dataframe): if input is dataframe it must contain as colnames: `Gene`, `Cell`, `Group`.
        - groupby (str): group cells by annotation (this must be store in `adata.obs`).
        - group (str): if colorize is a dataframe, use only markers label as `colorize['Group'] == xx`.
        - dot_min (float): Minimum fraction of cells expressing the gene.
        - dendodgram (bool): Whether to clusterize annotation selected in groupby.
    Output:
        - a scanpy-based doptplot.

    '''
    
    # -- Check if colorize is a dict or df and subset to group of interest
    if(type(colorize) is dict):
        colorize = colorize
    else:
        if ('All' in group):
            colorize = colorize
            colorize = colorize.groupby('Cell')['Gene'].apply(list).to_dict()
        else:
            idx = [ g in group for g in colorize['Group'] ]
            colorize = colorize[idx]

            # -- Create dictionary with genes of interst grouped by
            # -- Cell column
            colorize = colorize.groupby('Cell')['Gene'].apply(list).to_dict()

    # -- All genes present in the dictionary
    all_genes = sorted({x for v in colorize.values() for x in v})

    genes_to_plot = []
    cell_to_plot = []
    new_dict = dict()
    for key in colorize.keys():

        genes = colorize[key]

        genes_present = [ g for g in genes if g in adata.var.index ]
        new_dict[key] = genes_present

        genes_to_plot = genes_to_plot + genes_present


    sc.pl.dotplot(adata,
              var_names = new_dict,
              groupby = groupby,
              standard_scale = 'var',
              return_fig = True,
              dendrogram = dendrogram,
              use_raw = False).style(grid = True,
                      cmap = 'Reds',
                      dot_min = dot_min).show()
    
    # -- Print genes that are not present in the adata object
    genes_not_present = set(all_genes) - set(genes_to_plot)
    print(f'Genes not found in the object ({genes_not_present})')
