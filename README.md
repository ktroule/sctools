# sctools

Repository with some custom utilities for single-cell data analysis.

## Instructions to install

Check before installing dependencies that will be modified.

```bash
pip check git+https://github.com/ktroule/sctools.git
```
Install repository
```bash
pip install git+https://github.com/ktroule/sctools.git
```
## Modules

### Visualization
- **dotplot_genes**: creates a scanpy-based dotplot from a dictionary or dataframe.

### Metrics
- **soupy_ratio**: calculates the ratio of counts associated to empty droplets compared to all droplets (empty + cells) as called by cellranger.

### Tools
- **calculate_pca**: (to be tested). calculates the PCA on the adata without scaling the gene expression.
- **load_markers**: load set of expression markers.
- **pseudobulk_matrix**: creates a pseudobulk matrix by aggregating gene counts per user-defined group. 
