# sctools

Repository with utilities for single-cell data analysis.


## Modules

### Visualization
- **dotplot_genes**: creates a scanpy-based dotplot from a dictionary or dataframe.

### Metrics
- **soupy_ratio**: calculates the ratio of counts associated to empty droplets compared to all droplets (empty + cells) as called by cellranger.

### Tools
- **load_markers**: load set of expression markers.
- **pseudobulk_matrix**: creates a pseudobulk matrix by aggregating gene counts per user-defined group. 
