The correlation heatmap enables visualization of correlations in gene expression among samples or cells, as measured by the Pearson correlation coefficient (*R*).

1. Choose the number of genes and samples/cells to be included.

    2000 genes (or all if less than 2000) are selected by default, but this value can be changed.   
  
    Click **All** to use all genes.  
  
    Click **Default** to use the default number of genes.  
  
    If genes are downsampled, genes with higher mean expression levels across samples/cells are given priority. 
  
    A maximum of 500 cells or samples will be selected regardless of your selection. To downsample further, change **Samples** (for bulk data) or **Cells** (for single-cell data). To use the maximum allowed (also the default), click **Max**. 
  
    Samples/cells are selected randomly. To obtain the same selection of samples/cells, use the same **Seed** value. To obtain a new random selection of samples/cells, change the value.

2. Click **Build heatmap** to compute correlations and build the heatmap.
  
    Click on any point in the correlation heatmap to view a scatterplot (under **Correlation scatterplot**) comparing gene expression between the corresponding pair of samples/cells.
