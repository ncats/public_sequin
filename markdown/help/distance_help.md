The distance matrix heatmap enables visualization of differences in gene expression among samples as measured by Euclidean distance.

1. Choose the number of genes and samples to be included.

    2000 genes (or all if less than 2000) are selected by default, but this value can be changed.   
  
    Click **All** to use all genes.  
  
    Click **Default** to use the default number of genes.  
  
    If genes are downsampled, genes with higher mean expression levels across samples are given priority. 
  
    A maximum of 500 samples will be selected regardless of your selection. To downsample further, change **Samples**. To use the maximum allowed (also the default), click **Max**. 
  
    Samples are selected randomly. To obtain the same selection of samples, use the same **Seed** value. To obtain a new random selection of samples, change the value.

2. Click **Build heatmap** to compute distances and build the heatmap.
