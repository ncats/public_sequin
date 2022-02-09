This page includes a dotplot and table for exploring differentially expressed (DE) genes in a chosen cluster, based on the currently chosen resolution. 

A dendrogram above the dotplot depicts relationships among chosen genes, and a dendrogram on the left depicts relationships among cell clusters based on expression of the chosen genes. 

The size of each dot represents the detection rate of a gene in a given cell cluster. 

Each dot is colored based on the expression level of the corresponding gene in a given cluster.

The following parameters control the plot and table:

#### **Cluster**

The cluster to be explored in the dotplot and table.

#### **FDR threshold**

DE genes with a false discovery rate (FDR)-adjusted p-value less than or equal to this threshold are included.

#### **Dotplot genes**

* **DE vs rest**: DE genes in the chosen cluster vs. all other cells.
  
* **DE vs neighbor**: DE genes in the chosen cluster vs. the nearest neighbor cluster. 
  
#### **# genes per cluster**

The number of top DE genes to show in the dotplot.

### **DGE table**

Use **Abs. log2 fold-change** and **Adj. p-value** to include only genes that meet these thresholds, and click **Load table** to view the results.