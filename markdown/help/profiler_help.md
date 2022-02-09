The iPSC profiler includes visualizations based on module scores for gene sets (also referred to as modules or profiles). A module score is a measure of gene expression for a set of genes, relative to randomly selected control genes, in a sample or cell. The visualizations can be viewed in the sub-tabs and are described below.

### **Heatmap**

The heatmap displays module scores for each chosen module across all samples or cells. 

Choose a grouping factor for grouping samples or cells (columns) and a color palette for the heatmap.

By default, all gene modules are shown. To include one or more individual modules, un-check **Use all modules**, then select the modules(s) you want to display. To select all modules again, re-check **Use all modules**.

### **PCA/tSNE/UMAP**

Module scores for individual modules are projected onto dimensional reduction plots. 

Choose PCA, tSNE, or UMAP as the dimensional reduction method and select a module. 

The same plot is shown on the right with the samples or cells colored according to the chosen grouping factor to enable visualizing how module scores vary by group or cluster.

### **Violin plot**

The distributions of module scores for a chosen module are shown across samples or cells, grouped according to the chosen grouping factor.

### **Module info**

The module descriptions table shows information about each gene module, including a brief description and a publication source if available.