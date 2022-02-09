This feature enables interactive manual selection of groups of cells which can be explored and visualized using differential gene expression (DGE) analysis, gene set enrichment (GSE), volcano plots, and a heatmap. 

1. Select the cell embedding type (UMAP, tSNE, PCA) and the dimensions to display on the x and y axes.

2. Under **Metadata overlay and filtering**, select the factor or variable to be used for labeling and coloring the points (e.g., cell clusters). 

    If desired, click **+Add** and use the additional filters to select cells meeting certain criteria. 
    
    Click **-Remove** to remove the filter(s).

3. Select cells for set A using either Lasso Select (default) or Box Select. 

    **Lasso Select**: if not already selected, hover over the plot and click the lasso icon at the top right of the plot. Hold down your mouse button and draw a shape around the cells you wish to select. Then click **+Set A: Add Cells** to add the cells to Set A.
    
    **Box Select**: if not already selected, hover over the plot and click the dashed box icon at the top right of the plot. Hold down your mouse button and draw a box around the cells you wish to select. Then click **+Set A: Add Cells** to add the cells to Set A.

4. Repeat step 3 to select cells for set B.

5. Under **Short name for this comparison**, give the comparison a name.

6. Click **Calculate differential gene expression**.

The above steps can be repeated multiple times and all results can be further analyzed and visualized on the **Volcano plot**, **Gene set enrichment**, and **Heatmap** tabs.

DGE and GSE results for all comparisons can also be downloaded in the **Download results** section below the plot.