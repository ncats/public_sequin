This page enables manual interactive grouping of cells in a UMAP plot based on expression of one or more genes (e.g., cell type markers). The custom sets can be saved in a new metadata table for future analysis.

1. Select one or more genes under **Select genes**.

    For multiple genes, use **Summary metric** to choose whether the mean, median, or sum of their expression should be displayed on the UMAP plot.

2. Below the UMAP plot, click **+New set** as many times as desired to create sets for grouping cells. Use the default set names or create your own.

3. If desired, select a range of expression values to limit cell selection to only cells with expression values within the range. By default, the full range is selected.

    If **Inclusive** is selected for the lower or upper limit of the range, the filter will include the lower or upper limit, respectively. Otherwise, it will include only cells with gene expression below or above the limit.
    
4. If desired, select cells manually using either Lasso Select or Box Select. 

    **Lasso Select**: if not already selected, hover over the plot and click the lasso icon at the top right of the plot. Hold down your mouse button and draw a shape around the cells you wish to select.
    
    **Box Select**: if not already selected, hover over the plot and click the dashed box icon at the top right of the plot. Hold down your mouse button and draw a box around the cells you wish to select. 

5. Once you've made a selection, click **Add to \<set name\>** to add the cells to the desired set.

    You can add multiple cell selections to the same set. 
    
    The number of cells added to a set will be shown beside each **Add to \<set name\>** button.

6. Repeat steps 3-5 to assign cells to additional sets.

7. Provide a username, a short comment if desired, and give the updated table a name or accept the default name. 

8. Click **Save to database**. 

    An updated metadata table containing only the cells you assigned to sets will be saved. 
    
    The metadata table will contain a **Set** column containing the set name to which each cell was assigned. The next time you select the same dataset on the initial load page, the updated table will be available to use in place of the default metadata table.