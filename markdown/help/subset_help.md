### **Help**

This menu enables subsetting data to select specific samples or groups of samples to load. 

For bulk RNA-seq, PCA visualization may be used to identify and exclude outlier samples as described below.

To close the menu without subsetting, click **Close**.

#### **Subsetting data**
Metadata from each experiment selected is presented as a table. Each row corresponds to samples or cells from one level or combination of levels of the factor(s) chosen in **Select factor(s)**.

Click rows to include the corresponding samples or cells. Click a highlighted row to deselect it.

Click **Select all** or **Deselect all** to select all samples or deselect all samples from a dataset, respectively. 

**sample_id** represents a unique identifier for each sample in bulk RNA-seq data. Choose **sample_id** in **Select factor(s)** to view all individual samples and choose which samples to include.

Once you have subset your dataset(s), click **Submit** at the top of the page to return to the initial load page where you may load this subset by clicking **Load selected sample(s)**.

#### **Identify outliers**

1. Click **Identify outliers** at the top of the page.

2. Click **Run PCA** on the **Identify outliers by PCA** page.

    Each dataset is represented by a PCA plot. Samples in the plot are colored based on the factor chosen in **Grouping factor**. 
    
3. Inspect each PCA plot to decide if outlier samples are present.

4. To exclude a sample, click the corresponding point in the PCA plot. 

    The chosen sample will be shown just above the plot.

    You may zoom in on an area of the plot by holding down your mouse button and drawing a box around the area. Double-click anywhere in the plot to return to the original zoom level.
    
5. Click **+ Add sample** to add the chosen sample to the outliers list for a given dataset.

    The selected samples will be shown just above the plot.
    
    Click **-Clear samples** to clear all samples from the outliers list for a dataset.

6. Repeat steps 5 and 6 to select additional samples from any dataset to exclude.

7. Click **Submit outliers** at the top of the page.

    You will be returned to the **Subset data** page where all samples except the chosen outlier samples will be highlighted. 
    
    If no outlier samples were selected for a dataset, no samples will be highlighted. To include samples from a dataset in this case, select samples or groups of samples as described above in **Subsetting data**.
    
8. Click **Submit** at the top of the page as described above in **Subsetting data**.

