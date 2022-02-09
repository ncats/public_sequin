You may load one or more existing datasets, upload a custom dataset, or merge existing datasets with a custom dataset. Datasets may be loaded in full or subset to specific samples or experimental groups. Bulk RNA-seq data may be visualized by PCA to enable exclusion of outlier samples. 

First, choose between **Single-cell RNA-seq** and **Bulk RNA-seq**. 

### **Existing datasets**
On the **Existing datasets** tab, click rows in the table to choose one or more datasets. 

To undo a selection, click the corresponding row again, or to clear all selections, click **Clear selections**.

#### **Option: Select user-updated metadata**
For single-cell data, if an updated metadata table has been previously created for the chosen dataset using the **Merge clusters** or **Group cells by gene expression** tools, the updated metadata table may be selected and loaded in place of the default metadata dable.

### **Custom dataset**
On the **Custom dataset** tab, provide a name for the dataset and upload a counts matrix and metadata table.

Click **+ Add existing experiments** to combine existing experiments with the custom dataset. 

Click **- Clear existing experiments** to undo the existing dataset selection(s). 

### **Subset data**
Once existing datasets have been selected or a custom dataset has been uploaded, click **Subset data** to select specific samples or experimental groups based on the metadata table.

### **Single-cell RNA-seq data**
Choose whether to filter out pseudogenes and mitochondrial genes (recommended) and whether to use pre-computed clusters from metadata (only for a single dataset) for downstream analyses or to perform clustering over a range of resolutions using Seurat (required for merged or subset data).

#### **Downsample cells**
For single-cell data, a maximum of 10,000 cells per dataset will be loaded. 

To downsample, check **Downsample cells** and choose the maximum number of cells per dataset. 

Downsampling is applied to each dataset separately and prior to subsetting, so the number of resulting cells may be less than the maximum requested. 

To always generate the same selection of cells, keep the **Random seed** value the same. To generate a new random selection of cells, change the value.

### **Bulk RNA-seq data**
Select a read count cutoff, whether to use batch correction (if multiple datasets are selected) and which metadata variable to use, and choose a transformation method for counts.

### **Load data**
To include all cells or samples from each dataset, click **Load full dataset(s)**.

To include only cells or samples from your **Subset data** selections, click **Load selected sample(s)**.