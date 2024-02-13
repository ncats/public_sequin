SEQUIN is an R/Shiny application for analysis and visualization of bulk and single-cell RNA-seq data.

For complete details, please see our [manuscript](https://www.biorxiv.org/content/10.1101/2022.02.23.481646v1) or [tutorial.](https://htmlpreview.github.io/?https://github.com/ncats/public_sequin/blob/main/www/ncats-SEQUIN-tutorial.html)

For bug reports, please [submit an issue](https://github.com/ncats/public_sequin/issues) or email either [Andrew Weisman](mailto:andrew.weisman@nih.gov) or [Andrei Bombin](mailto:andrei.bombin@axleinfo.com). We are also very open to ideas for improvements. 

### Access SEQUIN in the web

Access the public app [here](https://sequin.ncats.io/).

### Install SEQUIN locally

To install SEQUIN locally, you will need the following tools installed:

* [git](https://git-scm.com/)
* [R](https://cloud.r-project.org/) (>= 4.3.2)
* [rtools4](https://cran.r-project.org/bin/windows/Rtools) (only needed in Windows and its version must be the same as the installed version of R)
* [RStudio](https://www.rstudio.com/) (optional but highly recommended)

Further, you must have cloned the public SEQUIN repository:

   ```bash
   git clone https://github.com/ncats/public_sequin.git
   ```

Installation can then be performed either in R or RStudio.

#### Instructions for R

1. Run the installation script to activate the SEQUIN project and install all required R packages (packages will be installed in the user's default directory for R libraries):

   ```r
   setwd("/path/to/public_sequin") # Use your local repo directory
   source("install_sequin.R")
   ```

1. After installation, restart R.

1. Launch SEQUIN using:

   ```r
   setwd("/path/to/public_sequin") # Use your local repo directory
   shiny::runApp(launch.browser = T) # Opens SEQUIN in browser
   ```

#### Instructions for RStudio

1. Run the installation scripts to activate the SEQUIN project and install all required R packages (packages will be installed in the `public_sequin` directory and will not interfere with user's default R environment):

   ```r
   setwd("/path/to/public_sequin") # Use your local repo directory
   source("install_1.R") # install renv and BioConductor
   rstudioapi::openProject('sctl-rshiny-complex.Rproj') # activates the project; be sure to save your current workspace when prompted if you have something to save
   source("install_2.R") # install the rest of dependencies
   # no need to restart Rstudio
   ```

1. Launch SEQUIN using:

   ```r
   setwd("/path/to/public_sequin") # Use your local repo directory
   rstudioapi::openProject('sctl-rshiny-complex.Rproj') # activates the project if needed
   shiny::runApp(launch.browser = T)
   ```

When you're finished using SEQUIN, deactivate the project to return to the default R environment.

```r
renv::deactivate()
```

To launch SEQUIN at a later time, restore the project and run the app.
```r
setwd("/path/to/public_sequin")
renv::restore()
shiny::runApp(launch.browser = T)
```

### Local data

Set the local directory where SEQUIN will read and write data by opening `app.R` and editing the first line.

```r
# Set local data directory for running as standalone app
options(localDir = "example_data")
```

Use the `example_data` directory as an example of how to format data. We recommend using a local directory outside of `public_sequin` for your real data to avoid git conflicts when you pull updated versions of SEQUIN. 

Each dataset must have the following:

1. **Experiment name**. Choose a unique experiment name with no spaces (e.g., `example_sc`). This name must be present in the `experiment_name` column in `data_info.csv` (see #4 below).
2. **Counts file**. The counts file must be in CSV format and named as follows: `<experiment name>_counts.csv` (e.g., `example_sc_counts.csv`).

    The counts file must have gene symbols in the first column and sample or cell names in subsequent columns. Counts data in each sample or cell column must consist of only integers. 
    
3. **Metadata file**. The metadata file must be in CSV format and named as follows: `<experiment name>_meta.csv` (e.g., `example_sc_meta.csv`).

    The metadata file must have sample or cell names in the first column and experimental variables in subsequent columns. At least one experimental variable is required. All experimental variables must be factors (e.g., treatment group, cell line, etc.). Some features of SEQUIN will interpret numeric variables as factors which may be undesirable.
    
4. **data_info.csv**. All datasets should have an entry in a CSV file called `data_info.csv` in your local data directory (e.g., `example_data/data_info.csv`). `data_info.csv` must have the columns listed below. All columns must be present but only those in **bold** must have an entry for each dataset.

* **`experiment_id`**: a unique integer ID value for each dataset.  
* **`experiment_name`**: a unique experiment name for each dataset (e.g., `example_sc`).
* `description`: a brief description of the dataset.
* `upload_date`: can be any value, but typically a date corresponding to the dataest.
* **`unique_table`**: must be identical to `experiment_name`. 
* `created_by`: can be any value, but typically a username or nickname for the person who created the data.
* **`type`**: `bulk` (for bulk RNA-seq) or `sc` (for single-cell RNA-seq).
* `publication`: hyperlink text to be displayed in the `Source` column on the initial load page in the app (e.g., `Walker et al., 2019`).
* `publication_link`: URL for the hyperlink in the `Source` column on the initial load page in the app (e.g., `https://www.nature.com/articles/s41598-019-56955-1`).  

#### User-updated metadata

SEQUIN includes two tools, **Merge clusters** and **Group cells by gene expression**, that enable the user to create updated metadata files that can be loaded instead of the default metadata.

SEQUIN will record information about user-updated metadata files in `sc_useradd.csv`.

To delete user-updated metadata files, delete the appropriate rows from `sc_useradd.csv` and delete the user-updated metadata files. The user-updated metadata files can be identified by looking at the `table_name` column in `sc_useradd.csv`. If you have deleted all user-updated metadata files, you may also delete `sc_useradd.csv`. 
