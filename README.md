# SEQUIN Exploration

## SEQUIN
The SEQUIN shiny web application builds off of the original IRIS-EDA (**I**nteractive **R**NA-seq analysis and **I**nterpretation using **S**hiny-**E**xpression **D**ata **A**nalysis) and the R library scClustViz. IRIS-EDA is a web-based tool for the analysis of bulk and single cell RNA-seq (scRNA-seq) count data. scClustViz is also an interactive R Shiny tool for visualizing single-cell RNA-seq clustering, select optimal clustering resolution and to annotate cell types. The SEQUIN shiny web application extends beyond IRIS-EDA and scClustViz by incorporating additional features that can be adjusted by the end-user. These additional features include:
1. Combining samples across multiple experiments
2. Utilizing pre-assigned clusters or calculating multiple resolutions for cluster assignment for scRNA-seq data
3. Performing differential gene expression (DGE) by cluster, any factor in the metadata or using lasso tool selection
4. Adding **G**ene**S**et **E**nrichment **A**nalysis (GSEA)
5. Adjusting cluster resolution by combining clusters or creating custom clustering via gene expression of gene(s)
6. Calculating module scores for all cells across clusters for a greater understanding of stem cell differentiation 

Two example datasets are included in the app: one bulk and one scRNA-Seq dataset.

## Running the SEQUIN web application

We recommend visiting the web application to explore data. You must be on the NIH VPN to access the app here: 
https://ipsceq-ci.ncats.io

## Running SEQUIN locally

Alternatively, if you choose to run SEQUIN locally in R, you will first need a github account. You can clone the repository by following the step below:
git clone https://github.com/ncats/sctl-rshiny-complex.git

SEQUIN is an R shiny application built in R and can be run on version >= 4.0.4. If you do not have R installed, see the two links below to install R and/or RStudio.

[Download (or upgrade) R here](https://cloud.r-project.org/)

[RStudio (a graphical user interface) can be downloaded here, but is not required](https://www.rstudio.com/products/rstudio/download3/)

## Running SEQUIN via Docker

The SEQUIN R shiny app can also be easily run using hte Docker container. To install Docker, go to the following link and follow the directions, below.

https://docs.docker.com/get-docker/

A Dockerfile has been created and can be found in this repository, which can be used for stable R version and package control.

On Mac, Windows or Linux, the end-user will need to create an Rprofile.site file. In a text editor or VIM, type the following:

      local({
      options(shiny.port = 3838, shiny.host = "0.0.0.0")
      })
      
This specifies the port in which to host the app locally.

Navigate to the cloned github repository. The Dockerfile and the install.R file should be located one directory above the directory named, app/. After Docker is installed locally on your machine, you must navigate to the location of the Dockerfile and build it by running the following command:

      docker build -t sequin_app .
      
The build process can take upwards of a couple of hours, which is why we recommmend accessing the app via website. After the build is complete, you can access the app by typing the following command in the terminal:

      docker run -p 3838:3838 sequin_app

You will need to navigate to your browser (we recommend chrome) and type the following:

      http://localhost:3838 

It will take a few minutes to load the app as there are several libraries loaded initially.

## Vignette
[A detailed tutorial can be found here](markdown/ncats-iPSC-tutorial.md)

## Running SEQUIN's user-friendly web app:
The SEQUIN R shiny application can be run directly in the R console by running the following code in your R console:

```{r}
library(shiny)
runAPP()
```
## Contact
If you encounter any problems running the software, have installation bugs or problems, please add an issue on the Issues tab or email ben.ernest@ranchobiosciences.com and/or marissa.hirst@ranchobiosciences.com. We are also very open to any comments or suggestions, including imporovements.
