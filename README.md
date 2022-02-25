### SEQUIN

SEQUIN is an R/Shiny application for analysis and visualization of bulk and single-cell RNA-seq data.

For complete details, please see our [manuscript](https://www.biorxiv.org/content/10.1101/2022.02.23.481646v1) or [tutorial.](www/ncats-SEQUIN-tutorial.html)

For bug reports, please [submit an issue](https://github.com/ncats/public_sequin/issues) or email either [Ben Ernest](mailto:ben.ernest@ranchobiosciences.com) or [Marissa Hirst](mailto:marissa.hirst@ranchobiosciences.com). We are also very open to ideas for improvements. 

### Access SEQUIN in the web

Access the public app [here](https://sequin.ncats.io/).

### Install SEQUIN locally

To install SEQUIN locally, you will need the following tools installed:

* [git](https://git-scm.com/)
* [R](https://cloud.r-project.org/) (>= 4.0.4)
* [rtools4](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) (only needed in Windows)
* [RStudio](https://www.rstudio.com/) (optional but highly recommended)

Installation steps:

1. In the terminal, clone the repo.
   ```bash
   git clone https://github.com/ncats/public_sequin.git
   ```

2. From R, run the installation script to activate the SEQUIN project and install all required R packages.

   ```r
   setwd("/path/to/public_sequin") # Use your local repo directory
   source("install_sequin.R")
   ```

3. From R, launch SEQUIN.

   ```r
   shiny::runApp(launch.browser = T) # Opens SEQUIN in browser
   ```

When you're finished using SEQUIN, deactivate the project to return to the default R environment.

```r
renv::deactivate()
```

To launch SEQUIN at a later time, restore the project and run the app.
```r
setwd("/path/to/public_sequin/")
renv::restore()
shiny::runApp(launch.browser = T)
```
