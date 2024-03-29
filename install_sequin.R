# Set a mirror to download packages from
options(repos = c(CRAN = "https://mirrors.nics.utk.edu/cran"))

# Install renv and BiocManager
install.packages("renv")
install.packages("BiocManager")

# Activate project
renv::activate()

# Install latest versions of packages from renv.lock
lockfile <- renv:::renv_lockfile_read("renv.lock")
biocPackages <- names(lockfile$Packages)[sapply(lockfile$Packages, function(pkg) pkg$Source == "Bioconductor")]
cranPackages <- names(lockfile$Packages)[sapply(lockfile$Packages, function(pkg) pkg$Source == "Repository")]
githubPackages <- names(lockfile$Packages)[sapply(lockfile$Packages, function(pkg) pkg$Source == "GitHub")]
githubRepos <- sapply(lockfile$Packages[githubPackages], function(pkg) paste0(pkg$RemoteUsername, "/", pkg$RemoteRepo))

options(repos = BiocManager::repositories())
renv::install(biocPackages)
renv::install(cranPackages)
renv::install(githubRepos)

# Update renv.lock with current package versions
renv::snapshot(prompt = F)