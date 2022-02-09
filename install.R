pkgs <- c(
    "Rhtslib",
    "geneplotter",
    "pheatmap",
    "DESeq2",
    "edgeR",
    "limma",
    "QUBIC",
    "impute",
    "preprocessCore",
    "AnnotationDbi",
    "RUVSeq",
    "EDASeq",
    "WGCNA",
    "multtest",
    "biomaRt",
    "xml2",
    "rvest",
    "ComplexHeatmap",
    "MAST")

versionSpecificPkgs <- c(
  "org.Hs.eg.db" = "3.12",
  "org.Mm.eg.db" = "3.12",
  "GO.db" = "3.12"
)

# Install specific versions of packages
sapply(names(versionSpecificPkgs), function(pkg) {
  pkgVersion <- versionSpecificPkgs[pkg]
  BiocManager::install(pkgs = pkg, update = FALSE, ask = FALSE, version = pkgVersion)
})

ap.db <- available.packages(contrib.url(BiocManager::repositories()))
ap <- rownames(ap.db)
fnd <- pkgs %in% ap
pkgs_to_install <- pkgs[fnd]

ok <- BiocManager::install(pkgs_to_install, update=FALSE, ask=FALSE) %in% rownames(installed.packages())

if (!all(fnd))
    message("Packages not found in a valid repository (skipped):\n  ",
            paste(pkgs[!fnd], collapse="  \n  "))
if (!all(ok))
    stop("Failed to install:\n  ",
         paste(pkgs_to_install[!ok], collapse="  \n  "))

