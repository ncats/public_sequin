#---------------------------------------------------------------------
# Title:         NCATS SEQUIN - Core Functions
# Author:        Marissa Hirst
# Author2:       Ben Ernest
# Last Modified: 2021-04-01
# --
# Created:       2018-01-26 11:29:39 CDT
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# CONTENTS
#---------------------------------------------------------------------
#
#   SEQUIN core function organization (most of this is bulk, from IRIS app)
#
#       01 HOUSE......... house-keeping functions
#       02 LIMMA-RF...... limma-voom return fit functions
#       03 EDGER-RF...... edgeR return fit functions
#       04 DESEQ-RF...... DESeq2 return functions
#       05 PLOTs......... all plotting functions
#       06 DGE........... roll-up/full rank check dge functions
#       07 GSE........... roll-up gse functions
#       08 SCCLUST-VIZ... scClustViz functions: Updated?
#       09 MYSQL-DB...... functions to push to SQL db
#
#---------------------------------------------------------------------
#
#   SEQUIN core function layout
#
#       |-  HOUSE
#           |- Specify seurat major version
#           |- Get normalized counts for selected genes
#               |- `getGenes()`
#           |- Get contrast tables from different object classes
#               |- `getContTable()`
#       |- LIMMA-RF
#           |- LIMMA functions for bulk RNA-Seq
#       |- EDGER-RF
#           |- EDGER functions for bulk RNA-Seq
#       |- DESEQ-RF
#           |- DESEQ2 functions for bulk RNA-Seq
#       |- PLOTS
#            |- QC plot fxns
#            |- PCA plot fxns
#            |- tSNE plot fxns
#            |- uMAP plot fxns
#            |- Heatmap plot fxns
#            |- Bulk DGE, MA, Volcano & scRNA-Seq DGE plot fxns
#            |- WGCNA plot fxns
#       |- DGE
#            |- DGE rollup functions & full rank check
#       |- GSE
#            |- GSE rollup functions
#       |- SCCLUST-VIZ
#            |- Updated scClustViz rollup/plotting functions
#       |- MYSQL-DB
#            |- Retrieva data and create a seurat object from MYSQL db
#
#---------------------------------------------------------------------

###################################################################
###################################################################
### SECTION 01 - HOUSE KEEPING FUNCTIONS (HOUSE)
###################################################################
###################################################################

# HOUSE - do not return class msg
options(getClass.msg = F)

# HOUSE - specify seurat major version
SeuratMajorVersion <- function() {
  vers <- as.character(packageVersion("Seurat"))
  as.numeric(strsplit(vers, split = "")[[1]][1])
}

# HOUSE - get size factor extraction (for normalization)
getSizeFact <- function(rc.data) {
  geomMean <- function(x) {
    prod(x)^(1 / length(x))
  }
  gm.mean <- apply(rc.data, 1, geomMean)
  gm.mean[gm.mean == 0] <- NA
  rc.data <- sweep(rc.data, 1, gm.mean, FUN = "/")
  sf <- apply(rc.data, 2, median, na.rm = TRUE)
  return(sf)
}

# HOUSE - get normalized counts from different object classes
getNormCounts <- function(rc.data, type) {
  if ("DESeqDataSet" %in% class(rc.data)) {
    nc.data <- BiocGenerics::counts(rc.data, normalize = TRUE)
    return(nc.data)
    # } else if (type %in% c("yes1", "yes2", "no")) {
  } else if (type %in% c("Bulk", "RASL")) {
    sf <- getSizeFact(rc.data)
    nc.data <- sweep(rc.data, 2, sf, FUN = "/")
    return(nc.data)
    # } else if (type %in% c("yes3", "scrna")) {
  } else if (type == "Single-cell") {
    nc.data <- NormalizeData(rc.data)
    return(nc.data)
  }
}

# HOUSE - get normalized counts for selected gene (heatmap interactivity)
getGenes <- function(rc.data, id, coldata, type) {
  nc.data <- getNormCounts(rc.data, type)
  nc.data <- as.data.frame(nc.data[rownames(nc.data) == id, ])
  names(nc.data) <- "counts"
  dat.l <- list(coldata, nc.data)
  nc.data <- Reduce(
    merge, lapply(dat.l, function(x) data.frame(x, Sample = row.names(x)))
  )
  return(nc.data)
}

# HOUSE - get contrast tables from different object classes
getContTable <- function(
  de.genes, coef, cts, expset, design, fact, fact5, fact6) {
  if (class(de.genes) == "MArrayLM") {
    tryCatch({
      de.genes2 <- topTable(
        fit = de.genes,
        coef = coef,
        number = nrow(cts)
      )
      de.genes2 <- de.genes2[order(rownames(de.genes2)), ]
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "avgexpr",
        "t",
        "pvalue",
        "padj",
        "B",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"
    }, error = function(e)
      NULL)
  } else if (class(de.genes) == "DGEGLM") {
    if (expset == "exp1" | expset == "exp2") {
      tryCatch({
        de.genes2 <- glmLRT(
          glmfit = de.genes,
          contrast = design[, coef]
        )
        de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2$baseMean <- rowMeans(cts)
        names(de.genes2) <- c(
          "log2FoldChange",
          "logCPM",
          "LR",
          "pvalue",
          "padj",
          "baseMean"
        )
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp3" | expset == "exp4") {
      tryCatch({
        de.genes2 <- glmLRT(
          glmfit = de.genes,
          coef = coef
        )
        de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2$baseMean <- rowMeans(cts)
        names(de.genes2) <- c(
          "log2FoldChange",
          "logCPM",
          "LR",
          "pvalue",
          "padj",
          "baseMean"
        )
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp5" | expset == "exp6") {
      tryCatch({
        de.genes2 <- glmLRT(
          glmfit = de.genes,
          contrast = design[, coef]
        )
        de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2$baseMean <- rowMeans(cts)
        names(de.genes2) <- c(
          "log2FoldChange",
          "logCPM",
          "LR",
          "pvalue",
          "padj",
          "baseMean"
        )
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp7") {
      tryCatch({
        de.genes2 <- glmLRT(
          glmfit = de.genes,
          coef = coef
        )
        de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2$baseMean <- rowMeans(cts)
        names(de.genes2) <- c(
          "log2FoldChange",
          "logCPM",
          "LR",
          "pvalue",
          "padj",
          "baseMean"
        )
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    }
  } else if (class(de.genes) == "DESeqDataSet") {
    if (expset == "exp1") {
      tryCatch({
        coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
        coef2 <- unlist(coef2)
        de.genes2 <- results(
          object = de.genes,
          contrast = c(fact, coef2[1], coef2[2])
        )
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp2") {
      tryCatch({
        coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
        coef2 <- unlist(coef2)
        de.genes2 <- results(
          object = de.genes,
          contrast = c("group", coef2[1], coef2[2])
        )
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"        
      }, error = function(e)
        NULL)
    } else if (expset == "exp3" | expset == "exp4") {
      tryCatch({
        names(mcols(de.genes))[grep(
          "log2 fold change", 
          mcols(mcols(de.genes))$description
        )] <- colnames(design)
        de.genes2 <- results(
          object = de.genes,
          contrast = list(coef)
        )
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp5") {
      tryCatch({
        coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
        coef2 <- unlist(coef2)
        de.genes2 <- results(
          object = de.genes,
          contrast = c(fact5, coef2[1], coef2[2])
        )
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp6") {
      tryCatch({
        coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
        coef2 <- unlist(coef2)
        de.genes2 <- results(
          object = de.genes,
          contrast = c(fact6, coef2[1], coef2[2])
        )
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    } else if (expset == "exp7") {
      tryCatch({
        de.genes2 <- results(
          object = de.genes,
          name = coef
        )
        de.genes2 <- as.data.frame(de.genes2)
        de.genes2 <- subset(de.genes2, baseMean != 0)
        de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
        names(de.genes2)[1] <- "id"
      }, error = function(e)
        NULL)
    }
  }
  if(!exists("de.genes2")) return()
  return(de.genes2)
}

# HOUSE - estimate housekeeping genes effect using RUVg
CalcW1 <- function(counts) {
  if(!is.matrix(counts)) counts <- as.matrix(counts)
  hkGenes <- c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 
               'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 
               'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 
               'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 
               'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 
               'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 
               'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 
               'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 
               'WDTC1', 'ZNF544','C1orf43', 'CHMP2A','EMC7','GPI',
               'PSMB2','PSMB4','RAB7A','REEP5', 'SNRPD3','VCP',
               'VPS29')
  hkGenes <- hkGenes[hkGenes %in% rownames(counts)]
  if(length(hkGenes) == 0) return()
  myResult <- RUVg(counts, cIdx = hkGenes, k = 1)
  w1 <- as.numeric(myResult$W[, "W_1"])
  names(w1) <- colnames(counts)
  return(w1)
}

###################################################################
###################################################################
### SECTION 02 - LIMMA RETURN FIT FUNCTIONS (LIMMA-RF)
###################################################################
###################################################################

# LIMMA - EXP 1 - two group comparisons
limma.exp1 <- function(fact, coldata, cts, perm.h, batchVar = NULL, w1 = NULL) {
  
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ 0 + coldata[, fact])
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, fact])
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(w1)) {
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(coldata[, fact])
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    coldata[, fact] <- relevel(coldata[, fact], ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, fact])
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(coldata[, fact])
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    coldata[, fact] <- relevel(coldata[, fact], ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact])
  }
  
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", replacement = "W_1", x = colnames(design))
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  colnames(design) <- designLevels
  
  v <- voom(cts, design)
  fit <- lmFit(v)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  return(list(fit.cont, cont))  
}

## LIMMA - EXP 2 - multiple factor comparisons
limma.exp2 <- function(fact1, fact2, coldata, cts, perm.h, batchVar = NULL,
                       w1 = NULL) {
  
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ 0 + group.c)
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + group.c)
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(w1)) {
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(group.c)
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    group.c <- relevel(group.c, ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, batchVar] + group.c)
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(group.c)
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    group.c <- relevel(group.c, ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + group.c)
  }

  colnames(design) <- gsub("group.c", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", replacement = "W_1", x = colnames(design))
  rownames(design) <- rownames(coldata)
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  colnames(design) <- designLevels
  
  v <- voom(cts, design)
  fit <- lmFit(v)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  return(list(fit.cont, cont)) 
}

## LIMMA - EXP 3 - classical interactions
limma.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl,
                       batchVar = NULL, w1 = NULL) {

  f1n <- length(levels(coldata[, fact1]))
  nBatchCols <- 0
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- 1
    f1n <- f1n + nBatchCols
  } else if(is.null(w1)) {
    design <- model.matrix(~ coldata[, batchVar] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
    f1n <- f1n + nBatchCols - 1
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
    f1n <- f1n + nBatchCols
  }
  
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[(nBatchCols + 2):ncol(tmp1)] <- paste0(colnames(tmp1)[(nBatchCols + 2):ncol(tmp1)], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)
  
  v <- voom(cts, design)
  fit <- lmFit(cts, design)
  fit.cont <- eBayes(fit)
  # if(nBatchCols > 0) design <- design[, -2:-(nBatchCols + 1), drop = FALSE]
  return(list(fit.cont, design))
}

## LIMMA - EXP 4 - added effects blocking and paired
limma.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl,
                       batchVar = NULL, w1 = NULL) {

  f1n <- length(levels(coldata[, fact1]))
  nBatchCols <- 0
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- 1
  } else if(is.null(w1)) {
    design <- model.matrix(~ coldata[, batchVar] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar])) - 1
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
  }

  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[(1 + nBatchCols + f1n):ncol(design)] <- paste0(
    colnames(design)[(1 + nBatchCols + f1n):ncol(design)],
    "_VS_", 
    fact2.rlvl
  )
  
  v <- voom(cts, design)
  fit <- lmFit(cts, design)
  fit.cont <- eBayes(fit)
  # if(nBatchCols > 0) design <- design[, -2:-(nBatchCols + 1), drop = FALSE]
  return(list(fit.cont, design))
}

## LIMMA - EXP 7 - user inputs
limma.exp7 <- function(cts, mod.matrix) {
  design <- mod.matrix
  fit <- lmFit(cts, design)
  fit <- eBayes(fit)
  return(list(fit, design))
}

###################################################################
###################################################################
### SECTION 03 - EDGER RETURN FIT FUNCTIONS (EDGER-RF)
###################################################################
###################################################################

## EDGER - EXP1 - two group comparisons
edger.exp1 <- function(fact, coldata, cts, perm.h, norm, batchVar = NULL,
                       w1 = NULL) {
  
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ 0 + coldata[, fact])
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, fact])
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(w1)) {
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(coldata[, fact])
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    coldata[, fact] <- relevel(coldata[, fact], ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, fact])
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(coldata[, fact])
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    coldata[, fact] <- relevel(coldata[, fact], ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact])
  }
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", replacement = "W_1", x = colnames(design))
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  colnames(design) <- designLevels

  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont)) 
}

## EDGER - EXP2 - multiple factor comparisons
edger.exp2 <- function(fact1, fact2, coldata, cts, perm.h, norm, 
                       batchVar = NULL, w1 = NULL) {
  
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ 0 + group.c)
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + group.c)
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(w1)) {
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(group.c)
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    group.c <- relevel(group.c, ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, batchVar] + group.c)
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(group.c)
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    group.c <- relevel(group.c, ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + group.c)
  }

  colnames(design) <- gsub("group.c", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", replacement = "W_1", x = colnames(design))
  rownames(design) <- rownames(coldata)
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  colnames(design) <- designLevels
  
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont)) 
}

## EDGER - EXP3 - classical interactions
edger.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl, 
                       norm, batchVar = NULL, w1 = NULL) {
  
  f1n <- length(levels(coldata[, fact1]))
  nBatchCols <- 0
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)

  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- 1
    f1n <- f1n + nBatchCols
  } else if(is.null(w1)) {
    design <- model.matrix(~ coldata[, batchVar] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
    f1n <- f1n + nBatchCols - 1
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
    f1n <- f1n + nBatchCols
  }
  
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[(nBatchCols + 2):ncol(tmp1)] <- paste0(colnames(tmp1)[(nBatchCols + 2):ncol(tmp1)], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, design))
}

## EDGER - EXP4 - added effects blocking and paired
edger.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl,
                       norm, batchVar = NULL, w1 = NULL) {
  
  f1n <- length(levels(coldata[, fact1]))
  nBatchCols <- 0
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)

  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- 1
  } else if(is.null(w1)) {
    design <- model.matrix(~ coldata[, batchVar] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar])) - 1
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
  }

  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[(1 + nBatchCols + f1n):ncol(design)] <- paste0(
    colnames(design)[(1 + nBatchCols + f1n):ncol(design)],
    "_VS_", 
    fact2.rlvl
  )
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, design))
}

## EDGER - EXP5 - main effects only
edger.exp5 <- function(fact, fact.levl, cts, coldata, perm.h, norm, 
                       batchVar = NULL, w1 = NULL) {
  
  coldata[, fact] <- relevel(x = coldata[, fact], ref = fact.levl)
  # Changed this from original because it was giving very different results 
  # compared to batch-corrected. Originally perm.c was just the comparison 1 level
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ 0 + coldata[, fact])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, fact])
  } else if(is.null(w1)) {
    design <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, fact])
    perm.c <- perm.h
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact])
    perm.c <- perm.h
  }
  
  colnames(design) <- gsub("coldata\\[, fact\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  
  rownames(design) <- rownames(coldata)
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  colnames(design) <- designLevels
  
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont))
}

## EDGER - EXP6 - main effects + grouping factor
edger.exp6 <- function(me.fact, me.levl, gp.fact, gp.levl, cts, coldata, 
                       perm.h, norm, batchVar = NULL, w1 = NULL) {
  
  cts <- cts[, which(coldata[, gp.fact] == gp.levl)]
  coldata <- coldata[which(coldata[, gp.fact] == gp.levl), ]
  coldata[] <- lapply(coldata, function(x) if(is.factor(x)) factor(x) else x)
  coldata[, me.fact] <- relevel(x = coldata[, me.fact], ref = me.levl)
  
  perm.c <- paste0(perm.h, "-", me.levl)
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ 0 + coldata[, me.fact])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, me.fact])
  } else if(is.null(w1)) {
    design <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, me.fact])
    perm.c <- perm.h
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, me.fact])
    perm.c <- perm.h
  }
  
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, me.fact\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  colnames(design) <- designLevels
  
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont))
}

## EDGER - EXP7 - user input
edger.exp7 <- function(cts, mod.matrix) {
  design <- mod.matrix
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  return(list(fit, design))
}

###################################################################
###################################################################
### SECTION 04 - DESEQ2 RETURN FUNCTIONS (DESEQ2-RF)
###################################################################
###################################################################

## DESeq2 - EXP1 - Two group comparisons
deseq.exp1 <- function(fact, coldata, cts, perm.h, batchVar = NULL, w1 = NULL) {
  
  if(is.null(batchVar) & is.null(w1)) {
    design0 <- model.matrix(~ 0 + coldata[, fact])
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    design <- paste0("~ ", fact)
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    design0 <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, fact])
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    design <- paste0("~ W_1 + ", fact)
  } else if(is.null(w1)) {
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(coldata[, fact])
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    coldata[, fact] <- relevel(coldata[, fact], ref = refLevel)
    design0 <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, fact])
    design <- paste0("~ ", batchVar, " + ", fact)
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(coldata[, fact])
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    coldata[, fact] <- relevel(coldata[, fact], ref = refLevel)
    design0 <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact])
    batchString <- paste0("W_1 + ", batchVar, " + ")
    design <- paste0("~ W_1 + ", batchVar, " + ", fact)
  }
  
  rownames(design0) <- rownames(coldata)
  colnames(design0) <- gsub("coldata\\[, fact\\]*", replacement = "", x = colnames(design0))
  colnames(design0) <- gsub("coldata\\[, batchVar\\]*", replacement = "", x = colnames(design0))
  colnames(design0) <- gsub("coldata\\[, \"W_1\"\\]*", replacement = "W_1", x = colnames(design0))
  
  # Create syntactically valid names for design0 and perm.c, then replace with original
  design0Levels <- colnames(design0)
  colnames(design0) <- paste0("X", colnames(design0))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h
  rownames(cont0) <- design0Levels
  
  design <- as.formula(design)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design
  )
  dds <- DESeq(dds)
  return(list(dds, cont0))
}

## DESeq2 - EXP2 - multiple factor comparisons
deseq.exp2 <- function(fact1, fact2, coldata, cts, perm.h, batchVar = NULL, w1 = NULL) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  if(is.null(batchVar) & is.null(w1)) {
    batchString <- ""
    design <- model.matrix(~ 0 + group.c)
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- "W_1 + "
    design <- model.matrix(~ 0 + coldata[, "W_1"] + group.c)
    perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  } else if(is.null(w1)) {
    batchString <- paste0(batchVar, " + ")
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(group.c)
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    group.c <- relevel(group.c, ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, batchVar] + group.c)
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- paste0("W_1 + ", batchVar, " + ")
    group1 <- strsplit(perm.h, split = "_VS_")[[1]][1]
    group2 <- strsplit(perm.h, split = "_VS_")[[1]][2]
    refLevel <- group2
    perm.c <- group1
    factLevels <- levels(group.c)
    if(length(factLevels) >= 3) {
      refLevel <- factLevels[!(factLevels %in% c(group1, group2))][1]
      perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
    }
    group.c <- relevel(group.c, ref = refLevel)
    design <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + group.c)
  }
  colnames(design) <- gsub("group.c", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", replacement = "", x = colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", replacement = "W_1", x = colnames(design))
  rownames(design) <- rownames(coldata)
  
  # Create syntactically valid names, then replace with original
  designLevels <- colnames(design)
  colnames(design) <- paste0("X", colnames(design))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  rownames(cont) <- designLevels
  
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ 1
  )
  dds$group <- group.c
  design(dds) <- as.formula(paste0("~ ", batchString, "group"))
  dds <- DESeq(dds)
  return(list(dds, cont))
}

## DESeq2 - EXP3 - classical interactions
deseq.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl, 
                       batchVar = NULL, w1 = NULL) {
  
  f1n <- length(levels(coldata[, fact1]))
  nBatchCols <- 0
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  batchString <- ""
  
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- "W_1 + "
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- 1
    f1n <- f1n + nBatchCols
  } else if(is.null(w1)) {
    batchString <- paste0(batchVar, " + ")
    design <- model.matrix(~ coldata[, batchVar] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
    f1n <- f1n + nBatchCols - 1
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- paste0("W_1 + ", batchVar, " + ")
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact1] * coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
    f1n <- f1n + nBatchCols
  }
  
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[(nBatchCols + 2):ncol(tmp1)] <- paste0(colnames(tmp1)[(nBatchCols + 2):ncol(tmp1)], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)

  design.dds <- paste0("~ ", batchString, fact1, " * ", fact2)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds)
  return(list(dds, design))
}

## DESeq2 - EXP4 - added effects blocking and paired
deseq.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl, 
                       batchVar = NULL, w1 = NULL) {
  
  f1n <- length(levels(coldata[, fact1]))
  nBatchCols <- 0
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  batchString <- ""
  
  if(is.null(batchVar) & is.null(w1)) {
    design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- "W_1 + "
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- 1
  } else if(is.null(w1)) {
    batchString <- paste0(batchVar, " + ")
    design <- model.matrix(~ coldata[, batchVar] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar])) - 1
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- paste0("W_1 + ", batchVar, " + ")
    design <- model.matrix(~ coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact1] + coldata[, fact2])
    nBatchCols <- length(levels(coldata[, batchVar]))
  }
  
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[(1 + nBatchCols + f1n):ncol(design)] <- paste0(
    colnames(design)[(1 + nBatchCols + f1n):ncol(design)],
    "_VS_", 
    fact2.rlvl
  )

  design.dds <- paste0("~ ", batchString, fact1, " + ", fact2)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds)
  # if(nBatchCols > 0) design <- design[, -2:-(nBatchCols + 1), drop = FALSE]
  return(list(dds, design))
}

## DESeq2 - EXP5 - main effects only
deseq.exp5 <- function(fact, fact.levl, cts, coldata, perm.h, batchVar = NULL,
                       w1 = NULL) {
  
  coldata[, fact] <- relevel(x = coldata[, fact], ref = fact.levl)
  # Changed this from original because it was giving very different results 
  # compared to batch-corrected. Originally perm.c was just the comparison 1 level
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  # perm.c <- paste0(perm.h, "-", fact.levl)
  batchString <- ""
  reducedModel <- as.formula("~1")
  
  if(is.null(batchVar) & is.null(w1)) {
    design0 <- model.matrix(~ 0 + coldata[, fact])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- "W_1 + "
    reducedModel <- as.formula("~ W_1")
    design0 <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, fact])
  } else if(is.null(w1)) {
    batchString <- paste0(batchVar, " + ")
    reducedModel <- as.formula(paste0("~ ", batchVar))
    design0 <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, fact])
    perm.c <- perm.h
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- paste0("W_1 + ", batchVar, " + ")
    reducedModel <- as.formula(paste0("~ W_1 + ", batchVar))
    design0 <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, fact])
    perm.c <- perm.h
  }
  
  colnames(design0) <- gsub("coldata\\[, fact\\]*", "", colnames(design0))
  colnames(design0) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design0))
  colnames(design0) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design0))
  
  # Create syntactically valid names, then replace with original
  design0Levels <- colnames(design0)
  colnames(design0) <- paste0("X", colnames(design0))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h
  rownames(cont0) <- design0Levels
  
  design.dds <- paste0("~ ", batchString,  fact)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds, test = "LRT", reduced = reducedModel)
  return(list(dds, cont0))
}

## DESeq2 - EXP6 - main effects only + grouping factor
deseq.exp6 <- function(
  me.fact, me.levl, gp.fact, gp.levl, cts, coldata, perm.h, batchVar = NULL,
  w1 = NULL) {
  
  cts <- cts[, which(coldata[, gp.fact] == gp.levl)]
  coldata <- coldata[which(coldata[, gp.fact] == gp.levl), ]
  coldata[] <- lapply(coldata, function(x) if(is.factor(x)) factor(x) else x)
  coldata[, me.fact] <- relevel(x = coldata[, me.fact], ref = me.levl)
  batchString <- ""
  reducedModel <- as.formula("~1")
  
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  if(is.null(batchVar) & is.null(w1)) {
    design0 <- model.matrix(~ 0 + coldata[, me.fact])
  } else if(is.null(batchVar)) {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- "W_1 + "
    reducedModel <- as.formula("~ W_1")
    design0 <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, me.fact])
  } else if(is.null(w1)) {
    batchString <- paste0(batchVar, " + ")
    reducedModel <- as.formula(paste0("~ ", batchVar))
    design0 <- model.matrix(~ 0 + coldata[, batchVar] + coldata[, me.fact])
    perm.c <- perm.h
  } else {
    coldata$W_1 <- w1[rownames(coldata)]
    batchString <- paste0("W_1 + ", batchVar, " + ")
    reducedModel <- as.formula(paste0("~ W_1 + ", batchVar))
    design0 <- model.matrix(~ 0 + coldata[, "W_1"] + coldata[, batchVar] + coldata[, me.fact])
    perm.c <- perm.h
  }
  
  rownames(design0) <- rownames(coldata)
  colnames(design0) <- gsub("coldata\\[, me.fact\\]*", "", colnames(design0))
  colnames(design0) <- gsub("coldata\\[, batchVar\\]*", "", colnames(design0))
  colnames(design0) <- gsub("coldata\\[, \"W_1\"\\]*", "W_1", colnames(design0))
  
  # Create syntactically valid names, then replace with original
  design0Levels <- colnames(design0)
  colnames(design0) <- paste0("X", colnames(design0))
  perm.c <- paste0("X", sub(pattern = "-", replacement = "-X", x = perm.c))
  
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h
  rownames(cont0) <- design0Levels

  design.dds <- paste0("~ ", batchString, me.fact)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds, test = "LRT", reduced = reducedModel)
  return(list(dds, cont0))
}

## DESeq2 - EXP7 - user input
deseq.exp7 <- function(cts, coldata, mod.matrix) {
  design <- mod.matrix
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ 1
  )
  dds <- DESeq(dds, full = design, modelMatrixType = "standard")
  return(list(dds, design))
}

###################################################################
###################################################################
### SECTION 05 - PLOTTING FUNCTIONS (PLOT)
###################################################################
###################################################################

#------------------------------------------------------------------
# Heatmapping
#------------------------------------------------------------------

# PLOT - heatmap function for download
qcHeatMap <- function(heat, color, rows, col, nam, pal, n = NULL, keyLabel = "Row z-score",
                      colorLabel = "Sample", colorbarLabel = F, type = "pdf") {
  if(is.null(n)) n = nrow(heat)
  n <- as.numeric(n)

  if(is.null(colorLabel) || colorLabel == "") colorLabel <- "Sample"
  row_ha <- HeatmapAnnotation(Sample = nam, col = list(Sample = c(pal)))
  row_ha@anno_list$Sample@name <- colorLabel
  row_ha@anno_list$Sample@label <- ifelse(test = (colorbarLabel), yes = colorLabel, no = "")
  row_ha@anno_list$Sample@color_mapping@levels <- rev(row_ha@anno_list$Sample@color_mapping@levels)
  names(row_ha) <- colorLabel
  
  showRowNames <- T
  if(nrow(heat) > 50) showRowNames <- F
  
  ht_opt(heatmap_column_title_gp = gpar(fontfamily = "noto-sans-jp", fontsize = 8),
         heatmap_row_title_gp = gpar(fontfamily = "noto-sans-jp", fontsize = 8),
         legend_title_gp = gpar(fontfamily = "noto-sans-jp", fontsize = 8),
         legend_labels_gp = gpar(fontfamily = "noto-sans-jp", fontsize = 8))

  if (type == "pdf") {
    p <- Heatmap(
      matrix = heat,
      col = color,
      name = keyLabel,
      cluster_rows = rows,
      cluster_columns = col,
      show_row_names = showRowNames,
      show_column_names = F,
      top_annotation = row_ha,
      row_names_gp = gpar(fontsize = 8)
    )
    
    myPadding <- unit(x = c(2, 2, 2, 2), units = "mm")
    if(nchar(colorLabel) > 20) {
      myPadding <- unit(x = c(2, 2, 2, 2 + (nchar(colorLabel) - 20) * 2), units = "mm")
    }
    
    draw(p, merge_legend = TRUE, padding = myPadding)
  } else {
    hmWidth <- unit(0.25, "cm")*nrow(heat)
    if(nrow(heat) < 50 || nrow(heat) > 125) hmWidth <- unit(1, "npc")
    p <- Heatmap(
      matrix = heat,
      col = color,
      name = keyLabel,
      cluster_rows = rows,
      cluster_columns = col,
      show_row_names = showRowNames,
      show_column_names = F,
      top_annotation = row_ha,
      row_names_gp = gpar(fontsize = 8),
      heatmap_width = hmWidth
    )
    
    myPadding <- unit(x = c(10, 10, 10, 10), units = "mm")
    if(nchar(colorLabel) > 20) {
      myPadding <- unit(x = c(10, 10, 10, 10 + (nchar(colorLabel) - 20) * 10), units = "mm")
    }
    draw(p, merge_legend = TRUE, padding = myPadding)
  }
}

# PLOT - heatmap count boxplots for download
qcHeatCount <- function(data, fact, var, xaxis, title, type) {
  if(!is.factor(fact)) fact <- factor(fact, levels = sort(unique(fact), decreasing = F))
  
  if (type == "bulk") {
    p <- ggplot(data, aes(x = fact, y = var, fill = fact)) +
      scale_fill_manual(values = hue_pal()(length(levels(fact)))) +
      geom_boxplot() +
      stat_boxplot(geom ='errorbar') + 
      geom_jitter(position = position_jitter(0.3)) +
      xlab(xaxis) +
      ylab("Normalized counts") +
      ggtitle(paste(title, "Counts")) +
      theme_light() +
      theme(text = element_text(family = "noto-sans-jp"),
            axis.text.x = element_text(size = 14, angle = -45, vjust = .55, hjust = 0),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.title = element_blank())
  } else {
    p <- ggplot(data, aes(x = fact, y = var, fill = fact)) +
      scale_fill_manual(values = rev(hue_pal()(length(levels(fact))))) +
      geom_boxplot() +
      stat_boxplot(geom ='errorbar') + 
      geom_jitter(position = position_jitter(0.3)) +
      xlab(xaxis) +
      ylab("Normalized counts") +
      ggtitle(paste(title, "Counts")) +
      theme_light() +
      theme(text = element_text(family = "noto-sans-jp"),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.title = element_blank())
  }
  print(p)
}

#------------------------------------------------------------------
# Differential gene expression (DGE)
#------------------------------------------------------------------

# PLOT - overview plot (bulk)
dgeOverPlot <- function(comp) {
  p <- ggplot(comp, aes(contrast, value)) +
    geom_bar(
      stat = "identity", 
      aes(fill = variable), 
      color = "#3d3d3d",
      position = position_dodge()
    ) +
    xlab("") +
    ylab("Number of IDs") +
    ggtitle("DGE Comparisons Overview") +
    theme_light() +
    scale_fill_manual(
      name = "Expression",
      values = hue_pal()(2)
    ) +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
}

# PLOT - ma plot (bulk)
dgeMAPlot <- function(dgeout2, p, l, cont) {
  group1 <- strsplit(cont, split = "_VS_")[[1]][1]
  group2 <- strsplit(cont, split = "_VS_")[[1]][2]
  if(!is.na(group2)) {
    plotTitle <- paste0(group1, " vs. ", group2)
  } else {
    plotTitle <- group1
  }
  
  dge <- dgeout2 %>%
    mutate(color = ifelse(log2FoldChange < 0, paste0("down in ", group1),
                          ifelse(log2FoldChange > 0, paste0("up in ", group1), NA)))
  dge$isPADJ <- ifelse(dge$padj <= p, TRUE, FALSE)
  dge$isPADJ[is.na(dge$isPADJ)] <- FALSE
  dge$isLFC <- ifelse(abs(dge$log2FoldChange) >= l, TRUE, FALSE)
  dge$isDGE <- "No differential expression"
  dge$isDGE[which(dge$isLFC & dge$isPADJ)] <- paste(
    "padj <=", p, "& LFC >=", l
  )
  dge$output <- ifelse(grepl("padj", dge$isDGE) & grepl("down", dge$color), "down",
                  ifelse(grepl("padj", dge$isDGE) & grepl("up", dge$color), "up", "No differential expression"))
  
  myColors <- c()
  if(sum(dge$output == "down") > 0) myColors <- c(myColors, "green")
  if(sum(dge$output == "No differential expression") > 0) myColors <- c(myColors, "gray")
  if(sum(dge$output == "up") > 0) myColors <- c(myColors, "red")
  
  p <- ggplot(dge, aes(log10(baseMean), log2FoldChange)) +
    geom_point(size = 0.7, aes(color = output)) +
    xlab(bquote("log"[10]*"(base mean)")) +
    ylab(bquote("log"[2]*"(fold change)")) +
    ggtitle(paste0("MA Plot: ", plotTitle)) +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(
      name = "",
      values = myColors
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    geom_hline(
      yintercept = 0, 
      color = "#3d3d3d", 
      linetype = "dashed"
    )
  print(p)
}


# PLOT - volcano plot (bulk)
dgeVolPlot <- function(dgeout2, p, l, cont) {
  group1 <- strsplit(cont, split = "_VS_")[[1]][1]
  group2 <- strsplit(cont, split = "_VS_")[[1]][2]
  if(!is.na(group2)) {
    plotTitle <- paste0(group1, " vs. ", group2)
  } else {
    plotTitle <- group1
  }
  
  dge <- dgeout2 %>%
    mutate(color = ifelse(log2FoldChange < 0, paste0("down in ", group1),
                          ifelse(log2FoldChange > 0, paste0("up in ", group1), NA)))
  dge$isPADJ <- ifelse(dge$padj <= p, TRUE, FALSE)
  dge$isPADJ[is.na(dge$isPADJ)] <- FALSE
  dge$isLFC <- ifelse(abs(dge$log2FoldChange) >= l, TRUE, FALSE)
  dge$isDGE <- "No differential expression"
  dge$isDGE[which(dge$isLFC & dge$isPADJ)] <- paste(
    "padj <=", p, "& LFC >=", l
  )
  dge$output <- ifelse(grepl("padj", dge$isDGE) & grepl("down", dge$color), "down",
                       ifelse(grepl("padj", dge$isDGE) & grepl("up", dge$color), "up", "No differential expression"))
  
  myColors <- c()
  if(sum(dge$output == "down") > 0) myColors <- c(myColors, "green")
  if(sum(dge$output == "No differential expression") > 0) myColors <- c(myColors, "gray")
  if(sum(dge$output == "up") > 0) myColors <- c(myColors, "red")
  
  p <- ggplot(dge, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(size = 0.6, aes(color = output)) +
    xlab(bquote("log"[2]*"(fold change)")) +
    ylab(bquote("-log"[10]*"(p-value)")) +
    ggtitle(paste0("Volcano Plot: ", plotTitle)) +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(
      name = "",
      # values = c("green", "gray", "red")
      values = myColors
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    geom_vline(
      xintercept = 0, 
      color = "#3d3d3d", 
      linetype = "dashed"
    )
  print(p)
}

#------------------------------------------------------------------
# WGCNA
#------------------------------------------------------------------

# PLOT - use ggplot to plot WGCNA sample dendrogram
GGDend <- function(clust, colorDF, leafLabels = T,
                   plotTitle = NULL, touchZero = T,
                   showBarNames = T, 
                   # blank = F,
                   axisLines = T, 
                   legendDF = NULL,
                   plotTitleSize = 20) {
  
  for(i in colnames(colorDF)) colorDF[, i] <- as.character(colorDF[, i])
  plotDF <- colorDF
  plotDF$sample_temp <- rownames(plotDF)
  p1_dendroData <- dendro_data(clust)
  
  if(!touchZero) {
    hangLength <- 0.01*(max(p1_dendroData$segments$yend)-min(p1_dendroData$segments$yend))
    p1_dendroData$segments$yend[p1_dendroData$segments$yend == 0] <- 
      p1_dendroData$segments$y[p1_dendroData$segments$yend == 0] - hangLength
    
  }
  yMin <- min(p1_dendroData$segments$yend)-(max(p1_dendroData$segments$yend)-min(p1_dendroData$segments$yend))*0.05
  yMax <- yMax_data <- max(p1_dendroData$segments$yend)
  
  if(!is.null(legendDF)) {
    legendDF_unique <- as.data.frame(matrix(NA, 
                                            nrow = max(apply(legendDF, 2, function(i) length(unique(i)))),
                                            ncol = ncol(legendDF)))
    colnames(legendDF_unique) <- colnames(legendDF)
    colorDF_unique <- legendDF_unique
    for(colname in colnames(legendDF)) {
      legendDF_unique[1:length(unique(legendDF[, colname])), colname] <- sort(unique(as.character(legendDF[, colname])))
      colorDF_unique[1:length(unique(legendDF[, colname])), colname] <- 
        colorDF[match(sort(unique(legendDF[, colname])),
                      table = legendDF[, colname]),
                colname]
    }
    colMat <- as.matrix(colorDF_unique)
    tt <- ttheme_default(base_size = 12, base_family = "noto-sans-jp", 
                         core = list(bg_params = list(fill = colMat)))
    legendDF_unique = replace(legendDF_unique, list = is.na(legendDF_unique), values = "")
    gTable <- tableGrob(legendDF_unique, rows = NULL, theme = tt)
    yMax <- yMax + (yMax - yMin) * 0.75
  }
  
  yStepChoices <- c(0.25)
  for(i in c(0.25,0.20,0.10,0.05,0.01)) {
    yStep <- i
    yBreaks <- seq(0,yMax_data, i)
    yBreaks <- yBreaks[yBreaks >= yMin & yBreaks <= yMax_data]
    if(length(yBreaks) >= 3) break
  }
  
  theme_set(theme_classic(base_family = "noto-sans-jp", 
                          base_size = 12))
  # Change the settings
  update_geom_defaults("text", list(family = theme_get()$text$family))
  
  p1_dendro <- ggdendrogram(p1_dendroData, labels = T, leaf_labels = leafLabels, theme_dendro = F) +
    coord_cartesian(xlim = c(-1, nrow(plotDF) + 1), clip = "off", 
                    expand = F) +
    scale_y_continuous("Height\n", breaks = yBreaks, limits = c(yMin,yMax), 
                       expand = c(0,0)) +
    theme(text = element_text(family = "noto-sans-jp"),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = unit(c(1,1,0,1), "cm"),
      plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      )
  if(!leafLabels) {
    theme_set(theme_classic(base_family = "noto-sans-jp", 
                            base_size = 12))
    # Change the settings
    update_geom_defaults("text", list(family = theme_get()$text$family))
    p1_dendro <- p1_dendro +
      theme(axis.text.x = element_blank())
  }
  if(axisLines) {
    theme_set(theme_classic(base_family = "noto-sans-jp", 
                            base_size = 12))
    # Change the settings
    update_geom_defaults("text", list(family = theme_get()$text$family))
    p1_dendro <- p1_dendro +
      theme(panel.border = element_blank(), 
            axis.line.x = element_line(size = 0.75),
            axis.line.y = element_blank()) +
      annotate("segment", x = -1, xend = -1, y = yMin, yend = yMax_data, size = 0.75) 
  }
  
  if(!is.null(legendDF)) {
    theme_set(theme_classic(base_family = "noto-sans-jp", 
                            base_size = 12))
    # Change the settings
    update_geom_defaults("text", list(family = theme_get()$text$family))
    p1_dendro <- p1_dendro +
      annotation_custom(gTable, 
                        xmin = 0, 
                        xmax =  nrow(legendDF) * 0.05 * ncol(legendDF),
                        ymin = yMax_data * 1.02)
  }
  
  plotList <- list(p1 = p1_dendro)
  for(colname in colnames(colorDF)) {
    yLab <- NULL
    if(showBarNames) yLab <- colname
    
    theme_set(theme_classic(base_family = "noto-sans-jp", 
                            base_size = 12))
    # Change the settings
    update_geom_defaults("text", list(family = theme_get()$text$family))
    
    p <- ggplot(plotDF, aes(sample_temp, y = 1)) + 
      geom_tile(aes(fill = sample_temp)) + theme_minimal() +
      scale_fill_manual(values = plotDF[as.character(p1_dendroData$labels$label), colname]) +
      coord_cartesian(xlim = c(-1, nrow(plotDF) + 1), expand = F) +
      ylab(yLab) +
      theme(
            text = element_text(family = "noto-sans-jp", 
                                size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            line = element_blank(),
            plot.margin = unit(c(0,1,0,1), "cm")) 
    plotList[[colname]] <- p
  }
  
  heights <- c(ncol(colorDF)*4, rep.int(1, times = ncol(colorDF)))
  plot_grid(plotlist = plotList, ncol = 1, align = "v", axis = "tblr",
            rel_heights = heights)
}

# PLOT - use ggplot2 to create TOM plot
GGTom <- function(clust, colors, distMat, plotTitle = NULL,
                  plotTitleSize = 20) {
  plotDF <- data.frame(module = colors, sample = 1:length(colors))
  for(i in 1:ncol(plotDF)) plotDF[, i] <- as.character(plotDF[, i])
  p1_dendroData <- dendro_data(clust)
  p3_dendroData <- p1_dendroData
  
  hangLength <- 0.01*(max(p1_dendroData$segments$yend)-min(p1_dendroData$segments$yend))
  p1_dendroData$segments$yend[p1_dendroData$segments$yend == 0] = 
    p1_dendroData$segments$y[p1_dendroData$segments$yend == 0] - hangLength
  
  yMin <- min(p1_dendroData$segments$yend)-(max(p1_dendroData$segments$yend)-min(p1_dendroData$segments$yend))*0.05
  yMax <- max(p1_dendroData$segments$yend)*1.05
  
  theme_set(theme_classic(base_family = "noto-sans-jp", 
                          base_size = 12))
  # Change the settings
  update_geom_defaults("text", list(family = theme_get()$text$family))
  
  p1_dendroX <- ggdendrogram(p1_dendroData, labels = F, leaf_labels = F, theme_dendro = T) +
    theme_void() +
    coord_cartesian(xlim = c(-1, nrow(plotDF) + 1),
                    ylim = c(yMin, yMax),
                    expand = F) +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          line = element_blank(),
          plot.margin=unit(c(0, 0, 0, 0), "cm")) 
  
  p2_colorBarX <- ggplot(plotDF,aes(sample,y=1)) + 
    theme_void() +
    geom_tile(aes(fill = sample)) + 
    scale_fill_manual(values = plotDF[as.character(p1_dendroData$labels$label), "module"]) +
    coord_cartesian(xlim = c(-1, nrow(plotDF) + 1), expand = F) +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  
  p3_dendroData$segments$yend[p3_dendroData$segments$yend == 0] = 
    p3_dendroData$segments$y[p3_dendroData$segments$yend == 0] - hangLength
  
  p3_dendroData$segments <- data.frame(x = max(p3_dendroData$segments$y) - p3_dendroData$segments$y,
                                       y = max(p3_dendroData$segments$x) - p3_dendroData$segments$x,
                                       xend = max(p3_dendroData$segments$yend) - p3_dendroData$segments$yend,
                                       yend = max(p3_dendroData$segments$xend) - p3_dendroData$segments$xend)
  xMin <- min(p3_dendroData$segments$xend)-(max(p3_dendroData$segments$xend)-min(p3_dendroData$segments$xend))*0.05
  xMax <- max(p3_dendroData$segments$xend)*1.05  
  p3_dendroY <- ggdendrogram(p3_dendroData, labels = F, leaf_labels = F, theme_dendro = T, rotate = F) +
    theme_void() +
    coord_cartesian(ylim = c(-1, nrow(plotDF) + 1),
                    xlim = c(xMin, xMax),
                    expand = F) +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p4_colorBarY <- ggplot(plotDF,aes(x=1,y=sample)) + 
    theme_void() +
    geom_tile(aes(fill = sample)) + 
    scale_fill_manual(values = rev(plotDF[as.character(p1_dendroData$labels$label), "module"])) +
    coord_cartesian(ylim = c(-1, nrow(plotDF) + 1), expand = F) +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  
  m <- reshape2::melt(distMat[clust$order, rev(clust$order)]^4)
  p5_heatmap <- ggplot(data = m, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    theme_void() +
    scale_fill_gradientn(colors = hcl.colors(12, "YlOrRd", rev = T)) +
    coord_cartesian(xlim = c(-1, nrow(plotDF) + 1),
                    ylim = c(-1, nrow(plotDF) + 1),
                    expand = F) +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  if(is.null(plotTitle)) {
    plot_grid(NULL, NULL, p1_dendroX, 
              NULL, NULL, p2_colorBarX, 
              p3_dendroY, p4_colorBarY, p5_heatmap, 
              ncol = 3,  
              align = "hv",
              axis = "l",
              rel_widths = c(1.25,0.35,3),
              rel_heights = c(1.25,0.35,3)
              )
    
  } else {
    title_gg <- ggplot() + 
      theme_void() + 
      coord_cartesian(xlim = c(xMin, xMax),
                      expand = F) +
      theme(plot.title = element_text(hjust = 0, 
                                      size = plotTitleSize),
            plot.margin = unit(c(0,0,0,0), "cm"))
    plot_grid(NULL, NULL, p1_dendroX, 
              NULL, NULL, p2_colorBarX, 
              p3_dendroY, p4_colorBarY, p5_heatmap, 
              ncol = 3, 
              align = "none",
              rel_widths = c(1.25,0.35,3),
              rel_heights = c(1.25,0.35,3)) +
      draw_label(plotTitle, y = 0.97, hjust = 0, 
                 size = plotTitleSize)
  }
}

# Helper function borrowed from scClustViz for setting up colors for 
# dendrogram
ColorDend <- function(dend, clusters, genes, cols) {
  dC <- dendrapply(dend, function(X) {
    if (is.leaf(X)) {
      attr(X,"edgePar") <- list(
        lwd=2,
        col=cols[which(attr(X,"label") == levels(clusters))]
      )
      attr(X,"nodePar") <- list(
        pch=NA,lab.font=2,lab.cex= 1,
        lab.col=cols[which(attr(X,"label") == levels(clusters))])
      if (attr(X,"label") != "Unselected") {
        if (attr(X,"label") %in% names(genes)) {
          attr(X,"label") <- paste0(attr(X,"label"),": ",
                                    length(genes[[attr(X,"label")]])," DE")
        } else {
          attr(X,"label") <- paste0(
            attr(X,"label"),": ",
            length(genes[[which(attr(X,"label") ==
                                    sapply(strsplit(names(genes),"-"),
                                           function(X) X[1]))]]),
            " DE")
        }
      }
    }
    return(X)
  })
}

# DGE - dge dotplot per cluster (sc only)
plot_deDotplot_ggplot <- function(
  scvData, deGenes, nGenes = 5, clust = 1, sizeFactor = 1, drData = NULL,
  mdgeData = NULL, hcData = NULL) {
  
  heatGenes <- names(deGenes[[as.character(clust)]])[!is.na(names(deGenes[[as.character(clust)]]))]
  heatGenes <- heatGenes[1:min(length(heatGenes), nGenes)]
  if(is.null(heatGenes) || length(heatGenes) == 0) return()
  
  # sizeFactor2 is for text labels. sizeFactor3 is for dots. These values are
  # scaled based on the number of genes/columns

  sizeFactor2 <- sizeFactor * (1 - (length(heatGenes)/5 - 1)/30)
  sizeFactor3 <- 2 /(length(heatGenes)/5)^0.6

  statsList <- ClustGeneStats(scvData)
  
  if(is.null(drData)) {
    drData <- as.data.frame(sapply(statsList, function(x) x$DR), row.names = rownames(statsList[[1]]))
  }
  if(is.null(mdgeData)) {
    mdgeData <- as.data.frame(sapply(statsList, function(x) x$MDGE), row.names = rownames(statsList[[1]]))
  }
  
  drData <- as.matrix(drData[heatGenes, , drop = F])
  mdgeData <- as.matrix(mdgeData[heatGenes, , drop = F])
  
  hG <- list(order = 1)
  if(nrow(drData) > 1) hG <- hclust(d = stats::dist(drData), method = "complete")
  
  hC <- hcData
  if(is.null(hC)) hC <- hclust(d = stats::as.dist(DEdist(scvData)), method = "single")
  
  myClusters <- Clusters(scvData)
  clustCols <- colorspace::qualitative_hcl(n = length(levels(myClusters)), palette = "Dark 3")
  clusterDend <- as.dendrogram(hC)
  
  dC <- ColorDend(dend = clusterDend, clusters = myClusters, genes = deGenes, cols = clustCols)
  dC_cols <- get_leaves_branches_col(dC)
  dC_labels <- labels(dC)
  
  drPlotDF <- as.data.frame(drData[hG$order, hC$order, drop = F]) %>% 
    rownames_to_column(var = "gene") %>%
    gather(key = cluster, value = dr, -gene)
  drPlotDF$cluster <- factor(drPlotDF$cluster, levels = hC$labels[hC$order])
  drPlotDF$dr <- drPlotDF$dr * 100
  
  mdgePlotDF <- as.data.frame(mdgeData[hG$order, hC$order, drop = F]) %>%
    rownames_to_column(var = "gene") %>%
    gather(key = cluster, value = "mdge", -gene)
  mdgePlotDF$cluster <- factor(mdgePlotDF$cluster, levels = hC$labels[hC$order])
  
  plotDF <- merge(drPlotDF, mdgePlotDF)
  
  geneNums <- 1:length(heatGenes)
  names(geneNums) <- heatGenes[hG$order]
  clustNums <- 1:length(levels(plotDF$cluster))
  names(clustNums) <- levels(plotDF$cluster)
  
  plotDF$x <- geneNums[plotDF$gene]
  plotDF$y <- clustNums[as.character(plotDF$cluster)]
  
  minMDGE <- min(plotDF$mdge)
  maxMDGE <- max(plotDF$mdge)
  
  geneLabAngle <- 0
  hJust <- 0.5
  if(length(heatGenes) > 8) {
    geneLabAngle <- 90
    hJust <- 1
  }
  
  dotPlot <- ggplot(data = plotDF) +
    geom_vline(xintercept = geneNums, color = "lightgrey", ) +
    geom_point(mapping = aes(x = x, y = y, size = dr, color = mdge)) +
    scale_x_continuous(breaks = geneNums, labels = names(geneNums)) +
    scale_y_continuous(
      breaks = clustNums,
      labels = dC_labels
    ) +
    scale_size_continuous(
      name = "Detection rate", 
      breaks = c(25, 50, 75),
      labels = c("25%", "50%", "75%"), 
      range = c(1, 6) * sizeFactor3,
      guide = guide_legend(title.position = "top", title.hjust = 0)
    ) +
    scale_color_continuous(
      name = "Mean detected expression",
      breaks = c(minMDGE, maxMDGE),
      labels = c(round(minMDGE, digits = 2), round(maxMDGE, digits = 2)),
      type = "viridis",
      guide = guide_colorbar(title.position = "top", title.hjust = 0, barwidth = 10),
      direction = -1
    ) +
    theme_nothing() +
    theme(
      legend.position = "right",
      legend.direction = "horizontal",
      legend.text = element_text(size = 15*sizeFactor),
      legend.title = element_text(size = 20*sizeFactor),
      axis.text.x = element_text(
        size = 15*sizeFactor2, margin = margin(t = 5), angle = geneLabAngle, hjust = hJust,
      ),
      axis.text.y = element_text(size = 15*sizeFactor, color = dC_cols, hjust = 0),
      text = element_text(family = "noto-sans-jp"),
      plot.margin = margin(c(0, 1, 1, 1), unit = "cm")
    )
  
  clusterDendPlot <- dC %>%
    set("branches_lwd", 0.5) %>%
    ggplot(horiz = T, labels = F, offset_labels = -2) +
    theme(plot.margin = margin(c(0,0,0,0), unit = "cm"))
  
  myWidths <- myHeights <- c(1, 8)
  if(nrow(drData) < 5) myWidths <- c(1, 1)
  
  if(nrow(drData) > 1) {
    geneDendPlot <- hG %>% as.dendrogram() %>% set("branches_lwd", 0.5) %>% 
      ggplot(labels = F)
    plot_spacer() + geneDendPlot + clusterDendPlot + dotPlot +
      plot_layout(widths = myWidths, heights = myHeights)
  } else {
    clusterDendPlot + dotPlot + plot_layout(widths = myWidths)
  }
}

###################################################################
###################################################################
### SECTION 06 - DIFFERENTIAL GENE EXPRESSION (DGE)
###################################################################
###################################################################

# DGE - get dge overview table (bulk)
dgeOverTbl <- function(cont.ls, lf, p) {
  p <- as.numeric(p)
  lf <- as.numeric(lf)
  for (i in 1:length(cont.ls)) {
    cont.ls[[i]] <- cont.ls[[i]][abs(cont.ls[[i]]$log2FoldChange) >= lf, ]
    cont.ls[[i]] <- cont.ls[[i]][cont.ls[[i]]$padj <= p, ]
  }
  cont.regup <- lapply(cont.ls, subset, log2FoldChange > 0)
  cont.regdn <- lapply(cont.ls, subset, log2FoldChange < 0)
  cont.regup <- lapply(cont.regup, nrow)
  cont.regup <- unlist(cont.regup)
  cont.regdn <- lapply(cont.regdn, nrow)
  cont.regdn <- unlist(cont.regdn)
  up.df <- data.frame(cont.regup)
  dn.df <- data.frame(cont.regdn)
  al.df <- merge(up.df, dn.df, by = 0, all = TRUE)
  colnames(al.df) <- c("contrast", "up", "down")
  al.df <- melt(al.df)
  al.df$contrast <- factor(al.df$contrast)
  al.df$variable <- factor(al.df$variable)
  return(al.df)
}

# DGE - check if factors are full rank (no pair of variables are collinear). 
# required for DESeq2; returns NULL if full rank and/or variables that are
# confounded
CheckFullRank <- function(df) {
  if(ncol(df) < 2) return(NULL)
  myList <- vector("list", length = ncol(df))
  names(myList) <- colnames(df)
  for(i in 1:(ncol(df)-1)) {
    for(j in (i+1):ncol(df)) {
      colname1 <- colnames(df)[i]
      colname2 <- colnames(df)[j]
      nUnique1 <- length(unique(df[,colname1]))
      nUnique2 <- length(unique(df[, colname2]))
      catCols <- paste(df[, colname1], df[, colname2], sep = "_")
      if(length(unique(catCols)) < (max(nUnique1,nUnique2) + 1)) myList[[colname1]] <- c(myList[[colname1]], colname2)
    }
  }
  myList <- myList[!(sapply(myList, is.null))]
  if(length(myList) > 0) {
    return(myList)
  } else {
    return(NULL)
  }
}

# Function to check whether there is at least 1 sample in each combination of levels of 
# 1 or more factors in dataframe. Used to check whether DGE can be run based on
# input selections.
CheckMultiFactorLevels <- function(df, factLevelsList, minCount = 1) {
  dfTmp <- df[, names(factLevelsList), drop = F]
  combs <- apply(dfTmp, MARGIN = 1, function(i) paste(i, collapse = "_"))
  requiredCombs <- apply(expand.grid(factLevelsList), MARGIN = 1, function(i) paste(i, collapse = "_"))
  countTable <- table(combs)
  return(all(requiredCombs %in% combs) && all(countTable[requiredCombs] >= minCount))
}

###################################################################
###################################################################
### SECTION 07 - GENESET ENRICHMENT (GSE)
###################################################################
###################################################################

# GSE - filter to return significant gse
getSigTerms <- function(enr, libs) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullSigRes <- data.frame(libName = character(), lib_rank = integer(), 
                           gene_count = integer(), term = character(), 
                           overlap = character(), pval = double(), 
                           adjPval = double(), oldPval = double(), 
                           oldAdjPval = double(), Z_score = double(),
                           score = double(), gene_list = character(), stringsAsFactors = FALSE)
  
  res <- lapply(1:length(enr), function(x) {
    res <- data.frame(enr[[x]])
    sig <- lapply(1:nrow(res), function(y) {
      if (res[y, 4] < 0.1) {
        currdf <- data.frame(c(libs[x], y, as.integer(unlist(strsplit(res[y, 2],"/"))[1]), res[y,]))
        names(currdf) <- col_names
        fullSigRes <- rbind(fullSigRes, currdf)
      } else {
        return(NULL)
      }
    })
    sig <- sig[sapply(sig, function(x) !is.null(x))]
    if (length(sig) > 1) {
      sig <- rbindlist(sig)
      return(sig)
    } else {
      return(NULL)
    }
  })
  return(res)
}

###################################################################
###################################################################
### SECTION 07 - SCCLUST-VIZ FUNCTIONS (SCCV)
###################################################################
###################################################################

# SCCV - pull genes from seurat_sc across clusters
geneSearch <- function(txt,st,CGS) {
  if (length(txt) < 1) { txt <- ""}
  geneNames <- rownames(CGS)
  names(geneNames) <- toupper(CGS$genes)
  temp <- switch(st,
                 comma={
                   temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                   temp_out <- geneNames[toupper(temp_in)]
                   names(temp_out) <- CGS[temp_out,"genes"]
                   temp_out
                 },
                 regex={
                   temp_in <- grep(txt,names(geneNames),ignore.case=T)
                   temp_out <- geneNames[temp_in]
                   names(temp_out) <- CGS[temp_out,"genes"]
                   temp_out
                 })
  temp <- temp[!is.na(temp)]
  if (length(temp) > 0) {
    return(temp)
  } else {
    return(switch(st,
                  comma={
                    temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                    return(geneNames[which(toupper(geneNames) %in% toupper(temp_in))])
                  },
                  regex=grep(txt,geneNames,value=T,ignore.case=T)))
  }
}

# SCCV - set up data in an sCVdata format
sCVdata <- setClass(Class="sCVdata",
                    slots=c(Clusters="factor", #need to strip '-' from cluster names
                            ClustGeneStats="list", #list of length(getClusters(sCVdata)) containing data frames 
                            DEvsRest="list", #list of length(getClusters(sCVdata)) containing data frames
                            DEcombn="list", #list of length(getClusters(sCVdata)) containing lists of length(getClusters(sCVdata))-1 containing data frames 
                            Silhouette="ANY",
                            params="sCVparams"))

# IS THIS NEEDED? CAN'T FIND THIS IN SERVER SIDE
# SCCV - testing validity for scClustViz data
setValidity("sCVdata", function(object) {
  if (length(object@Clusters) > 0) {
    if (is.null(names(object@Clusters))) {
      return("Clusters should be a named factor where names are cells (colnames) of the input data.")
    }
    if (any(grepl("-",levels(object@Clusters)))) {
      return("Cluster names cannot contain '-'.")
    }
    if (any(table(object@Clusters) <= 1)) {
      return("All clusters must contain more than one cell.")
    }
  }
  if (length(object@ClustGeneStats) > 0) {
    if (!identical(names(object@ClustGeneStats),levels(object@Clusters)) |
        !all(sapply(object@ClustGeneStats,is.data.frame))) {
      return(paste("ClustGeneStats should be built using function 'calcClustGeneStats'.",
                   "Expected format is a named list with a data frame for each level in @Clusters."))
    }
  }
  if (length(object@DEvsRest) > 0) {
    if (!identical(names(object@DEvsRest),levels(object@Clusters)) |
        !all(sapply(object@DEvsRest,is.data.frame))) {
      return(paste("DEvsRest should be a named list with data frames for each level in @Clusters."))
    }
    if (!all(sapply(object@DEvsRest,
                    function(X) 
                      all(c("logGER","FDR") %in% names(X))))) {
      return("All DEvsRest data frames must contain variables 'logGER' and 'FDR'.")
    }
  }
  if (length(object@DEcombn) > 0) {
    if (length(object@DEcombn) != choose(length(levels(object@Clusters)),2) |
        !all(sapply(object@DEcombn,is.data.frame))) {
      return(paste("DEcombn should be a named list with data frames",
                   "for pairwise combinations of levels in @Clusters."))
    }
    tempNames <- apply(combn(levels(object@Clusters),2),2,
                       function(X) c(paste(X,collapse="-"),
                                     paste(X[2:1],collapse="-")))
    if (!all(sapply(seq_along(object@DEcombn),function(X) 
      names(object@DEcombn)[X] %in% tempNames[,X]))) {
      return(paste("Each entry name in DEcombn should be the names of the pair",
                   "of clusters (levels of @Clusters) separated by '-'.",
                   "See ?CalcDEcombn for example.",sep="\n"))
    }
    if (!all(sapply(object@DEcombn,
                    function(X) 
                      all(c("logGER","dDR","FDR") %in% names(X))))) {
      return("All DEcombn data frames must contain variables 'logGER', 'dDR', and 'FDR'.")
    }
  }
  if (!is.null(object@Silhouette)) {
    if (is(object@Silhouette)[1] != "silhouette") {
      return("Silhouette should be of the class 'silhouette' as returned by cluster::silhouette.")
    }
  }
  validObject(object@params)
})

# SCCV - calc gene detection rate/cluster; mean detected gene exp/cluster;
# mean gene exp/cluster using scClustViz
fx_calcCGS <- function(nge, cl, exponent, pseudocount) {
  message("-- Calculating gene detection rate per cluster --")
  DR <- pbapply::pbsapply(sapply(levels(cl),function(i) nge[,cl %in% i,drop=F],simplify=F),
                          function(X) apply(X,1,function(Y) sum(Y > 0)/length(Y)),simplify=F)
  message("-- Calculating mean detected gene expression per cluster --")
  MDGE <- pbapply::pbsapply(sapply(levels(cl),function(i) nge[,cl %in% i,drop=F],simplify=F),
                            function(X) apply(X,1,function(Y) {
                              temp <- meanLogX(Y[Y > 0],
                                               ncell=ncol(nge),
                                               ex=exponent,
                                               pc=pseudocount)
                              if (is.na(temp)) { temp <- 0 }
                              return(temp)
                            }),simplify=F)
  message("-- Calculating mean gene expression per cluster --")
  MGE <- pbapply::pbsapply(sapply(levels(cl),function(i) nge[,cl %in% i,drop=F],simplify=F),
                           function(X) apply(X,1,function(Y)
                             meanLogX(Y,
                                      ncell=ncol(nge),
                                      ex=exponent,
                                      pc=pseudocount)),simplify=F)
  return(sapply(levels(cl),function(X) 
    data.frame(DR=DR[[X]],MDGE=MDGE[[X]],MGE=MGE[[X]]),simplify=F))
}

# SCCV - calc diff gene exp cluster vs rest effect size scClustViz
fx_calcESvsRest <- function(nge,cl,CGS,exponent,pseudocount,DRthresh) {
  message("-- Calculating differential expression cluster vs rest effect size --")
  return(pbapply::pbsapply(levels(cl),function(i) {
    data.frame(
      logGER=CGS[[i]][CGS[[i]]$DR >= DRthresh,"MGE"] - 
        apply(nge[CGS[[i]]$DR >= DRthresh,(!cl %in% i | is.na(cl))],
              MARGIN=1,
              FUN=meanLogX,
              ncell=ncol(nge),
              ex=exponent,
              pc=pseudocount
        ),
      Wstat=NA,
      pVal=NA,
      FDR=NA,
      row.names=rownames(CGS[[i]])[CGS[[i]]$DR >= DRthresh]
    )
  },simplify=F))
}

# SCCV - dge cluster vs rest effect size scClustViz
fx_calcDEvsRest <- function(nge,cl,deTes) {
  message("-- Testing differential expression cluster vs rest --")
  deT_pVal <- presto::wilcoxauc(X=nge,y=cl)
  for (i in names(deTes)) {
    tempRows <- deT_pVal$feature %in% rownames(deTes[[i]]) & deT_pVal$group == i
    deTes[[i]][deT_pVal[tempRows,"feature"],"Wstat"] <- deT_pVal[tempRows,"statistic"]
    deTes[[i]][deT_pVal[tempRows,"feature"],"pVal"] <- deT_pVal[tempRows,"pval"]
    deTes[[i]][deT_pVal[tempRows,"feature"],"FDR"] <- p.adjust(deT_pVal[tempRows,"pval"],"fdr")
  } 
  return(deTes)
}

# SCCV - effect size for cluster A vs B scClustViz
fx_calcEScombn <- function(cl,CGS,DRthresh) {
  combos <- combn(levels(cl),2)
  colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
  return(apply(combos,2,function(i) {
    l <- CGS[[i[1]]]$DR >= DRthresh | CGS[[i[2]]]$DR >= DRthresh
    data.frame(
      logGER=CGS[[i[1]]]$MGE[l] - CGS[[i[2]]]$MGE[l],
      dDR=CGS[[i[1]]]$DR[l] - CGS[[i[2]]]$DR[l],
      row.names=rownames(CGS[[i[1]]])[l]
    )
  }))
}

# SCCV - dge cluster A vs B scClustViz
fx_calcDEcombn <- function(nge,cl,deMes) {
  combosL <- strsplit(names(deMes),"-")
  message("-- Testing differential expression between clusters --")
  deM_pVal <- pbapply::pbsapply(combosL,function(G) {
    temp <- presto::wilcoxauc(X=nge,y=cl,groups_use=G)
    temp <- temp[temp$group == G[1],]
    return(temp)
  },simplify=F)
  names(deM_pVal) <- names(deMes)
  for (i in names(deMes)) {
    tempRows <- deM_pVal[[i]]$feature %in% rownames(deMes[[i]])
    deMes[[i]][deM_pVal[[i]][tempRows,"feature"],"Wstat"] <- deM_pVal[[i]][tempRows,"statistic"]
    deMes[[i]][deM_pVal[[i]][tempRows,"feature"],"pVal"] <- deM_pVal[[i]][tempRows,"pval"]
    deMes[[i]][deM_pVal[[i]][tempRows,"feature"],"FDR"] <- p.adjust(deM_pVal[[i]][tempRows,"pval"],"fdr")
  } 
  return(deMes)
}

# SCCV - scClustViz fxn to plot character metadata
plot_mdBoxplotX <- function(MD,sel_clust,md_log,size) {
  temp_par <- par(no.readonly=T)
  
  par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
  #par(mar=c(3,3,2,1),mgp=2:0,cex=1.1)
  if (any(MD$sel_cells)) {
    temp1 <- tapply(MD[!MD$sel_cells,2],as.factor(MD[!MD$sel_cells,1]),c)
    temp2 <- tapply(MD[MD$sel_cells,2],as.factor(MD[MD$sel_cells,1]),c)
    plot(x=NULL,y=NULL,ylim=range(MD[,2]),
         xlim=c(0,length(levels(as.factor(MD[,1]))) * 3),
         log=md_log,xaxt="n",
         xlab=names(MD)[1],ylab=names(MD)[2])
    boxplot(temp1,add=T,xaxt="n",
            at=seq(1,length(levels(as.factor(MD[,1]))) * 3,by=3))
    boxplot(temp2,add=T,xaxt="n",border="red",
            at=seq(2,length(levels(as.factor(MD[,1]))) * 3,by=3))
    axis(side=1,at=seq(1.5,length(levels(as.factor(MD[,1]))) * 3,by=3),
         labels=names(temp1))
    legend(x=par("usr")[1],y=switch(sub("x","",md_log),"y"=10^par("usr")[4],par("usr")[4]),
           xjust=0,yjust=0.2,xpd=NA,bty="n",pch=0,col="red",pt.bg=alpha("red",0.5),
           legend=paste("Cluster",names(sel_clust),sel_clust))
    
  } else {
    boxplot(tapply(MD[,2],as.factor(MD[,1]),c),log=md_log,
            xlab=names(MD)[1],ylab=names(MD)[2])
  }
  par(temp_par)
}

# SCCV - scClustViz fxn to plot character metadata
plot_mdBoxplotY <- function(MD,sel_clust,md_log,size) {
  temp_par <- par(no.readonly=T)
  
  par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
  #par(mar=c(3,3,2,1),mgp=2:0,cex=1.1)
  if (any(MD$sel_cells)) {
    temp1 <- tapply(MD[!MD$sel_cells,1],as.factor(MD[!MD$sel_cells,2]),c)
    temp2 <- tapply(MD[MD$sel_cells,1],as.factor(MD[MD$sel_cells,2]),c)
    plot(x=NULL,y=NULL,xlim=range(MD[,1]),
         ylim=c(0,length(levels(as.factor(MD[,2]))) * 3),
         log=md_log,yaxt="n",
         xlab=names(MD)[1],ylab=names(MD)[2])
    boxplot(temp1,add=T,yaxt="n",horizontal=T,
            at=seq(1,length(levels(as.factor(MD[,2]))) * 3,by=3))
    boxplot(temp2,add=T,yaxt="n",border="red",horizontal=T,
            at=seq(2,length(levels(as.factor(MD[,2]))) * 3,by=3))
    axis(side=2,at=seq(1.5,length(levels(as.factor(MD[,2]))) * 3,by=3),
         labels=names(temp1))
    legend(x=switch(sub("y","",md_log),"x"=10^par("usr")[1],par("usr")[1]),y=par("usr")[4],
           xjust=0,yjust=0.2,xpd=NA,bty="n",pch=0,col="red",pt.bg=alpha("red",0.5),
           legend=paste("Cluster",names(sel_clust),sel_clust))
    
  } else {
    boxplot(tapply(MD[,1],as.factor(MD[,2]),c),log=md_log,
            horizontal=T,xlab=names(MD)[1],ylab=names(MD)[2])
  }
  par(temp_par)
}

# SCCV - scClustViz fxn to plot non-character metdata
plot_mdScatter <- function(MD,sel_clust,md_log,size) {
  if (nrow(MD) > 1e4) {
    temp_pch <- "."
    temp_cex <- 2
  } else {
    temp_pch <- 21
    temp_cex <- 1
  }
  temp_par <- par(no.readonly=T)
  layout(matrix(c(2,1,0,3),2),c(5,1),c(1,5))
  
  par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
  #par(mar=c(3,3,0,0),mgp=2:0,cex=1.1)
  plot(MD[!MD$sel_cells,1:2],log=md_log,xlim=range(MD[,1]),ylim=range(MD[,2]),
       pch=temp_pch,cex=temp_cex,col=alpha("black",0.2),bg=alpha("black",0.1))
  points(MD[MD$sel_cells,1:2],pch=temp_pch,cex=temp_cex + .5,
         col=alpha("red",0.4),bg=alpha("red",0.2))
  par(mar=c(0,3,1,0))
  boxplot(tapply(MD[,1],MD$sel_cells,c),log=sub("y","",md_log),
          horizontal=T,xaxt="n",yaxt="n",border=c("black","red"))
  if (any(MD$sel_cells)) {
    legend(x=switch(sub("y","",md_log),"x"=10^par("usr")[1],par("usr")[1]),y=par("usr")[4],
           xjust=0,yjust=0.2,xpd=NA,bty="n",pch=21,col="red",pt.bg=alpha("red",0.5),
           legend=paste("Cluster",names(sel_clust),sel_clust))
  }
  par(mar=c(3,0,0,1))
  boxplot(tapply(MD[,2],MD$sel_cells,c),log=sub("x","",md_log),
          horizontal=F,xaxt="n",yaxt="n",border=c("black","red"))
  par(temp_par)
}

# SCCV - scClustViz fxn to plot stacked barplots of character metdata
plot_mdBarplot <- function(MD,opt,size) {
  len <- length(levels(MD$cl))
  temp_par <- par(no.readonly=T)
  id0 <- as.factor(MD[,1])
  id <- switch(opt,
               "relative"=tapply(id0,MD$cl,function(X) table(X) / length(X)),
               "absolute"=tapply(id0,MD$cl,table))
  if (is.list(id)) { id <- do.call(cbind,id) }
  idylab <- switch(opt,
                   "relative"="Proportion of cells per cluster",
                   "absolute"="Number of cells per cluster")
  len2 <- length(levels(id0))
  par(mar=c(4,5,3,5) +.3, oma = c(0,0,0,5), family = "noto-sans-jp", cex = size)
  barplot(id, col = hue_pal()(len2),
          ylab = idylab, xlab = "Clusters",
          yaxt = "n", mgp = c(2,0,0))
  axis(2)
  abline(h=0)
  legend(x = par("usr")[2], y = par("usr")[4], 
         legend = levels(id0), cex = 1, xpd = NA,
         fill = hue_pal()(len2), bty = "n",
         inset = c(-.05, 0)
  )
}

# SCCV - scClustViz fxn to plot boxplots of non-character metdata
plot_mdBoxplot <- function(MD,opt,inp,size) {
  if (inp == "seurat_res") {
    len <- length(unique(MD$cl))
    temp_par <- par(no.readonly=T)
    
    par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
    boxplot(MD[,1] ~ MD$cl,log=opt,
            ylab=names(MD)[1],xlab="Clusters",
            border = hue_pal()(len),
            col = hue_pal()(len))
    par(temp_par)
  } else {
    len <- length(unique(MD$cl))
    temp_par <- par(no.readonly=T)
    
    par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
    boxplot(tapply(MD[,1],MD$cl,c),log=opt,
            ylab=names(MD)[1],xlab="Clusters",
            border = hue_pal()(len),
            col = hue_pal()(len))
    par(temp_par)
  }
}

# Adapted scClustViz functions for metadata relationship plot to ggplot2
plot_mdCompare_ggplot_scatter <- function(
  df, x, y, logX = F, logY = F, xRotate = F, sizeFactor = 1) {
  
  if(!is.numeric(df[[x]])) logX <- F
  if(!is.numeric(df[[y]])) logY <- F
  
  df <- data.frame(x = df[[x]], y = df[[y]])
  
  xLab <- x
  yLab <- y
  
  if(logX) {
    df$x <- log(df$x)
    xLab <- paste0(x, " (log scale)")
  }
  
  if(logY) {
    df$y <- log(df$y)
    yLab <- paste0(y, " (log scale)")
  }
  
  tickAngle <- 0
  vJust <- 1
  if(xRotate) {
    tickAngle <- 45
    vJust <- 0.7
  }
  
  p1 <- ggplot(df, aes(x = x, y = y)) + 
    geom_point(
      shape = 21, 
      fill = alpha("black", alpha = 0.1), 
      color = alpha("black", alpha = 0.2)
    ) +
    expand_limits(y = c(min(df$y) - 0.1 * diff(range(df$y)),
                        max(df$y) + 0.1 * diff(range(df$y)))) +
    expand_limits(x = c(min(df$x) - 0.1 * diff(range(df$x)),
                        max(df$x) + 0.1 * diff(range(df$x)))) +
    theme_nothing() +
    xlab(xLab) +
    ylab(yLab) +
    theme(
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.text.x = element_text(size = 15*sizeFactor, angle = tickAngle, vjust = vJust, margin = margin(t = 5)),
      axis.text.y = element_text(size = 15*sizeFactor, margin = margin(r = 5)),
      axis.title.x = element_text(size = 20*sizeFactor, margin = margin(t = 10)),
      axis.title.y = element_text(size = 20*sizeFactor, margin = margin(r = 10), angle = 90),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(x = 0.2, units = "cm"),
      text = element_text(family = "noto-sans-jp")
    ) 
  
  p2 <- ggplot(df, aes(x = factor(1), y = x)) + 
    stat_boxplot(geom = "errorbar") +
    geom_boxplot(outlier.shape = 1, color = "black", fill = "lightgray") + 
    expand_limits(y = c(min(df$x) - 0.1 * diff(range(df$x)),
                        max(df$x) + 0.1 * diff(range(df$x)))) +
    coord_flip() + 
    theme_nothing() +
    theme(
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  p3 <- ggplot(df, aes(x = factor(1), y = y)) + 
    stat_boxplot(geom = "errorbar") +
    geom_boxplot(outlier.shape = 1, color = "black", fill = "lightgray") + 
    expand_limits(y = c(min(df$y) - 0.1 * diff(range(df$y)),
                        max(df$y) + 0.1 * diff(range(df$y)))) +
    theme_nothing() +
    theme(
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  p2 + plot_spacer() + p1 + p3 +
    plot_layout(widths = c(6,1), heights = c(1,6), ncol = 2)
  
  
}

plot_mdCompare_ggplot_box <- function(df, x, y, logX = F, logY = F, xRotate = F,
                                      sizeFactor = 1, color = F) {
  
  if(!is.numeric(df[[x]])) logX <- F
  if(!is.numeric(df[[y]])) logY <- F
  
  df <- data.frame(x = df[[x]], y = df[[y]])
  
  xLab <- x
  yLab <- y
  
  if(logX) {
    df$x <- log(df$x)
    xLab <- paste0(x, " (log scale)")
  }
  
  if(logY) {
    df$y <- log(df$y)
    yLab <- paste0(y, " (log scale)")
  }
  
  tickAngle <- 0
  vJust <- 1
  if(xRotate) {
    tickAngle <- 45
    vJust <- 0.7
  }
  
  geom_boxplotValue <- geom_boxplot(outlier.shape = 1, color = "black", fill = "lightgray")
  stat_boxplotValue <- stat_boxplot(geom = "errorbar", color = "black")
  if(color) {
    geom_boxplotValue <- geom_boxplot(outlier.shape = 1)
    stat_boxplotValue <- stat_boxplot(geom = "errorbar")
  }
  
  ggplot(df, aes(x = x, y = y, fill = x, color = x)) + 
    stat_boxplotValue +
    geom_boxplotValue +
    theme_nothing() +
    xlab(xLab) +
    ylab(yLab) +
    theme(
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.text.x = element_text(size = 15*sizeFactor, angle = tickAngle, vjust = vJust, margin = margin(t = 5)),
      axis.text.y = element_text(size = 15*sizeFactor, margin = margin(r = 5)),
      axis.title.x = element_text(size = 20*sizeFactor, margin = margin(t = 10)),
      axis.title.y = element_text(size = 20*sizeFactor, margin = margin(r = 10), angle = 90),
      plot.margin = unit(x = c(1,1,1,1), units = "cm"),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(x = 0.2, units = "cm"),
      text = element_text(family = "noto-sans-jp")
    ) 
}

plot_mdCompare_ggplot_stackedbar <- function(
  df, x, y, xRotate = F, sizeFactor = 1, plotType = "absolute") {
  
  xLab <- x
  yLab <- y
  
  tickAngle <- 0
  vJust <- 1
  if(xRotate) {
    tickAngle <- 45
    vJust <- 0.7
  }
  
  df <- data.frame(x = df[[x]], y = df[[y]])
  
  geom_barValue <- geom_bar()
  yAxisLabel <- "Number of cells"
  if(plotType == "relative") {
    geom_barValue <- geom_bar(position = "fill")
    yAxisLabel <- "Proportion of cells per cluster"
  }
  
  ggplot(df, aes(x = x, fill = y)) +
    geom_barValue +
    scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
    theme_nothing() +
    xlab(xLab) +
    ylab(yAxisLabel) +
    guides(
      fill = guide_legend(
        title = yLab, 
        title.theme = element_text(size = 15*sizeFactor, family = "noto-sans-jp")
      )
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 15*sizeFactor),
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.text.x = element_text(
        size = 15*sizeFactor, angle = tickAngle, vjust = vJust, 
        margin = margin(t = 5)
      ),
      axis.text.y = element_text(size = 15*sizeFactor, margin = margin(r = 5)),
      axis.title.x = element_text(size = 20*sizeFactor, margin = margin(t = 10)),
      axis.title.y = element_text(size = 20*sizeFactor, margin = margin(r = 10), angle = 90),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(x = 0.2, units = "cm"),
      text = element_text(family = "noto-sans-jp"),
    ) 
  
}


plot_mdCompare_ggplot <- function(df, x, y, logX = F, logY = F, xRotate = F, sizeFactor = 1) {
  if(is.numeric(df[[x]]) && is.numeric(df[[y]])) {
    plot_mdCompare_ggplot_scatter(
      df = df, x = x, y = y, logX = logX, logY = logY, xRotate = xRotate,
      sizeFactor = sizeFactor
    )
  } else if(is.numeric(df[[x]]) || is.numeric(df[[y]])) {
    plot_mdCompare_ggplot_box(
      df = df, x = x, y = y, logX = logX, logY = logY, xRotate = xRotate,
      sizeFactor = sizeFactor
    )
  } else {
    plot_mdCompare_ggplot_stackedbar(
      df = df, x = x, y = y, xRotate = xRotate, sizeFactor = sizeFactor
    )
  }
}

# scClustViz gene expression by cluster plot using ggplot
# If supplied, hcData should be the result of running hclust.
# If supplied, drData should be a dataframe with columns "Cluster" and "DR"
# containing cluster numbers and detection rates for a gene
plot_GEboxplot_ggplot <- function(seuratData, scvData, gene, hcData = NULL, 
                                  plotJitter = T, drData = NULL, sizeFactor = 1) {
  
  if(is.null(hcData)) hcData <- hclust(d = stats::as.dist(m = DEdist(scvData)), method = "single")
  
  axisTitleSize <- 18 * sizeFactor
  axisTextSize <- 12 * sizeFactor
  
  myExp <- Param(scvData, param = "exponent")
  yLab <- paste0(gene, " normalized gene expression (log", myExp, " scale)")
  if(myExp == exp(1)) yLab <- paste0(gene, " normalized gene expression (natural log scale)")
  
  if(!(gene %in% rownames(seuratData))) return()
  expData <- FetchData(seuratData, vars = gene)
  myClusters <- Clusters(scvData)
  
  pos <- 1
  if(length(levels(myClusters)) > 1) pos <- levels(myClusters)[hcData$order]
  
  len <- length(levels(myClusters))

  levelNums <- as.numeric(1:len)
  names(levelNums) <- as.character(pos)
  group <- factor(myClusters[rownames(expData)], levels = pos)
  df <- data.frame(
    x = levelNums[as.character(myClusters)], 
    y = expData[[1]], 
    group = group
  )

  drValue <- NULL
  rightYAxisTitle <- rightYAxisText <- element_blank()
  if(!is.null(drData)) {
    drData <- data.frame(x = levelNums[drData$Cluster], y = drData$DR * max(df$y))
    drValue <- geom_point(
      mapping = aes(x = x, y = y), shape = "-", size = 20*sizeFactor, 
      data = drData
    )
    rightYAxisTitle <- element_text(size = axisTitleSize, margin = margin(l = 10), angle = 90)
    rightYAxisText <- element_text(size = axisTextSize, margin = margin(l = 5))
  }
  
  set.seed(10)
  
  p1 <- ggdendrogram(hcData) +
    coord_cartesian(xlim = c(0.5, len + 0.5)) +
    theme_nothing() +
    theme(
      text = element_text(family = "noto-sans-jp", angle = 90),
      plot.margin = unit(x = c(0,1,0,1), units = "cm")
    ) 
  
  p2 <- ggplot(df, mapping = aes(x = x, y = y)) +
    scale_x_continuous(
      name = "Clusters, ordered by heatmap dendrogram", 
      breaks = levelNums, 
      labels = names(levelNums)
    ) +
    coord_cartesian(xlim = c(0.5, len + 0.5))

  if(plotJitter) p2 <- p2 + geom_jitter(position = position_jitter(0.2), color = "black")
  
  p2 <- p2 +
    stat_boxplot(
      geom = "errorbar", mapping = aes(x = x, y = y, group = group),
      color = "black"
    ) +
    geom_boxplot(
      outlier.shape = NA, 
      mapping = aes(x = x, y = y, fill = group), 
      color = "black"
    ) +
    drValue
  
  if(!is.null(drData)) {
    p2 <- p2 +
      scale_y_continuous(
        sec.axis = sec_axis(
          trans = ~ . + 0,
          name = "- Gene detection rate per cluster",
          breaks = c(0, 0.25, 0.5, 0.75, 1) * max(df$y),
          labels = c("0%", "25%", "50%", "75%", "100%")
        )
      )
  }
  
  p2 <- p2 +
    theme_nothing() +
    xlab("Clusters, ordered by heatmap dendrogram") +
    ylab(yLab) +
    theme(
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.text.x = element_text(size = axisTextSize, margin = margin(t = 5)),
      axis.title.x = element_text(size = axisTitleSize, margin = margin(t = 10)),
      axis.text.y.left = element_text(size = axisTextSize, margin = margin(r = 5)),
      axis.title.y.left = element_text(size = axisTitleSize, margin = margin(r = 10), angle = 90),
      axis.text.y.right = rightYAxisText,
      axis.title.y.right = rightYAxisTitle,
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(x = 0.2, units = "cm"),
      text = element_text(family = "noto-sans-jp"),
      plot.margin = unit(x = c(0,1,1,1), units = "cm")
    ) 
  
  p1 + p2 + plot_layout(ncol = 1, heights = c(1,5))
}

# SCCV - Gene expression by cluster plot using ggplot
plot_goi_ggplot <- function(
  seuratData, gene, method = "umap", clusterLabels = T, plotType = "clust", 
  sizeFactor = 1, dims = c(1, 2)) {
  Idents(seuratData) <- "seurat_clusters"
  
  if(is.null(dims) || length(dims) != 2) {
    dims <- c(1, 2)
  }
  xLab <- paste0("PC ", dims[1])
  yLab <- paste0("PC ", dims[2])
  if(method == "umap") {
    xLab <- paste0("UMAP ", dims[1])
    yLab <- paste0("UMAP ", dims[2])
  }
  if(method == "tsne") {
    xLab <- paste0("tSNE ", dims[1])
    yLab <- paste0("tSNE ", dims[2])
  }
  
  if(plotType == "goi") {
    expData <- FetchData(seuratData, vars = gene)
    myMax <- max(expData)
    myMin <- min(expData)
    p <- FeaturePlot(
      object = seuratData, reduction = method, dims = dims, features = gene, 
      label = clusterLabels, label.size = 6*sizeFactor
    ) +
      scale_color_gradientn(
        colors = c("lightgrey", "blue"),
        labels = c(round(myMin, digits = 2), round(myMax, digits = 2)),
        breaks = c(myMin, myMax),
        guide = guide_colorbar(barwidth = 15)
      ) 
    myTitle <- gene
    hJust <- vJust <- NULL
    legendPos <- "top"
    legendJust <- "center"
    
  } else {
    p <- DimPlot(
      object = seuratData, reduction = method, dims = dims,
      group.by = "seurat_clusters", label = clusterLabels, 
      label.size = 6*sizeFactor
    ) 
    
    myTitle <- "Clusters"
    hJust <- 0
    vJust <- 2
    legendPos <- c(1,1)
    legendJust <- c(1,0)
    
  }
  
  p +
    ggtitle(myTitle) +
    theme_nothing() +
    ggtitle(myTitle) +
    xlab(xLab) +
    ylab(yLab) +
    theme(
      plot.title = element_text(size = 25*sizeFactor, face = "bold", hjust = hJust, vjust = vJust),
      legend.position = legendPos,
      legend.direction = "horizontal",
      legend.justification = legendJust,
      legend.text = element_text(size = 15*sizeFactor),
      panel.border = element_rect(fill = NA, size = 0.5),
      axis.text.x = element_text(size = 15*sizeFactor, margin = margin(t = 5)),
      axis.text.y = element_text(size = 15*sizeFactor, margin = margin(r = 5)),
      axis.title.x = element_text(size = 20*sizeFactor, margin = margin(t = 10)),
      axis.title.y = element_text(size = 20*sizeFactor, margin = margin(r = 10), angle = 90),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(x = 0.2, units = "cm"),
      plot.margin = margin(c(1.4,1.4,1.4,1.4), unit = "cm"),
      text = element_text(family = "noto-sans-jp")
    )    
  
}

# SCCV - scClustViz fxn modified to create cluster by gene
# overlay in gene expression tab
plot_tsne_adj <- function (cell_coord, md, md_title, md_cols = NULL, md_log = F, 
                           label = NULL, sel_cells, sel_cells_A, sel_cells_B,size) {
  if (is.null(md_title)) {
    id <- as.factor(md)
    if (is.null(md_cols)) {
      idcol <- colorspace::qualitative_hcl(length(levels(id)), 
                                           palette = "Dark 3")
    }
    else {
      idcol <- md_cols
    }
    if (any(is.na(id))) {
      levels(id) <- c(levels(id), "Unselected")
      id[is.na(id)] <- "Unselected"
      idcol <- c(idcol, "grey80")
    }
    par(mar=c(4, 4, 2.5, 1) +.3, family = "noto-sans-jp", 
        cex = size)
  }
  else if (is.factor(md) | is.character(md)) {
    id <- as.factor(md)
    idcol <- colorspace::qualitative_hcl(length(levels(id)), 
                                         palette = "Dark 3")
    par(mar=c(4, 4, ceiling(length(levels(id))/4) + 1, 
              1) +.3, family = "noto-sans-jp", 
        cex = size)
  }
  else if (any(md < 0)) {
    if (md_log) {
      warning("Can't log-scale md because it contains negative values.")
    }
    temp_down <- cut(c(0, md[md <= 0]), 50, labels = F)[-1]
    temp_up <- cut(c(0, md[md > 0]), 50, labels = F)[-1]
    id <- rep(NA, length(md))
    id[md <= 0] <- temp_down
    id[md > 0] <- temp_up + 50
    idcol <- colorspace::diverge_hcl(100, palette = "Blue-Red")
    #par(mar = c(3, 3, 2.5, 1), mgp = 2:0)
    par(mar=c(4, 4, 2.5, 1) +.3, family = "noto-sans-jp", 
        cex = size)
  }
  else {
    if (md_log) {
      id <- cut(log10(md), 100)
    }
    else {
      id <- cut(md, 100)
    }
    idcol <- colorspace::sequential_hcl(100, palette = "Viridis", 
                                        rev = T)
    par(mar=c(4, 4, 2.5, 1) +.3, family = "noto-sans-jp", 
        cex = size)
  }
  if (missing(sel_cells)) {
    sel_cells <- character()
  }
  if (nrow(cell_coord) > 10000) {
    temp_pch <- "."
    temp_cex <- 2
  }
  else {
    temp_pch <- 21
    temp_cex <- 1
  }
  plot(x = NULL, y = NULL, xlab = colnames(cell_coord)[1], 
       ylab = colnames(cell_coord)[2], xlim = range(cell_coord[,1]), 
       ylim = range(cell_coord[, 2]))
  if (length(sel_cells) > 0) {
    points(cell_coord[!rownames(cell_coord) %in% sel_cells, , drop = F], 
           pch = temp_pch, cex = temp_cex, col = alpha(idcol, 0.6)[id[!rownames(cell_coord) %in% sel_cells]], 
           bg = alpha(idcol, 0.3)[id[!rownames(cell_coord) %in% sel_cells]])
    points(cell_coord[sel_cells, , drop = F], pch = temp_pch, 
           cex = temp_cex + 0.5, 
           col = alpha(idcol, 1)[id[rownames(cell_coord) %in% sel_cells]], 
           bg = alpha(idcol, 0.6)[id[rownames(cell_coord) %in% sel_cells]])
  }
  else {
    points(cell_coord, pch = temp_pch, cex = temp_cex, 
           col = alpha(idcol, 0.8)[id], bg = alpha(idcol, 0.4)[id])
  }
  if (!missing(sel_cells_A) & !missing(sel_cells_B)) {
    points(x = cell_coord[sel_cells_A, 1], 
           y = cell_coord[sel_cells_A, 2], 
           pch = 19, col = "#a50026")
    points(x = cell_coord[sel_cells_B, 1], 
           y = cell_coord[sel_cells_B, 2], 
           pch = 19, col = "#313695")
    points(x = cell_coord[intersect(sel_cells_A, sel_cells_B), 1], 
           y = cell_coord[intersect(sel_cells_A, sel_cells_B), 2], 
           pch = 19, col = "#ffffbf")
    points(x = cell_coord[intersect(sel_cells_A, sel_cells_B), 1], 
           y = cell_coord[intersect(sel_cells_A, sel_cells_B), 2], 
           pch = 4, col = "red")
  }
  if (!is.null(label)) {
    text(label, labels = rownames(label), font = 2, cex = size)
  }
  if (is.null(md_title)) {
  }
  else if (is.factor(md) | is.character(md)) {
    legend(x = par("usr")[2], y = par("usr")[4], xjust = 1, 
           yjust = 0.2, xpd = NA, bty = "n", 
           ncol = switch(as.character(length(levels(id)) < 4), `TRUE` = length(levels(id)), `FALSE` = 4), 
           legend = levels(id), pch = 21, col = idcol, pt.bg = alpha(idcol, 0.5), cex = size)
    mtext(md_title, side = 3, adj = 0, font = 2, 
          line = ceiling(length(levels(id))/4) - 1, 
          cex = size)
  }
  else if (any(md < 0)) {
    temp_x <- c(seq(from = par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.15, 
                    to = par("usr")[1] + (par("usr")[2] - par("usr")[1])/2 - strwidth("0"), length.out = 51), 
                seq(from = par("usr")[2] - (par("usr")[2] - par("usr")[1])/2 + strwidth("0"), 
                    to = par("usr")[2] - (par("usr")[2] - par("usr")[1]) * 0.15, length.out = 51))
    for (i in 1:50) {
      rect(xleft = temp_x[i], xright = temp_x[i + 1], 
           ybottom = par("usr")[4] + (par("usr")[4] - par("usr")[3]) * 0.001, 
           ytop = par("usr")[4] + strheight(md_title), 
           col = idcol[i], border = NA, xpd = NA)
    }
    for (i in 52:102) {
      rect(xleft = temp_x[i], xright = temp_x[i + 1], 
           ybottom = par("usr")[4] + (par("usr")[4] - par("usr")[3]) * 
             0.001, ytop = par("usr")[4] + strheight(md_title), 
           col = idcol[i - 1], border = NA, xpd = NA)
    }
    mtext(round(min(md), 2), side = 3, line = 0, at = temp_x[1], 
          adj = 1.1, cex = size)
    mtext(round(max(md), 2), side = 3, line = 0, at = temp_x[102], 
          adj = -0.1, cex = size)
    mtext(0, side = 3, line = 0, adj = 0.5, 
          at = par("usr")[1] + (par("usr")[2] - par("usr")[1])/2)
    mtext(md_title, side = 3, line = 1, adj = 0.5, font = 2, 
          cex = size, at = par("usr")[1] + (par("usr")[2] - par("usr")[1])/2)
  }
  else {
    if (md_log) {
      md_title <- paste(md_title, "(log scale)")
    }
    temp_x <- seq(from = par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.15, 
                  to = par("usr")[2] - (par("usr")[2] - par("usr")[1]) * 0.15, length.out = 101)
    for (i in seq_along(idcol)) {
      rect(xleft = temp_x[i], xright = temp_x[i + 1], 
           ybottom = par("usr")[4] + (par("usr")[4] - par("usr")[3]) * 0.001, 
           ytop = par("usr")[4] + strheight(md_title), 
           col = idcol[i], border = NA, xpd = NA)
    }
    mtext(round(min(md), 2), side = 3, line = 0, at = temp_x[1], 
          adj = 1.1, cex = size)
    mtext(round(max(md), 2), side = 3, line = 0, at = temp_x[101], 
          adj = -0.1, cex = size)
    mtext(md_title, side = 3, line = 1, at = temp_x[51], 
          adj = 0.5, font = 2, cex = size)
  }
}

###################################################################
###################################################################
### SECTION 08 - MYSQL-DB FUNCTIONS (SQL)
###################################################################
###################################################################

# SQL - retrieve data from database
DBTableNRows <- function(db, myTab) {
  as.integer(DBI::dbGetQuery(db, statement = paste0("SELECT COUNT(*) FROM ", myTab))[1,1])
}

# SQL - down sample counts to a chosen number of cells
DownsampleCells <- function(db, expName, nCells = 2000, mySeed = NULL) {
  myCells <- DBI::dbReadTable(db, name = paste0(expName, "_meta"))$cell_id
  if(nCells <= length(myCells)) {
    if(!is.null(mySeed)) set.seed(mySeed)
    myCells <- sample(myCells, size = nCells)
  }
  myCells <- sub(pattern = "-", replacement = ".", x = myCells)
  cellString <- paste0("('", paste(myCells, collapse = "','"), "')")
  myStatement <- paste0("SELECT cell_id, gene, counts FROM ", expName, "_counts WHERE counts > 0 AND cell_id IN ", cellString)
  DBI::dbGetQuery(db, statement = myStatement)
}

# SQL - generate a seurat object from SQL db count/metadata (sc only)
SeuratProcess <- function(dataset, samples = NULL, nCells = NULL, mySeed = NULL) {
  # connect to mysql database
  mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                    dbname = scdb, host = ec_host, port = p)
  # read the table containing info on the scRNA-seq exp
  # in the mysql db
  dat <- dbReadTable(mydb, "sc")
  dat <- dat %>%
    dplyr::filter(experiment_name %in% dataset) %>%
    dplyr::pull(unique_table)
  
  metadata <- dbReadTable(mydb, name = paste(dat, "_meta", sep = ""))
  if(!is.null(nCells)) {
    if(nCells < nrow(metadata)) {
      counts <- DownsampleCells(db = mydb, expName = dat, nCells = nCells, mySeed = mySeed)
    } else {
      counts <- dbGetQuery(mydb, statement = paste0("SELECT cell_id, gene, counts FROM ", dat, "_counts WHERE counts > 0"))
    }
  } else {
    counts <- dbGetQuery(mydb, statement = paste0("SELECT cell_id, gene, counts FROM ", dat, "_counts WHERE counts > 0"))
  }
  dbDisconnect(mydb)
  
  if(!is.null(samples)) {
    rowsToKeep <- SubsetMulti(fullDF = metadata, subDF = samples)
    pull_cell_id <- metadata$cell_id[rowsToKeep]
    counts <- counts %>%
      filter(cell_id %in% pull_cell_id) %>%
      spread(cell_id, counts, fill = 0, convert = T) %>%
      column_to_rownames("gene")
    
    metadata <- metadata %>%
      filter(cell_id %in% pull_cell_id)
    metadata <- metadata %>%
      column_to_rownames("cell_id")
  } else {
    counts <- counts %>%
      spread(cell_id, counts, fill = 0, convert = T) %>%
      column_to_rownames("gene")
    
    metadata <- metadata %>%
      column_to_rownames("cell_id")
  }
  
  # Replace hyphens with periods
  rownames(metadata) <- sub(pattern = "-", replacement = ".", x = rownames(metadata))
  colnames(counts) <- sub(pattern = "-", replacement = ".", x = colnames(counts))
  
  seuratData <- CreateSeuratObject(
    counts = counts, 
    assay = "RNA", 
    meta.data = metadata[colnames(counts), , drop = F]
  )
  return(seuratData)
}

# Takes a full dataframe (e.g., metadata) and a second dataframe containing
# levels of each column to keep, and returns corresponding row numbers
SubsetMulti <- function(fullDF, subDF) {
  colNames <- colnames(subDF)
  if(length(colNames) == 1) {
    return(which(fullDF[[colNames]] %in% subDF[[colNames]]))
  }
  which(
    apply(as.matrix(fullDF[, colNames, drop = F]), 1, function(myVals) paste(myVals, collapse = "_AND_")) 
    %in%
    apply(as.matrix(subDF), 1, function(myVals) paste(myVals, collapse = "_AND_"))
  )
}

# Generate DESeqTransform object with original values (used for RASL-seq) from DESeqDataset
RawDDSTrans <- function(dds) {
  se <- SummarizedExperiment(
    assays = counts(dds, normalized = F), 
    colData = colData(dds), 
    rowRanges = rowRanges(dds), 
    metadata = metadata(dds)
  )
  DESeqTransform(se)
}

# Adapted PCA.DESeqTransform function from DESeq2 to handle data not in 
# DESeqDataset or DESeqTransform format (e.g., RASL-seq)
plotPCA_matrix <- function(mat, meta, intgroup = "condition", ntop = 500, returnData = FALSE) {
  mat <- na.omit(mat)
  meta <- meta[colnames(mat), ]
  # calculate the variance for each gene
  rv <- rowVars(mat)
  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(mat[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(intgroup %in% names(meta))) {
    stop("the argument 'intgroup' should specify columns of 'meta'")
  }
  
  intgroup.df <- as.data.frame(meta[, intgroup, drop = FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    meta[[intgroup]]
  }
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(mat))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() + theme(text = element_text(family = "noto-sans-jp"))
}

# Function to read in bulk dataset and store as DESeqDataset
RawToDDS <- function(exp) {
  mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                    dbname = bdb, host = ec_host, port = p)
  countsMat <- as.matrix(
    dbReadTable(mydb, name = paste(exp, "_counts", sep = "")) %>%
      spread(sample_id, counts, convert = T) %>%
      column_to_rownames("gene")
  )
  metaDF <- dbReadTable(mydb, name = paste(exp, "_meta", sep = "")) %>%
    column_to_rownames("sample_id")
  dbDisconnect(mydb)
  commonSamples <- intersect(colnames(countsMat), rownames(metaDF))
  countsMat <- countsMat[, commonSamples, drop = F]
  metaDF <- metaDF[commonSamples, , drop = F]
  metaDF <- metaDF[, !(apply(metaDF, 2, function(myCol) length(unique(myCol)) == 1)), drop = F]
  
  for(colname in colnames(metaDF)) {
    if(is.character(metaDF[, colname]) | is.factor(metaDF[, colname])) {
      metaDF[, colname] <- factor(gsub("-", replacement = "_", x = metaDF[, colname]))
    } else {
      metaDF[, colname] <- as.numeric(metaDF[, colname])
    }
  }
  filtGenes <- scan("data/filter_genes.txt", what = "character", quiet = T)
  countsMat <- countsMat[!(rownames(countsMat) %in% filtGenes), ]
  countsMat <- countsMat[rowSums(countsMat) >= 10, ]
  dds <- DESeqDataSetFromMatrix(
    countData = countsMat,
    colData = metaDF,
    design = ~ 1
  )
  return(dds)
}

# Reorders columns of a matrix or dataframe ('x'). Columns are grouped by 'fact'
# (a factor with names corresponding to column names of 'x'), and within each 
# group are sorted by colSums. Default is decreasing order. 
ReorderColsByFactor <- function(x, fact, decreasing = T) {
  myList <- lapply(levels(fact), function(myLevel) {
    colNames <- names(fact)[fact == myLevel]
    tmp <- x[, colNames, drop = F]
    return(colNames[order(colSums(tmp, na.rm = T), decreasing = decreasing)])
  })
  return(x[, unlist(myList), drop = F])
}

# Returns T if x is an integer or can be coerced to one, otherwise F
IsInteger <- function(x) {
  if(is.null(x) || !is.finite(x)) return(F)
  !suppressWarnings(is.na(as.integer(na.omit(x))))
}

# Returns T if x is a numeric or can be coerced to one, otherwise F
IsNumeric <- function(x) {
  if(is.null(x) || !is.finite(x)) return(F)
  !suppressWarnings(is.na(as.numeric(na.omit(x))))
}


# Takes a Seurat object that has slot "scale.data" and has had PCA run,
# and returns the number of PCs for which the summed fraction of total variance
# explained is as close as possible to 'varFrac' but no higher
PCADims <- function(seuratObject, varFrac = 1) {
  mat <- GetAssayData(object = seuratObject, slot = "scale.data")
  totalVariance <- sum(matrixStats::rowVars(mat))
  eigenValues <- (seuratObject[["pca"]]@stdev)^2
  varianceExplained <- eigenValues / totalVariance
  length(varianceExplained[cumsum(varianceExplained) < varFrac])
}

SeuratDGE <- function(seuratData, fact, group1, group2 = NULL,
                      logfc.threshold = 0.25, test.use = "wilcox",
                      min.pct = 0.1, min.diff.pct = -Inf, ...) {
  group1Levels <- group1
  group2_alias <- group2
  if(group1 == "All") {
    group1Levels <- sort(unique(as.character(seuratData@meta.data[[fact]])))
    group1Levels <- group1Levels[group1Levels != group2]
  }
  if(group2 == "Rest") group2_alias <- NULL
  
  seuratMajorVersion <- SeuratMajorVersion()
  if(test.use == "roc") {
    if(seuratMajorVersion <= 3) {
      colOrder <- c("Gene", "avg_logFC", "myAUC", "avg_diff", "power", "pct.1", "pct.2")
      colNames <- c("Gene", "Log fold-change", "AUC", "Avg diff", "Power",
                    paste0("Pct. ", group1), paste0("Pct. ", group2))
    } else {
      colOrder <- c("Gene", "avg_log2FC", "myAUC", "avg_diff", "power", "pct.1", "pct.2")
      colNames <- c("Gene", "Log2 fold-change", "AUC", "Avg diff", "Power",
                    paste0("Pct. ", group1), paste0("Pct. ", group2))
    }
  } else {
    if(seuratMajorVersion <= 3) {
      colOrder <- c("Gene", "avg_logFC", "p_val", "p_val_adj", "pct.1", "pct.2")
      colNames <- c("Gene", "Log fold-change", "P-value", "P-value (adj)",
                    paste0("Pct. ", group1), paste0("Pct. ", group2))
    } else {
      colOrder <- c("Gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
      colNames <- c("Gene", "Log2 fold-change", "P-value", "P-value (adj)",
                    paste0("Pct. ", group1), paste0("Pct. ", group2))
    }
  }
  
  myList <- lapply(group1Levels, function(myGroup) {
    tryCatch(
      suppressWarnings(
        FindMarkers(
          object = seuratData,
          ident.1 = myGroup,
          ident.2 = group2_alias,
          group.by = fact,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          verbose = F,
          ...
        )
      ),
      error = function(e) data.frame()
    )
  })
  names(myList) <- paste0(group1Levels, "-", group2)
  for(i in 1:length(myList)) {
    colNamesTmp <- colNames
    colNamesTmp[length(colNamesTmp) - 1] <- paste0("Pct. ", group1Levels[i])
    if(nrow(myList[[i]]) == 0) {
      myList[[i]] <- data.frame(matrix(nrow = 0, ncol = length(colNamesTmp)))
      colnames(myList[[i]]) <- colNamesTmp
    } else {
      myList[[i]]$Gene <- as.character(rownames(myList[[i]]))
      myList[[i]] <- myList[[i]][, colOrder]
      colnames(myList[[i]]) <- colNamesTmp
    }
  }
  return(myList)
}