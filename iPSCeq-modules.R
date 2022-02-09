#---------------------------------------------------------------------
# Title:         NCATS SEQUIN - MODULES
# Author:        Marissa Hirst
# Author2:       Ben Ernest
# Last Modified: 2021-09-14
# --
# Created:       2018-01-26 11:29:39 CDT
#---------------------------------------------------------------------


#' Layout for showing help buttons at top of each tab
#'
#' This function defines the layout for help info at the top of each tab
#'
#' @param id a callModule shiny output
#'   \code{\link{callModule}}.
#'
#' @examples
#' \dontrun{
#' helpButtonUI("overview_help")
#' }
#'
#' @export

## Help button module
helpButtonUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(1, uiOutput(ns("showhidebutton"))),
    column(9, uiOutput(ns("message")), offset = 1)
  )
}

#' Show help buttons at top of each tab
#'
#' This function toggles on/off help info for each tab
#'
#' @param id a callModule shiny output
#'   \code{\link{callModule}}.
#'
#' @examples
#' \dontrun{
#' helpButton
#' }
#'
#' @export

helpButton <- function(input, output, session, ...) {
  showhidebutton_label <- reactiveVal("Show help")
  
  output$showhidebutton <- renderUI({
    actionButton(session$ns("showhidebutton"),
                 label = showhidebutton_label())
  })
  
  observeEvent(input$showhidebutton, {
    if (showhidebutton_label() == "Show help") {
      showhidebutton_label("Hide help")
    } else {
      showhidebutton_label("Show help")
    }
  })
  
  output$message <- renderUI({
    if(showhidebutton_label() == "Show help") return(NULL)
    fluidPage(class = "myRow",
              ...,
              tags$head(tags$style(".myRow{background-color: white;}"))
    )
  })
  
}

# Cluster Solution DE boxplots -------------

#' scClustViz plot: Cluster separation boxplots
#'
#' This function plots metrics of cluster solution cohesion or overfitting as a
#' function of the number of clusters found.
#'
#' @param sCVdL A named list of sCVdata objects, output of
#'   \code{\link{CalcAllSCV}}.
#' @param DEtype One of "DEneighb", "DEmarker", or "silWidth". "DEneighb" shows
#'   number of significantly differentially expressed genes between nearest
#'   neighbouring clusters. "DEmarker" shows number of marker genes per cluster,
#'   significantly positively differentially expressed genes in all pairwise
#'   comparisons with other clusters. "silWidth" shows silhouette widths with
#'   average silhouette width as a trace across all clustering solutions. (see
#'   \code{\link[cluster]{silhouette}}).
#' @param FDRthresh Default=0.05. The false discovery rate threshold for
#'   determining significance of differential gene expression.
#' @param res Optional. Name of cluster resolution to highlight. Must be one of
#'   \code{names(sCVdL)}.
#' @param Xlim Optional. Passed to
#'   \code{\link[graphics]{plot.default}(xlim=Xlim)}.
#' @param Ylim Optional. Passed to
#'   \code{\link[graphics]{plot.default}(ylim=Ylim)}.
#'
#' @examples
#' \dontrun{
#' plot_clustSep(sCVdL,DEtype="DEneighb",FDRthresh=0.05,res="res.0.8")
#' }
#'
#' @export

options(getClass.msg = F)

plot_clustSep <- function(sCVdL,DEtype,FDRthresh=0.05,res,Xlim,Ylim,size) {
  if (missing(Xlim)) { Xlim <- NULL }
  if (missing(Ylim)) { Ylim <- NULL }
  if (missing(res)) { res <- "" }
  if (!res %in% c(names(sCVdL),"")) {
    warning(paste(paste0("res = '",res,"' not found in cluster resolutions."),
                  "Cluster resolutions are names(sCVdL):",
                  paste(names(sCVdL),collapse=", "),sep="\n  "))
  }
  if (!DEtype %in% c("DEneighb","DEmarker","silWidth")) {
    stop('DEtype must be one of "DEneighb", "DEmarker", or "silWidth".')
  }
  numClust <- sapply(sCVdL,function(X) length(levels(Clusters(X))))
  for (X in unique(numClust[duplicated(numClust)])) {
    numClust[numClust == X] <- seq(X-.25,X+.25,length.out=sum(numClust == X))
  }
  
  if (is.null(Xlim)) { Xlim <- range(numClust) }
  bpData <- sapply(sCVdL,function(X) switch(DEtype,
                                            DEneighb=sapply(DEneighb(X,FDRthresh),nrow),
                                            DEmarker=sapply(DEmarker(X,FDRthresh),nrow),
                                            silWidth=Silhouette(X)[,"sil_width"]),
                   simplify=F)
  if (is.null(Ylim)) { Ylim <- range(unlist(bpData)) }
  
  if (grepl("^Comp",res)) {
    par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA,
         family = "noto-sans-jp", 
         cex = size)
    text(.5,.5,paste("Press 'View clusters at this resolution'",
                     "to view the comparison",
                     sub("Comp.","",res,fixed=T),sep="\n"))
  } else {
    par(mar=c(4,4,3,1) +.3, family = "noto-sans-jp", cex = size)
    if (DEtype == "silWidth") {
      plot(x=NA,y=NA,xlim=Xlim + c(-.5,.5),ylim=Ylim,xaxt="n",
           xlab="Number of clusters",ylab="Silhouette width per cluster")
    } else {
      plot(x=numClust,y=sapply(bpData,median),type="l",xaxt="n",
           xlim=Xlim + c(-.5,.5),ylim=Ylim,xlab="Number of clusters",
           ylab=switch(DEtype,
                       scoreDE="Distance between clusters by differential expression score",
                       DEmarker="Positive DE genes per cluster to all other clusters",
                       DEneighb="Positive DE genes per cluster to nearest cluster"))
    }
    axis(side=3,at=seq(round(min(numClust)) - 0.5,round(max(numClust)) + 0.5,by=1),
         labels=F,tick=T,pos=par("usr")[3])
    axis(side=1,at=seq(round(min(numClust)) - 0.5,round(max(numClust)) + 0.5,by=1),
         labels=F,tick=T,pos=par("usr")[3])
    axis(side=1,at=seq(round(min(numClust)),round(max(numClust)),by=1),labels=T,tick=F)
    
    abline(h=seq(0,max(unlist(bpData)),switch(as.character(diff(Ylim) > 1000),
                                              "FALSE"=10,"TRUE"=100)),
           lty=3,col=alpha(1,0.3))
    for (i in names(bpData)[names(bpData) != res]) {
      boxplot(bpData[[i]],add=T,at=numClust[i],yaxt="n",col=alpha("white",.5))
    }
    if (any(names(bpData) == res)) {
      if (DEtype == "silWidth") {
        boxplot(bpData[[res]],add=T,at=numClust[res],border="red")
      } else {
        boxplot(bpData[[res]],add=T,at=numClust[res],border="red",outline=F)
        points(jitter(rep(numClust[res],length(bpData[[res]])),amount=.2),
               bpData[[res]],col=alpha("red",.5),pch=20)
      }
    }
    if (DEtype == "silWidth") {
      temp_avSil <- sapply(bpData,mean)
      lines(numClust,y=temp_avSil,type="b",col="darkred",pch=16)
      points(numClust[res],temp_avSil[res],col="red",pch=16)
      legend(x=par("usr")[2],y=par("usr")[4],
             xjust=1,yjust=0.2,xpd=NA,bty="n",horiz=T,
             legend=c("Average silhouette width",
                      paste("Selected resolution:",res)),
             col=c("darkred","red"),pch=c(16,0),lty=c(1,NA))
    } else {
      legend(x=par("usr")[2],y=par("usr")[4],
             xjust=1,yjust=0.2,xpd=NA,bty="n",
             legend=paste("Selected resolution:",res),
             col="red",pch=0)
    }
  }
}

# Silhouette plot ------

#' scClustViz plot: Silhouette plot
#'
#' This function is a wrapper to \code{plot(silhouette(x))}.
#'
#' @param sCVd An \code{\link{sCVdata}} object with a non-null \code{Silhouette}
#'   slot.
#'
#' @export

plot_sil <- function(sCVd, size) {
  len <- length(levels(Clusters(sCVd)))
  n <- nrow(Silhouette(sCVd))
  x <- sortSilhouette(Silhouette(sCVd))
  s <- rev(x[, "sil_width"])
  space <- c(0, rev(diff(cli <- x[, "cluster"])))
  space[space != 0] <- 0.5
  axisnames <- (n < 40)
  main <- NULL
  
  if (is.null(main)) {
    main <- "Silhouette plot"
    if (!is.null(cll <- attr(x, "call"))) {
      if (!is.na(charmatch("silhouette", deparse(cll[[1]])))) 
        cll[[1]] <- as.name("FF")
      main <- paste(main, "of", sub("^FF", "", deparse(cll)))
    }
  }
  smry <- summary(x)
  k <- length(nj <- smry$clus.sizes)
  sub <- paste("Average silhouette width : ", round(smry$avg.width, 
                                                    digits = 2))
  
  par(mar=c(5.5,2,3.5,1.5) +.3, family = "noto-sans-jp",
      cex = size)
  plot(Silhouette(sCVd),
       nmax.lab = 0,
       beside = T,
       border = NA,
       main = NA,
       col = hue_pal()(len),
       do.n.k = F,
       do.clus.stat = F) +
    mtext(paste("n =", n), adj = 0, cex = 1.5) +
    mtext(substitute(k ~ ~"clusters" ~ ~C[j], list(k = k)), 
        adj = 1, cex = 1.5)
}

# tsnePlot -------------------

#' scClustViz plot element: Cluster names on cluster centroid.
#'
#' See \code{\link{plot_tsne}} for application.
#'
#' @param sCVd An sCVdata object.
#' @param cell_coord A numeric matrix where named rows are cells, and two
#'   columns are the x and y dimensions of the cell embedding.
#' @param lab_type One of "ClusterNames", "ClusterNamesAll", or "Clusters".
#'   "ClusterNames" places cluster names (added to sCVdata object by
#'   \code{\link{labelCellTypes}}) at the centroid of all points sharing that
#'   cluster name (can span clusters). "ClusterNamesAll" places cluster names at
#'   the centroid of each cluster. "Clusters" places cluster ID
#'   (\code{levels(Clusters(sCVd))}) at the centroid of each cluster.
#'   
#' @export

tsne_labels <- function(sCVd,cell_coord,lab_type) {
  if (!lab_type %in% c("ClusterNames","ClusterNamesAll","Clusters")) {
    stop('lab_type must be one of "ClusterNames","ClusterNamesAll","Clusters"')
  }
  if (lab_type == "ClusterNames") {
    temp_labelNames <- sapply(unique(attr(Clusters(sCVd),"ClusterNames")),function(X) 
      names(which(attr(Clusters(sCVd),"ClusterNames") == X)),simplify=F)
    temp_labels <- apply(cell_coord,2,function(Y) 
      tapply(Y,apply(sapply(temp_labelNames,function(X) Clusters(sCVd) %in% X),1,which),mean))
    if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
    rownames(temp_labels) <- names(temp_labelNames)
  } else if (lab_type == "ClusterNamesAll") {
    temp_labels <- apply(cell_coord,2,function(X) tapply(X,Clusters(sCVd),mean))
    if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
    rownames(temp_labels) <- attr(Clusters(sCVd),"ClusterNames")
  } else if (lab_type == "Clusters") {
    temp_labels <- apply(cell_coord,2,function(X) tapply(X,Clusters(sCVd),mean))
    if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
    rownames(temp_labels) <- levels(Clusters(sCVd))
  } else {
    stop("lab_type should be one of 'ClusterNames', 'ClusterNamesAll', or 'Clusters'.")
  }
  return(temp_labels)
}

#' scClustViz plot: Plot cell embedding in 2D
#'
#' This function plots cells in two dimensions, with various overlays.
#'
#' @param cell_coord A numeric matrix where named rows are cells, and two
#'   columns are the x and y dimensions of the cell embedding.
#' @param md The overlay information. Either a factor or numeric vector matching
#'   the rows (cells) of the \code{cell_coord} matrix. If this is a factor, the
#'   cells will be coloured by the factor levels. If a positive numeric vector,
#'   the cells will be coloured using the Viridis sequential colourscale
#'   implemented in \code{\link[colorspace]{sequential_hcl}}. Otherwise a
#'   diverging red-blue colourscale from \code{\link[colorspace]{diverging_hcl}}
#'   will be used.
#' @param md_title NULL or a character vector of one. If NULL, \code{md} is
#'   assumed to be cluster assignments. Otherwise this should be the title of
#'   the overlay represented by \code{md}.
#' @param md_log Default=FALSE. Logical vector of length one indicating whether
#'   \code{md} should be log-transformed. Only to be used when \code{md} is
#'   numeric.
#' @param label Default=NULL. The output of \code{\link{tsne_labels}} to have
#'   cluster names overlaid on the plot.
#' @param sel_cells Optional. A character vector of cell names (rownames of
#'   \code{cell_coord}) to highlight in the plot.
#' @param sel_cells_A Optional. Alternative highlighting method to sel_cells,
#'   can be used in conjunction. Meant for indicating a selected set of cells
#'   when building manual cell set comparisons, in conjunction with
#'   \code{sel_cells_B}.
#' @param sel_cells_B Optional. See \code{sel_cells_A}.
#'
#' @examples
#' \dontrun{
#' # Cluster overlay:
#' plot_tsne(cell_coord=getEmb(input_data_obj,"tsne"),
#'           md=Clusters(sCVdata),
#'           md_title=NULL,
#'           label=tsne_labels(sCVd=sCVdata,
#'                             cell_coord=getEmb(input_data_obj,"tsne"),
#'                             lab_type="ClusterNames"))
#'
#' # Metadata overlay:
#' plot_tsne(cell_coord=getEmb(input_data_obj,"tsne"),
#'           md=getMD(input_data_obj)$total_counts,
#'           md_title="Library Size",
#'           md_log=TRUE,
#'           label=tsne_labels(sCVd=sCVdata,
#'                             cell_coord=getEmb(input_data_obj,"tsne"),
#'                             lab_type="ClusterNames"))
#'
#' # Gene expression overlay:
#' plot_tsne(cell_coord=getEmb(input_data_obj,"tsne"),
#'           md=getExpr(input_data_obj,Param(sCVdata,"assayType"))["Actb",],
#'           md_title="Actb")
#' }
#'
#' @export
#' 

# SCCV - updated scClustViz fxn to allow choice of color palette
# for gene expression overlay
plot_tsne2 <- function(cell_coord,md,md_title,md_log=F,label=NULL,
                       sel_cells,sel_cells_A,sel_cells_B, idcol = NULL,size) {
  if (is.null(md_title)) {
    id <- as.factor(md)
    idcol <- colorspace::qualitative_hcl(length(levels(id)),
                                         palette="Dark 3")
    if (any(is.na(id))) {
      levels(id) <- c(levels(id),"Unselected")
      id[is.na(id)] <- "Unselected"
      idcol <- c(idcol,"grey80")
    }
    
    par(mar=c(4,4,2.5,1) +.3, family = "noto-sans-jp", 
        cex = size)
  } else if (is.factor(md) | is.character(md)) {
    id <- as.factor(md)
    idcol <- colorspace::qualitative_hcl(length(levels(id)),palette="Dark 3")
    
    par(mar=c(4,4,ceiling(length(levels(id))/4)+1,1) +.3, family = "noto-sans-jp", 
        cex = size)
  } else if (any(md < 0)) {
    if (md_log) {
      warning("Can't log-scale md because it contains negative values.")
    }
    temp_down <- cut(c(0,md[md <= 0]),50,labels=F)[-1]
    temp_up <- cut(c(0,md[md > 0]),50,labels=F)[-1]
    id <- rep(NA,length(md))
    id[md <= 0] <- temp_down
    id[md > 0] <- temp_up + 50
    if(is.null(idcol)) {
      idcol <- colorspace::diverge_hcl(100,palette="Blue-Red")
    }
    
    par(mar=c(4,4,2.5,1) +.3, family = "noto-sans-jp", 
        cex = size)
  } else{
    if (md_log) {
      id <- cut(log10(md),100)
    } else {
      id <- cut(md,100)
    }
    if(is.null(idcol)) {
      idcol <- colorspace::sequential_hcl(100,palette="Viridis",rev=T)
    }
    
    par(mar=c(4,4,2.5,1) +.3, family = "noto-sans-jp", 
        cex = size)
  }
  if (missing(sel_cells)) { sel_cells <- character() }
  if (nrow(cell_coord) > 1e4) {
    temp_pch <- "."
    temp_cex <- 2
  } else {
    temp_pch <- 21
    temp_cex <- 1
  }
  
  plot(x=NULL,y=NULL,xlab=gsub("_", " ", colnames(cell_coord)[1]),ylab=gsub("_", " ", colnames(cell_coord)[2]),
       xlim=range(cell_coord[,1]),ylim=range(cell_coord[,2]))
  if (length(sel_cells) > 0) {
    points(cell_coord[!rownames(cell_coord) %in% sel_cells,],pch=temp_pch,cex=temp_cex,
           col=alpha(idcol,.6)[id[!rownames(cell_coord) %in% sel_cells]],
           bg=alpha(idcol,0.3)[id[!rownames(cell_coord) %in% sel_cells]])
    points(cell_coord[sel_cells,],pch=temp_pch,cex=temp_cex + .5,
           col=alpha(idcol,1)[id[rownames(cell_coord) %in% sel_cells]],
           bg=alpha(idcol,0.6)[id[rownames(cell_coord) %in% sel_cells]])
  } else {
    points(cell_coord,pch=temp_pch,cex=temp_cex,col=alpha(idcol,.8)[id],bg=alpha(idcol,0.4)[id])
  }
  
  if (!missing(sel_cells_A) & !missing(sel_cells_B)) {
    points(x=cell_coord[sel_cells_A,1],
           y=cell_coord[sel_cells_A,2],
           pch=19,col="#a50026")
    points(x=cell_coord[sel_cells_B,1],
           y=cell_coord[sel_cells_B,2],
           pch=19,col="#313695")
    points(x=cell_coord[intersect(sel_cells_A,sel_cells_B),1],
           y=cell_coord[intersect(sel_cells_A,sel_cells_B),2],
           pch=19,col="#ffffbf")
    points(x=cell_coord[intersect(sel_cells_A,sel_cells_B),1],
           y=cell_coord[intersect(sel_cells_A,sel_cells_B),2],
           pch=4,col="red")
  }
  if (!is.null(label)) {
    text(label,labels=rownames(label),font=2,cex=size)
  }
  if (is.null(md_title)) {
  } else if (is.factor(md) | is.character(md)) {
    legend(x=par("usr")[2],y=par("usr")[4],
           xjust=1,yjust=0.2,xpd=NA,bty="n",
           ncol=switch(as.character(length(levels(id)) < 4),"TRUE"=length(levels(id)),"FALSE"=4),
           legend=levels(id),pch=21,col=idcol,pt.bg=alpha(idcol,0.5))
    mtext(md_title,side=3,adj=0,font=2,line=ceiling(length(levels(id))/4)-1,cex=size)
  } else if (any(md < 0)) {
    temp_x <- c(
      seq(from=par("usr")[1] + (par("usr")[2] - par("usr")[1]) * .15,
          to=par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2 - strwidth("0"),
          length.out=51),
      seq(from=par("usr")[2] - (par("usr")[2] - par("usr")[1]) / 2 + strwidth("0"),
          to=par("usr")[2] - (par("usr")[2] - par("usr")[1]) * .15,
          length.out=51)
    )
    for (i in 1:50) {
      rect(xleft=temp_x[i],xright=temp_x[i+1],
           ybottom=par("usr")[4] + (par("usr")[4] - par("usr")[3]) * .001,
           ytop=par("usr")[4] + strheight(md_title),
           col=idcol[i],border=NA,xpd=NA)
    }
    for (i in 52:102) {
      rect(xleft=temp_x[i],xright=temp_x[i+1],
           ybottom=par("usr")[4] + (par("usr")[4] - par("usr")[3]) * .001,
           ytop=par("usr")[4] + strheight(md_title),
           col=idcol[i-1],border=NA,xpd=NA)
    }
    mtext(round(min(md),2),side=3,line=0,at=temp_x[1],adj=1.1, cex = size)
    mtext(round(max(md),2),side=3,line=0,at=temp_x[102],adj=-0.1, cex = size)
    mtext(0,side=3,line=0,adj=.5,
          at=par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2)
    mtext(md_title,side=3,line=1,adj=.5,font=2, cex = size,
          at=par("usr")[1] + (par("usr")[2] - par("usr")[1]) / 2)
  } else {
    if (md_log) { md_title <- paste(md_title,"(log scale)") } 
    temp_x <- seq(from=par("usr")[1] + (par("usr")[2] - par("usr")[1]) * .15,
                  to=par("usr")[2] - (par("usr")[2] - par("usr")[1]) * .15,
                  length.out=101)
    for (i in seq_along(idcol)) {
      rect(xleft=temp_x[i],xright=temp_x[i+1],
           ybottom=par("usr")[4] + (par("usr")[4] - par("usr")[3]) * .001,
           ytop=par("usr")[4] + strheight(md_title),
           col=idcol[i],border=NA,xpd=NA)
    }
    mtext(round(min(md),2),side=3,line=0,at=temp_x[1],adj=1.1, cex = size)
    mtext(round(max(md),2),side=3,line=0,at=temp_x[101],adj=-0.1, cex = size)
    mtext(md_title,side=3,line=1,at=temp_x[51],adj=.5,font=2,cex=size)
  }
}

#' scClustViz plot: Plot an MA plot
#'
#' This function plots an MA plot for two different cluster comparisons.
#'
#' @param sCVd scClustViz object
#' @param clA cluster A for the first comparison
#' @param clB cluster B for the second comparison
#' @param dataType For MA-style plots comparing difference and mean of gene
#'   summary statistics, one of: \code{"DR"} (detection rate); \code{"MGE"}
#'   (mean gene expression); \code{"MDGE"} (mean detected gene expression). For
#'   volcano plots, the effect size measure can be one of: \code{"dDR"}
#'   (difference in detection rate); \code{"logGER"} (log gene expression
#'   ratio). To compare relationship between difference in detection rate and
#'   log gene expression ratio, use \code{"GERvDDR"}.
#' @param labType Default="de". A character vector indicating which genes to
#'   highlight. One of \code{"de"} (most statistically significant genes),
#'   \code{"diff"} (most different by dataType shown), or \code{"search"}
#'   (specified genes).
#' @param labGenes Only required if \code{labType="search"}. Gene names to
#'   highlight.
#' @param labNum Default=5. Number of genes to highlight per side.
#'
#' @examples
#' \dontrun{
#' # Cluster overlay:
#' plot_compareClusts_MAplot2(sCVd = seurat_sc[[1]], clA = 1, clB = 2,
#'                            dataType = "DR", labType = "de",
#'                            labNum = 5)
#'
#' }
#'
#' @export
#' 

plot_compareClusts_MAplot2 <- function(sCVd, clA, clB, labType, labNum, labGenes,
                                       sizeFactor = 1) {
  
  CGS <- compareClusts_DF(sCVd = sCVd, clA = clA, clB = clB, dataType = "MGE")
  temp_label <- paste0("mean normalized gene expression (log", Param(sCVd,"exponent")," scale)")
  if(Param(sCVd,"exponent") == exp(1)) {
    temp_label <- paste0("mean normalized gene expression (natural log scale)")
  }
  
  gnA <- rownames(head(CGS[order(CGS$x_diff, decreasing = T), ], labNum))
  gnB <- rownames(tail(CGS[order(CGS$x_diff, decreasing = T), ], labNum))
  if(labType == "de") {
    ts <- order(CGS$FDR, na.last=T)
    gnA <- rownames(CGS)[ts[CGS[ts, "dir"] == clA][1:labNum]]
    gnB <- rownames(CGS)[ts[CGS[ts, "dir"] == clB][1:labNum]]
  }
  
  dat <- data.frame(x = CGS$x_diff, y = CGS$y_mean, gene = rownames(CGS))
  rownames(dat) <- dat$gene
  dat$col <- alpha("black", 0.3)
  clustCols <- scales::hue_pal()(2)
  dat[gnA, "col"] <- alpha(clustCols,alpha = 0.8)[1]
  dat[gnB, "col"] <- alpha(clustCols, alpha = 0.8)[2]
  col <- as.character(dat$col)
  names(col) <- as.character(dat$col)
  datLab <- dat[c(gnA, gnB), ]
  datLab[gnA, "col"] <- alpha(clustCols, alpha = 0.8)[1]
  datLab[gnB, "col"] <- alpha(clustCols, alpha = 0.8)[2]
  
  xLab <- paste0("Difference in ", temp_label)
  yLab <- paste0("Avg ", temp_label)
  mainLab <- paste0("Modified MA plot of mean gene expression ", clA, " vs. ", clB)
  subLab <- paste(
    paste("Cosine similarity:",
          round(scClustViz:::cosineSim(ClustGeneStats(sCVd)[[clA]][, "MGE"], 
                                       ClustGeneStats(sCVd)[[clB]][, "MGE"]), 2)),
    paste("Spearman's Rho:",
          round(cor(x = ClustGeneStats(sCVd)[[clA]][, "MGE"], 
                    y = ClustGeneStats(sCVd)[[clB]][, "MGE"],
                    method = "spearman"), 2)),
    sep = ", ")
  
  theme_set(theme_classic(base_family = "noto-sans-jp", 
                          base_size = 12))
  # Change the settings
  update_geom_defaults("text", list(family = theme_get()$text$family))
  
  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(aes(color = col)) +
    scale_color_manual(values = col) +
    annotate("text", label = paste0("Higher in ", clA), size = 8 * sizeFactor,
             color = clustCols[1],
             x = Inf, y = -Inf, hjust = 1.03, vjust = -0.5,
             family= theme_get()$text[["family"]]) +
    annotate("text", label = paste0("Higher in ", clB), size = 8 * sizeFactor,
             color = clustCols[2],
             x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
             family= theme_get()$text[["family"]]) +
    labs(title = mainLab, subtitle = subLab) +
    xlab(xLab) +
    ylab(yLab) +
    geom_text_repel(data = datLab, mapping = aes(x, y, label = gene, color = col),
                    size = 5 * sizeFactor, fontface = "bold", family = "noto-sans-jp",
                    max.overlaps = 30
    ) +
    geom_vline(xintercept = 0, color = "gray50") +
    theme(
      text = element_text(family = "noto-sans-jp"),
      axis.text = element_text(size = 15 * sizeFactor),
      axis.title = element_text(size = 20 * sizeFactor), 
      plot.subtitle = element_text(size = 15 * sizeFactor),
      plot.title = element_text(size = 25 * sizeFactor),
      legend.position = "none"
    )
  print(p)
} 


#' scClustViz plot: Plot a scatterplot
#'
#' This function plots the relationship between difference in detection rate and
#' log gene expression ratio
#'
#' @param sCVd scClustViz object
#' @param clA cluster A for the first comparison
#' @param clB cluster B for the second comparison
#' @param dataType For MA-style plots comparing difference and mean of gene
#'   summary statistics, one of: \code{"DR"} (detection rate); \code{"MGE"}
#'   (mean gene expression); \code{"MDGE"} (mean detected gene expression). For
#'   volcano plots, the effect size measure can be one of: \code{"dDR"}
#'   (difference in detection rate); \code{"logGER"} (log gene expression
#'   ratio). To compare relationship between difference in detection rate and
#'   log gene expression ratio, use \code{"GERvDDR"}.
#' @param labType Default="de". A character vector indicating which genes to
#'   highlight. One of \code{"de"} (most statistically significant genes),
#'   \code{"diff"} (most different by dataType shown), or \code{"search"}
#'   (specified genes).
#' @param labGenes Only required if \code{labType="search"}. Gene names to
#'   highlight.
#' @param labNum Default=5. Number of genes to highlight per side.
#' @param labTypeDiff Default="logGER". Only required if
#'   \code{dataType="GERvDDR"} and \code{labType="diff"}. Which axis to use for
#'   difference calculation. One of \code{"dDR"} (difference in detection rate)
#'   or \code{"logGER"} (log gene expression ratio).

#'
#' @examples
#' \dontrun{
#' # Cluster overlay:
#' plot_compareClusts_DEscatter2(sCVd = seurat_sc[[1]], clA = 1, clB = 2,
#'                            dataType = "DR", labType = "diff",
#'                            labNum = 5, labTypeDiff = "logGER")
#'
#' }
#'
#' @export
#' 

plot_compareClusts_DEscatter2 <- function(sCVd,clA,clB,dataType,labType,
                                          labTypeDiff,labNum,labGenes) {
  CGS <- compareClusts_DF(sCVd,clA,clB,dataType)
  temp_exp <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                     "TRUE"="(natural log scale)",
                     "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
  if (labType == "diff") {
    gnA <- rownames(head(CGS[order(CGS[[labTypeDiff]],decreasing=T),],labNum))
    gnB <- rownames(tail(CGS[order(CGS[[labTypeDiff]],decreasing=T),],labNum))
  } else if (labType == "de") {
    ts <- order(CGS$FDR,na.last=T)
    gnA <- rownames(CGS)[ts[CGS[ts,"dir"] == clA][1:labNum]]
    gnB <- rownames(CGS)[ts[CGS[ts,"dir"] == clB][1:labNum]]
  }
  clustCols <- scales::hue_pal()(2)
  
  dat <- data.frame(x = CGS$dDR, y = CGS$logGER, gene = rownames(CGS))
  rownames(dat) <- dat$gene
  dat$col <- alpha("black",0.3)
  dat[gnA, "col"] <- alpha(clustCols,alpha=.8)[1]
  dat[gnB, "col"] <- alpha(clustCols,alpha=.8)[2]
  col <- as.character(dat$col)
  names(col) <- as.character(dat$col)
  datLab <- dat[c(gnA,gnB), ]
  datLab[gnA, "col"] <- alpha(clustCols,alpha=.8)[1]
  datLab[gnB, "col"] <- alpha(clustCols,alpha=.8)[2]
  
  xLab <- "Difference in detection rate"
  yLab <- paste0("Gene expression ratio ", temp_exp)
  mainLab <- paste0("Expression difference effect sizes (",clA," vs. ",clB,")")
  
  theme_set(theme_classic(base_family = "noto-sans-jp", 
                          base_size = 12))
  # Change the settings
  update_geom_defaults("text", list(family = theme_get()$text$family))
  
  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(aes(color = col)) +
    scale_color_manual(values = col) +
    annotate("text", label = paste0("Higher in ", clA), size = 4,
             color = clustCols[1],
             x = Inf, y = -Inf, hjust = 1.03, vjust = -0.5,
             family= theme_get()$text[["family"]]) +
    annotate("text", label = paste0("Higher in ", clB), size = 4,
             color = clustCols[2],
             x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
             family= theme_get()$text[["family"]]) +
    theme(text = element_text(family = "noto-sans-jp"),
          legend.position = "none") +
    geom_text_repel(data = datLab, mapping = aes(x, y, label = gene, color = col),
                    size = 4, fontface = "bold", family = "noto-sans-jp"
    ) +
    labs(title = mainLab) +
    xlab(xLab) +
    ylab(yLab) +
    geom_vline(xintercept = 0, color = "gray50")  
  print(p)
}

#' scClustViz plot: Plot volcano plot
#'
#' This function plots a volcano plot comparing cluster A to B.
#'
#' @param sCVd scClustViz object
#' @param clA cluster A for the first comparison
#' @param clB cluster B for the second comparison
#' @param dataType For MA-style plots comparing difference and mean of gene
#'   summary statistics, one of: \code{"DR"} (detection rate); \code{"MGE"}
#'   (mean gene expression); \code{"MDGE"} (mean detected gene expression). For
#'   volcano plots, the effect size measure can be one of: \code{"dDR"}
#'   (difference in detection rate); \code{"logGER"} (log gene expression
#'   ratio). To compare relationship between difference in detection rate and
#'   log gene expression ratio, use \code{"GERvDDR"}.
#' @param labType Default="de". A character vector indicating which genes to
#'   highlight. One of \code{"de"} (most statistically significant genes),
#'   \code{"diff"} (most different by dataType shown), or \code{"search"}
#'   (specified genes).
#' @param labGenes Only required if \code{labType="search"}. Gene names to
#'   highlight.
#' @param labNum Default=5. Number of genes to highlight per side.
#' @param labTypeDiff Default="logGER". Only required if
#'   \code{dataType="GERvDDR"} and \code{labType="diff"}. Which axis to use for
#'   difference calculation. One of \code{"dDR"} (difference in detection rate)
#'   or \code{"logGER"} (log gene expression ratio).

#'
#' @examples
#' \dontrun{
#' # Cluster overlay:
#' plot_compareClusts_DEscatter2(sCVd = seurat_sc[[1]], clA = 1, clB = 2,
#'                            dataType = "DR", labType = "diff",
#'                            labNum = 5, labTypeDiff = "logGER")
#'
#' }
#'
#' @export

plot_compareClusts_volcano2 <- function(sCVd, clA, clB, dataType, labType, 
                                        labNum, labGenes, sizeFactor = 1) {
  CGS <- compareClusts_DF(sCVd = sCVd, clA = clA, clB = clB, dataType = dataType)
  CGS <- CGS[!is.na(CGS$FDR), ]
  CGS$FDR <- -log10(CGS$FDR)
  CGS$FDR[CGS$FDR == Inf] <- 300 # Min positive non-zero value is 1e-300
  
  temp_exp <- paste0("(log", Param(sCVd,"exponent"), " scale)")
  if(Param(sCVd,"exponent") == exp(1)) {
    temp_exp <- "(natural log scale)"
  }
  
  gnA <- rownames(head(CGS[order(CGS[[dataType]], decreasing = T), ], labNum))
  gnB <- rownames(tail(CGS[order(CGS[[dataType]], decreasing = T), ], labNum))
  if(labType == "de") {
    ts <- order(CGS$FDR, decreasing = T, na.last = T)
    gnA <- rownames(CGS)[ts[CGS[ts, "dir"] == clA][1:labNum]]
    gnB <- rownames(CGS)[ts[CGS[ts, "dir"] == clB][1:labNum]]
  }
  
  clustCols <- scales::hue_pal()(2)
  dat <- data.frame(x = CGS[[dataType]], y = CGS$FDR, gene = rownames(CGS))
  rownames(dat) <- dat$gene
  dat$col <- alpha("black", 0.3)
  dat[gnA, "col"] <- alpha(clustCols, alpha = 0.8)[1]
  dat[gnB, "col"] <- alpha(clustCols, alpha = 0.8)[2]
  col <- as.character(dat$col)
  names(col) <- as.character(dat$col)
  datLab <- dat[c(gnA, gnB), ]
  xLab <- "Difference in detection rate"
  if(dataType == "logGER") xLab <- paste0("Gene expression ratio ", temp_exp)
  yLab <- "-log10 FDR-adjusted p-value"
  mainLab <- paste0("Volcano plot of DE genes ", clA, " vs. ", clB)
  
  theme_set(theme_classic(base_family = "noto-sans-jp", 
                          base_size = 12))
  # Change the settings
  update_geom_defaults("text", list(family = theme_get()$text$family))
  
  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(aes(color = col)) +
    scale_color_manual(values = col) +
    annotate("text", label = paste0("Higher in ", clA), size = 8 * sizeFactor,
             color = clustCols[1],
             x = Inf, y = -Inf, hjust = 1.03, vjust = -0.5,
             family= theme_get()$text[["family"]]) +
    annotate("text", label = paste0("Higher in ", clB), size = 8 * sizeFactor,
             color = clustCols[2],
             x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
             family= theme_get()$text[["family"]]) +
    geom_text_repel(
      data = datLab, mapping = aes(x, y, label = gene, color = col),
      size = 5 * sizeFactor, fontface = "bold", family = "noto-sans-jp",
      max.overlaps = 30
    ) +
    labs(title = mainLab) +
    xlab(xLab) +
    ylab(yLab) +
    geom_vline(xintercept = 0, color = "gray50") +
    theme(
      text = element_text(family = "noto-sans-jp"),
      axis.text = element_text(size = 15 * sizeFactor),
      axis.title = element_text(size = 20 * sizeFactor),
      plot.title = element_text(size = 25 * sizeFactor),
      legend.position = "none"
    )
  print(p)
}

#' scClustViz plot: Plot to compare cell metadata
#'
#' This function makes scatter/boxplots comparing cellular metadata.
#'
#' @param MD A dataframe of cellular metadata. See \code{\link{getMD}}.
#' @param mdX A character vector of one refering to the variable name from
#'   \code{MD} to plot on the x-axis.
#' @param mdY A character vector of one refering to the variable name from
#'   \code{MD} to plot on the y-axis.
#' @param sel_cells Optional. A character vector of cell names (rownames of
#'   \code{MD}) to highlight in the plot.
#' @param sel_clust Optional. The name of the selected cluster
#'   (\code{sel_cells}) to include in the legend. If
#'   \code{\link{labelCellTypes}} has been run, pass the appropriate element of
#'   \code{attr(Clusters(sCV),"ClusterNames")} to this argument to show both
#'   cluster number and cell type label in the legend.
#' @param md_log Optional. A character vector indicating which axes should be
#'   log scaled. \code{c("x","y")} to log-scale both axes.
#'
#' @examples
#' \dontrun{
#' plot_mdCompare(MD=getMD(input_data_obj),
#'                mdX="total_counts",
#'                mdY="total_features",
#'                sel_cells=names(Clusters(sCVdata))[Clusters(sCVdata) == "1"],
#'                sel_clust="1",
#'                md_log="x")
#' }
#'
#' @export

plot_mdCompare <- function(MD,mdX,mdY,sel_cells,sel_clust,md_log,size) {
  if (missing(sel_cells)) { sel_cells <- "" }
  if (missing(sel_clust)) { sel_clust <- "" }
  if (missing(md_log)) { md_log <- "" }
  MD <- data.frame(MD[,c(mdX,mdY)])
  MD$sel_cells <- rownames(MD) %in% sel_cells
  if ("x" %in% md_log & !(is.factor(MD[,1]) | is.character(MD[,1]))) {
    tempLX <- "x"
    if (any(MD[,1] <= 0)) {
      names(MD)[1] <- paste(names(MD)[1],
                            paste0("(log scale: ",sum(MD[,1] <= 0),
                                   " values <= 0 omitted)"))
      MD <- MD[MD[,1] > 0,]
    } else {
      names(MD)[1] <- paste(names(MD)[1],"(log scale)")
    }
  } else {
    tempLX <- ""
  }
  if ("y" %in% md_log & !(is.factor(MD[,2]) | is.character(MD[,2]))) { 
    tempLY <- "y" 
    if (any(MD[,2] <= 0)) {
      names(MD)[2] <- paste(names(MD)[2],
                            paste0("(log scale: ",sum(MD[,2] <= 0),
                                   " values <= 0 omitted)"))
      MD <- MD[MD[,2] > 0,]
    } else {
      names(MD)[2] <- paste(names(MD)[2],"(log scale)")
    } 
  } else {
    tempLY <- ""
  }
  md_log <- paste(c(tempLX,tempLY),collapse="")
  
  if ((is.factor(MD[,1]) | is.character(MD[,1])) &
      (is.factor(MD[,2]) | is.character(MD[,2]))) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,"This figure is not designed to compare to categorical variables.")
  } else if (is.factor(MD[,1]) | is.character(MD[,1])) {
    plot_mdBoxplotX(MD,sel_clust,md_log,size)
  } else if (is.factor(MD[,2]) | is.character(MD[,2])) {
    plot_mdBoxplotY(MD,sel_clust,md_log,size)
  } else {
    plot_mdScatter(MD,sel_clust,md_log,size)
  }
  
}

#' scClustViz plot: Plot to view cellular metadata by cluster
#'
#' This function makes boxplots / stacked barplots of cellular metadata
#' separated by cluster.
#'
#' @param MD A dataframe of cellular metadata. See \code{\link{getMD}}.
#' @param sel A character vector of one refering to the variable name from
#'   \code{MD} to plot.
#' @param cl A factor of cluster assignments. See \code{\link{Cluster}}.
#' @param opt Default="absolute". A character vector of plotting options. One of
#'   \code{"absolute"}, \code{"relative"}, or \code{"y"}. \code{"y"} sets
#'   log-scales the data for postive numerical metadata. For categorical
#'   metadata, \code{"absolute"} plots a stacked barplot of raw counts, whereas
#'   \code{"relative"} plots the proportion of each cluster represented by each
#'   category.
#'
#' @examples
#' \dontrun{
#' plot_mdPerClust(MD=getMD(input_data_obj),
#'                 sel="cyclonePhases",
#'                 cl=Clusters(sCVdata),
#'                 opt="relative")
#' }
#'
#' @export

plot_mdPerClust <- function(MD,sel,cl,opt="absolute",inp,size) {
  MD <- MD[sel]
  MD$cl <- as.factor(cl)
  if ("y" %in% opt & !(is.factor(MD[,1]) | is.character(MD[,1]))) { 
    if (any(MD[,1] <= 0)) {
      names(MD)[1] <- paste(names(MD)[1],
                            paste0("(log scale: ",sum(MD[,1] <= 0),
                                   " values <= 0 omitted)"))
      MD <- MD[MD[,1] > 0,]
    } else {
      names(MD)[1] <- paste(names(MD)[1],"(log scale)")
    } 
  }
  if (is.factor(MD[,1]) | is.character(MD[,1])) {
    plot_mdBarplot(MD,opt,size)
  } else {
    plot_mdBoxplot(MD,opt,inp,size)
  }
}

# DE gene dotplot

#' scClustViz plot helper function: Return DE genes per cluster
#'
#' This function returns a named numeric vector of FDR-corrected p-values for
#' statistically significant differentially expressed genes for a set comparison
#' type and FDR threshold. For \code{"DEmarker"}, the returned value is the max
#' of all comparisons.
#'
#' @param sCVd The sCVdata object.
#' @param DEtype One of: \code{"DEvsRest"} - see \code{\link{DEvsRest}};
#'   \code{"DEneighb"} - see \code{\link{DEneighb}}; \code{"DEmarker"} - see
#'   \code{\link{DEmarker}}.
#' @param FDRthresh A numeric vector of length 1 setting a false discovery rate
#'   threshold for statistical significance.
#'
#' @examples
#' \dontrun{
#' dotplotDEgenes(sCVdata,
#'                DEtype="DEneighb",
#'                FDRthresh=0.01)
#' }
#'
#' @export

dotplotDEgenes <- function(sCVd,DEtype,FDRthresh) {
  if (missing(FDRthresh)) { FDRthresh <- 1 }
  if (DEtype == "DEvsRest") {
    return(lapply(DEvsRest(sCVd),function(X) {
      temp <- X[which(X$FDR <= FDRthresh),"FDR",drop=F]
      out <- unlist(temp,use.names=F)
      names(out) <- rownames(temp)
      return(sort(out))
    }))
  } else if (DEtype == "DEneighb") {
    outL <- lapply(DEneighb(sCVd,FDRthresh), function(X) {
      if (nrow(X) < 1) { return(numeric(0)) }
      out <- X[,grep("^FDR_",names(X))]
      names(out) <- rownames(X)
      return(sort(out))
    })
    names(outL) <- levels(Clusters(sCVd))
    return(outL)
  } else if (DEtype == "DEmarker") {
    outL <- lapply(DEmarker(sCVd,FDRthresh), function(X) {
      if (nrow(X) < 1) { return(numeric(0)) }
      out <- apply(X[,grep("^FDR_",names(X)),drop=F],1,max)
      return(sort(out))
    })
    return(outL)
  }
}

#' scClustViz plot: Plot gene expression dotplots.
#'
#' This function makes dotplots (a heatmap analogue) showing gene expression for
#' a set of genes across all clusters.
#'
#' When generated in an interactive context (i.e. RStudio), this can sometimes
#' result in a \code{figure margins too large} error. See example for suggested
#' dimensions of the graphic device.
#'
#' @param sCVd The sCVdata object.
#' @param DEgenes The output of \code{\link{dotplotDEgenes}}.
#' @param DEnum Single integer representing the maximum number of DE genes per
#'   cluster to include in the dotplot.
#'
#' @examples
#' \dontrun{
#' pdf("filepath.pdf",width=11,height=7)
#' plot_deDotplot(sCVd=sCVdata,
#'                DEgenes=dotplotDEgenes(sCVdata,
#'                                       DEtype="DEneighb",
#'                                       FDRthresh=0.01)
#'                DEnum=5)
#' dev.off()
#' }
#'
#' @export

plot_deDotplot <- function(sCVd,DEgenes,DEnum,size) {
  # ^ Setup ----
  heatGenes <- unique(unlist(lapply(DEgenes,function(X) names(X)[1:DEnum])))
  heatGenes <- heatGenes[!is.na(heatGenes)]
  
  if (is.null(heatGenes)) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,
         "No genes were statistically significant at the current false discovery rate")
    return(invisible())
  }
  
  temp_DR <- sapply(ClustGeneStats(sCVd),function(X) X[heatGenes,"DR"])
  if (is.vector(temp_DR)) {
    temp_DR <- matrix(temp_DR,1,dimnames=list(NULL,names(temp_DR))) 
  }
  temp_MDGE <- sapply(ClustGeneStats(sCVd),function(X) X[heatGenes,"MDGE"])
  if (is.vector(temp_MDGE)) {
    temp_MDGE <- matrix(temp_MDGE,1,dimnames=list(NULL,names(temp_MDGE))) 
  }
  rownames(temp_DR) <- rownames(temp_MDGE) <- heatGenes
  
  if (nrow(temp_DR) > 1) {
    hG <- hclust(stats::dist(temp_DR),"complete")
  } else {
    hG <- list(order=1)
  }
  
  hC <- hclust(stats::as.dist(DEdist(sCVd)),"single")
  
  clustCols <- colorspace::qualitative_hcl(length(levels(Clusters(sCVd))),palette="Dark 3")
  
  dC <- dendrapply(as.dendrogram(hC),function(X) {
    if (is.leaf(X)) {
      attr(X,"edgePar") <- list(
        lwd=2,
        col=clustCols[which(attr(X,"label") == levels(Clusters(sCVd)))]
      )
      attr(X,"nodePar") <- list(
        pch=NA,lab.font=2,lab.cex=size,
        lab.col=clustCols[which(attr(X,"label") == levels(Clusters(sCVd)))])
      if (attr(X,"label") != "Unselected") {
        if (attr(X,"label") %in% names(DEgenes)) {
          attr(X,"label") <- paste0(attr(X,"label"),": ",
                                    length(DEgenes[[attr(X,"label")]])," DE")
        } else {
          attr(X,"label") <- paste0(
            attr(X,"label"),": ",
            length(DEgenes[[which(attr(X,"label") ==
                                    sapply(strsplit(names(DEgenes),"-"),
                                           function(X) X[1]))]]),
            " DE")
        }
      }
    }
    return(X)
  })
  
  if ("genes" %in% names(ClustGeneStats(sCVd)[[1]])) {
    tempLabCol <- ClustGeneStats(sCVd)[[1]][heatGenes,"genes"]
  } else {
    tempLabCol <- rownames(ClustGeneStats(sCVd)[[1]][heatGenes,])
  }
  DR <- temp_DR[hG$order,hC$order,drop=F]
  temp <- range(sapply(ClustGeneStats(sCVd),function(X) X[,"MDGE"]))
  temp <- seq(temp[1],temp[2],length.out=101)
  MDGE <- findInterval(as.vector(temp_MDGE[hG$order,hC$order]),
                       vec=temp,all.inside=T)
  
  # ^ Plot dotplot ----
  temp_par <- par(no.readonly=T)
  if (length(levels(Clusters(sCVd))) <= 1) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Heatmap cannot be computed",
                     "with less than two clusters.",sep="\n"))
  } else if (length(heatGenes) < 1) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("There are no differentially expressed genes at",
                     "false discovery rate threshold."))
  } else {
    layout(matrix(c(0,2,3,1),2),widths=c(1,5),heights=c(1,5))
    par(mar=c(9,0,0,.5))
    plot(x=NULL,y=NULL,xlim=c(0.5,nrow(DR)+.5),ylim=c(0.5,ncol(DR)+.5),
         xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
    abline(v=1:nrow(DR),col="grey90")
    symbols(x=rep(1:nrow(DR),ncol(DR)),
            y=as.vector(sapply(1:ncol(DR),function(X) rep(X,nrow(DR)))),
            circles=as.vector(DR)/2,inches=F,add=T,xpd=NA,
            fg=colorspace::sequential_hcl(100,palette="Viridis",rev=T)[MDGE],
            bg=colorspace::sequential_hcl(100,palette="Viridis",rev=T)[MDGE])
    axis(side=1,at=1:nrow(DR),lwd=0,labels=tempLabCol[hG$order],las=2,cex.axis=size)
    
    # Legend:
    tx0 <- par("usr")[1]
    tx <- (par("usr")[2] - par("usr")[1])
    ty0 <- par("usr")[3]
    ty <- par("usr")[4] - par("usr")[3]
    segments(x0=tx0 - seq(.15,.03,length.out=1000) * tx,
             y0=ty0 - 0.02 * ty,y1=ty0 - 0.05 * ty,
             col=colorspace::sequential_hcl(1000,palette="Viridis",rev=T),xpd=NA)
    text(x=tx0 - c(.15,.09,.03) * tx,
         y=ty0 - c(0.035,0.02,0.035) * ty,
         labels=c(round(min(temp_MDGE),2),
                  "Mean detected expression",
                  round(max(temp_MDGE),2)),pos=2:4,xpd=NA)
    symbols(x=tx0 - c(.15,.09,.03) * tx,
            y=ty0 - rep(.14,3) * ty,add=T,xpd=NA,
            circles=c(0.25,0.5,0.75)/2,inches=F,bg="black")
    text(x=tx0 - c(.149,.089,0.029,.09) * tx,
         y=ty0 - c(rep(.23,3),.26) * ty,xpd=NA,
         labels=c("25%","50%","75%","Detection Rate"))
    
    
    par(mar=c(9,0,0,7))
    plot(dC,horiz=T,xpd=NA,
         ylim=c(0.5,length(hC$order)+.5),yaxs="i",yaxt="n")
    
    par(mar=c(0,0,0,.5))
    if (class(hG) == "hclust") {
      plot(as.dendrogram(hG),leaflab="none",
           xlim=c(0.5,length(hG$order)+.5),xaxs="i",yaxt="n")
    }
  }
  par(temp_par)
}

#' geneSearch: 
#'
#' This function generates boxplots comparing normalized gene abundance across
#' all clusters.
#'
#' @param txt The gene of interest.
#' @param st The search type: For this app, only gene list not regex.
#' @param CGS The gene stats by cluster from seurat_sc.
#'
#' @examples
#' \dontrun{
#' geneSearch(txt = "A2M", 
#'            st = "",
#'            CGS = ClustGeneStats(d$SCV[[d$res]])[[1]])
#' }
#'
#' @export

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

# Gene expression boxplots

#' scClustViz plot: Compare gene expression across clusters
#'
#' This function generates boxplots comparing normalized gene abundance across
#' all clusters.
#'
#' @param nge The gene expression matrix, see \code{\link{getExprs}}.
#' @param sCVd The sCVdata object.
#' @param gene The gene to display.
#' @param geneName Optional. A named character vector of length one. The element
#'   is the full gene name, and the name is the gene symbol.
#' @param opts Default=\code{c("sct","dr")}. A character vector with plotting
#'   options. If it includes \code{"sct"}, data points will be overlaid as a
#'   jitter over the boxplot. If it includes \code{"dr"}, detection rate per
#'   cluster will be plotted as a small black bar over each boxplot, with the
#'   corresponding axis on the right.
#'
#' @examples
#' \dontrun{
#' plot_GEboxplot(getExpr(input_data_obj),
#'                sCVd=sCVdata,
#'                gene="Actb")
#' }
#'
#' @export

plot_GEboxplot <- function(nge,sCVd,gene,geneName,opts=c("sct","dr"),size) {
  if (gene == "") {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Select a gene by either clicking on the plot above",
                     "or searching for genes of interest in the search bar above,",
                     "then pick the gene from the list just above this figure",
                     "to see a comparison of that gene's expression across all clusters.",
                     sep="\n"))
  } else {
    # ^ setup -----
    hC <- hclust(stats::as.dist(DEdist(sCVd)),"single")
    temp_ylab <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                        "TRUE"="(natural log scale)",
                        "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
    temp_pos <- switch(as.character(length(levels(Clusters(sCVd))) > 1),
                       "TRUE"=hC$order,"FALSE"=1)
    if ("sct" %in% opts) {
      len <- length(levels(Clusters(sCVd)))
      bxpCol <- hue_pal()(len)
    } else {
      len <- length(levels(Clusters(sCVd)))
      bxpCol <- hue_pal()(len)
    }
    
    # ^ plot boxplot -----
    temp_par <- par(no.readonly=T)
    layout(matrix(2:1,nrow=2),heights=c(1,4))
    par(mar=c(6.5,4.5,0,4), family = "noto-sans-jp", cex = size)
    suppressWarnings(boxplot(
      vector("list",length(levels(Clusters(sCVd))[levels(Clusters(sCVd)) != "Unselected"])),
      ylim=range(nge[gene,]),xaxt="n",xlab=NA,
      ylab=paste(gene,"normalized gene expression", temp_ylab),
      cex.lab = 0.75
    ))
    mtext(levels(Clusters(sCVd))[temp_pos], side=1, line=1, at=seq_along(temp_pos), cex = size)
    mtext("Clusters, ordered by heatmap dendrogram", side=1, line=3, cex = size)
    if (missing(geneName)) { geneName <- NULL }
    if (is.null(geneName)) { 
      mtext(paste(gene,collapse="\n"),side = 1,line = 5, cex = size) 
    } else {
      mtext(paste(paste0(names(geneName),": ",geneName), collapse="\n"),
            side=1,line=2,cex = size) 
    }
    for (i in temp_pos) {
      if ("sct" %in% opts) {
        points(jitter(rep(which(temp_pos == i),
                          sum(Clusters(sCVd) %in% levels(Clusters(sCVd))[i])),
                      amount=.2),
               nge[gene,Clusters(sCVd) %in% levels(Clusters(sCVd))[i]],
               pch=".",cex=3,
               col= "black")
      }
      boxplot(nge[gene,Clusters(sCVd) %in% levels(Clusters(sCVd))[i]],add=T,
              at=which(temp_pos == i),col=bxpCol[i],outline=F)
    }
    if ("dr" %in% opts) {
      points(x=seq_along(ClustGeneStats(sCVd)),
             y=sapply(ClustGeneStats(sCVd)[temp_pos],function(X) X[gene,"DR"]) * 
               max(nge[gene,]) + min(nge[gene,]),
             pch="-",cex=2)
      axis(side=4,at=seq(0,1,.25) * max(nge[gene,]) + min(nge[gene,]),
           labels=paste0(seq(0,1,.25) * 100,"%"), cex.axis = 1)
      mtext(side=4,line=3, text="- Gene detection rate per cluster", cex = size)
    }
    if (length(temp_pos) > 1) { 
      par(mar=c(0,4.5,1,3) +.3, family = "noto-sans-jp", cex = size)
      plot(as.dendrogram(hC),leaflab="none") 
    }
    par(temp_par)
  }
}

#' scClustViz plot: Compare gene expression across clusters
#'
#' This function extracts cluster comparison data for selected
#' clusters A and B based on data type.
#'
#' @param sCVd scClustViz object
#' @param clA cluster A for the first comparison
#' @param clB cluster B for the second comparison
#' @param dataType For MA-style plots comparing difference and mean of gene
#'   summary statistics, one of: \code{"DR"} (detection rate); \code{"MGE"}
#'   (mean gene expression); \code{"MDGE"} (mean detected gene expression). For
#'   volcano plots, the effect size measure can be one of: \code{"dDR"}
#'   (difference in detection rate); \code{"logGER"} (log gene expression
#'   ratio). To compare relationship between difference in detection rate and
#'   log gene expression ratio, use \code{"GERvDDR"}.
#'
#' @examples
#' \dontrun{
#' compareClusts_DF(sCVd = seurat_sc[[1]], clA = 1, clB = 1, dataType = "DR")
#' }
#'
#' @export

compareClusts_DF <- function(sCVd,clA,clB,dataType) {
  if (dataType %in% c("MGE","MDGE","DR")) {
    loc1 <- c(paste(clA,clB,sep="-"),paste(clB,clA,sep="-"))
    loc <- loc1[loc1 %in% names(DEcombn(sCVd))]
    loc1 <- which(loc1 %in% names(DEcombn(sCVd)))
    if (loc1 == 2) { loc1 <- -1 }
    if ("Wstat" %in% colnames(DEcombn(sCVd)[[loc]])) {
      tempW <- DEcombn(sCVd)[[loc]]$Wstat - 
        DEcombn(sCVd)[[loc]]$Wstat[which.max(DEcombn(sCVd)[[loc]]$pVal)]
    } else {
      tempW <- DEcombn(sCVd)[[loc]]$logGER
    }
    temp <- data.frame(x_diff=ClustGeneStats(sCVd)[[clA]][,dataType] - 
                         ClustGeneStats(sCVd)[[clB]][,dataType],
                       y_mean=rowMeans(cbind(ClustGeneStats(sCVd)[[clA]][,dataType],
                                             ClustGeneStats(sCVd)[[clB]][,dataType])),
                       logGER=NA,FDR=NA,dir=NA)
    rownames(temp) <- rownames(ClustGeneStats(sCVd)[[clA]])
    temp[rownames(DEcombn(sCVd)[[loc]]),"logGER"] <- DEcombn(sCVd)[[loc]]$logGER
    temp[rownames(DEcombn(sCVd)[[loc]]),"FDR"] <- DEcombn(sCVd)[[loc]]$FDR
    temp[rownames(DEcombn(sCVd)[[loc]]),"dir"] <- c(clB,clA)[(tempW * loc1 > 0) + 1]
    return(temp)
  } else if (dataType %in% c("GERvDDR","logGER","dDR")) {
    loc1 <- which(c(paste(clA,clB,sep="-"),paste(clB,clA,sep="-")) %in% names(DEcombn(sCVd)))
    if (loc1 == 2) { loc1 <- -1 }
    loc <- which(names(DEcombn(sCVd)) %in% c(paste(clA,clB,sep="-"),paste(clB,clA,sep="-")))
    temp <- DEcombn(sCVd)[[loc]][,c("logGER","dDR","FDR")]
    temp <- as.data.frame(mapply("*",temp,c(loc1,loc1,1))) 
    rownames(temp) <- rownames(DEcombn(sCVd)[[loc]])
    if ("Wstat" %in% colnames(DEcombn(sCVd)[[loc]])) {
      tempW <- DEcombn(sCVd)[[loc]]$Wstat - 
        DEcombn(sCVd)[[loc]]$Wstat[which.max(DEcombn(sCVd)[[loc]]$pVal)]
    } else {
      tempW <- DEcombn(sCVd)[[loc]]$logGER
    }
    temp$dir <- c(clB,clA)[(tempW * loc1 > 0) + 1]
    return(temp)
  } 
}

#' scClustViz plot: Volcano and MA-style plots to compare clusters
#'
#' This function generates scatterplots inspired by volcano and MA plots for
#' comparing gene expression between pairs of clusters.
#'
#' @param sCVd The sCVdata object.
#' @param clA Cluster identifier for side A of the comparison.
#' @param clB Cluster identifier for side B of the comparison.
#' @param dataType For MA-style plots comparing difference and mean of gene
#'   summary statistics, one of: \code{"DR"} (detection rate); \code{"MGE"}
#'   (mean gene expression); \code{"MDGE"} (mean detected gene expression). For
#'   volcano plots, the effect size measure can be one of: \code{"dDR"}
#'   (difference in detection rate); \code{"logGER"} (log gene expression
#'   ratio). To compare relationship between difference in detection rate and
#'   log gene expression ratio, use \code{"GERvDDR"}.
#' @param labType Default="de". A character vector indicating which genes to
#'   highlight. One of \code{"de"} (most statistically significant genes),
#'   \code{"diff"} (most different by dataType shown), or \code{"search"}
#'   (specified genes).
#' @param labGenes Only required if \code{labType="search"}. Gene names to
#'   highlight.
#' @param labNum Default=5. Number of genes to highlight per side.
#' @param labTypeDiff Default="logGER". Only required if
#'   \code{dataType="GERvDDR"} and \code{labType="diff"}. Which axis to use for
#'   difference calculation. One of \code{"dDR"} (difference in detection rate)
#'   or \code{"logGER"} (log gene expression ratio).
#'
#' @examples
#' \dontrun{
#' plot_compareClusts(sCVdata,
#'                    clA="1",
#'                    clB="2",
#'                    dataType="GERvDDR",
#'                    labType="search",
#'                    labGenes="Actb")
#' }
#'
#' @export

plot_compareClusts2 <- function(sCVd, clA, clB, dataType,
                                labType = "de", labGenes,
                                labNum = 5, labTypeDiff = "logGER",
                                sizeFactor = 1) {
    if(dataType == "MGE") {
      plot_compareClusts_MAplot2(
        sCVd = sCVd, 
        clA = clA, 
        clB = clB, 
        labType = labType, 
        labNum = labNum, 
        labGenes = labGenes,
        sizeFactor = sizeFactor
      )
    } else {
      plot_compareClusts_volcano2(
        sCVd = sCVd, 
        clA = clA, 
        clB = clB, 
        dataType = dataType, 
        labType = labType, 
        labNum = labNum, 
        labGenes = labGenes,
        sizeFactor = sizeFactor
      )
    }
}

#################################################################################################

## Profiler module

# Profiler UI
profilerUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 2,
      uiOutput(ns("showHeatmap")),
      uiOutput(ns("showDimRed")),
      uiOutput(ns("showViolin")),
      uiOutput(ns("showModuleDescriptions"))
    ),
    mainPanel(
      width = 10,
      div(
        style = "display: inline-block; float: right;",
        actionLink(ns("profiler_help_button"), label = "Help")
      ),
      tabsetPanel(
        id = ns("profiler_tabsetPanel"),
        type = "hidden",
        tabPanelBody(
          value = "profiler_heatmap",
          fluidPage(
            fluidRow(uiOutput(ns("heatmap_moduleScore_header"))),
            fluidRow(
              div(
                style = "display: inline-block; width: 200px; vertical-align: top;",
                uiOutput(ns("heatmap_profiler_factor")),
              ),
              div(
                style = "display: inline-block; width: 240px; margin-left: 20px; vertical-align: top;",
                uiOutput(ns("heatmapcolors_profiler"))
              ),
              div(
                style = "display: inline-block; width: 150px; margin-left: 20px; margin-top: 30px; vertical-align: top;",
                uiOutput(ns("select_all_profiler"))
              ),
              div(
                style = "display: inline-block; width: 280px; margin-left: 20px;",
                uiOutput(ns("profile_select"))
              )
            ),
            fluidRow(
              div(
                style = "margin-top: 0px; vertical-align: top;",
                withSpinner(
                  plotOutput(ns("heatmap_moduleScore"), width = 800, height = 800),
                  type = 8, size = 1, color = "black", proxy.height = "300px"
                )
              )
            ),
            fluidRow(
              div(
                style = "display: inline-block;",
                uiOutput(ns("heatmap_moduleScore_download_pdf"))
              ),
              div(
                style = "display: inline-block;",
                uiOutput(ns("heatmap_moduleScore_download_png"))
              )
            ),
            br()
          )
        ),
        tabPanelBody(
          value = "profiler_dimred",
          fluidPage(
            fluidRow(uiOutput(ns("dimred_profiler_header"))),
            fluidRow(
              column(
                width = 6,
                fluidRow(
                  div(
                    style = "display: inline-block; width: 100px;",
                    uiOutput(ns("dimredmethod_profiler"))
                  ),
                  div(
                    style = "display: inline-block; width: 200px; margin-left: 10px;",
                    uiOutput(ns("dimredmodule_profiler"))
                  )
                ),
                fluidRow(
                  plotOutput(ns("dimred_profiler"), width = "100%", height = "400px")
                ),
                br(),
                fluidRow(
                  div(
                    style = "display: inline-block;",
                    uiOutput(ns("dimred_profiler_download_pdf"))
                  ),
                  div(
                    style = "display: inline-block;",
                    uiOutput(ns("dimred_profiler_download_png"))
                  )
                )
              ),
              column(
                width = 6,
                fluidRow(
                  div(
                    style = "width: 200px;",
                    uiOutput(ns("dimred_profiler_factor"))
                  )
                ),
                fluidRow(
                  plotOutput(ns("dimred_profiler_2"), width = "100%", height = "400px")
                ),
                br(),
                fluidRow(
                  div(
                    style = "display: inline-block;",
                    uiOutput(ns("dimred_profiler_2_download_pdf"))
                  ),
                  div(
                    style = "display: inline-block;",
                    uiOutput(ns("dimred_profiler_2_download_png"))
                  )
                )
              )
            ),
            br()
          )
        ),
        tabPanelBody(
          value = "profiler_violin",
          fluidPage(
            fluidRow(uiOutput(ns("violinplot_profiler_header"))),
            fluidRow(
              div(
                style = "display: inline-block; width: 200px;",
                uiOutput(ns("violinplot_profiler_factor"))
              ),
              div(
                style = "display: inline-block; width: 200px; margin-left: 10px;",
                uiOutput(ns("violinplotmodule_profiler"))
              )
            ),
            fluidRow(plotOutput(ns("violinplot_profiler"), width = 700, height = 500)),
            fluidRow(
              div(
                style = "display: inline-block;",
                uiOutput(ns("violinplot_profiler_download_pdf"))
              ),
              div(
                style = "display: inline-block; margin-left: 0px;",
                uiOutput(ns("violinplot_profiler_download_png"))
              )                  
            ),
            br()
          )
        ),
        tabPanelBody(
          value = "profiler_moduleDescriptions",
          fluidPage(
            fluidRow(
              div(
                style = "display: inline-block;",
                uiOutput(ns("module_descriptions_header"))
              ),
              div(
                style = "display: inline-block; float: right; margin-right: 20px;",
                uiOutput(ns("gene_list"))
              )
            ),
            fluidRow(DT::dataTableOutput(ns("moduleDescriptionsTable"), width = 1000)),
            br()
          )
        )
      )
    )
  )
}

# Plot module score heatmap
PlotModuleScoreHeatmap <- function(moduleScoreSeurat, modules, fact, plotColors) {
  labAngle <- 45
  hjust <- 0
  # If max number of characters in group labels is 2 or less 
  # (e.g., cluster numbers), print labels horizontally
  maxChars <- max(nchar(as.character(unique(moduleScoreSeurat@meta.data[, fact]))))
  if(maxChars <= 2) {
    labAngle <- 0
    hjust <- 0.5
  }
  
  # If > 10 group labels, use legend, else use put group labels above colorbar
  rmLegend <- NULL
  plotMargins <- c(1,1,1,1)
  colLabs <- F
  if(length(unique(moduleScoreSeurat@meta.data[, fact])) <= 10) {
    colLabs <- T
    rmLegend <- guides(color = F)
    plotMargins <- c(maxChars/6,1,1,1)
  }
  if(maxChars <= 2) {
    labAngle <- 0
    hjust <- 0.5
  }
  
  theme_set(theme_classic(base_family = "noto-sans-jp", base_size = 12))
  # Change the settings
  update_geom_defaults("text", list(family = theme_get()$text$family))
  
  g <- DoHeatmap(
    object = moduleScoreSeurat, 
    slot = "counts",
    features = modules,
    group.by = fact,
    label = colLabs,
    angle = labAngle,
    draw.lines = F,
    hjust = hjust,
    size = 4
  ) +
    theme(plot.margin = grid::unit(x = plotMargins, units = "cm"), legend.key = element_blank()) +
    labs(fill = "Module score") +
    scale_fill_gradientn(colors = plotColors) +
    rmLegend
  
  if(!colLabs) {
    g$labels$colour <- fact
    g$guides$colour$order <- 1
  }
  
  myLevels <- gsub("-", replacement = " ", x = levels(g$data$Feature))
  myFeatures <- gsub("-", replacement = " ", x = g$data$Feature)
  g$data$Feature <- factor(myFeatures, levels = myLevels)
  
  return(g)
}

# Plot module score dimred
ModuleDimRedPlot <- function(
  seuratModuleScore, markerList, method = "pca", module = "Cluster1", 
  pctVars = NULL, sizeFactor = 1) {
  
  if(method == "pca") {
    if(is.null(pctVars)) pctVars <- SeuratVarExplained(seuratModuleScore)
    xLab <- paste0("PC 1: ", round(pctVars[1]), " % variance")
    yLab <- paste0("PC 2: ", round(pctVars[2]), " % variance")
  } else if(method == "umap") {
    xLab <- "UMAP 1"
    yLab <- "UMAP 2"
  } else if(method == "tsne") {
    xLab <- "tSNE 1"
    yLab <- "tSNE 2"
  }
  
  myTitle <- names(markerList)[paste0("Cluster", 1:length(markerList)) == module]
  myTitle <- gsub(pattern = "-", replacement = " ", x = myTitle)
  
  g <- FeaturePlot(
    object = seuratModuleScore, 
    reduction = method, 
    features = module
  ) +
    scale_colour_gradientn(
      colours = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
    ) +
    xlab(xLab) + ylab(yLab) +
    labs(title = myTitle, color = "Module score") +
    theme(
      text = element_text(family = "noto-sans-jp"),
      axis.text = element_text(size = 15 * sizeFactor),
      axis.title = element_text(size = 20 * sizeFactor),
      plot.title = element_text(size = 25 * sizeFactor),
      legend.text = element_text(size = 15 * sizeFactor),
      legend.title = element_text(size = 20 * sizeFactor),
      aspect.ratio = 1,
      plot.margin = unit(c(0,0,0,0), units = "cm")
    )
  return(g)
}

# Feature plot with same format as module dimred plot
FactorDimRedPlot <- function(
  seuratData, fact, method = "pca", pctVars = NULL, sizeFactor = 1) {
  if(method == "pca") {
    if(is.null(pctVars)) pctVars <- SeuratVarExplained(seuratModuleScore)
    xLab <- paste0("PC 1: ", round(pctVars[1]), " % variance")
    yLab <- paste0("PC 2: ", round(pctVars[2]), " % variance")
  } else if(method == "umap") {
    xLab <- "UMAP 1"
    yLab <- "UMAP 2"
  } else if(method == "tsne") {
    xLab <- "tSNE 1"
    yLab <- "tSNE 2"
  }
  
  g <- DimPlot(
    object = seuratData, 
    reduction = method, 
    group.by = fact
  ) +
    xlab(xLab) + ylab(yLab) +
    labs(title = fact) +
    theme(
      text = element_text(family = "noto-sans-jp"),
      axis.text = element_text(size = 15 * sizeFactor),
      axis.title = element_text(size = 20 * sizeFactor),
      plot.title = element_text(size = 25 * sizeFactor, hjust = 0.5),
      legend.text = element_text(size = 15 * sizeFactor),
      legend.title = element_text(size = 20 * sizeFactor),
      aspect.ratio = 1,
      plot.margin = unit(c(0,0,0,0), units = "cm")
    )
  
  return(g)
}

# Profiler server
profiler <- function(input, output, session, seuratData) {
  
  output$showHeatmap <- renderUI({
    req(input$profiler_tabsetPanel)
    myLab <- "Heatmap"
    if(input$profiler_tabsetPanel == "profiler_heatmap") myLab <- strong(myLab)
    actionLink(
      inputId = session$ns("showHeatmap"), 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  }) 
  
  output$showDimRed <- renderUI({
    req(input$profiler_tabsetPanel)
    myLab <- "PCA/tSNE/UMAP"
    if(input$profiler_tabsetPanel == "profiler_dimred") myLab <- strong(myLab)
    actionLink(
      inputId = session$ns("showDimRed"), 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  }) 
  
  output$showViolin <- renderUI({
    req(input$profiler_tabsetPanel)
    myLab <- "Violin plot"
    if(input$profiler_tabsetPanel == "profiler_violin") myLab <- strong(myLab)
    actionLink(
      inputId = session$ns("showViolin"), 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  }) 
  
  output$showModuleDescriptions <- renderUI({
    req(input$profiler_tabsetPanel)
    myLab <- "Module info"
    if(input$profiler_tabsetPanel == "profiler_moduleDescriptions") myLab <- strong(myLab)
    actionLink(
      inputId = session$ns("showModuleDescriptions"), 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })         
  
  observeEvent(input$showHeatmap, {
    updateTabsetPanel(inputId = "profiler_tabsetPanel", selected = "profiler_heatmap")
  })
  
  observeEvent(input$showDimRed, {
    updateTabsetPanel(inputId = "profiler_tabsetPanel", selected = "profiler_dimred")
  })
  
  observeEvent(input$showViolin, {
    updateTabsetPanel(inputId = "profiler_tabsetPanel", selected = "profiler_violin")
  })
  
  observeEvent(input$showModuleDescriptions, {
    updateTabsetPanel(inputId = "profiler_tabsetPanel", selected = "profiler_moduleDescriptions")
  })
  
  # Returns Seurat input data. If seuratData is a list, run CreateSeuratObject
  # using counts and metadata
  GetSeurat <- reactive({
    req(seuratData)
    if(is.list(seuratData)) {
      dat <- CreateSeuratObject(counts = seuratData[["counts"]], 
                                meta.data = seuratData[["metadata"]])
      dat <- NormalizeData(dat)
      dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
      dat <- ScaleData(dat)
      nPCs <- min(50, nrow(dat)-1, ncol(dat)-1)
      dat <- RunPCA(dat, features = VariableFeatures(object = dat), npcs = nPCs)
      dat <- RunUMAP(dat, dims = 1:min(10, nPCs), n.neighbors = min(30, nrow(dat), ncol(dat)), verbose = F)
      dat <- RunTSNE(dat, perplexity = min(30, (ncol(dat)-1)/3), check_duplicates = F)
      return(dat)
    }
    return(seuratData)
  })
  
  # Returns the maximum valid value of "ctrl" parameter for AddModuleScore.
  # We then use either this value or 100 (default for AddModuleScore) for 
  # "ctrl", whichever is smaller, in GetSeuratWithModuleScore().
  # AddModuleScore splits all genes into 24 bins based on average expression 
  # across cells/samples, so the "ctrl" value cannot exceed the size of the 
  # smallest bin. The code used here is from the Seurat AddModuleScore function.
  MaxCtrlVal <- reactive({
    req(GetSeurat())
    dataAvg <- rowMeans(GetAssayData(GetSeurat()))
    dataCut <- cut_number(x = dataAvg + rnorm(n = length(dataAvg))/1e30, n = 24, 
                          labels = FALSE, right = FALSE)
    nGenesPerBin <- as.integer(table(dataCut))
    return(min(nGenesPerBin))
  })
  
  # Returns seurat object with module score column
  GetSeuratWithModuleScore <- reactive({
    req(GetSeurat(), MarkerList(), MaxCtrlVal())
    ctrlVal <- min(MaxCtrlVal(), 100)
    AddModuleScore(GetSeurat(), features = MarkerList(), ctrl = ctrlVal)
  })
  
  # Returns module scores stored as Seurat object
  ModuleScoreSeurat <- reactive({
    req(GetSeurat())
    moduleNamesRefVec <- paste0("Cluster", 1:length(MarkerList()))
    names(moduleNamesRefVec) <- names(MarkerList())
    names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
    moduleScoreMat <- t(GetSeuratWithModuleScore()@meta.data[, moduleNamesRefVec])
    rownames(moduleScoreMat) <- names(MarkerList())
    seuratData_moduleScore <- CreateSeuratObject(moduleScoreMat)
    seuratData_moduleScore@meta.data[, colnames(GetSeurat()@meta.data)] <-
      GetSeurat()@meta.data
    return(seuratData_moduleScore)
  })
  
  # Help for iPSC profiler
  observeEvent(input$profiler_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("iPSC profiler")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/profiler_help.md")
    ))
  })
  
  # Returns list of vectors of markers genes with module names as names
  MarkerList <- reactive({
    req(GetSeurat())
    markerList <- as.list(readWorkbook("./data/iPSC_Profiler_curated_markers.xlsx", sheet = 1))
    markerList <- lapply(markerList, function(markers) markers[!is.na(markers)])
    markerList <- lapply(markerList, function(markers) markers[markers %in% rownames(GetSeurat())])
    markerList <- markerList[sapply(markerList, length) > 0]
    names(markerList) <- gsub("_", replacement = "-", x = names(markerList))
    
    return(markerList)
  })
  
  # Returns dataframe containing module names and descriptions
  ModuleDescriptions <- reactive({
    req(MarkerList())
    moduleDescriptionsDF <- readWorkbook("./data/iPSC_Profiler_curated_markers.xlsx", sheet = 2)
    rownames(moduleDescriptionsDF) <- moduleDescriptionsDF$Module
    moduleDescriptionsDF$Module <- gsub("_", replacement = "-", x = moduleDescriptionsDF$Module)
    moduleDescriptionsDF <- moduleDescriptionsDF[moduleDescriptionsDF$Module %in% names(MarkerList()), ]
    moduleDescriptionsDF$Module <- gsub("-", replacement = " ", x = moduleDescriptionsDF$Module)
    return(moduleDescriptionsDF)
  })
  
  # Check box to use all modules. If checked, all modules are used in heatmap.
  # If not, user selects individual modules
  output$select_all_profiler <- renderUI({
    req(MarkerList())
    prettyCheckbox(session$ns("select_all_profiler"), label = "Use all modules", 
                   value = T, status = "default", icon = icon("check"))
    
  })
  
  # If select_all_profiler unchecked, user selects individual modules
  output$profile_select <- renderUI({
    req(MarkerList())
    if(is.null(input$select_all_profiler)) return()
    if(input$select_all_profiler) return()
    choices <- names(MarkerList())
    names(choices) <- gsub("-", replacement = " ", x = choices)
    selectInput(session$ns("profile_select"), label = "Select modules", 
                choices = choices, 
                selected = NULL,
                multiple = TRUE)
  })
  
  # Returns vector of module names selected by user (either all or individual)
  SelectedModules <- reactive({
    if(is.null(input$select_all_profiler)) return()
    if(input$select_all_profiler) return(names(MarkerList()))
    return(input$profile_select)
  })
  
  # Updates module selection when select_all_profiler is unchecked
  observeEvent(input$select_all_profiler, {
    if(!(input$select_all_profiler)) {
      choices <- names(MarkerList())
      names(choices) <- gsub("-", replacement = " ", x = choices)
      updateSelectInput(session, inputId = "profile_select", label = "Select profiles",
                        choices = choices, selected = NULL)
    } 
  })
  
  # Heatmap header
  output$heatmap_moduleScore_header <- renderUI({
    h3("Heatmap of module scores", style = "margin-top: 0px;")      
  })
  
  # Select factor for heatmap colorbar
  output$heatmap_profiler_factor <- renderUI({
    req(GetSeurat())
    choices <- colnames(GetSeurat()@meta.data)
    choices <- choices[!(choices %in% c("nCount_RNA","nFeature_RNA"))]
    choices <- choices[!(apply(GetSeurat()@meta.data[, choices, drop = FALSE], 2, is.numeric))]
    if(is.list(seuratData)) choices <- choices[choices != "orig.ident"]
    selected <- choices[1]
    selectInput(session$ns("heatmap_profiler_factor"), label = "Grouping factor",
                choices = choices, selected = selected)
  })
  
  # Colors for heatmap
  output$heatmapcolors_profiler <- renderUI({
    choices <- c("Blue-white-red (default)" = "Default", 
                 "Red-white-blue" = "Red-white-blue", 
                 "Viridis" = "Viridis", 
                 "Green-yellow-red" = "Green-yellow-red"
    )
    selectInput(
      inputId = session$ns("heatmapcolors_profiler"),
      label = "Color palette",
      choices = choices,
      selected = "Default"
    )
  })
  
  # IPSC - download handler for ipsc profiler heatmap (PDF)
  output$gene_list <- renderUI({
    req(ModuleScoreSeurat(), SelectedModules(), MarkerList())
    downloadButton(session$ns("gene_list_download"), label = "Download iPSC profiler gene lists")
  })
  
  output$gene_list_download <- downloadHandler(
    filename = function() {
      "iPSC_profiler_gene_list.csv"
    }, content = function(file) {
      if(is.null(GetSeurat())) return()
      markerList <- rbindlist(lapply(1:length(MarkerList()), function(x) {
        y <- data.frame("Module" = names(MarkerList())[x],
                        "Genes" = MarkerList()[[x]])
      }))
      write.csv(markerList, file = file, row.names = F, quote = F)
    }
  )
  
  # Returns color palette based on input
  heatCols_profiler <-  reactive({
    if(is.null(input$heatmapcolors_profiler)) return()
    if (input$heatmapcolors_profiler == "Red-white-blue") {
      cols <- colorRampPalette(c("red", "white", "blue"))(100)
    } else if (input$heatmapcolors_profiler == "Viridis") {
      cols <- viridis(100)
    } else if (input$heatmapcolors_profiler == "Green-yellow-red") {
      cols <- rev(brewer.pal(n=11, name = "RdYlGn"))
    } else {
      cols <- rev(brewer.pal(n = 11, name = "RdBu"))
    }
    return(cols)
  })
  
  # Render heatmap
  output$heatmap_moduleScore <- renderPlot({
    req(ModuleScoreSeurat(), input$heatmap_profiler_factor)
    if(!(input$heatmap_profiler_factor %in% colnames(ModuleScoreSeurat()@meta.data))) return()
    if(length(SelectedModules()) == 0 || is.null(heatCols_profiler())) {
      ggplot() + xlim(0,1) + ylim(0,1) + theme_void() + 
        annotate("text", x = 0.5, y = 1, size = 6,
                 label = "Please select at least one module") +
        theme(text = element_text(family = "noto-sans-jp"))
    } else {
      PlotModuleScoreHeatmap(
        moduleScoreSeurat = ModuleScoreSeurat(),
        modules = SelectedModules(),
        fact = input$heatmap_profiler_factor,
        plotColors = heatCols_profiler()
      )
      
    }
  })
  
  # Download heatmap as pdf
  output$heatmap_moduleScore_download_pdf <- renderUI({
    req(ModuleScoreSeurat(), input$heatmap_profiler_factor)
    if(!(input$heatmap_profiler_factor %in% colnames(ModuleScoreSeurat()@meta.data))) return()
    if(length(SelectedModules()) == 0 || is.null(heatCols_profiler())) return()
      downloadButton(session$ns("heatmap_moduleScore_pdf_img"), "Download plot (PDF)")
  })
  output$heatmap_moduleScore_pdf_img <- downloadHandler(
    filename = "Heatmap_ModuleScores.pdf",
    content = function(file) {
      pdf(file, width = 8, height = 11)
      g <- PlotModuleScoreHeatmap(moduleScoreSeurat = ModuleScoreSeurat(),
                                  modules = SelectedModules(),
                                  fact = input$heatmap_profiler_factor,
                                  plotColors = heatCols_profiler())
      print(g)
      dev.off()
    }
  )
  
  # Download heatmap as png
  output$heatmap_moduleScore_download_png <- renderUI({
    req(ModuleScoreSeurat(), input$heatmap_profiler_factor)
    if(!(input$heatmap_profiler_factor %in% colnames(ModuleScoreSeurat()@meta.data))) return()
    if(length(SelectedModules()) == 0 || is.null(heatCols_profiler())) return()
    downloadButton(session$ns("heatmap_moduleScore_png_img"), "Download plot (PNG)")
  })
  output$heatmap_moduleScore_png_img <- downloadHandler(
    filename = "Heatmap_ModuleScores.png",
    content = function(file) {
      png(file, width = 600, height = 800)
      g <- PlotModuleScoreHeatmap(moduleScoreSeurat = ModuleScoreSeurat(),
                                  modules = SelectedModules(),
                                  fact = input$heatmap_profiler_factor,
                                  plotColors = heatCols_profiler())
      print(g)
      dev.off()
    }
  )
  
  # Dimred header
  output$dimred_profiler_header <- renderUI({
    h3("Dimensional reduction", style = "margin-top: 0px;")
  })
  
  # Dimred method
  output$dimredmethod_profiler <- renderUI({
    selectInput(
      inputId = session$ns("dimredmethod_profiler"),
      label = "Method",
      choices = c("PCA" = "pca", "tSNE" = "tsne", "UMAP" = "umap"),
      selected = "pca"
    )
  })
  
  # Module for dimred
  output$dimredmodule_profiler <- renderUI({
    choices <- paste0("Cluster", 1:length(MarkerList()))
    names(choices) <- names(MarkerList())
    names(choices) <- gsub("-", replacement = " ", x = names(choices))
    selectInput(
      inputId = session$ns("dimredmodule_profiler"),
      label = "Module",
      choices = choices
    )
  })
  
  # Returns % variance explained by principal components
  GetPctVarianceExplained <- reactive({
    req(GetSeuratWithModuleScore())
    SeuratVarExplained(GetSeuratWithModuleScore())
  })
  
  # Render dimred plot
  output$dimred_profiler <- renderPlot({
    req(GetSeuratWithModuleScore(), GetPctVarianceExplained(), input$dimredmethod_profiler, 
        input$dimredmodule_profiler)

    if(!(input$dimredmodule_profiler %in% colnames(GetSeuratWithModuleScore()@meta.data))) return()
    
    ModuleDimRedPlot(
      seuratModuleScore = GetSeuratWithModuleScore(),
      markerList = MarkerList(),
      method = input$dimredmethod_profiler,
      module = input$dimredmodule_profiler,
      pctVars = GetPctVarianceExplained()
    )
  })
  
  # Download dimred plot as pdf
  output$dimred_profiler_download_pdf <- renderUI({
    req(GetSeuratWithModuleScore(), GetPctVarianceExplained(), input$dimredmethod_profiler, 
        input$dimredmodule_profiler)
    downloadButton(session$ns("dimred_profiler_download_pdf_img"), "Download plot (PDF)")
  })
  output$dimred_profiler_download_pdf_img <- downloadHandler(
    filename = paste0(toupper(input$dimredmethod_profiler), "_Profiler.pdf"),
    content = function(file) {
      p <- ModuleDimRedPlot(
        seuratModuleScore = GetSeuratWithModuleScore(),
        markerList = MarkerList(),
        method = input$dimredmethod_profiler,
        module = input$dimredmodule_profiler,
        pctVars = GetPctVarianceExplained(),
        sizeFactor = 0.7
      )
      
      ggsave(file, plot = p, width = 6, height = 6)
    }
  )   
  
  # Download dimred plot as png
  output$dimred_profiler_download_png <- renderUI({
    req(GetSeuratWithModuleScore(), GetPctVarianceExplained(), input$dimredmethod_profiler, 
        input$dimredmodule_profiler)
    downloadButton(session$ns("dimred_profiler_download_png_img"), "Download plot (PNG)")
  })
  output$dimred_profiler_download_png_img <- downloadHandler(
    filename = paste0(toupper(input$dimredmethod_profiler), "_Profiler.png"),
    content = function(file) {
      
      p <- ModuleDimRedPlot(
        seuratModuleScore = GetSeuratWithModuleScore(),
        markerList = MarkerList(),
        method = input$dimredmethod_profiler,
        module = input$dimredmodule_profiler,
        pctVars = GetPctVarianceExplained(),
        sizeFactor = 2
      )
      
      ggsave(file, plot = p, width = 6, height = 6, bg = "white")
    }
  ) 
  
  # Dimred factor
  output$dimred_profiler_factor <- renderUI({
    req(GetSeurat())
    choices <- colnames(GetSeurat()@meta.data)
    choices <- choices[!(choices %in% c("nCount_RNA","nFeature_RNA"))]
    choices <- choices[!(apply(GetSeurat()@meta.data[, choices, drop = F], 2, is.numeric))]
    if(is.list(seuratData)) choices <- choices[choices != "orig.ident"]
    selected <- choices[1]
    selectInput(session$ns("dimred_profiler_factor"), label = "Grouping factor",
                choices = choices, selected = selected)
  })
  
  # Second dimred plot
  output$dimred_profiler_2 <- renderPlot({
    req(GetSeurat(), GetSeuratWithModuleScore(), GetPctVarianceExplained(), 
        input$dimredmethod_profiler, input$dimred_profiler_factor)
    
    FactorDimRedPlot(
      seuratData = GetSeurat(),
      fact = input$dimred_profiler_factor,
      method = input$dimredmethod_profiler,
      pctVars = GetPctVarianceExplained()
    )
    
  })
  
  # Download second dimred plot as pdf
  output$dimred_profiler_2_download_pdf <- renderUI({
    req(GetSeuratWithModuleScore(), GetPctVarianceExplained(), input$dimredmethod_profiler, 
        input$dimredmodule_profiler)
    downloadButton(session$ns("dimred_profiler_2_download_pdf_img"), "Download plot (PDF)")
  })
  output$dimred_profiler_2_download_pdf_img <- downloadHandler(
    filename = paste0(toupper(input$dimredmethod_profiler), "_Profiler_Grouped.pdf"),
    content = function(file) {
      
      p <- FactorDimRedPlot(
        seuratData = GetSeurat(),
        fact = input$dimred_profiler_factor,
        method = input$dimredmethod_profiler,
        pctVars = GetPctVarianceExplained(),
        sizeFactor = 0.7
      )
      
      ggsave(file, plot = p, width = 6, height = 6)
    }
  )   
  
  # Download second dimred plot as png
  output$dimred_profiler_2_download_png <- renderUI({
    req(GetSeuratWithModuleScore(), GetPctVarianceExplained(), input$dimredmethod_profiler, 
        input$dimredmodule_profiler)
    downloadButton(session$ns("dimred_profiler_2_download_png_img"), "Download plot (PNG)")
  })
  output$dimred_profiler_2_download_png_img <- downloadHandler(
    filename = paste0(toupper(input$dimredmethod_profiler), "_Profiler_Grouped.png"),
    content = function(file) {
      
      p <- FactorDimRedPlot(
        seuratData = GetSeurat(),
        fact = input$dimred_profiler_factor,
        method = input$dimredmethod_profiler,
        pctVars = GetPctVarianceExplained(),
        sizeFactor = 2
      )
      
      ggsave(file, plot = p, width = 6, height = 6)
    }
  ) 
  
  # Violin plot header
  output$violinplot_profiler_header <- renderUI({
    req(GetSeuratWithModuleScore())
    h3("Violin plot", style = "margin-top: 0px;")
  })
  
  # Violin plot factor
  output$violinplot_profiler_factor <- renderUI({
    req(GetSeurat())
    choices <- colnames(GetSeurat()@meta.data)
    choices <- choices[!(choices %in% c("nCount_RNA","nFeature_RNA"))]
    choices <- choices[!(apply(GetSeurat()@meta.data[, choices], 2, is.numeric))]
    if(is.list(seuratData)) choices <- choices[choices != "orig.ident"]
    selected <- choices[1]
    selectInput(session$ns("violinplot_profiler_factor"), label = "Grouping factor",
                choices = choices, selected = selected)
    
  })
  
  # Module for violin plot
  output$violinplotmodule_profiler <- renderUI({
    req(MarkerList())
    moduleNamesRefVec <- paste0("Cluster", 1:length(MarkerList()))
    names(moduleNamesRefVec) <- names(MarkerList())
    names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
    choices <- moduleNamesRefVec
    selectInput(
      inputId = session$ns("violinplotmodule_profiler"),
      label = "Module",
      choices = choices
    )
  })
  
  # Render violin plot
  output$violinplot_profiler <- renderPlot({
    req(GetSeuratWithModuleScore(), input$violinplot_profiler_factor, 
        input$violinplotmodule_profiler)
    labAngle <- 45
    hjust <- NULL
    if(max(nchar(as.character(unique(GetSeuratWithModuleScore()@meta.data[, input$violinplot_profiler_factor])))) <= 2) {
      labAngle <- 0
      hjust <- 0.5
    }
    
    moduleNamesRefVec <- paste0("Cluster", 1:length(MarkerList()))
    names(moduleNamesRefVec) <- names(MarkerList())
    names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
    
    VlnPlot(GetSeuratWithModuleScore(), features = input$violinplotmodule_profiler, 
            group.by = input$violinplot_profiler_factor) +
      theme(text = element_text(family = "noto-sans-jp", 
                                size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text.x = element_text(angle = labAngle, hjust = hjust, size = 12)) +
      NoLegend() +
      xlab(input$violinplot_profiler_factor) +
      ylab("Module score") +
      labs(title = names(moduleNamesRefVec)[moduleNamesRefVec == input$violinplotmodule_profiler])
  })
  
  # Download violin plot as pdf
  output$violinplot_profiler_download_pdf <- renderUI({
    req(GetSeuratWithModuleScore(), input$violinplot_profiler_factor, 
        input$violinplotmodule_profiler)
    downloadButton(session$ns("violinplot_profiler_pdf_img"), "Download plot (PDF)")
  })
  output$violinplot_profiler_pdf_img <- downloadHandler(
    filename = "ViolinPlot_Profiler.pdf",
    content = function(file) {
      labAngle <- 45
      hjust <- NULL
      if(max(nchar(as.character(unique(GetSeuratWithModuleScore()@meta.data[, input$violinplotmodule_profiler])))) <= 2) {
        labAngle <- 0
        hjust <- 0.5
      }
      
      moduleNamesRefVec <- paste0("Cluster", 1:length(MarkerList()))
      names(moduleNamesRefVec) <- names(MarkerList())
      names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
      
      pdf(file, width = 6, height = 4)
      g <- VlnPlot(GetSeuratWithModuleScore(), features = input$violinplotmodule_profiler, 
                   group.by = input$violinplot_profiler_factor) +
        theme(text = element_text(family = "noto-sans-jp", 
                                  size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(angle = labAngle, hjust = hjust, size = 12)) +
        NoLegend() +
        xlab(input$violinplot_profiler_factor) +
        ylab("Module score") +
        labs(title = names(moduleNamesRefVec)[moduleNamesRefVec == input$violinplotmodule_profiler])
      print(g)
      dev.off()
    }
  )
  
  # Download violin plot as png
  output$violinplot_profiler_download_png <- renderUI({
    req(GetSeuratWithModuleScore(), input$violinplot_profiler_factor, 
        input$violinplotmodule_profiler)
    downloadButton(session$ns("violinplot_profiler_png_img"), "Download plot (PNG)")
  })
  output$violinplot_profiler_png_img <- downloadHandler(
    filename = "ViolinPlot_Profiler.png",
    content = function(file) {
      labAngle <- 45
      hjust <- NULL
      if(max(nchar(as.character(unique(GetSeuratWithModuleScore()@meta.data[, input$violinplotmodule_profiler])))) <= 2) {
        labAngle <- 0
        hjust <- 0.5
      }
      
      moduleNamesRefVec <- paste0("Cluster", 1:length(MarkerList()))
      names(moduleNamesRefVec) <- names(MarkerList())
      names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
      
      png(file, width = 600, height = 400)
      g <- VlnPlot(GetSeuratWithModuleScore(), features = input$violinplotmodule_profiler, 
                   group.by = input$violinplot_profiler_factor) +
        theme(text = element_text(family = "noto-sans-jp", 
                                  size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(angle = labAngle, hjust = hjust, size = 12)) +
        NoLegend() +
        xlab(input$violinplot_profiler_factor) +
        ylab("Module score") +
        labs(title = names(moduleNamesRefVec)[moduleNamesRefVec == input$violinplotmodule_profiler])
      print(g)
      dev.off()
    }
  )
  
  # Module descriptions table header
  output$module_descriptions_header <- renderUI({
    req(GetSeuratWithModuleScore())
    h3("Module descriptions", style = "margin-top: 0px;")
  })
  
  # Module descriptions table
  output$moduleDescriptionsTable = DT::renderDataTable({
    req(ModuleDescriptions())
    df <- ModuleDescriptions()
    df$Source[!is.na(df$Source)] <- 
      paste0("<a href='", df$Source[!is.na(df$Source)], "' target='_blank'>", df$Text[!is.na(df$Source)],"</a>")
    DT::datatable(df[, c("Module", "Description", "Source")], rownames = F, escape = F,
                  selection = "none",
                  options = list(lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
                                 ordering = F, 
                                 scrollCollapse = T)) %>%
      DT::formatStyle("Source","white-space"="nowrap")
  })
  
}

#################################################################################################

## Dimensional reduction module

dimredUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      div(
        style = "display: inline-block; width: 100px;",
        uiOutput(ns("dimredmethod"))
      ),
      div(
        style = "display: inline-block; width: 200px; margin-left: 10px;",
        uiOutput(ns("dimredfact"))
      ),
      div(
        style = "display: inline-block; vertical-align: top; margin-top: 23px;",
        uiOutput(ns("addSample_button")),
      ),
      div(
        style = "display: inline-block; vertical-align: top; margin-top: 23px;",
        uiOutput(ns("clearSamples_button")),
      ),
      div(
        style = "display: inline-block; float: right;",
        uiOutput(ns("datasetLabel"))
      )
    ),
    fluidRow(uiOutput(ns("chosenSample"))),
    fluidRow(uiOutput(ns("selectedSamples"))),
    fluidRow(
      div(style = "display: inline-block;", uiOutput(ns("dimredPlot_ui")))
    ),
    br(),
    fluidRow(
      div(style = "display: inline-block;", uiOutput(ns("dlqcdimredpdf"))),
      div(style = "display:inline-block", uiOutput(ns("dlqcdimredpng")))
    ),
    br()
  )
}

#------------------------------------------------------------------
# PCA plot functions
#------------------------------------------------------------------

# PLOT - returns a vector with length equal to the number of PCs; 
# percent of the total variance explained by each PC
SeuratVarExplained = function(seuratData) {
  totalVariance = sum(rowVars(seuratData@assays$RNA@scale.data))
  eigValues = (seuratData@reductions$pca@stdev)^2
  pctVarExplained = eigValues / totalVariance * 100
  return(pctVarExplained)
}

# PLOT - PCA plot (bulk)
PCABulk <- function(tmp, fact, sourceName = "A") {
  pca.lab <- plotPCA(
    object = tmp,
    intgroup = fact
  )
  pca <- plotPCA(
    object = tmp,
    intgroup = fact,
    returnData = TRUE
  )
  len <- length(unique(pca$group))
  tooltips <- paste0(
    "<b>Sample:</b> ", rownames(pca), "<br />",
    "<b>PC 1:</b> ", round(pca$PC1, 3), "<br />",
    "<b>PC 2:</b> ", round(pca$PC2, 3)
  )
  
  font <- list(
    family = "Noto Sans JP",
    size = 12,
    color = "white"
  )
  label <- list(
    bgcolor = "transparent",
    bordercolor = "transparent",
    font = font
  )
  
  plot_ly(
    data = pca,
    type = "scatter",
    mode = "markers",
    x = ~PC1,
    y = ~PC2,
    marker = list(size = 9),
    text = tooltips,
    hoverinfo = "text",
    hoverlabel = label,
    color = ~group,
    colors = hue_pal()(len),
    source = sourceName
  ) %>%
    layout(
      xaxis = list(title = gsub("PC1", "PC 1", pca.lab$labels$x)),
      yaxis = list(title = gsub("PC2", "PC 2", pca.lab$labels$y)),
      font = list(family = "Noto Sans JP")
    )
}

# PLOT - PCA plot (rasl)
PCARASL <- function(tmp, meta, fact) {
  pca.lab <- plotPCA_matrix(
    mat = tmp,
    meta = meta,
    intgroup = fact
  )
  pca <- plotPCA_matrix(
    mat = tmp,
    meta = meta,
    intgroup = fact,
    returnData = TRUE
  )
  len <- length(unique(pca$group))
  tooltips <- paste0(
    "<b>Sample:</b> ", rownames(pca), "<br />",
    "<b>PC 1:</b> ", round(pca$PC1, 3), "<br />",
    "<b>PC 2:</b> ", round(pca$PC2, 3)
  )
  
  font <- list(
    family = "Noto Sans JP",
    size = 12,
    color = "white"
  )
  label <- list(
    bgcolor = "transparent",
    bordercolor = "transparent",
    font = font
  )
  
  plot_ly(
    data = pca,
    type = "scatter",
    mode = "markers",
    x = ~PC1,
    y = ~PC2,
    marker = list(size = 9),
    text = tooltips,
    hoverinfo = "text",
    hoverlabel = label,
    color = ~group,
    colors = hue_pal()(len)
  ) %>%
    layout(
      xaxis = list(title = gsub("PC1", "PC 1", pca.lab$labels$x)),
      yaxis = list(title = gsub("PC2", "PC 2", pca.lab$labels$y)),
      font = list(family = "Noto Sans JP")
    )
}

# PLOT - PCA plot (sc only)
PCASeurat <- function(tmp, fact) {
  varExplained = SeuratVarExplained(tmp)
  pca <- DimPlot(
    tmp,
    reduction = "pca",
    group.by = fact,
    label = T
  )
  
  pca <- pca + xlab(paste0("PC 1: ", round(varExplained[1]), "% variance")) +
    ylab(paste0("PC 2: ", round(varExplained[2]), "% variance")) +
    theme_light()
  
  font <- list(
    family = "Noto Sans JP",
    size = 12,
    color = "white"
  )
  label <- list(
    bgcolor = "transparent",
    bordercolor = "transparent",
    font = font
  )
  
  gg <- ggplotly(pca, tooltip = c(fact, "x", "y"), hoverinfo = "text") %>%
    layout(font = list(family = "Noto Sans JP")) %>%
    style(hoverlabel = label)
  return(gg)
}

# PLOT - plot bulk RNA-seq PC1 vs PC2
qcPCAPlotBulk <- function(tmp, fact) {
  pca.lab <- plotPCA(
    tmp,
    intgroup = fact
  )
  pca <- plotPCA(
    tmp,
    intgroup = fact,
    returnData = TRUE
  )
  p <- ggplot(pca, aes(PC1, PC2)) +
    geom_point(aes(color = group), size = 2) +
    xlab(pca.lab$labels$x) +
    ylab(pca.lab$labels$y) +
    ggtitle("Principle Component Analysis") +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  print(p)
}

# PLOT - plot rasl RNA-seq PC1 vs PC2
qcPCAPlotRASL <- function(mat, meta, fact) {
  pca.lab <- plotPCA_matrix(
    mat = mat,
    meta = meta,
    intgroup = fact
  )
  pca <- plotPCA_matrix(
    mat = mat,
    meta = meta,
    intgroup = fact,
    returnData = TRUE
  )
  p <- ggplot(pca, aes(PC1, PC2)) +
    geom_point(aes(color = group), size = 2) +
    xlab(pca.lab$labels$x) +
    ylab(pca.lab$labels$y) +
    ggtitle("Principle Component Analysis") +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  print(p)
}

# PLOT - plot scRNA-seq PC1 vs PC2
qcPCAPlotSeurat <- function(tmp, fact) {
  varExplained = SeuratVarExplained(tmp)
  p <- DimPlot(tmp, reduction = "pca", label = T, group.by = fact) +
    ggtitle("Principle Component Analysis") +
    xlab(paste0("PC1: ", round(varExplained[1]), "% variance")) +
    ylab(paste0("PC2: ", round(varExplained[2]), "% variance")) +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
}

#------------------------------------------------------------------
# tSNE plot functions (tSNE)
#------------------------------------------------------------------

# PLOT - tSNE plot (sc only)
TSNESeurat <- function(tmp, fact) {
  tsne <- DimPlot(
    tmp,
    reduction = "tsne",
    group.by = fact,
    label = T
  ) 
  tsne <- tsne + xlab("tSNE 1") +
    ylab("tSNE 2") +
    theme_light()
  
  font <- list(
    family = "Noto Sans JP",
    size = 12,
    color = "white"
  )
  label <- list(
    bgcolor = "transparent",
    bordercolor = "transparent",
    font = font
  )
  
  gg <- ggplotly(tsne, tooltip = c(fact, "x", "y"), hoverinfo = "text") %>%
    layout(font = list(family = "Noto Sans JP")) %>%
    style(hoverlabel = label)
  return(gg)
}

# PLOT - tSNE plot (sc only)
qcTSNEPlotSeurat <- function(tmp, fact) {
  tsne <- DimPlot(
    tmp,
    reduction = "tsne",
    group.by = fact,
    label = T
  ) 
  tsne <- tsne + xlab("tSNE coordinate 1") +
    ggtitle("tSNE Analysis") +
    ylab("tSNE coordinate 2") +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(tsne)
}

#------------------------------------------------------------------
# uMAP plot functions
#------------------------------------------------------------------

# PLOT - UMAP plot (sc)
UMAPSeurat <- function(tmp, fact) {
  umapPlot <- DimPlot(
    tmp,
    reduction = "umap",
    group.by = fact,
    label = T
  ) 
  umapPlot <- umapPlot + xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_light()
  
  font <- list(
    family = "Noto Sans JP",
    size = 12,
    color = "white"
  )
  label <- list(
    bgcolor = "transparent",
    bordercolor = "transparent",
    font = font
  )
  
  gg <- ggplotly(umapPlot, tooltip = c(fact, "x", "y"), hoverinfo = "text") %>%
    layout(font = list(family = "Noto Sans JP")) %>%
    style(hoverlabel = label)
  return(gg)
}

# PLOT - UMAP plot (Seurat)
qcUMAPPlotSeurat <- function(tmp, fact) {
  umapPlot <- DimPlot(
    tmp,
    reduction = "umap",
    group.by = fact,
    label = T
  ) 
  umapPlot <- umapPlot + xlab("UMAP 1") +
    ggtitle("UMAP Analysis") +
    ylab("UMAP 1") +
    theme_light() +
    theme(text = element_text(family = "noto-sans-jp"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(umapPlot)
}

dimred <- function(input, output, session, dat, dimredMethods = "pca", 
                   plotWidth = 800, plotHeight = 600, datasetLabel = NULL,
                   returnData = F) {
  
  d <- reactiveValues(
    plotClick = NULL,
    selectedSample = NULL,
    selectedSamples = NULL
  )
  
  output$datasetLabel <- renderUI({
    req(datasetLabel)
    h4(datasetLabel, style = "margin-top: 0px;")
  })
  
  # Return metadata
  GetMetadata <- reactive({
    req(dat)
    if(is.list(dat)) return(dat$metadata)
    return(dat@meta.data)
  })
  
  # Returns named vector of methods to allow
  GetMethods <- reactive({
    myMethods = tolower(dimredMethods)
    names(myMethods)[myMethods == "pca"] <- "PCA"
    names(myMethods)[myMethods == "tsne"] <- "tSNE"
    names(myMethods)[myMethods == "umap"] <- "UMAP"
    return(myMethods)
  })
  
  # DDA - choose dimensionality reduction method
  output$dimredmethod <- renderUI({
    req(GetMethods())
    selectInput(
      inputId = session$ns("dimredmethod"),
      label = "Method",
      choices = GetMethods(),
      selected = GetMethods()[1]
    )
  })
  
  # DDA - select input - choose factor - all dimensionality reduction methods
  output$dimredfact <- renderUI({
    req(GetMetadata())
    choices <- colnames(GetMetadata())
    choices <- choices[!(apply(GetMetadata(), 2, is.numeric) | choices %in% c("nCount_RNA","nFeature_RNA"))]

    selectInput(
      inputId = session$ns("dimredfact"),
      label = "Grouping factor",
      choices = choices
    )
  })
  
  # Button to add selected sample to list
  output$addSample_button <- renderUI({
    if(is.null(returnData) || !returnData || 
       is.null(d$selectedSample)) return()
    div(
      style = "margin-left: 20px;",
      actionButton(
        inputId = session$ns("addSample_button"),
        label = "+ Add sample"
      )
    )
  })
  
  observeEvent(input$addSample_button, {
    d$selectedSamples <- unique(c(d$selectedSamples, d$selectedSample))
  })
  
  # Button to clear selected samples
  output$clearSamples_button <- renderUI({
    if(is.null(returnData) || !returnData || 
       is.null(d$selectedSamples)) return()
    div(
      style = "margin-left: 20px;",
      actionButton(
        inputId = session$ns("clearSamples_button"),
        label = "- Clear samples"
      )
    )
  })
  observeEvent(input$clearSamples_button, {
    d$selectedSamples <- NULL
  })
  
  output$chosenSample <- renderUI({
    if(is.null(returnData) || !returnData || is.null(d$selectedSample)) return()
    p(strong("Current selection:"), d$selectedSample)
  })
  
  output$selectedSamples <- renderUI({
    if(is.null(returnData) || !returnData || is.null(d$selectedSamples)) return()
    txt <- paste(d$selectedSamples, collapse = ", ")
    p(strong("Selected samples:"), txt)
  })
  
  # DDA - dimensionality plot: PCA, tSNE (scRNA-seq only), 
  # uMAP (scRNA-seq only)
  output$dimred <- renderPlotly({
    req(dat, input$dimredfact, input$dimredmethod)
    isolate({
      if(is.list(dat)) {
        if(!(input$dimredfact %in% colnames(dat$metadata))) return()
        if("DESeqTransform" %in% class(dat$data)) {
          sourceName <- "A"
          if(!is.null(datasetLabel)) sourceName <- datasetLabel
          PCABulk(dat$data, fact = input$dimredfact, sourceName = sourceName)
        } else {
          PCARASL(dat$data, meta = dat$metadata, fact = input$dimredfact)
        }
      } else {
        if(!(input$dimredfact %in% colnames(dat@meta.data))) return()
        if(input$dimredmethod == "pca") DimRedSeurat = PCASeurat
        if(input$dimredmethod == "tsne") DimRedSeurat = TSNESeurat
        if(input$dimredmethod == "umap") DimRedSeurat = UMAPSeurat
        DimRedSeurat(dat, fact = input$dimredfact)
      }
    })
  })
  
  output$dimredPlot_ui <- renderUI({
    req(dat, input$dimredfact, input$dimredmethod, plotWidth, plotHeight)
    plotlyOutput(session$ns("dimred"), width = plotWidth, height = plotHeight)
  })
  
  # Returns PCA data so we can use plotly event_data to get sample ID
  PCAData <- reactive({
    req(input$dimredfact)
    if(is.null(returnData) || !returnData || is.null(dimredMethods) ||
       !("pca" %in% dimredMethods) || is.null(dat) || !is.list(dat) || 
       is.null(dat$data)) {
      return()
    }
    plotPCA(dat$data, intgroup = input$dimredfact, returnData = T)
  })
  
  # Observer for pca plot click
  observe({
    req(datasetLabel)
    d$plotClick <- event_data("plotly_click", source = datasetLabel)
  })
  
  observeEvent(d$plotClick, {
    plotData <- PCAData()
    plotData <- PCAData()[PCAData()$group == levels(PCAData()$group)[d$plotClick$curveNumber + 1], ]
    d$selectedSample <- rownames(plotData)[d$plotClick$pointNumber + 1]
  })
  
  
  
  # DDA - download button dim reduction plot (PDF)
  output$dlqcdimredpdf <- renderUI({
    downloadButton(session$ns("dlqcdimredpdfimg"), label = "Download static plot (PDF)")
  })
  
  # DDA - download file dim reduction plot (PDF)
  output$dlqcdimredpdfimg <- downloadHandler(
    filename =  function() {
      if(input$dimredmethod == "pca") return("qc-pca.pdf")
      if(input$dimredmethod == "tsne") return("qc-tsne.pdf")
      if(input$dimredmethod == "umap") return("qc-umap.pdf")
    },
    content = function(file) {
      pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
      if(is.list(dat)) {
        if("DESeqTransform" %in% class(dat$data)) {
          qcPCAPlotBulk(
            tmp = dat$data,
            fact = input$dimredfact
          )
        } else {
          qcPCAPlotRASL(
            mat = dat$data,
            meta = dat$metadata,
            fact = input$dimredfact
          )
        }
      } else {
        if(input$dimredmethod == "pca") DimRedFunc = qcPCAPlotSeurat
        if(input$dimredmethod == "tsne") DimRedFunc = qcTSNEPlotSeurat
        if(input$dimredmethod == "umap") DimRedFunc = qcUMAPPlotSeurat
        DimRedFunc(dat, fact = input$dimredfact)
      }
      dev.off()
    }
  )
  
  # DDA - download button dim reduction plot (PNG)
  output$dlqcdimredpng <- renderUI({
    downloadButton(session$ns("dlqcdimredpngimg"), "Download static plot (PNG)")
  })
  
  # DDA - download file dim reduction plot (PNG)
  output$dlqcdimredpngimg <- downloadHandler(
    filename =  function() {
      if(input$dimredmethod == "pca") return("qc-pca.png")
      if(input$dimredmethod == "tsne") return("qc-tsne.png")
      if(input$dimredmethod == "umap") return("qc-umap.png")
    },
    content = function(file) {
      png(file, width = 800, height = 750)
      if(is.list(dat)) {
        if("DESeqTransform" %in% class(dat$data)) {
          qcPCAPlotBulk(
            tmp = dat$data,
            fact = input$dimredfact
          )
        } else {
          qcPCAPlotRASL(
            mat = dat$data,
            meta = dat$metadata,
            fact = input$dimredfact
          )
        }
      } else {
        if(input$dimredmethod == "pca") DimRedFunc = qcPCAPlotSeurat
        if(input$dimredmethod == "tsne") DimRedFunc = qcTSNEPlotSeurat
        if(input$dimredmethod == "umap") DimRedFunc = qcUMAPPlotSeurat
        DimRedFunc(dat, fact = input$dimredfact)
      }
      dev.off()
    }
  )
  
  if(returnData) return(d)
}

#################################################################################################

## Choose samples module

choosesamplesUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidPage(
      div(style = "display: inline-block;", uiOutput(ns("close"))),
      div(style = "display: inline-block;", uiOutput(ns("showFactorChoice"))),
      div(style = "display: inline-block;", uiOutput(ns("showPCAOutliers"))),
      div(style = "display: inline-block;", uiOutput(ns("submit"))),
      div(
        style = "display: inline-block; float: right;",
        actionLink(ns("choosesamples_help_button"), label = "Help")
      )
    ),
    br(),
    fluidPage(uiOutput(ns("nSubsetSamples"))),
    # br(),
    tabsetPanel(
      id = ns("subsetData_tabsetPanel"),
      type = "hidden",
      selected = "factorChoice",
      tabPanelBody(
        value = "factorChoice",
        uiOutput(ns("select_samples_factorchoice"))
      ),
      tabPanelBody(
        value = "pcaOutliers",
        fluidPage(
          hr(),
          h3("Identify outliers by PCA", style = "margin-top: 0px;"),
          div(
            div(
              style = "display: inline-block;",
              actionButton(
                inputId = ns("runPCA"), 
                label = "Run PCA",
                icon = icon("space-shuttle")
              ),
            ),
            div(
              style = "display: inline-block; margin-left: 30px;",
              uiOutput(ns("submitOutliers_button"))
            )
          ),
          br(), br(),
          uiOutput(ns("pcaPlots"))
        )
      ),
      tabPanelBody(
        value = "subset_help",
        fluidPage(
          div(
            style = "border: 1px solid #eee; padding-left: 10px;",
            includeMarkdown("markdown/help/subset_help.md")
          )
        )
      )
    )
  )
}

# Takes dataframe and column names and returns a (possibly smaller) dataframe
# with only the unique combinations of levels of each corresponding column
UniqueDF <- function(df, colNames = colnames(df)) {
  if(length(colNames) == 1) {
    newDF <- data.frame(sort(unique(df[, colNames])))
    colnames(newDF) <- colNames
    return(newDF)
  }
  concatVals <- apply(df[, colNames, drop = F], 1, function(myRow) {
    paste(myRow, collapse = "_AND_")
  })
  concatVals <- unique(concatVals)
  newMat <- matrix(nrow = length(concatVals), ncol = length(colNames))
  colnames(newMat) <- colNames
  newMat <- t(sapply(concatVals, function(concatVal) strsplit(concatVal, split = "_AND_")[[1]]))
  colnames(newMat) <- colNames
  rownames(newMat) <- NULL
  newDF <- as.data.frame(newMat)
  return(newDF)
}

# Returns names of list 'listX' for which there either isn't a corresponding name in 
# list 'listY' or the contents are different
ListDiff <- function(listX, listY) {
  if(isTRUE(all.equal(target = listX, current = listY))) return()
  extraNames <- setdiff(x = names(listX), y = names(listY))
  commonNames <- intersect(x = names(listX), y = names(listY))
  if(length(commonNames) == 0) return(names(listX))
  commonNamesNotEqual <- sapply(commonNames, function(commonName) {
    !isTRUE(all.equal(target = listX[[commonName]], current = listY[[commonName]]))
  })
  commonNamesNotEqual <- commonNames[commonNamesNotEqual]
  diffNames <- c(extraNames, commonNamesNotEqual)
  if(length(diffNames) == 0) return()
  return(diffNames)
}

choosesamples <- function(input, 
                          output, 
                          session, 
                          selectSamplesTableVarsList,
                          allowMulti = F,
                          usePCA = F,
                          data_type = NULL) {
  observe({
    if(!is.null(usePCA) && usePCA && !is.null(d$pcaNames) &&
       (is.null(input$runPCA) || input$runPCA == 0)) {
      d$pcaNames <- NULL
      d$pcaOutput <- list()
      updateTabsetPanel(inputId = "subsetData_tabsetPanel", selected = "factorChoice")
    }
  })
  
  d <- reactiveValues(
    inputFactors = NULL,
    inputFactorsTmp = NULL,
    selectSamplesTableList = NULL,
    submit = 0,
    close = 0,
    datasetsRowsSelected = NULL,
    datasetsRowsSelectedTmp = NULL,
    pcaNames = NULL,
    pcaOutput = list(),
    updateTables = T
  )
  
  observeEvent(input$submit, {
    # Create list with datasets as names which store row numbers of selected 
    # sample groups as vectors 
    datasetsRowsSelected <- list()
    for(dataset in names(d$selectSamplesTableList)) {
      if(length(input[[paste0(dataset, "_rows_selected")]]) >= 1) {
        datasetsRowsSelected[[dataset]] <- input[[paste0(dataset, "_rows_selected")]]
      }
    }
    d$datasetsRowsSelected <- datasetsRowsSelected
    d$submit <- d$submit + 1
  }, ignoreInit = T)
  
  observeEvent(input$close, {
    d$close <- d$close + 1
    d$pcaNames <- NULL
    d$pcaOutput <- list()
  }, ignoreInit = T)
  
  observeEvent(input$choosesamples_help_button, {
    updateTabsetPanel(inputId = "subsetData_tabsetPanel", selected = "subset_help")
  })
  
  observe({
    req(selectSamplesTableVarsList)
    myLabel <- "Select factor"
    if(allowMulti) myLabel <- "Select factor(s)"
    for(i in 1:length(selectSamplesTableVarsList)) {
      local({
        dataset <- names(selectSamplesTableVarsList)[i]
        factorData <- selectSamplesTableVarsList[[dataset]]
        myFactor <- colnames(factorData)[1]
        samples <- mixedsort(unique(factorData[, myFactor]))
        dfTmp <- as.data.frame(matrix(nrow = length(samples), ncol = 3))
        colnames(dfTmp) <- c("Dataset", myFactor, "Dataset2")
        dfTmp[1,1] <- dataset
        dfTmp[, 3] <- dataset
        dfTmp[, 2] <- samples
        isolate({
          d$selectSamplesTableList[[dataset]] = dfTmp
        })
        colnames(dfTmp) <- c("Experiment", myFactor, "Dataset2")
        output[[dataset]] <- DT::renderDataTable({
          DT::datatable(
            data.frame(Checkbox = "", dfTmp[, c("Experiment", myFactor)]), 
            colnames = c("", "Experiment", myFactor),
            rownames = F, escape = F,
            selection = "none",
            extensions = c("Select", "Buttons"),
            options = list(
              lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
              scrollX = 700,
              ordering = F, 
              dom = "Bt",
              select = list(style = "multi"),
              buttons = c("selectAll", "selectNone"),
              columnDefs = list(list(className = "select-checkbox", targets = 0))
            )
          )
        }, server = F)
        output[[paste0(dataset, "_factor")]] <- renderUI({
          selectInput(session$ns(paste0(dataset, "_factor")),
                      label = myLabel,
                      choices = mixedsort(colnames(factorData)),
                      selected = myFactor,
                      multiple = allowMulti)
        })
      })
    }    
  }, autoDestroy = T)
  
  # Submit button for select samples. Appears once samples are selected
  output$submit <- renderUI({
    if(is.null(SubsetNSamples()) || sum(SubsetNSamples()) < 6) return()
    div(
      style = "margin-left: 20px;",
      actionButton(session$ns("submit"), label = "Submit")
    )
  })
  
  
  # Close button for select samples
  output$close <- renderUI({
    actionButton(session$ns("close"), label = "Close")
  })
  
  # Switch to factor choice tab
  output$showFactorChoice <- renderUI({
    req(input$subsetData_tabsetPanel, usePCA)
    if(input$subsetData_tabsetPanel == "factorChoice") return()
    div(
      style = "margin-left: 20px;",
      actionButton(
        inputId = session$ns("showFactorChoice"),
        label = "Back",
        icon = icon("arrow-left")
      )
    )
  })
  
  observeEvent(input$showFactorChoice, {
    updateTabsetPanel(inputId = "subsetData_tabsetPanel", selected = "factorChoice")
  })
  
  # Switch to PCA outliers tab
  output$showPCAOutliers <- renderUI({
    req(input$subsetData_tabsetPanel, usePCA)
    if(input$subsetData_tabsetPanel == "pcaOutliers") return()
    div(
      style = "margin-left: 20px;",
      actionButton(
        inputId = session$ns("showPCAOutliers"),
        label = "Identify outliers"
      )
    )
  })
  
  observeEvent(input$showPCAOutliers, {
    updateTabsetPanel(inputId = "subsetData_tabsetPanel", selected = "pcaOutliers")
  })
  
  # DATA - creates a ui element for factor and levels (DT)
  # of that factor for chosen samples from given dataset(s)
  output$select_samples_factorchoice = renderUI({
    req(selectSamplesTableVarsList)
    lapply(names(selectSamplesTableVarsList), function(i) {
      if(is.null(selectSamplesTableVarsList[[i]])) return(NULL)
      fluidPage(
        div(
          style = "border: 1px solid #eee; width: 740px; padding: 20px;",
          div(
            div(style = "display: inline-block;", h4(i)),
            div(style = "display: inline-block; float: right; width: 250px;", uiOutput(session$ns(paste0(i, "_factor"))))
          ),
          div(
            div(style = "display: inline-block; width: 700px;", DT::dataTableOutput(session$ns(i)))
          )
        ),
        br()
      )
    })
  })
  
  # DATA - specify d$inputFactors for chosen samples for selected dataset(s)
  observe({
    req(selectSamplesTableVarsList)
    datasets <- names(selectSamplesTableVarsList)
    if(sum(names(input) %in% paste0(datasets, "_factor")) == 0) {
      d$inputFactors <- d$inputFactorsTmp <- NULL
      return()
    }
    inputFactors <- lapply(datasets, function(dataset) input[[paste0(dataset, "_factor")]])
    names(inputFactors) <- datasets
    
    diffNames <- ListDiff(listX = inputFactors, listY = d$inputFactors)
    d$inputFactors <- inputFactors
    if(length(diffNames) > 0) {
      d$inputFactorsTmp <- inputFactors[diffNames]
    }
  })
  
  # DATA - select factor for selected dataset(s) under choose samples
  observeEvent(d$inputFactorsTmp, {
    if(is.null(d$inputFactorsTmp)) return()
    if(!d$updateTables) {
      d$updateTables <- T
      return()
    }
    datasets <- names(d$inputFactorsTmp)
    for(dataset in datasets) {
      local({
        factorData <- selectSamplesTableVarsList[[dataset]]
        myFactors <- d$inputFactorsTmp[[dataset]]
        if(is.null(myFactors)) {
          dfTmp <- data.frame(Experiment = c(dataset, "Please select at least one factor."))
          d$selectSamplesTableList[[dataset]] <- NULL
          output[[dataset]] <- DT::renderDataTable({
            DT::datatable(
              dfTmp,
              rownames = F, escape = F, selection = "none",
              options = list(lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
                             ordering = F)
            )
          })
        } else {
          samplesDF <- UniqueDF(df = factorData, colNames = myFactors)
          samplesDF <- samplesDF[do.call(order, samplesDF), , drop = F]
          dfTmp <- as.data.frame(matrix(nrow = nrow(samplesDF), ncol = length(myFactors) + 2))
          colnames(dfTmp) <- c("Dataset", myFactors, "Dataset2")
          dfTmp[1,1] <- dataset
          dfTmp[, ncol(dfTmp)] <- dataset
          dfTmp[, 2:(ncol(dfTmp)-1)] <- samplesDF
          d$selectSamplesTableList[[dataset]] <- dfTmp
          colnames(dfTmp) <- c("Experiment", myFactors, "Dataset2")
          output[[dataset]] <- DT::renderDataTable({
            DT::datatable(
              data.frame(Checkbox = "", dfTmp[, c("Experiment", myFactors)]), 
              colnames = c("", "Experiment", myFactors),
              extensions = c("Select", "Buttons"),
              rownames = F, escape = F,
              selection = "none",
              options = list(
                lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
                scrollX = 700,
                ordering = F, 
                dom = "Bt",
                select = list(style = "multi"),
                buttons = c("selectAll", "selectNone"),
                columnDefs = list(list(className = "select-checkbox", targets = 0))
              )
            )
          }, server = F)
        }
      })
    }
  }, ignoreNULL = F, ignoreInit = T)
  
  observeEvent(input$runPCA, {
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                      dbname = bdb, host = ec_host, port = p)
    datasets <- names(selectSamplesTableVarsList)
    names(datasets) <- datasets
    datasets <- datasets[paste0(datasets, "_meta") %in% dbListTables(mydb)]
    dbDisconnect(mydb)
    for(i in 1:length(datasets)) {
      local({
        dataset <- datasets[i]
        dds <- RawToDDS(dataset)
        meta <- as.data.frame(colData(dds))
        myTran <- normTransform(dds)
        d$pcaOutput[[dataset]] <- callModule(
          module = dimred, 
          id = paste0(dataset, "_pca"), 
          dat = list(data = myTran, metadata = meta),
          plotWidth = 650, plotHeight = 425,
          datasetLabel = dataset,
          returnData = T
        ) 
      })
    }
    d$pcaNames <- datasets
  })
  
  output$pcaPlots <- renderUI({
    req(input$runPCA, d$pcaNames)
    lapply(d$pcaNames, function(dataset) {
      div(
        div(
          style = "border: 1px solid #eee; padding: 20px;",
          dimredUI(session$ns(paste0(dataset, "_pca")))
        ),
        br()
      )
    })
    
  })
  
  output$submitOutliers_button <- renderUI({
    if(is.null(d$pcaOutput) || length(d$pcaOutput) == 0) return()
    sampleList <- lapply(d$pcaOutput, function(dataset) dataset[["selectedSamples"]])
    if(length(unlist(sampleList)) == 0) return()
    actionButton(
      inputId = session$ns("submitOutliers_button"),
      label = "Submit outliers"
    )
  })
  
  observeEvent(input$submitOutliers_button, {
    if(is.null(d$pcaOutput) || length(d$pcaOutput) == 0) return()
    sampleList <- lapply(d$pcaOutput, function(dataset) dataset[["selectedSamples"]])
    if(length(unlist(sampleList)) == 0) return()
    
    updateTabsetPanel(inputId = "subsetData_tabsetPanel", selected = "factorChoice")
    datasets <- names(d$pcaOutput)
    for(i in 1:length(datasets)) {
      local({
        dataset <- datasets[i]
        if(!is.null(d$pcaOutput[[dataset]]) && !is.null(d$pcaOutput[[dataset]]$selectedSamples)) {
          factorData <- selectSamplesTableVarsList[[dataset]]
          factorChoices <- colnames(factorData)
          inputID <- paste0(dataset, "_factor")
          if(!is.null(input[[inputID]]) && "sample_id" %in% factorChoices) {
            d$updateTables <- F
            updateSelectInput(inputId = inputID, selected = "sample_id")
            factorData <- factorData[order(factorData$sample_id), , drop = F]
            dfTmp <- as.data.frame(matrix(nrow = nrow(factorData), ncol = 3))
            colnames(dfTmp) <- c("Dataset", "sample_id", "Dataset2")
            dfTmp[1,1] <- dataset
            dfTmp[, ncol(dfTmp)] <- dataset
            dfTmp[, 2] <- factorData$sample_id
            d$selectSamplesTableList[[dataset]] <- dfTmp
            colnames(dfTmp) <- c("Experiment", "sample_id", "Dataset2")
            rowsToKeep <- which(!(dfTmp$sample_id %in% d$pcaOutput[[dataset]]$selectedSamples))
            output[[dataset]] <- DT::renderDataTable({
              DT::datatable(
                data.frame(Checkbox = "", dfTmp[, c("Experiment", "sample_id")]),
                colnames = c("", "Experiment", "sample_id"),
                extensions = c("Select", "Buttons"),
                rownames = F, escape = F,
                selection = list(mode = "multiple", selected = rowsToKeep),
                options = list(
                  lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
                  scrollX = 700,
                  ordering = F, 
                  dom = "Bt",
                  select = list(style = "multi"),
                  buttons = c("selectAll", "selectNone"),
                  columnDefs = list(list(className = "select-checkbox", targets = 0))
                )
              )
            }, server = F)
          }
        }
      })
    }
  })
  
  # Returns vector with subset experiment names as names and numbers of samples 
  # from each subset dataset as values
  SubsetNSamples <- reactive({
    req(selectSamplesTableVarsList, d$selectSamplesTableList, d$inputFactors)
    expNames <- names(selectSamplesTableVarsList)
    nSamples <- vector(mode = "integer", length = length(expNames))
    names(nSamples) <- expNames
    for(expName in expNames) {
      myRows <- input[[paste0(expName, "_rows_selected")]]
      if(length(myRows) == 0) next
      subDF <- d$selectSamplesTableList[[expName]][myRows, d$inputFactors[[expName]], drop = F]
      rowsIncluded <- SubsetMulti(
        fullDF = selectSamplesTableVarsList[[expName]], 
        subDF = subDF
      )
      nSamples[expName] <- length(rowsIncluded)
    }
    return(nSamples)
  })
  
  # Show number of samples or cells currently selected in subset data
  output$nSubsetSamples <- renderUI({
    req(data_type)
    if(is.null(SubsetNSamples()) || sum(SubsetNSamples()) == 0) {
      return(
        tagList(
          br()
        )
      )
    }
    nSamples <- sum(SubsetNSamples())
    
    sampleTypeMulti <- "cells"
    sampleTypeSingle <- "cell"
    if(data_type == "Bulk") {
      sampleTypeMulti <- "samples"
      sampleTypeSingle <- "sample"
    }

    txt <- paste0(nSamples, " ", sampleTypeMulti, " selected.")
    if(nSamples == 1) {
      txt <- paste0(nSamples, " ", sampleTypeSingle, " selected. Please select at least 6 ", sampleTypeMulti, ".")
    } else if(nSamples < 6) {
      txt <- paste0(nSamples, " ", sampleTypeMulti, " selected. Please select at least 6 ", sampleTypeMulti, ".")
    }
    tagList(
      # br(),
      HTML(txt)
    )
  })
  
  return(d)
}

# Load data module
LoadDataUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      div(
        style = "display: inline-block; vertical-align: top; float: right;",
        actionLink(ns("load_data_help_button"), label = "Help")
      ),
      tabsetPanel(
        id = "load_data_tabset",
        tabPanel(
          title = "Existing datasets",
          fluidPage(
            fluidRow(
              style = "margin-top: 10px;",
              radioGroupButtons(
                inputId = ns("data_type"),
                choices = c(
                  "Single-cell RNA-seq" = "Single-cell", 
                  "Bulk RNA-seq" = "Bulk"
                )                  
              )
            ),
            fluidRow(
              column(
                width = 10,
                fluidRow(h4("Select dataset(s) from existing experiment")),
                fluidRow(DT::dataTableOutput(ns("select_datasets_table"))),
                fluidRow(uiOutput(ns("clearExisting_existing"))),
                fluidRow(textOutput(ns("end_warn")), style = "color:red"),
                br(),
                fluidRow(uiOutput(ns("user_made_meta_ui")))
              )
            )
          )
        ),
        tabPanel(
          title = "Custom dataset",
          fluidPage(
            fluidRow(
              style = "margin-top: 10px;",
              radioGroupButtons(
                inputId = ns("data_type_custom"),
                choices = c(
                  "Single-cell RNA-seq" = "Single-cell", 
                  "Bulk RNA-seq" = "Bulk"
                )                  
              )
            ),
            fluidRow(
              textInput(
                inputId = ns("caption"), 
                label = h4("Name uploaded dataset"), 
                value = NULL
              )
            ),            
            br(),
            fluidRow(
              column(
                width = 5,
                style = "background-color: #eee;",
                fluidRow(
                  column(width = 12, align = "right", actionLink(ns("upload_cts_help_button"), label = "Help"))
                ),
                fluidRow(
                  column(width = 10, uiOutput(ns("end_count"))),
                  column(
                    width = 1,
                    style = "padding-top: 30px; margin-left: -10px;",
                    actionLink(
                      inputId = ns("end_count_delete_button"),
                      label = "",
                      icon = icon("trash")
                    )
                  )
                )
              ),
              column(
                width = 5,
                style = "background-color: #eee; margin-left: 20px;",
                fluidRow(
                  column(width = 12, align = "right", actionLink(ns("upload_meta_help_button"), label = "Help"))
                ),
                fluidRow(
                  column(width = 10, uiOutput(ns("end_meta"))),
                  column(
                    width = 1,
                    style = "padding-top: 30px; margin-left: -10px;",
                    actionLink(
                      inputId = ns("end_meta_delete_button"),
                      label = "",
                      icon = icon("trash")
                    )
                  )
                )
              )
            ),
            fluidRow(uiOutput(ns("addClearExisting_ui"))),
          )
        )
      )
    ),
    fluidRow(hr()),
    fluidRow(uiOutput(ns("choose_samples_ui"))),
    fluidRow(uiOutput(ns("filtAndTransformOpts_ui"))),
    fluidRow(
      div(
        style = "display: inline-block; vertical-align: bottom;",
        uiOutput(ns("res"))
      ),
      div(
        style = "display: inline-block; vertical-align: bottom; width: 80px; margin-left: 10px;", 
        uiOutput(ns("resMin"))
      ),
      div(
        style = "display: inline-block; vertical-align: bottom; width: 80px; margin-left: 10px;", 
        uiOutput(ns("resMax"))
      ),
      div(
        style = paste0(
          "display: inline-block; vertical-align: bottom; width: 80px;",
          " margin-left: 10px;"
        ), 
        uiOutput(ns("resStep"))
      ),
      div(
        style = "display: inline-block; margin-bottom: 15px; margin-left: 10px; margin-right: 10px;",
        uiOutput(ns("resDefault")),
      ),
      div(
        style = paste0(
          "display: inline-block; vertical-align: bottom; margin-bottom: 25px;",
          " font-size: 16px;"
        ),
        uiOutput(ns("resValues"))
      )
    ),
    fluidRow(uiOutput(ns("filt_genes"))),
    fluidRow(uiOutput(ns("downsample_select"))),
    fluidRow(uiOutput(ns("downsample_options_ui"))),
    fluidRow(uiOutput(ns("showAdvancedOptions"))),
    fluidRow(uiOutput(ns("advanced_options_ui"))),
    fluidRow(uiOutput(ns("sc_largedata_msg"))),
    fluidRow(uiOutput(ns("goqc_warning_msg"))),
    fluidRow(
      div(style = "display: inline-block;", uiOutput(ns("load_data"))),
      div(style = "display: inline-block; margin-left: 20px;", uiOutput(ns("save_load_selected")))
    ),
    br()
  )
}

LoadData <- function(input, output, session, maxSamples = 10000) {
  
  # Set input$data_type and input$data_type_custom to be equal
  observeEvent(input$data_type, {
    updateRadioGroupButtons(
      session = session, 
      inputId = "data_type_custom", 
      selected = input$data_type
    )
  })
  
  observeEvent(input$data_type_custom, {
    updateRadioGroupButtons(
      session = session, 
      inputId = "data_type", 
      selected = input$data_type_custom
    )
  })
  
  
  # Reactive expression to load list of bulk/scRNA-seq datasets on app load
  ExpData <- reactive({
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                      dbname = scdb, host = ec_host, port = p)
    sc <- dbReadTable(mydb, "sc") %>%
      filter(paste0(unique_table, "_counts") %in% dbListTables(mydb)) %>%
      filter(paste0(unique_table, "_meta") %in% dbListTables(mydb))
    sc <- sc[c(which(sc$experiment_name == "example_sc"), which(sc$experiment_name != "example_sc")), ]
    
    names(sc)[1] <- "experiment_id"
    sc$Type <- "sc"
    # Count rows in metadata
    nSamples <- sapply(sc$unique_table, function(expName) DBTableNRows(mydb, myTab = paste0(expName, "_meta")))
    sc <- data.frame(sc, NSamples = nSamples)
    dbDisconnect(mydb)
    
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                      dbname = bdb, host = ec_host, port = p)
    b <- dbReadTable(mydb, "bulk") %>%
      filter(paste0(unique_table, "_counts") %in% dbListTables(mydb)) %>%
      filter(paste0(unique_table, "_meta") %in% dbListTables(mydb))
    b <- b[c(which(b$experiment_name == "example_bulk"), which(b$experiment_name != "example_bulk")), ]
    
    b$Type <- "bulk"
    # Count rows in metadata
    nSamples <- sapply(b$unique_table, function(expName) DBTableNRows(mydb, myTab = paste0(expName, "_meta")))
    b <- data.frame(b, NSamples = nSamples)
    
    dbDisconnect(mydb)
    
    df <- rbind(b, sc)
    
    # Add publication link as "Source" column
    df$Source <- NA
    df$Source[!is.na(df$publication_link)] <- 
      paste0(
        "<a href='", 
        df$publication_link[!is.na(df$publication_link)], 
        "' target='_blank'>", 
        df$publication[!is.na(df$publication_link)],
        "</a>"
      )

    return(df)
  })
  
  d <- reactiveValues(
    selectedDatasetsTable = NULL,
    selectedData = NULL,
    selectSamplesTableList = NULL,
    datasetsRowsSelected = NULL,
    cts_out = NULL,
    selectSamplesTable = NULL,
    selectSamplesTableVarsList = NULL,
    selectSamplesTableSelectedVar = NULL,
    datasetNames = NULL,
    inputFactors = NULL,
    allowPrecomputedClusters = F,
    goqc_warning = NULL,
    newExpSetup = T,
    drt = 0.1,
    pcaVarianceFraction = 75,
    data_type = NULL,
    ddsone = NULL,
    ddstran = NULL,
    resType = NULL,
    goqc = 0,
    select_datasets_table_rows_selected = NULL,
    select_datasets_table_addExisting_rows_selected = NULL,
    newDataLoaded = 0,
    end_count = NULL,
    end_meta = NULL,
    end_ready = F,
    end_partial = F,
    datasets_table = NULL,
    count_content = NULL,
    meta_content = NULL,
    resValues = seq(from = 0.4, to = 2.8, by = 0.4),
    showAdvancedOptions = F,
    tooFewSamplesMessage = "Minimum 6 cells required. Please include more cells."
  )
  
  # Observes whether user-uploaded data is ready and whether it is 
  # partially ready
  observe({
    sumInputs <- sum(
      !is.null(count_content()), 
      !is.null(meta_content()), 
      !(is.null(input$caption) || input$caption == "")
    )
    d$end_ready <- (sumInputs == 3)
    d$end_partial <- (sumInputs > 0 & sumInputs < 3)
  })
  
  observeEvent(input$load_data_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Loading data")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html#selection-of-data",
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/load_data_help.md")
    ))
  })
  
  # DATA - reactive expression to subset dataset into bulk or scRNA-Seq
  # depending on what the user selects
  load_datasets_table <- reactive({
    req(input$data_type, ExpData())
    df <- ExpData()
    if(input$data_type == "Single-cell") {
      df <- df[df$Type == "sc", ]
    } else if(input$data_type == "Bulk") {
      df <- df[df$Type == "bulk", ]
    } else {
      df <- df[df$Type == "rasl", ]
    }
    if(nrow(df) == 0) return()
    return(df)
  })
  
  observe({
    d$datasets_table <- load_datasets_table()
  })
  
  # DATA - reactive expression to pull bulk and scRNA-seq
  # data set from mysql db
  data_info <- reactive({
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                      dbname = scdb, host = ec_host, port = p)
    sc <- dbReadTable(mydb, "sc") %>%
      filter(paste0(unique_table, "_counts") %in% dbListTables(mydb)) %>%
      filter(paste0(unique_table, "_meta") %in% dbListTables(mydb))
    
    names(sc)[1] <- "experiment_id"
    sc$Type <- "sc"
    nSamples <- sapply(sc$unique_table, function(expName) DBTableNRows(mydb, myTab = paste0(expName, "_meta")))
    sc <- data.frame(sc, NSamples = nSamples)
    dbDisconnect(mydb)
    
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                      dbname = bdb, host = ec_host, port = p)
    b <- dbReadTable(mydb, "bulk") %>%
      filter(paste0(unique_table, "_counts") %in% dbListTables(mydb)) %>%
      filter(paste0(unique_table, "_meta") %in% dbListTables(mydb))
    
    b$Type <- "bulk"
    nSamples <- sapply(b$unique_table, function(expName) DBTableNRows(mydb, myTab = paste0(expName, "_meta")))
    b <- data.frame(b, NSamples = nSamples)
    dbDisconnect(mydb)
    
    df <- rbind(b, sc)
    return(df)
  })
  
  output$submit_addExisting <- renderUI({
    req(input$select_datasets_table_addExisting_rows_selected)
    actionButton(session$ns("submit_addExisting"), label = "Submit")
  })
  
  observeEvent(input$combine_existing_button, {
    showModal(modalDialog(
      title = "Add existing experiments",
      size = "l",
      footer = NULL,
      fluidPage(
        fluidRow(
          div(style = "display: inline-block;", actionButton(session$ns("close_addExisting"), label = "Close")),
          div(style = "display: inline-block;", uiOutput(session$ns("submit_addExisting")))
        ),
        br(),
        DT::dataTableOutput(session$ns("select_datasets_table_addExisting"))
      )
    ))
  })
  
  output$addExisting <- renderUI({
    if(is.null(load_datasets_table()) || nrow(load_datasets_table()) == 0) return()
    actionLink(
      inputId = session$ns("combine_existing_button"), 
      label = "+ Add existing experiments"
    )
  })
  
  # Button to clear selection from popup datasets table
  output$clearExisting <- renderUI({
    if((is.null(d$select_datasets_table_addExisting_rows_selected) || 
        length(d$select_datasets_table_addExisting_rows_selected) == 0) &&
       (is.null(d$select_datasets_table_rows_selected) ||
        length(d$select_datasets_table_rows_selected) == 0)) return()
    
    actionLink(
      inputId = session$ns("clearExisting"),
      label = "- Clear existing experiments"
    )
  })
  
  output$addClearExisting_ui <- renderUI({
    req(input$data_type)
    if(input$data_type == "RASL") return()
    tagList(
      br(),
      div(
        style = "display: inline-block;",
        uiOutput(session$ns("addExisting")),
      ),
      div(
        style = "display: inline-block; margin-left: 30px;",
        uiOutput(session$ns("clearExisting"))
      )
    )
  })
  
  # DATA - this is a DT for selecting bulk or sc data
  output$select_datasets_table <- DT::renderDataTable({
    req(input$data_type)
    
    df <- load_datasets_table()
    if(input$data_type == "Single-cell") {
      colNames <- c("", "Experiment", "Cells", "Description", "Source")
    } else {
      colNames <- c("", "Experiment", "Samples", "Description", "Source")
    }

    if(!is.null(df) && nrow(df) > 0) {
      df <- data.frame(Checkbox = "", df)
    } else {
      df <- as.data.frame(matrix(nrow = 0, ncol = 5))
      colnames(df) <- c("Checkbox", "experiment_name", "NSamples", "description", "Source")
    }
    
    DT::datatable(
      df[, c("Checkbox", "experiment_name", "NSamples", "description", "Source")],
      colnames = colNames, 
      rownames = F, escape = F,
      extensions = "Select",
      selection = "none",
      options = list(
        lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
        ordering = F, scrollY = 200, scrollCollapse = T,
        columnDefs = list(
          list(className = "select-checkbox", targets = 0),
          list(width = "500px", targets = 3),
          list(width = "200px", targets = 4)
        ),
        select = list(style = "multi"),
        language = list(zeroRecords = "No experiments available")
      )
    )
  }, server = F)

  
  # Allow clearing datasets selection
  output$clearExisting_existing <- renderUI({
    if((is.null(d$select_datasets_table_addExisting_rows_selected) || 
       length(d$select_datasets_table_addExisting_rows_selected) == 0) &&
       (is.null(d$select_datasets_table_rows_selected) ||
        length(d$select_datasets_table_rows_selected) == 0)) return()

    tagList(
      br(),
      actionLink(inputId = session$ns("clearExisting_existing"), label = "Clear selections")
    )
  })
  
  observeEvent(input$clearExisting_existing, {
    d$select_datasets_table_addExisting_rows_selected <- NULL
    d$select_datasets_table_rows_selected <- NULL
    output$select_datasets_table <- DT::renderDataTable({
      req(input$data_type, load_datasets_table())
      df <- load_datasets_table()
      if(input$data_type == "Single-cell") {
        colNames <- c("", "Experiment", "Cells", "Description", "Source")
      } else {
        colNames <- c("", "Experiment", "Samples", "Description", "Source")
      }
      df <- data.frame(Checkbox = "", df)
      DT::datatable(
        df[, c("Checkbox", "experiment_name", "NSamples", "description", "Source")],
        colnames = colNames,
        rownames = F, escape = F,
        extensions = "Select",
        selection = "none",
        options = list(
          lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
          ordering = F, scrollY = 200, scrollCollapse = T,
          columnDefs = list(list(className = "select-checkbox", targets = 0)),
          select = list(style = "multi")
        )
      )
    }, server = F)
  })
  
  # DATA - this is a DT for selecting bulk or sc data
  output$select_datasets_table_addExisting <- DT::renderDataTable({
    req(input$data_type, load_datasets_table())
    df <- load_datasets_table()
    if(input$data_type == "Single-cell") {
      colNames <- c("", "Experiment", "Cells", "Description", "Source")
    } else {
      colNames <- c("", "Experiment", "Samples", "Description", "Source")
    }
    df <- data.frame(Checkbox = "", df)
    DT::datatable(
      df[, c("Checkbox", "experiment_name", "NSamples", "description", "Source")],
      colnames = colNames, 
      rownames = F, escape = F,
      extensions = "Select",
      selection = "none",
      options = list(
        lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
        ordering = F, scrollY = 200, scrollCollapse = T,
        columnDefs = list(list(className = "select-checkbox", targets = 0)),
        select = list(style = "multi")
      )
    )
  }, server = F)
  
  # DATA - reactive expression to return user-created metadata table
  user_made_meta_table <- reactive({
    req(load_datasets_table(), d$select_datasets_table_rows_selected)
    
    if(length(d$select_datasets_table_rows_selected) == 0) return()
    
    selectedDatasetsTable <- load_datasets_table()[d$select_datasets_table_rows_selected, ]
    names(selectedDatasetsTable)[1] <- gsub("^X\\.", "", names(selectedDatasetsTable)[1])
    selectedDatasetsTable <- selectedDatasetsTable %>%
      dplyr::select(unique_table)
    
    if(nrow(selectedDatasetsTable) > 1) {
      selectedDatasetsTable$unique_table <- paste(selectedDatasetsTable$unique_table, collapse = "_")
    }
    
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                      dbname = scdb, host = ec_host, port = p)
    
    sc <- dbReadTable(conn = mydb, name = "sc_useradd") %>%
      filter(unique_table == selectedDatasetsTable$unique_table) %>%
      filter(table_name %in% dbListTables(mydb))
    
    dbDisconnect(mydb)
    
    if(nrow(sc) == 0) return()
    
    sc <- sc %>%
      dplyr::select(c(unique_table, username, uploaded_date, table_name, note))
    
    return(sc)
  })
  
  # DATA - actual user-created DT metadata table
  output$user_made_meta <- DT::renderDataTable({
    req(user_made_meta_table())
    if(nrow(user_made_meta_table()) == 0) return()
    df <- user_made_meta_table()[, c("unique_table", "username", "uploaded_date", "table_name", "note")]
    df <- data.frame(Checkbox = "", df) 
    DT::datatable(
      df, 
      rownames = F,
      extensions = "Select",
      selection = "none",
      colnames = c("", "Experiment", "Username", "Date", "Table name", "Note"),
      escape = F, 
      options = list(
        lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
        ordering = F, scrollX = T,  scrollY = 240, scrollCollapse = T,
        columnDefs = list(list(className = "select-checkbox", targets = 0)),
        select = list(style = "single")
      )
    )
  }, server = F)
  
  # User renderUI to avoid blank space if header and table aren't loaded
  output$user_made_meta_ui <- renderUI({
    req(user_made_meta_table())
    if(nrow(user_made_meta_table()) == 0) return()
    tagList(
      h4("Option: Select user-updated metadata"),
      DT::dataTableOutput(session$ns("user_made_meta")),
      br()
    )
  })
  
  # DATA - action button to load sample selection
  output$choose_samples <- renderUI({
    req(input$data_type)
    if(input$data_type == "RASL") return()
    actionButton(session$ns("choose_samples"), label = "Subset data")
  })
  output$choose_samples_ui <- renderUI({
    req(input$data_type)
    if(input$data_type == "RASL") return()
    if (is.null(d$select_datasets_table_rows_selected) & !d$end_ready) return()
    if(input$data_type == "Single-cell" & !is.null(input$user_made_meta_rows_selected)) return()
    tagList(
      uiOutput(session$ns("choose_samples"))
    )
  })
  
  # DATA - warning message for user-selected count/metadata
  output$end_warn <- renderText({
    req(input$choose_samples, meta_content(), count_content())
    m <- meta_content()
    c <- count_content()
    
    txt <- NULL
    if(is.null(input$caption)) {
      if(!all(colnames(c) %in% rownames(m))) {
        if(is.null(rownames(c))) {
          txt <- paste("Please name your uploaded data.",
                       "The column names in the count matrix are not the same as the row names in the metadata.",
                       "The row names names in the count matrix are null. They need to contain gene symbols.", sep = "\n")
        } else {
          txt <- paste("Please name your uploaded data.",
                       "The column names in the count matrix are not the same as the row names in the metadata.", sep = "\n")
        }
      } else {
        txt <- paste("Please name your uploaded data.",
                     "The row names names in the count matrix are null. They need to contain gene symbols.", sep = "\n")
      }
    } else if(!all(colnames(c) %in% rownames(m))) {
      if(is.null(rownames(c))) {
        txt <- paste("The column names in the count matrix are not the same as the row names in the metadata.",
                     "The row names names in the count matrix are null. They need to contain gene symbols.", sep = "\n")
      } else {
        txt <- "The column names in the count matrix are not the same as the row names in the metadata."
      }
    } else if(is.null(rownames(c))) {
      txt <- "The row names names in the count matrix are null. They need to contain gene symbols."
    } else if(input$caption == "") {
      txt <- "Please name your uploaded data."
    }
    if(is.null(txt)) return()
    HTML(txt)
  })
  
  # DATA - set d$selectSamplesTableList back to NULL when data
  # type changes & add a tab for ipsc profiler for scRNA-Seq
  observeEvent(input$data_type, {
    d$selectSamplesTableList <- d$datasetsRowsSelected <- d$cts_out <- 
      d$selectedData <- d$end_count <- d$end_meta <- 
      d$select_datasets_table_rows_selected <- d$selectSamplesTableVarsList <- NULL
    d$data_type <- input$data_type
    updateTextInput(
      session = session,
      inputId = "caption",
      value = ""
    )
  })
  
  # Store selected rows from popup when user clicks Submit
  observeEvent(input$submit_addExisting, {
    d$select_datasets_table_addExisting_rows_selected <- input$select_datasets_table_addExisting_rows_selected
    removeModal()
  })
  
  # Clear selected rows from popup when user clicks Close
  observeEvent(input$close_addExisting, {
    d$select_datasets_table_addExisting_rows_selected <- NULL
    removeModal()
  })
  
  # Clear selected rows from popup when user clicks Clear existing. 
  observeEvent(input$clearExisting, {
    d$select_datasets_table_addExisting_rows_selected <- NULL
    d$select_datasets_table_rows_selected <- NULL
    output$select_datasets_table <- DT::renderDataTable({
      req(input$data_type, load_datasets_table())
      df <- load_datasets_table()
      if(input$data_type == "Single-cell") {
        colNames <- c("", "Experiment", "Cells", "Description", "Source")
      } else {
        colNames <- c("", "Experiment", "Samples", "Description", "Source")
      }
      df <- data.frame(Checkbox = "", df)
      DT::datatable(
        df[, c("Checkbox", "experiment_name", "NSamples", "description", "Source")],
        colnames = colNames,
        rownames = F, escape = F,
        extensions = "Select",
        selection = "none",
        options = list(
          lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
          ordering = F, scrollY = 200, scrollCollapse = T,
          columnDefs = list(list(className = "select-checkbox", targets = 0)),
          select = list(style = "multi")
        )
      )
    }, server = F)
  })
  
  # Combine selected rows from popup table with any selected in Existing datasets table
  observe({
    d$select_datasets_table_rows_selected <- unique(
      c(input$select_datasets_table_rows_selected, d$select_datasets_table_addExisting_rows_selected)
    )
  })
  
  # Observer to reset when select datasets table is changed
  observeEvent(d$select_datasets_table_rows_selected, {
    d$selectSamplesTableList <- d$datasetsRowsSelected <- d$cts_out <- 
      d$selectedData <- NULL
    if(is.null(count_content()) & is.null(meta_content())) {
      updateTextInput(
        session = session,
        inputId = "caption",
        value = ""
      )
    }
  }, ignoreNULL = F, ignoreInit = T)
  
  # DATA - after selecting choose samples, pull factors for each
  # selected dataset
  observeEvent(input$choose_samples, {
    
    d$selectedDatasetsTable <- NULL
    d$selectSamplesTable <- NULL
    d$selectSamplesTableList <- NULL
    d$selectSamplesTableVarsList <- NULL
    d$selectSamplesTableSelectedVar <- NULL
    d$selectedData <- NULL
    d$datasetsRowsSelected <- NULL
    d$datasetNames <- NULL
    
    data <- list()
    excludeCols <- c("V1", "cell_id", "orig.ident", "nFeature_RNA", "nCount_RNA")
    if(d$end_ready) {
      meta <- meta_content()
      factorData <- meta[, !names(meta) %in% excludeCols, drop = F]
      data[["meta"]][[input$caption]] <- factorData
      data[["factor"]][[input$caption]] <- colnames(factorData)[1]
    }
    
    if(!is.null(d$select_datasets_table_rows_selected)) {
      selectedDatasetsTable = load_datasets_table()[d$select_datasets_table_rows_selected, ]
      if(input$data_type == "Single-cell") {
        mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                          dbname = scdb, host = ec_host, port = p)
        tablesDF <- dbReadTable(mydb, "sc")
      } else {
        mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                          dbname = bdb, host = ec_host, port = p)
        tablesDF <- dbReadTable(mydb, "bulk")
      }
      for(i in 1:nrow(selectedDatasetsTable)) {
        exp <- selectedDatasetsTable[i,] %>% dplyr::pull(experiment_name)
        dat <- tablesDF %>%
          dplyr::filter(experiment_name %in% exp) %>%
          dplyr::pull(unique_table)
        meta <- dbReadTable(mydb, name = paste(dat, "_meta", sep = ""))
        # ignore the following columns
        factorData <- meta[, !names(meta) %in% excludeCols, drop = F]
        data[["meta"]][[selectedDatasetsTable$experiment_name[i]]] <- factorData
        data[["factor"]][[selectedDatasetsTable$experiment_name[i]]] <- colnames(factorData)[1]
      }
      dbDisconnect(mydb)
    }
    
    d$selectSamplesTableVarsList <- data[["meta"]]
    d$selectSamplesTableSelectedVar <- data[["factor"]]
    
  }, priority = 1)
  
  # Returns list from reactiveValues "d" object from choosesamples module.
  # Contains inputFactors, selectSamplesTableList, submit, close, and
  # datasetsRowsSelected.
  ChooseSamplesData <- reactive({
    req(input$choose_samples, input$data_type)
    if(input$data_type == "RASL") return()
    if(is.null(d$selectSamplesTableVarsList)) return()
    if(input$data_type == "Bulk") {
      callModule(
        choosesamples,
        id = "choose_samples_modal_bulk",
        selectSamplesTableVarsList = d$selectSamplesTableVarsList,
        allowMulti = T,
        usePCA = T,
        data_type = "Bulk"
      )
    } else {
      callModule(
        choosesamples,
        id = "choose_samples_modal_sc",
        selectSamplesTableVarsList = d$selectSamplesTableVarsList,
        allowMulti = F,
        usePCA = F,
        data_type = "Single-cell"
      )
    }
  })
  
  # When user clicks "Submit" inside Choose samples modal, return
  # "selectedSamplesTableList", "datasetsRowsSelected", and "inputFactors",
  # and close modal dialog
  observeEvent(ChooseSamplesData()$submit, {
    if(is.null(ChooseSamplesData()$submit)) return()
    if(ChooseSamplesData()$submit == 0) return()
    d$selectSamplesTableList <- ChooseSamplesData()$selectSamplesTableList
    d$datasetsRowsSelected <- ChooseSamplesData()$datasetsRowsSelected
    d$inputFactors <- ChooseSamplesData()$inputFactors
    removeModal()
  })
  
  # Set d objects returned by ChooseSamplesData() to NULL and
  # close modal dialog if user clicks "Close"
  observeEvent(ChooseSamplesData()$close, {
    if(is.null(ChooseSamplesData()$close) || ChooseSamplesData()$close == 0) {
      d$selectSamplesTableList <- d$datasetsRowsSelected <- d$inputFactors <- 
        d$inputFactorsTmp <- NULL
      return()
    }
    removeModal()
  }, ignoreNULL = T)
  
  observeEvent(input$choose_samples, {
    d$selectSamplesTableList <- d$datasetsRowsSelected <- d$inputFactors <- d$inputFactorsTmp <- NULL
    id <- "choose_samples_modal_bulk"
    if(input$data_type == "Single-cell") id <- "choose_samples_modal_sc"
    showModal(modalDialog(
      title = "Subset data",
      footer = NULL,
      size = "l",
      choosesamplesUI(session$ns(id))
    ))
  }, priority = 0)
  
  # DATA - action button for loading selected samples
  output$save_load_selected <- renderUI({
    req(d$selectSamplesTableList, d$datasetsRowsSelected, d$inputFactors)
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    actionButton(session$ns("save_load_selected"), label = "Load selected sample(s)")
  })
  
  # DATA - help for uploading data
  observeEvent(input$upload_cts_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Upload count data")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/upload_count_data_help.md")
    ))
  })
  
  # DATA - end-user to provide count matrix
  output$end_count <- renderUI({
    req(input$data_type)
    input$end_count_delete_button
    input$data_type
    myLabel <- "Upload a count matrix "
    if(input$data_type == "RASL") myLabel <- "Upload gene-level data "
    fileInput(
      inputId = session$ns("end_count"), 
      label = div(
        h4(myLabel, style = "display: inline;"),
        h5("(txt/csv)\n", style = "display: inline;")
      ),
      width = "600px",
      multiple = FALSE,
      accept = c(
        "text/csv",
        "text/comma-separated-values,text/plain",
        ".csv")
    )
  })
  
  count_content <- reactive({
    if(is.null(d$end_count)) return()
    fread(d$end_count$datapath, header = T) %>%
      mutate_all(funs(replace_na(.,0))) %>% 
      column_to_rownames(names(.)[1])
  })
  observe({
    d$count_content <- count_content()
  })
  
  
  # DATA - help for uploading data
  observeEvent(input$upload_meta_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Upload metadata")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/upload_metadata_help.md")
    ))
  })
  
  # DATA - end-user to provide metadata matrix
  output$end_meta <- renderUI({
    input$end_meta_delete_button
    input$data_type
    myLabel <- "Upload metadata "
    if(input$data_type == "RASL") myLabel <- "Upload treatment-level data "
    fileInput(
      inputId = session$ns("end_meta"),
      label = div(
        h4(myLabel, style = "display: inline;"),
        h5("(txt/csv)\n", style = "display: inline;"),
      ),
      width = "600px",
      multiple = FALSE,
      accept = c(
        "text/csv",
        "text/comma-separated-values,text/plain",
        ".csv")
    )
  })
  
  meta_content <- reactive({
    if(is.null(d$end_meta)) return()
    fread(d$end_meta$datapath, header = T) %>%
      column_to_rownames(names(.)[1])
  })
  observe({
    d$meta_content <- meta_content()
  })
  
  # Use d object for end_count and end_meta so we can remove uploaded data
  observeEvent(input$end_count, {
    d$end_count <- input$end_count
  })
  
  observeEvent(input$end_meta, {
    d$end_meta <- input$end_meta
  })
  
  observeEvent(input$end_count_delete_button, {
    d$end_count <- NULL
  })
  
  observeEvent(input$end_meta_delete_button, {
    d$end_meta <- NULL
  })
  
  # DATA - checkbox choice of whether to downsample
  output$downsample_select <- renderUI({
    req(input$data_type, SelectedDataNCells())
    if(input$data_type %in% c("Bulk", "RASL")) return()
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    
    if(!is.null(input$user_made_meta_rows_selected)) return()
    defaultValue <- T
    if(max(SelectedDataNCells()) < 10000) defaultValue <- F
    
    prettyCheckbox(
      inputId = session$ns("downsample_select"), 
      label = "Downsample cells",
      value = defaultValue, 
      status = "default", 
      icon = icon("check")
    )
  })
  
  # DATA - reactive expression to return a vector of numbers
  # of cells in datasets for down sampling
  SelectedDataNCells <- reactive({
    req(input$data_type, load_datasets_table())
    if(input$data_type != "Single-cell") return()
    if (is.null(d$select_datasets_table_rows_selected)) {
      if (!is.null(count_content()) & is.null(meta_content())) return()
      else if (is.null(count_content()) & !is.null(meta_content())) return()
      else if (!is.null(count_content()) & !is.null(meta_content())) {
        meta <- meta_content()
        nCells <- nrow(meta)
        return(nCells)
      }
    } else {
      expNames <- load_datasets_table()$unique_table[d$select_datasets_table_rows_selected]
      mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                        dbname = scdb, host = ec_host, port = p)
      nCells <- sapply(expNames, function(expName) DBTableNRows(mydb, myTab = paste0(expName, "_meta")))
      dbDisconnect(mydb)
      if (!is.null(count_content()) & is.null(meta_content())) return()
      else if (is.null(count_content()) & !is.null(meta_content())) return()
      else if (!is.null(count_content()) & !is.null(meta_content())) {
        meta <- meta_content()
        nCells2 <- nrow(meta)
        nCells <- c(nCells, nCells2)
        return(nCells)
      }
      return(nCells)
    }
  })
  
  # DATA - select total number of cells to down sample
  output$downsample_nCells <- renderUI({
    req(input$data_type, input$downsample_select, SelectedDataNCells())
    
    if(input$data_type %in% c("Bulk", "RASL")) return()
    if(!is.null(input$user_made_meta_rows_selected)) return()
    nCellsMax <- max(SelectedDataNCells())
    nCellsDefault <- min(nCellsMax, 10000)
    
    numericInput(
      inputId = session$ns("downsample_nCells"), 
      label = "Max cells", 
      value = nCellsDefault, 
      min = 30, 
      max = nCellsMax, 
      step = 1
    )
  })
  
  # DATA - selected sample seed to start down sampling at
  output$downsample_seed <- renderUI({
    req(input$data_type, input$downsample_select, SelectedDataNCells())
    if(input$data_type %in% c("Bulk", "RASL")) return()
    if(!is.null(input$user_made_meta_rows_selected)) return()
    
    numericInput(
      inputId = session$ns("downsample_seed"), 
      label = "Random seed", 
      value = 10, 
      step = 1
    )
  })
  
  output$downsample_options_ui <- renderUI({
    req(input$data_type, input$downsample_select, SelectedDataNCells())
    if(input$data_type != "Single-cell") return()
    if(!is.null(input$user_made_meta_rows_selected)) return()
    
    tagList(
      div(
        style = "display: inline-block; width: 150px;", 
        uiOutput(session$ns("downsample_nCells"))
      ),
      div(
        style = "display: inline-block; width: 150px; margin-left: 10px;", 
        uiOutput(session$ns("downsample_seed"))
      ),
    )
  })
  
  # Action link to show or hide advanced options
  output$showAdvancedOptions <- renderUI({
    req(input$data_type)
    if(input$data_type != "Single-cell") return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    tagList(
      actionLink(inputId = session$ns("showAdvancedOptions"), label = "Advanced options"),
      br(), br()
    )
  })
  
  # Toggle showing advanced options when "Advanced options" action link is clicked
  observeEvent(input$showAdvancedOptions, {
    d$showAdvancedOptions <- ifelse(test = d$showAdvancedOptions, yes = F, no = T)
  })
  
  observe({
    if(d$showAdvancedOptions) {
      shinyjs::show("advanced_options_ui")
    } else {
      shinyjs::hide("advanced_options_ui")
    }
  })
  
  # Single-cell advanced options
  output$advanced_options_ui <- renderUI({
    req(input$data_type)
    if(input$data_type != "Single-cell" || is.null(input$res) ||
       is.null(d$end_ready) || is.null(d$end_partial)) return()
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    myStyle <- "border: 1px solid #eee; width: 220px; padding: 10px;"
    if(input$res == "seurat_res") {
      myStyle <- "border: 1px solid #eee; width: 440px; padding: 10px;"
    }

    tagList(
      div(
        style = myStyle,
        div(
          style = "display: inline-block; width: 180px; vertical-align: top;",
          uiOutput(session$ns("drt"))
        ),
        div(
          style = "display: inline-block;",
          uiOutput(session$ns("pcaVarianceFraction"))
        )
      ),
      br()
    )
  })
  
  #DATA - header warning for large datasets if user doesn't select downsampling
  output$sc_largedata_msg <- renderUI({
    if(is.null(input$data_type)) return()
    if(is.null(input$downsample_select)) return()
    if(!input$downsample_select) return()
    if(is.null(SelectedDataNCells())) return()
    if(input$data_type %in% c("Bulk", "RASL")) return()
    if(!is.null(input$user_made_meta_rows_selected)) return()
    
    if(max(SelectedDataNCells()) <= maxSamples) return()
    if(!is.null(input$downsample_nCells)) {
      if(input$downsample_select && is.finite(input$downsample_nCells) && 
         input$downsample_nCells <= maxSamples) return()
    }
    maxCellsFormatted <- formatC(maxSamples, big.mark = ",")
    txt <- paste0("Dataset(s) will be downsampled to a maximum of ", maxCellsFormatted, " cells.")
    h4(txt, style = "color: red;")
  })
  
  # DATA - warning if data can't be processed
  output$goqc_warning_msg <- renderUI({
    if(is.null(d$goqc_warning)) return()
    tagList(
      div(style = "color: red;", d$goqc_warning),
      br()
    )
  })
  
  # DATA - action button to load sample data for selected datasets
  output$load_data <- renderUI({
    req(input$data_type)
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    myLabel <- "Load full dataset(s)"
    if(input$data_type == "RASL") myLabel <- "Load data"
    actionButton(session$ns("load_data"), label = myLabel)
  })
  
  # DATA - warning message for user-selected count/metadata for load full datasets
  output$end_warn_loadFull <- renderText({
    req(input$load_data, meta_content(), count_content())
    m <- meta_content()
    c <- count_content()
    
    txt <- NULL
    if (is.null(input$caption)) {
      if (!all(colnames(c) %in% rownames(m))) {
        if (is.null(rownames(c))) {
          txt <- paste("Please name your uploaded data.",
                       "The column names in the count matrix are not the same as the row names in the metadata.",
                       "The row names names in the count matrix are null. They need to contain gene symbols.", sep = "\n")
        } else {
          txt <- paste("Please name your uploaded data.",
                       "The column names in the count matrix are not the same as the row names in the metadata.", sep = "\n")
        }
      } else {
        txt <- paste("Please name your uploaded data.",
                     "The row names names in the count matrix are null. They need to contain gene symbols.", sep = "\n")
      }
    } else if (!all(colnames(c) %in% rownames(m))) {
      if (is.null(rownames(c))) {
        txt <- paste("The column names in the count matrix are not the same as the row names in the metadata.",
                     "The row names names in the count matrix are null. They need to contain gene symbols.", sep = "\n")
      } else {
        txt <- "The column names in the count matrix are not the same as the row names in the metadata."
      }
    } else if (is.null(rownames(c))) {
      txt <- "The row names names in the count matrix are null. They need to contain gene symbols."
    } else if (input$caption == "") {
      txt <- "Please name your uploaded data."
    }
    if(is.null(txt)) return()
    HTML(txt)
  })
  
  # DATA - create a new merged dataset when multiple datasets
  # and samples in each dataset are selected
  observeEvent(input$save_load_selected, {
    loadMessage <- "Loading data, this may take a few seconds..."
    if(!is.null(SelectedDataNCells()) && sum(SelectedDataNCells()) > 1000) {
      loadMessage <- "Loading data, this may take a few minutes..."
    }
    withProgress(message = loadMessage, value = 0, {
      incProgress(1/3)
      
      if(is.null(d$select_datasets_table_rows_selected) & input$data_type == "Single-cell") {
        if(!is.null(count_content()) & is.null(meta_content())) return()
        if(is.null(count_content()) & !is.null(meta_content())) return()
        if(!is.null(count_content()) & !is.null(meta_content())) {
          datasetList <- list()
          cts <- count_content()
          meta <- meta_content()
          
          # Downsample if selected or if number of cells is greater than 10000
          if(input$downsample_select || nrow(meta) > 10000) {
            nCells <- min(input$downsample_nCells, 10000)
            set.seed(seed = input$downsample_seed)
            myCells <- sample(x = rownames(meta), size = nCells)
            cts <- cts[, myCells, drop = F]
            meta <- meta[myCells, , drop = F]
          }
          
          nam <- input$caption
          
          # Subset cells
          if(!(nam %in% names(d$selectSamplesTableList))) return()
          if(!(nam %in% names(d$datasetsRowsSelected))) return()
          if(!(nam %in% names(d$inputFactors))) return()
          selectedSamplesTable <- d$selectSamplesTableList[[nam]][d$datasetsRowsSelected[[nam]], ]
          myFactors <- d$inputFactors[[nam]]
          samplesDF <- selectedSamplesTable[, myFactors, drop = F]
          rowsToKeep <- SubsetMulti(fullDF = meta, subDF = samplesDF)
          meta <- meta[rowsToKeep, , drop = F]
          cts <- cts[, rownames(meta)]
          
          upl_seurat <- CreateSeuratObject(cts, assay = "RNA", meta.data = meta)
          
          datasetList[[nam]] <- upl_seurat
          seuratMerged <- datasetList[[1]]
          dataMerged <- list(counts = t(FetchData(datasetList[[1]], vars = rownames(datasetList[[1]]),
                                                  slot = "counts")),
                             metadata = datasetList[[1]]@meta.data)
        }
      } else if (is.null(d$select_datasets_table_rows_selected) & input$data_type == "Bulk") {
        if (!is.null(count_content()) & is.null(meta_content())) return()
        if (is.null(count_content()) & !is.null(meta_content())) return()
        if (!is.null(count_content()) & !is.null(meta_content())) {
          datasetList <- list()
          cts <- count_content()
          meta <- meta_content()
          nam <- input$caption
          
          # Subset samples
          if(!(nam %in% names(d$selectSamplesTableList))) return()
          if(!(nam %in% names(d$datasetsRowsSelected))) return()
          if(!(nam %in% names(d$inputFactors))) return()
          selectedSamplesTable <- d$selectSamplesTableList[[nam]][d$datasetsRowsSelected[[nam]], ]
          myFactors <- d$inputFactors[[nam]]
          samplesDF <- selectedSamplesTable[, myFactors, drop = F]
          rowsToKeep <- SubsetMulti(fullDF = meta, subDF = samplesDF)
          meta <- meta[rowsToKeep, , drop = F]
          cts <- cts[, rownames(meta)]
          
          datasetList[[nam]][["counts"]] <- cts
          datasetList[[nam]][["metadata"]] <- meta
          
          dataMerged <- datasetList[[1]]
        }
      } else {
        datasets <- names(d$datasetsRowsSelected)
        selectedDatasets <- c()
        for(dataset in datasets) {
          if(length(d$datasetsRowsSelected[[dataset]]) >= 1) {
            selectedDatasets <- c(selectedDatasets, dataset)
          }
        }
        if(length(selectedDatasets) == 0) return()
        
        selectedDatasetsTable <- load_datasets_table()[d$select_datasets_table_rows_selected, ]
        
        datasetList <- list()
        if(input$data_type == "Single-cell") {
          maxCells <- mySeed <- NULL
          if(input$downsample_select) {
            maxCells <- input$downsample_nCells
            mySeed <- input$downsample_seed
          }
          for(dataset in names(d$selectSamplesTableList)) {
            if(length(d$datasetsRowsSelected[[dataset]]) < 1) next
            selectedSamplesTable <- d$selectSamplesTableList[[dataset]][d$datasetsRowsSelected[[dataset]], ]
            myFactors <- d$inputFactors[[dataset]]
            samplesDF <- selectedSamplesTable[, myFactors, drop = F]
            
            if(!is.null(input$caption)) {
              if(dataset == input$caption) {
                if (!is.null(count_content()) & !is.null(meta_content())) {
                  cts <- count_content()
                  meta <- meta_content()
                  
                  # Downsample if selected or if number of cells is greater than 10000
                  if(input$downsample_select || nrow(meta) > 10000) {
                    nCells <- min(input$downsample_nCells, 10000)
                    set.seed(seed = input$downsample_seed)
                    myCells <- sample(x = rownames(meta), size = nCells)
                    cts <- cts[, myCells, drop = F]
                    meta <- meta[myCells, , drop = F]
                  }
                  
                  # Subset cells
                  if(!(dataset %in% names(d$selectSamplesTableList))) next
                  if(!(dataset %in% names(d$datasetsRowsSelected))) next
                  if(!(dataset %in% names(d$inputFactors))) next
                  selectedSamplesTable <- d$selectSamplesTableList[[dataset]][d$datasetsRowsSelected[[dataset]], ]
                  myFactors <- d$inputFactors[[dataset]]
                  samplesDF <- selectedSamplesTable[, myFactors, drop = F]
                  rowsToKeep <- SubsetMulti(fullDF = meta, subDF = samplesDF)
                  meta <- meta[rowsToKeep, , drop = F]
                  cts <- cts[, rownames(meta)]
                  
                  upl_seurat <- CreateSeuratObject(cts, assay = "RNA", meta.data = meta)
                  datasetList[[dataset]] <- upl_seurat
                  
                  # substitute sample column for orig.ident if orig.ident has 1 level
                  # and sample has more
                  if(length(unique(datasetList[[dataset]]@meta.data$orig.ident)) == 1) {
                    if(length(unique(datasetList[[dataset]]@meta.data$sample)) > 1) {
                      datasetList[[dataset]]@meta.data$orig.ident <- datasetList[[dataset]]@meta.data$sample
                    }
                  }
                  
                  next
                }
              }
            }
            
            datasetList[[dataset]] <-
              SeuratProcess(dataset,
                            samples = samplesDF,
                            nCells = maxCells,
                            mySeed = mySeed)
            
            # substitute sample column for orig.ident if orig.ident has 1 level
            # and sample has more
            if(length(unique(datasetList[[dataset]]@meta.data$orig.ident)) == 1) {
              if(length(unique(datasetList[[dataset]]@meta.data$sample)) > 1) {
                datasetList[[dataset]]@meta.data$orig.ident <- datasetList[[dataset]]@meta.data$sample
              }
            }
          }
          
          if(length(datasetList) > 1) {
            for(i in 1:length(datasetList)) {
              dataset <- names(datasetList)[i]
              myFactors <- d$inputFactors[[dataset]]
              datasetList[[dataset]] <- RenameCells(datasetList[[dataset]],
                                                    new.names = paste0(colnames(datasetList[[dataset]]), ".", dataset))
              for(colname in colnames(datasetList[[dataset]]@meta.data)[!(colnames(datasetList[[dataset]]@meta.data) %in% c("nCount_RNA","nFeature_RNA"))]) {
                datasetList[[dataset]]@meta.data[, colname] <- paste0(datasetList[[dataset]]@meta.data[, colname], ".", dataset)
              }
              datasetList[[dataset]]@meta.data$Dataset <- dataset
            }
            seuratMerged <- datasetList[[1]]
            dataMerged <- list(counts = t(FetchData(datasetList[[1]], vars = rownames(datasetList[[1]]),
                                                    slot = "counts")),
                               metadata = datasetList[[1]]@meta.data)
            for(i in 2:length(datasetList)) {
              dataset <- names(datasetList)[i]
              first <- dataMerged[["counts"]]
              nxt <- list(counts = t(FetchData(datasetList[[i]], vars = rownames(datasetList[[i]]),
                                               slot = "counts")),
                          metadata = datasetList[[i]]@meta.data)
              nxt <- nxt[["counts"]]
              genes <- rownames(first)[rownames(first) %in% rownames(nxt)]
              first <- first %>%
                data.frame(.) %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              nxt <- nxt %>%
                data.frame(.) %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              dataMerged[["counts"]] <- first %>%
                left_join(., nxt, by = "gene") %>%
                column_to_rownames("gene")
              colnamesMergedOnly <- colnames(dataMerged[["metadata"]])[!(colnames(dataMerged[["metadata"]]) %in% colnames(datasetList[[i]]@meta.data))]
              if(length(colnamesMergedOnly) > 0) {
                for(colname in colnamesMergedOnly) datasetList[[i]]@meta.data[, colname] <- dataset
              }
              colnamesNewOnly <- colnames(datasetList[[i]]@meta.data)[!(colnames(datasetList[[i]]@meta.data) %in% colnames(dataMerged[["metadata"]]))]
              if(length(colnamesNewOnly) > 0) {
                for(colname in colnamesNewOnly) dataMerged[["metadata"]][, colname] <- dataMerged[["metadata"]][, "Dataset"]
              }
              dataMerged[["metadata"]] <- rbind(dataMerged[["metadata"]],
                                                datasetList[[i]]@meta.data[, colnames(dataMerged[["metadata"]])])
              seuratMerged <- merge(seuratMerged, datasetList[[i]])
            }
          } else {
            dataMerged <- list(counts = t(FetchData(datasetList[[1]], vars = rownames(datasetList[[1]]),
                                                    slot = "counts")),
                               metadata = datasetList[[1]]@meta.data)
            seuratMerged <- datasetList[[1]]
          }
        } else {
          selectedDatasetsTable <- load_datasets_table()[d$select_datasets_table_rows_selected, ]
          for(dataset in names(d$selectSamplesTableList)) {
            if(length(d$datasetsRowsSelected[[dataset]]) < 1) next
            selectedSamplesTable <- d$selectSamplesTableList[[dataset]][d$datasetsRowsSelected[[dataset]], ]
            
            if(!is.null(input$caption)) {
              if(dataset == input$caption) {
                if (!is.null(count_content()) & !is.null(meta_content())) {
                  cts <- count_content()
                  meta <- meta_content()
                  
                  # Subset samples
                  if(!(dataset %in% names(d$selectSamplesTableList))) next
                  if(!(dataset %in% names(d$datasetsRowsSelected))) next
                  if(!(dataset %in% names(d$inputFactors))) next
                  selectedSamplesTable <- d$selectSamplesTableList[[dataset]][d$datasetsRowsSelected[[dataset]], ]
                  myFactors <- d$inputFactors[[dataset]]
                  samplesDF <- selectedSamplesTable[, myFactors, drop = F]
                  rowsToKeep <- SubsetMulti(fullDF = meta, subDF = samplesDF)
                  meta <- meta[rowsToKeep, , drop = F]
                  cts <- cts[, rownames(meta)]
                  
                  # "-" in counts sample names gets converted to "." in data.frame()
                  # which causes errors. Converting "-" in sample names to "_" in
                  # counts and metadata
                  for(colname in colnames(meta)) {
                    if(is.character(meta[, colname]) | is.factor(meta[, colname])) {
                      meta[, colname] <- gsub("-", replacement = "_", x = meta[, colname])
                    }
                  }
                  
                  commonSamples <- intersect(colnames(cts), rownames(meta))
                  cts <- cts[, commonSamples, drop = F]
                  meta <- meta[commonSamples, , drop = F]
                  
                  datasetList[[dataset]][["counts"]] <- cts
                  datasetList[[dataset]][["metadata"]] <- meta
                  next
                }
              }
            } 
            
            # read the table containing info on the scRNA-seq exp
            # in the mysql db
            mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                              dbname = bdb, host = ec_host, port = p)
            dat <- dbReadTable(mydb, "bulk")
            dat <- dat %>%
              dplyr::filter(experiment_name %in% dataset) %>%
              dplyr::pull(unique_table)
            metaDF <- dbReadTable(mydb, name = paste(dat, "_meta", sep = ""))
            countsMatrix <- dbReadTable(mydb, name = paste(dat, "_counts", sep = ""))
            dbDisconnect(mydb)
            
            myFactors <- d$inputFactors[[dataset]]
            if(!is.null(myFactors)) {
              samplesDF <- selectedSamplesTable[, myFactors, drop = F]
              rowsToKeep <- SubsetMulti(fullDF = metaDF, subDF = samplesDF)
              pull_sample_id <- metaDF$sample_id[rowsToKeep]
              countsMatrix <- countsMatrix %>%
                filter(sample_id %in% pull_sample_id)
              metaDF <- metaDF %>%
                filter(sample_id %in% pull_sample_id)
            }
            countsMatrix <- countsMatrix %>%
              spread(sample_id, counts, convert = T) %>%
              column_to_rownames("gene")
            metaDF <- metaDF %>%
              column_to_rownames("sample_id")
            
            # "-" in counts sample names gets converted to "." in data.frame()
            # which causes errors. Converting "-" in sample names to "_" in
            # counts and metadata
            for(colname in colnames(metaDF)) {
              if(is.character(metaDF[, colname]) | is.factor(metaDF[, colname])) {
                metaDF[, colname] <- gsub("-", replacement = "_", x = metaDF[, colname])
              }
            }
            
            commonSamples <- intersect(colnames(countsMatrix), rownames(metaDF))
            countsMatrix <- countsMatrix[, commonSamples, drop = F]
            metaDF <- metaDF[commonSamples, , drop = F]
            
            datasetList[[dataset]][["counts"]] <- countsMatrix
            datasetList[[dataset]][["metadata"]] <- metaDF
          }
          
          if(length(datasetList) > 1) {
            for(i in 1:length(datasetList)) {
              dataset <- names(datasetList)[i]
              colnames(datasetList[[dataset]][["counts"]]) <-
                paste0(colnames(datasetList[[dataset]][["counts"]]), ".", dataset)
              rownames(datasetList[[dataset]][["metadata"]]) <-
                paste0(rownames(datasetList[[dataset]][["metadata"]]), ".", dataset)
              datasetList[[dataset]][["metadata"]]$Dataset <- dataset
            }
            dataMerged <- datasetList[[1]]
            for(i in 2:length(datasetList)) {
              dataset <- names(datasetList)[i]
              first <- dataMerged[["counts"]]
              nxt <- datasetList[[i]][["counts"]]
              genes <- rownames(first)[rownames(first) %in% rownames(nxt)]
              first <- first %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              nxt <- nxt %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              dataMerged[["counts"]] <- first %>%
                left_join(., nxt, by = "gene") %>%
                column_to_rownames("gene")
              colnamesMergedOnly <- colnames(dataMerged[["metadata"]])[!(colnames(dataMerged[["metadata"]]) %in% colnames(datasetList[[i]][["metadata"]]))]
              if(length(colnamesMergedOnly) > 0) {
                for(colname in colnamesMergedOnly) datasetList[[i]][["metadata"]][, colname] <- dataset
              }
              colnamesNewOnly <- colnames(datasetList[[i]][["metadata"]])[!(colnames(datasetList[[i]][["metadata"]]) %in% colnames(dataMerged[["metadata"]]))]
              if(length(colnamesNewOnly) > 0) {
                for(colname in colnamesNewOnly) dataMerged[["metadata"]][, colname] = dataMerged[["metadata"]][, "Dataset"]
              }
              dataMerged[["metadata"]] <- rbind(dataMerged[["metadata"]],
                                                datasetList[[i]][["metadata"]][, colnames(dataMerged[["metadata"]])])
            }
          } else {
            dataMerged <- datasetList[[1]]
          }
        }
      }
      
      incProgress(2/3)
      d$selectedData <- dataMerged
      d$cts_out <- NULL
      d$datasetNames <- names(datasetList)
    })
    d$newDataLoaded <- d$newDataLoaded + 1
  }, priority = 100)
  
  # DATA - radio button to filter MT, pseudogenes, etc
  output$filt_genes <- renderUI({
    req(input$data_type)
    if(input$data_type == "RASL") return()
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    prettyCheckbox(
      inputId = session$ns("filt_genes"), 
      label = "Exclude RNA/MT/pseudo genes",
      value = T, 
      status = "default", 
      icon = icon("check")
    )
    
  })
  
  # DATA - transform data for bulk RNA-Seq
  output$transformSelection <- renderUI({
    req(input$data_type)
    
    if(input$data_type != "Bulk") return()
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    tempChoices = c(
      "Normal log: log2(count + 1)" = "log",
      "Regularized log: rlog(n)" = "rlog",
      "Variance stabilizing transform: vst(n)" = "vst",
      "No transformation" = "raw"
    )
    
    selectInput(
      inputId = session$ns("transform"),
      label = "Count transformation method",
      choices = tempChoices
    )
  })
  
  output$filtAndTransformOpts_ui <- renderUI({
    req(input$data_type)
    if(input$data_type != "Bulk") return()
    tagList(
      br(),
      div(style = "display: inline-block; width: 100px; vertical-align: middle;", uiOutput(session$ns("filter_choice"))),
      div(style = "display: inline-block; margin-left: 100px; vertical-align: middle;", uiOutput(session$ns("transformSelection")))
    )
  })
  
  # DATA - lets user know if batch effect factor only has 1 level
  output$goqc_msg <- renderText({
    req(d$selectedData)
    req(input$data_type, d$select_datasets_table_rows_selected, d$selectedDatasetsTable)
    if(length(d$select_datasets_table_rows_selected) < 1) return()
    if(all(is.null(c(input$load_data, input$save_load_selected)))) return()
    if(all(c(input$load_data, input$save_load_selected) == 0)) return()
    if(input$data_type == "Single-cell") return()
    
    # Exit if either no data selected user-uploaded data incomplete
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    d$goqc_warning
  })
  
  # DATA - filter choice for bulk RNA-seq data
  output$filter_choice <- renderUI({
    req(input$data_type)
    if(input$data_type != "Bulk") return()
    
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    textInput(
      inputId = session$ns("prefilt"),
      label = div(style = "width: 200px;", "Min. total counts per gene"),
      value = as.numeric(10)
    )
  })
  
  # DATA - reactive - load and add data to DESeqDataSet class
  ddsout <- eventReactive(d$newDataLoaded, {
    cts <- as.matrix(d$selectedData[["counts"]])
    coldata <- d$selectedData[["metadata"]]
    
    # remove metadata columns that have only 1 level
    # oneLevelVars <- colnames(coldata)[apply(coldata, 2, function(col) length(unique(col)) == 1)]
    # coldata <- coldata[, !(colnames(coldata) %in% oneLevelVars), drop = F]
    if(input$data_type == "Single-cell") {
      if("nCount_RNA" %in% colnames(coldata)) {
        coldata[, "nCount_RNA"] <- as.numeric(coldata[, "nCount_RNA"])
      }
      if("nFeature_RNA" %in% colnames(coldata)) {
        coldata[, "nFeature_RNA"] <- as.numeric(coldata[, "nFeature_RNA"])
      }
      if(any(grepl("silWidth", x = colnames(coldata)))) {
        for(i in grep("silWidth", x = colnames(coldata))) {
          coldata[, i] = as.numeric(coldata[, i])
        }
      }
    }
    cts <- cts[sort(rownames(cts)), ]
    cols <- colnames(coldata)
    coldata[cols] <- lapply(coldata[cols], function(x) {
      if(is.numeric(x)) {
        as.numeric(x)
      } else {
        factor(x)
      }
    })
    
    cts.filt <- cts
    if(input$data_type != "RASL" && !is.null(input$filt_genes) && input$filt_genes) {
      gen_to_filt <- fread("./data/filter_genes.txt", header = F)
      # Use check.names = FALSE to avoid converting "-" to "."
      cts.filt <- data.frame(cts.filt, check.names = FALSE)
      # ensures coldata rows and count file are in the same order
      cts.filt <- cts.filt[, rownames(coldata)]
      cts.filt <- cts.filt %>%
        mutate(genes = rownames(.)) %>%
        filter(!genes %in% gen_to_filt$V1)
      rownames(cts.filt) <- cts.filt$genes
      cts.filt$genes <- NULL
      cts.filt <- as.matrix(cts.filt)
    }
    if(input$data_type == "Bulk" && !is.null(input$prefilt) && is.finite(input$prefilt)) {
      cts.filt <- cts.filt[rowSums(cts.filt) >= as.numeric(input$prefilt), ]
    }
    dds <- NULL
    if(input$data_type != "RASL") {
      dds <- DESeqDataSetFromMatrix(
        countData = cts.filt,
        colData = coldata,
        design = ~ 1
      )
    }
    # Return objects for downstream analysis
    return(list(dds, coldata, cts.filt, cts))
  })
  
  # If one full dataset is loaded and no user-uploaded data is loaded,
  # allow use of pre-computed clusters from metadata
  observe({
    isolate({
      d$allowPrecomputedClusters <- F
    })
    if(is.null(d$select_datasets_table_rows_selected)) return()
    if(length(d$select_datasets_table_rows_selected) == 0) return()
    selectedDatasetsTable <- load_datasets_table()[d$select_datasets_table_rows_selected, ]
    if(length(d$select_datasets_table_rows_selected) == 1 & is.null(count_content()) &
       is.null(meta_content()) & is.null(d$datasetsRowsSelected)) {
      isolate({
        d$allowPrecomputedClusters <- T
      })
    } else if(length(d$select_datasets_table_rows_selected) == 1 & is.null(count_content()) &
              is.null(meta_content())) {
      datasets <- names(d$datasetsRowsSelected)
      selectedDatasets <- c()
      for(dataset in datasets) {
        if(length(d$datasetsRowsSelected[[dataset]]) >= 1) {
          selectedDatasets <- c(selectedDatasets, dataset)
        }
      }
      if(length(selectedDatasets) == 1) {
        if(length(d$datasetsRowsSelected[[selectedDatasets]]) == nrow(d$selectSamplesTableList[[selectedDatasets]])) {
          isolate({
            d$allowPrecomputedClusters <- T
          })
        }
      }
    }
  })
  
  # DATA - observeEvent for input$load_data. Stores merged data in d$selectedData
  # as list with "metadata" and "counts" dataframes
  observeEvent(input$load_data, {
    
    # Exit if either no data selected or no name provided for user-uploaded data
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    loadMessage <- "Loading data, this may take a few seconds..."
    if(!is.null(SelectedDataNCells()) && sum(SelectedDataNCells()) > 1000) {
      loadMessage <- "Loading data, this may take a few minutes..."
    }
    withProgress(message = loadMessage, value = 0, {
      incProgress(1/3)
      
      if(is.null(d$select_datasets_table_rows_selected) && input$data_type == "Single-cell") {
        if(is.null(count_content()) || is.null(meta_content())) return()
        datasetList <- list()
        cts <- count_content()
        meta <- meta_content()
        
        # Downsample if selected or if number of cells is greater than 10000
        if(input$downsample_select || nrow(meta) > 10000) {
          nCells <- min(input$downsample_nCells, 10000)
          set.seed(seed = input$downsample_seed)
          myCells <- sample(x = rownames(meta), size = nCells)
          cts <- cts[, myCells, drop = F]
          meta <- meta[myCells, , drop = F]
        }
        
        nam <- input$caption
        upl_seurat <- CreateSeuratObject(cts, assay = "RNA", meta.data = meta[colnames(cts), , drop = F])
        
        datasetList[[nam]] <- upl_seurat
        seuratMerged <- datasetList[[1]]
        dataMerged <- list(counts = t(FetchData(datasetList[[1]], vars = rownames(datasetList[[1]]),
                                                slot = "counts")),
                           metadata = datasetList[[1]]@meta.data)
        selectedDatasetsTable <- NULL
      } else if (is.null(d$select_datasets_table_rows_selected) & input$data_type %in% c("Bulk", "RASL")) {
        if(is.null(count_content()) | is.null(meta_content())) return()
        datasetList <- list()
        cts <- count_content()
        meta <- meta_content()
        nam <- input$caption
        datasetList[[nam]][["counts"]] <- cts
        datasetList[[nam]][["metadata"]] <- meta
        
        dataMerged <- datasetList[[1]]
        selectedDatasetsTable <- NULL
      } else {
        # Return selected dataset(s) from "Select data" table into dataframe
        selectedDatasetsTable <- load_datasets_table()[d$select_datasets_table_rows_selected, ]
        
        datasetList = list()
        if(input$data_type == "Single-cell") {
          maxCells <- mySeed <- NULL
          if(input$downsample_select) {
            maxCells <- input$downsample_nCells
            mySeed <- input$downsample_seed
          }
          for(i in 1:nrow(selectedDatasetsTable)) {
            dataset <- selectedDatasetsTable$experiment_name[i]
            datasetList[[dataset]] <-
              SeuratProcess(dataset,
                            samples = NULL,
                            nCells = maxCells,
                            mySeed = mySeed)
            
            # Substitute sample column for orig.ident if orig.ident has 1 level
            # and sample has more
            if(length(unique(datasetList[[dataset]]@meta.data$orig.ident)) == 1) {
              if(length(unique(datasetList[[dataset]]@meta.data$sample)) > 1) {
                datasetList[[dataset]]@meta.data$orig.ident <- datasetList[[dataset]]@meta.data$sample
              }
            }
          }
          
          if (!is.null(count_content()) && !is.null(meta_content())) {
            cts <- count_content()
            meta <- meta_content()
            
            # Downsample if selected or if number of cells is greater than 10000
            if(input$downsample_select || nrow(meta) > 10000) {
              nCells <- min(input$downsample_nCells, 10000)
              set.seed(seed = input$downsample_seed)
              myCells <- sample(x = rownames(meta), size = nCells)
              cts <- cts[, myCells, drop = F]
              meta <- meta[myCells, , drop = F]
            }
            
            nam <- input$caption
            upl_seurat <- CreateSeuratObject(cts, assay = "RNA", meta.data = meta[colnames(cts), , drop = F])
            
            datasetList[[nam]] <- upl_seurat
          }
          
          
          if(length(datasetList) > 1) {
            for(i in 1:length(datasetList)) {
              dataset <- names(datasetList)[i]
              datasetList[[dataset]] <- RenameCells(datasetList[[dataset]],
                                                    new.names = paste0(colnames(datasetList[[dataset]]), ".", dataset))
              for(colname in colnames(datasetList[[dataset]]@meta.data)[!(colnames(datasetList[[dataset]]@meta.data) %in% c("nCount_RNA","nFeature_RNA"))]) {
                datasetList[[dataset]]@meta.data[, colname] <- paste0(datasetList[[dataset]]@meta.data[, colname], ".", dataset)
              }
              datasetList[[dataset]]@meta.data$Dataset <- dataset
            }
            dataMerged <- list(counts = t(FetchData(datasetList[[1]], vars = rownames(datasetList[[1]]),
                                                    slot = "counts")),
                               metadata = datasetList[[1]]@meta.data)
            for(i in 2:length(datasetList)) {
              dataset <- names(datasetList)[i]
              dataset <- names(datasetList)[i]
              first <- dataMerged[["counts"]]
              nxt <- list(counts = t(FetchData(datasetList[[i]], vars = rownames(datasetList[[i]]),
                                               slot = "counts")),
                          metadata = datasetList[[i]]@meta.data)
              nxt <- nxt[["counts"]]
              genes <- rownames(first)[rownames(first) %in% rownames(nxt)]
              first <- first %>%
                data.frame(.) %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              nxt <- nxt %>%
                data.frame(.) %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              
              dataMerged[["counts"]] = first %>%
                left_join(., nxt, by = "gene") %>%
                column_to_rownames("gene")
              colnamesMergedOnly = colnames(dataMerged[["metadata"]])[!(colnames(dataMerged[["metadata"]]) %in% colnames(datasetList[[i]]@meta.data))]
              if(length(colnamesMergedOnly) > 0) {
                for(colname in colnamesMergedOnly) datasetList[[i]]@meta.data[, colname] = dataset
              }
              colnamesNewOnly = colnames(datasetList[[i]]@meta.data)[!(colnames(datasetList[[i]]@meta.data) %in% colnames(dataMerged[["metadata"]]))]
              if(length(colnamesNewOnly) > 0) {
                for(colname in colnamesNewOnly) dataMerged[["metadata"]][, colname] = dataMerged[["metadata"]][, "Dataset"]
              }
              dataMerged[["metadata"]] = rbind(dataMerged[["metadata"]],
                                               datasetList[[i]]@meta.data[, colnames(dataMerged[["metadata"]])])
            }
          } else {
            dataMerged = list(counts = t(FetchData(datasetList[[1]], vars = rownames(datasetList[[1]]),
                                                   slot = "counts")),
                              metadata = datasetList[[1]]@meta.data)
          }
          if(!is.null(user_made_meta_table()) & !is.null(input$user_made_meta_rows_selected)) {
            mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                              dbname = scdb, host = ec_host, port = p)
            tableName <- user_made_meta_table()$table_name[input$user_made_meta_rows_selected]
            
            if(tableName %in% dbListTables(mydb)) {
              userMadeMetaDF <- dbReadTable(mydb, tableName)
              if("sample_DUP" %in% colnames(userMadeMetaDF) && !("sample" %in% colnames(userMadeMetaDF))) {
                colnames(userMadeMetaDF)[colnames(userMadeMetaDF) == "sample_DUP"] = "sample"
              }
              rownames(userMadeMetaDF) <- userMadeMetaDF$Sample
              dataMerged$counts <- dataMerged$counts[, rownames(userMadeMetaDF)]
              dataMerged$metadata <- dataMerged$metadata[rownames(userMadeMetaDF), ]
              if("clus_comb" %in% colnames(userMadeMetaDF)) {
                dataMerged$metadata$seurat_clusters <- userMadeMetaDF$clus_comb
              }
              if("Set" %in% colnames(userMadeMetaDF)) {
                dataMerged$metadata$Set <- userMadeMetaDF$Set
              }
            }
            dbDisconnect(mydb)
          }
          dataMerged$metadata$seurat_clusters <- as.factor(dataMerged$metadata$seurat_clusters)
          
        } else {
          for(i in 1:nrow(selectedDatasetsTable)) {
            dataset <- selectedDatasetsTable$experiment_name[i]
            
            # read the table containing info on the scRNA-seq exp
            # in the mysql db
            mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_bulk, password = pwd_bulk,
                              dbname = bdb, host = ec_host, port = p)
            dat <- dbReadTable(mydb, "bulk")
            dat <- dat %>%
              dplyr::filter(experiment_name %in% dataset) %>%
              dplyr::pull(unique_table)
            metaDF <- dbReadTable(mydb, name = paste(dat, "_meta", sep = ""))
            countsMatrix <- dbReadTable(mydb, name = paste(dat, "_counts", sep = ""))
            dbDisconnect(mydb)
            
            countsMatrix <- countsMatrix %>%
              spread(sample_id, counts, convert = T) %>%
              column_to_rownames("gene")
            
            metaDF <- metaDF %>%
              column_to_rownames("sample_id")
            # "-" in counts sample names gets converted to "." in data.frame()
            # which causes errors. Converting "-" in sample names to "_" in
            # counts and metadata
            for(colname in colnames(metaDF)) {
              if(is.character(metaDF[, colname]) | is.factor(metaDF[, colname])) {
                metaDF[, colname] <- gsub("-", replacement = "_", x = metaDF[, colname])
              }
            }
            commonSamples <- intersect(colnames(countsMatrix), rownames(metaDF))
            countsMatrix <- countsMatrix[, commonSamples, drop = F]
            metaDF <- metaDF[commonSamples, , drop = F]
            
            datasetList[[dataset]][["counts"]] = countsMatrix
            datasetList[[dataset]][["metadata"]] = metaDF
          }
          
          if (!is.null(count_content()) & !is.null(meta_content())) {
            cts <- count_content()
            meta <- meta_content()
            nam <- input$caption
            datasetList[[nam]][["counts"]] <- cts
            datasetList[[nam]][["metadata"]] <- meta
          }
          
          if(length(datasetList) > 1) {
            for(i in 1:length(datasetList)) {
              dataset = names(datasetList)[i]
              colnames(datasetList[[dataset]][["counts"]]) =
                paste0(colnames(datasetList[[dataset]][["counts"]]), ".", dataset)
              rownames(datasetList[[dataset]][["metadata"]]) =
                paste0(rownames(datasetList[[dataset]][["metadata"]]), ".", dataset)
              datasetList[[dataset]][["metadata"]]$Dataset = dataset
            }
            dataMerged = datasetList[[1]]
            for(i in 2:length(datasetList)) {
              dataset = names(datasetList)[i]
              
              first <- dataMerged[["counts"]]
              nxt <- datasetList[[i]][["counts"]]
              genes <- rownames(first)[rownames(first) %in% rownames(nxt)]
              first <- first %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              nxt <- nxt %>%
                rownames_to_column("gene") %>%
                filter(gene %in% genes)
              
              dataMerged[["counts"]] = first %>%
                left_join(., nxt, by = "gene") %>%
                column_to_rownames("gene")
              colnamesMergedOnly = colnames(dataMerged[["metadata"]])[!(colnames(dataMerged[["metadata"]]) %in% colnames(datasetList[[i]][["metadata"]]))]
              if(length(colnamesMergedOnly) > 0) {
                for(colname in colnamesMergedOnly) datasetList[[i]][["metadata"]][, colname] = dataset
              }
              colnamesNewOnly = colnames(datasetList[[i]][["metadata"]])[!(colnames(datasetList[[i]][["metadata"]]) %in% colnames(dataMerged[["metadata"]]))]
              if(length(colnamesNewOnly) > 0) {
                for(colname in colnamesNewOnly) dataMerged[["metadata"]][, colname] = dataMerged[["metadata"]][, "Dataset"]
              }
              dataMerged[["metadata"]] = rbind(dataMerged[["metadata"]],
                                               datasetList[[i]][["metadata"]][, colnames(dataMerged[["metadata"]])])
            }
          } else {
            dataMerged = datasetList[[1]]
          }
          
        }
      }
      incProgress(2/3)
      
      # oneLevelVars <- colnames(dataMerged[["metadata"]])[apply(dataMerged[["metadata"]], 2, function(col) length(unique(col)) == 1)]
      # if(length(oneLevelVars) >= 1) {
      #   dataMerged[["metadata"]] <- dataMerged[["metadata"]][, !(colnames(dataMerged[["metadata"]]) %in% oneLevelVars), drop = F]
      #   d$goqc_warning <- paste("The following metadata variables have only 1 level and have been excluded:",
      #                           paste(oneLevelVars, collapse = ", "))
      # } else {
      #   d$goqc_warning <- NULL
      # }
      d$selectedData <- dataMerged
      d$selectedDatasetsTable <- selectedDatasetsTable
      d$cts_out <- NULL
      d$datasetNames <- names(datasetList)
    })
    d$newDataLoaded <- d$newDataLoaded + 1
  }, priority = 100)
  
  ## DATA - select resolution for scRNA-seq
  output$res <- renderUI({
    req(input$data_type)
    if(input$data_type != "Single-cell") return()
    if(is.null(d$allowPrecomputedClusters)) return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    choices <- c("Run multiple resolutions using seurat" = "seurat_res")
    if(d$allowPrecomputedClusters) {
      choices <- c("Use pre-assigned clusters in metadata" = "clus_meta",
                   choices)
    }
    
    tagList(
      br(),
      selectInput(
        inputId = session$ns("res"),
        label = "Type of clustering",
        choices = choices
      )
    )
    
  })
  
  # DATA - fraction of total variance to use for choosing number of principal 
  # components in Seurat processing
  output$pcaVarianceFraction <- renderUI({
    req(input$data_type)
    if(input$data_type != "Single-cell" || is.null(input$res) || 
       input$res != "seurat_res" || is.null(d$end_ready) || is.null(d$end_partial)) return()
    
    # Exit if either no data selected or no name provided for user-uploaded data
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    div(
      style = "width: 200px; margin-left: 20px;",
      numericInput(
        inputId = session$ns("pcaVarianceFraction"),
        label = div(style = "width: 200px;", "PCA cumulative variance %"),
        value = 75,
        min = 0,
        max = 100,
        step = 1,
        width = "100px"
      )
    )
  })
  
  # DATA - update drt value when user changes it
  observeEvent(input$pcaVarianceFraction, {
    if(is.null(input$pcaVarianceFraction)) return()
    if(IsNumeric(input$pcaVarianceFraction)) {
      if(input$pcaVarianceFraction > 0 && input$drt <= 100) {
        d$newExpSetup <- T
        d$pcaVarianceFraction <- as.numeric(input$pcaVarianceFraction)
        return()
      }
    }
    updateNumericInput(
      session = session,
      inputId = "pcaVarianceFraction",
      value = 75,
    )
  })
  
  # DATA - specify detection rate threshold for scRNA-Seq data
  output$drt <- renderUI({
    req(input$data_type)
    if(input$data_type != "Single-cell") return()
    
    # Exit if either no data selected or no name provided for user-uploaded data
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    numericInput(
      inputId = session$ns("drt"),
      label = div(style = "width: 180px;", "Detection rate threshold"),      
      value = 0.1,
      min = 0.01,
      max = 1,
      step = 0.01,
      width = "100px"
    )
  })
  
  # Seurat resolution values
  output$resMin <- renderUI({
    req(input$data_type, input$res)
    if(input$data_type != "Single-cell" || input$res != "seurat_res") return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    numericInput(
      inputId = session$ns("resMin"),
      label = "From",
      value = 0.4,
      min = 0.1,
      step = 0.1
    )
  })
  
  output$resMax <- renderUI({
    req(input$data_type, input$res)
    if(input$data_type != "Single-cell" || input$res != "seurat_res") return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    numericInput(
      inputId = session$ns("resMax"),
      label = "To",
      value = 2.8,
      min = 0.1,
      max = 6,
      step = 0.1
    )
  })
  
  output$resStep <- renderUI({
    req(input$data_type, input$res)
    if(input$data_type != "Single-cell" || input$res != "seurat_res") return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    numericInput(
      inputId = session$ns("resStep"),
      label = "By",
      value = 0.4,
      min = 0,
      step = 0.05
    )
  })
  
  # Button to restore resolution range defaults
  output$resDefault <- renderUI({
    req(input$data_type, input$res)
    if(input$data_type != "Single-cell" || input$res != "seurat_res") return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    actionButton(
      inputId = session$ns("resDefault"),
      label = "Default",
      icon = icon("redo")
    )
  })
  
  observeEvent(input$resDefault, {
    updateNumericInput(session, inputId = "resMin", value = 0.4)
    updateNumericInput(session, inputId = "resMax", value = 2.8)
    updateNumericInput(session, inputId = "resStep", value = 0.4)
  })
  
  # Get d$resValues from min, max, and step
  observe({
    if(is.null(input$data_type) || input$data_type != "Single-cell" ||
       is.null(input$res) || input$res != "seurat_res" || 
       !is.numeric(input$resMin) || !is.numeric(input$resMax) ||
       !is.numeric(input$resStep)) return()

    # If max is less than min and min > 0
    if(input$resMax < input$resMin && input$resMin > 0) {
      if(input$resMin > 6) {
        updateNumericInput(session, inputId = "resMin", value = input$resMax)
      } else {
        updateNumericInput(session, inputId = "resMax", value = input$resMin)
      }
      return()
    } 
    
    # If min < 0
    if(input$resMin < 0) {
      updateNumericInput(session, inputId = "resMin", value = 0)
      return()
    }
    
    # If max < 0
    if(input$resMax < 0) {
      updateNumericInput(session, inputId = "resMax", value = 0)
      return()
    }
    
    # If max and min are both 0
    if(input$resMin == 0 && input$resMax == 0) {
      updateNumericInput(session, inputId = "resMax", value = 0.1)
      return()
    } 
    
    # Limit max to 6 or less
    if(input$resMax > 6) {
      updateNumericInput(session, inputId = "resMax", value = 6)
      return()
    }
    
    if(input$resStep == 0) {
      resValues <- c(input$resMin, input$resMax)
    } else {
      resValues <- seq(from = input$resMin, to = input$resMax, by = input$resStep)
    }
    resValues <- sort(unique(c(input$resMin, input$resMax, resValues)))
    resValues <- resValues[resValues > 0]
    d$resValues <- resValues
  })
  
  # Show user Seurat resolutions chosen
  output$resValues <- renderUI({
    req(input$data_type, input$res, d$resValues)
    if(input$data_type != "Single-cell" || input$res != "seurat_res" ||
       is.null(d$allowPrecomputedClusters)) return()
    
    # Exit if either no data selected or no user-uploaded data provided
    if(is.null(d$select_datasets_table_rows_selected) & !d$end_ready) {
      return()
    }
    if(d$end_partial) return()
    
    txt <- paste0(
      strong("Resolutions: "),
      paste(d$resValues, collapse = ", ")
    )
    HTML(txt)
  })
  
  # DATA - update drt value when user changes it
  observeEvent(input$drt, {
    if(is.null(input$drt)) return()
    if(IsNumeric(input$drt)) {
      if(input$drt >= 0.01 & input$drt <= 1) {
        d$newExpSetup <- T
        d$drt <- as.numeric(input$drt)
        return()
      }
    }
    updateNumericInput(
      session = session,
      inputId = "drt",
      label = "Detection rate threshold",
      value = 0.1,
      min = 0.01, max = 1,
      step = 0.01
    )
  })
  
  # DATA - reactive data transformation
  ddstran <- eventReactive(d$newDataLoaded, {
    if(input$data_type == "RASL") return()
    dds <- ddsout()[[1]]
    if(input$data_type == "Bulk") {
      if(input$transform == "log") {
        tmp <- normTransform(dds)
        lab <- "log<sub>2</sub>(counts + 1)"
      } else if(input$transform == "rlog") {
        tmp <- rlog(dds)
        lab <- "rlog(counts)"
      } else if(input$transform == "vst") {
        tmp <- vst(dds)
        lab <- "vst(counts)"
      } else if(input$transform == "raw") {
        tmp <- RawDDSTrans(dds)
        lab <- "Raw Counts"
      }
    } else {
      tmp <- dds
      lab <- "Normalized Counts"
    }
    return(list(tmp, lab))
  })
  
  # SC-DGE-OVER - create seurat object for downstream analysis
  seurat_only <- reactive({
    # req(SubmitData$data_type, d$resType)
    if(is.null(input$data_type) || input$data_type != "Single-cell" || 
       is.null(input$res) || is.null(ddsout())) return()
    withProgress(message = "Building Seurat analysis...", value = 0, {
      incProgress(1/3)
      
      dat <- CreateSeuratObject(ddsout()[[3]], assay = "RNA",
                                meta.data = data.frame(colData(ddsout()[[1]])))
      # Use number of principal components based on input fraction of variance explained
      # and use verbose = F to avoid error from RunUMAP step small number of cells
      dat <- NormalizeData(dat, verbose = F)
      dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000, verbose = F)
      dat <- ScaleData(dat, verbose = F)
      nPCs <- min(50, nrow(dat)-1, ncol(dat)-1)
      dat <- RunPCA(dat, features = VariableFeatures(object = dat), npcs = nPCs, verbose = F)
      nPCs_targetFrac <- PCADims(seuratObject = dat, varFrac = d$pcaVarianceFraction/100)
      if(nPCs_targetFrac > 1) nPCs <- nPCs_targetFrac
      dat <- RunUMAP(dat, dims = 1:nPCs, n.neighbors = min(30, nrow(dat), ncol(dat)), verbose = F)
      dat <- RunTSNE(dat, perplexity = min(30, (ncol(dat)-1)/3), check_duplicates = F)
      dat <- FindNeighbors(dat, dims = 1:nPCs, k.param = min(20, ncol(dat) - 1))
      
      incProgress(2/3)
      # for data pulled from mysql, rename orig.ident
      if (unique(dat@meta.data$orig.ident == "SeuratProject")) {
        dat$orig.ident <- gsub(".*\\.", "", rownames(dat@meta.data))
      }
      # for several seurat resolution types
      if (input$res == "seurat_res") {
        seur_obj <- list()
        # for (i in seq(.4, 2.8, by = .4)) {
        for (i in d$resValues) {
          seurat_resolution <- 0 + i
          # ^ iteratively incrementing resolution parameter
          obj <- tryCatch(Seurat::FindClusters(dat, resolution = seurat_resolution),
                          error = function(e) NULL)
          if(is.null(obj)) next
          if(length(unique(obj$seurat_clusters)) == 1) next
          seur_obj[[paste0("res.", seurat_resolution)]] <- obj
        }
        return(seur_obj)
        # for selected metadata file
      } else {
        # Remove cells that are the only cell in a cluster (prevents errors with scClustViz)
        clustCol <- grep("^res|^Cluster|^seurat", x = colnames(dat@meta.data), value = T)[1]
        clustCellCounts <- table(dat@meta.data[, clustCol])
        singleClusters <- names(clustCellCounts[clustCellCounts == 1])
        if(length(singleClusters) > 0) {
          dat <- dat[, !(dat@meta.data[, clustCol] %in% singleClusters)]
          dat@meta.data[, clustCol] <- droplevels(dat@meta.data[, clustCol])
        }
        return(dat)
      }
    })
  })
  
  # SC-DGE-OVER - convert seurat object to a sCVdata object
  seurat_sc <- reactive({
    if(is.null(input$data_type) || input$data_type != "Single-cell" || 
       is.null(input$res) || is.null(seurat_only())) return()
    if(input$res == "seurat_res" && length(seurat_only()) == 0) return()
    withProgress(message = "Building scClustViz analysis...", value = 0, {
      incProgress(1/3)
      if (input$res == "seurat_res") {
        DE_bw_clust <- TRUE
        sCVdata_list <- list()
        # for (i in seq(.4, 2.8, by = .4)) {
        for (i in d$resValues) {
          seurat_resolution <- 0 + i
          seuratObjName <- paste0("res.", i)
          if(is.null(seurat_only()[[seuratObjName]])) next
          if(length(unique(seurat_only()[[seuratObjName]]$seurat_clusters)) == 1) next
          cluster_cols <- paste0("RNA_snn_res.", seurat_resolution)
          cluster_res <- getMD(seurat_only()[[seuratObjName]])[cluster_cols]
          curr_sCVdata <- CalcAllSCV(inD = seurat_only()[[seuratObjName]],
                                     assayType = "RNA",
                                     clusterDF = cluster_res,
                                     DRforClust = "pca",
                                     exponent = exp(1),
                                     pseudocount = 1,
                                     DRthresh = d$drt,
                                     testAll = T,
                                     calcSil = T,
                                     calcDEvsRest = T,
                                     calcDEcombn = T)
          DE_bw_NN <- sapply(DEneighb(curr_sCVdata[[1]], 0.05), length)
          # ^ counts # of DE genes between neighboring clusters at 5% FDR
          if (min(DE_bw_NN) < 1) { DE_bw_clust <- FALSE }
          # ^ If no DE genes between nearest neighbors, don't loop again.
          sCVdata_list[[paste0("res.", seurat_resolution)]] <- curr_sCVdata
        }
        incProgress(2/3)
        sCVdata_list <- unlist(sCVdata_list, recursive = F)
        names(sCVdata_list) <- gsub("^res\\.[0-9]{1}\\.[0-9]*\\.*", "", names(sCVdata_list))
        return(sCVdata_list)
      } else {
        
        cols <- names(seurat_only()@meta.data)
        # may need to update this!
        cols <- cols[grepl(pattern = "^res|^Cluster|^seurat", x = cols)][1]
        meta <- seurat_only()@meta.data
        seurat_clusters <- meta[, cols]
        dat <- CalcAllSCV(inD = seurat_only(),
                          assayType = "RNA",
                          clusterDF = seurat_clusters,
                          DRforClust = "pca",
                          exponent = exp(1),
                          pseudocount = 1,
                          DRthresh = d$drt,
                          testAll = T,
                          calcSil = T,
                          calcDEvsRest = T,
                          calcDEcombn = T)
        return(dat)
      }
    })
  })
  
  observeEvent(d$newDataLoaded, {
    withProgress(message = "Filtering data...", value = 0, {
      incProgress(1/3)
      d$data_type <- input$data_type
      d$ddsout <- ddsout()
      incProgress(1/3)
      d$ddstran <- ddstran()
      d$resType <- input$res
      if(input$data_type == "Single-cell" && input$res == "seurat_res" && length(seurat_only()) == 0) {
        d$goqc_warning <- "Clustering could not be performed. Please try different data or settings."
      } else {
        d$goqc_warning <- NULL
      }
      d$seurat_only <- seurat_only()
      d$seurat_sc <- seurat_sc()
      d$goqc <- d$goqc + 1
      incProgress(1/3)
    })
  }, ignoreInit = T)
  
  return(d)
}
