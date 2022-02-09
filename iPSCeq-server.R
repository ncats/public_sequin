#---------------------------------------------------------------------
# Title:         NCATS SEQUIN - Server Script
# Author:        Marissa Hirst
# Author2:       Ben Ernest
# Last Modified: 2022-02-08
# --
# Created:       2018-01-26 11:29:39 CDT
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# CONTENTS
#---------------------------------------------------------------------
#
#   NCATS server function organization
#
#       01 SHINY OPTIONS
#       02 LAYOUT OF APP... ui elements for bulk/sc layout
#       03 DATA & QC....... data loading & quality control for preliminary analysis
#       04 DISCOVERY....... correlation by cells/samples, distances, pca/tsne/umap
#       05 DGE............. differential gene expression analysis: bulk or scRNA-seq
#       06 iPSC PROFILER... module score profiler for scRNA-seq
#       07 MORE............ tutorial, faq, about us, session info
#
#---------------------------------------------------------------------
#
#   NCATS SEQUIN application layout
#
#       SEQUIN
#       |-  Shiny options & bulk/scRNA-seq dataset list load from SQL
#       |-  Layout of app for bulk/scRNA-Seq side
#           |-  d: Store reactive values in the app here
#           |-  Defaults: Image type (png), annotation type (HS) and cell markers
#           |-  Toggle tab panel elements on bulk/scRNA-Seq
#           |-  ui element layout for bulk/scRNA-Seq
#       |-  Submit and quality control
#           |-  Select dataset: bulk, scRNA-Seq
#           |-  File summary
#           |-  Count summary
#       |-  Discovery-driven analysis
#           |-  Correlation & scatterplot
#           |-  Distance matrix (bulk RNA-Seq only)
#           |-  Clustering: PCA (bulk only), tSNE, uMAP
#       |-  DGE analysis:
#           |- Bulk RNA-seq:
#             |-  Overview: BULK-DGE-OVER
#             |-  Volcano plots: BULK-DGE-VOL
#             |-  Gene set enrichment & heatmap: BULK-DGE-GSE
#             |-  Clustering: BULK-DGE-CLUS
#           |- scRNA-seq:
#             |-  Overview: SC-DGE-OVER
#             |-  Gene expression: SC-DGE-GE
#             |-  DGE analysis: SC-DGE-DGE
#             |-  Custom DGE analysis: SC-DGE-CUST
#             |-  Volcano plots: SC-DGE-VOL
#             |-  Gene set enrichment & heatmap: SC-DGE-GSE
#             |-  Manually select cells for DGE: SC-DGE-MSC
#             |-  Clustering: SC-DGE-CLUS
#             |-  Combine clusters or cells: SC-DGE-CCC
#       |   iPSC profiler for scRNA-Seq
#       |-  More
#           |-  Tutorial
#           |-  FAQ
#           |-  About us
#           |-  Session info

# Change file upload size to 200 MB and sanitize errors
options(
    shiny.maxRequestSize = 2000 * 1024^2,
    shiny.sanitize.errors = TRUE,
    getClass.msg = F,
    spinner.color = "black", 
    spinner.color.background = "white"
)

# Server function
iPSCeqServer <- function(input, output, session) {
  
  # Hide tabs initially
  hideTab(inputId = "tab_structure", target = "val4")
  hideTab(inputId = "tab_structure", target = "val3")
  hideTab(inputId = "tab_structure", target = "val2")
  hideTab(inputId = "tab_structure", target = "summary")
  
  # Initial popup
  showModal(modalDialog(
    size = "m",
    title = fluidPage(
      div(
        style = "height: 80px",
        img(
          src = "ncats-logo.png", height = 60, 
          style = "display: inline-block; vertical-align: top; margin-top: 12px; padding-left: 10px;"
        ),
        div(
          style = paste0(
            "display: inline-block; height: 60px; color: #672e6c; font-size: 60px;",
            " font-weight: 350; float: right; padding-right: 10px; margin-top: -5px;"
          ),
          "SEQUIN"
        )
      ),
      div(
        style = "font-size: 16px; margin-top: 5px; padding-left: 10px; color: #636569;",
        "Stem Cell Translation Laboratory"
      )
    ),
    footer = div(style = "float: right;", actionButton(inputId = "closeSplash", label = "Enter")),
    includeMarkdown(path = "markdown/splash.md")
  ))
  
  observeEvent(input$closeSplash, {
    removeModal()
    updateTabsetPanel(inputId = "main_tabsetPanel", selected = "mainApp")
  })

  output$placeholderDesc <- renderUI({
    if(!is.null(d$cts_out)) return()
    div(
      style = "width: 800px; margin-left: 50px; font-size: 14px; font-weight: 300;",
      HTML(paste0(
        "Load data to get started. For help, see the tutorial under ",
        strong("More"), " or click the ", strong("Help"), " link to ",
        "the right."
      ))
    )
  })
  
  # Store all reactive values used in app
  d <- reactiveValues(
    MD = NULL, 
    SCV = NULL, 
    inD = NULL, 
    res = NULL,
    a = NULL,
    b = NULL,
    compNames = NULL,
    cts_out = NULL,
    dgeContrast = NULL,
    newExpSetup = T,
    newGSE = T,
    newHeat = T,
    newCluster = T,
    dge_warning_msg = "",
    x = NULL, # for lasso selection tool for inter-cluster DE boxplots
    y = NULL, # for lasso selection tool for inter-cluster DE boxplots
    genes_debug = NULL,
    n_genes = 100,
    n_genes_heatmap = 100,
    n_genes_sc = 100,
    n_genes_heatmap_sc = 100,
    n_genes_manuallyselectcells_sc = 100,
    fea_all_genes = F, # whether to use all genes in GSE
    fea_all_genes_heatmap = F,
    sc_mytabTemp = NULL,
    filts = NULL,
    gse_warning_msg = "",
    heatmap_warning_msg = "",
    calcTextVal = "",
    dge_info = NULL,
    dge_info_heatmap = NULL,
    dge_info_sc = NULL,
    dge_info_heatmap_sc = NULL,
    dge_info_manuallyselectcells_sc = NULL,
    clust_time = 0,
    bulkExpFactorInfo = NULL,
    drt = 0.1,
    pcaVarianceFraction = 1,
    fdr = 0.05,
    goqc_or_relaunch = 0,
    sc_meta = NULL,
    clustList = NULL,
    numClust = NULL,
    seurat_only = NULL,
    seurat_sc = NULL,
    test = "none",
    selectedDatasetsTable = NULL,
    selectSamplesTable = NULL,
    selectSamplesTableList = NULL,
    selectSamplesTableVarsList = NULL,
    selectSamplesTableSelectedVar = NULL,
    selectedData = NULL,
    upload_warn = NULL,
    inD_orig = NULL,
    setNames = NULL,
    cellSets = list(),
    addCellsClicks = list(),
    observeCount = 0, 
    addCellsOverlap = list(),
    data_info = NULL,
    corplotClick = NULL,
    heatmap1Click = NULL,
    heatmap1Hover = NULL,
    legendList_sc = NULL,
    moduleScoreSeurat = NULL, 
    inD_withModuleScore = NULL,
    moduleNamesRefVec = NULL,
    selected_modules = NULL,
    inD_goi = NULL,
    goqc_warning = NULL,
    wgcna_warning = NULL,
    wgcna_minModuleSizeTried = NULL,
    wgcna_warning_sc = NULL,
    wgcna_minModuleSizeTried_sc = NULL,
    inputFactors = NULL,
    timeStamps = list(),
    newCorplot = T,
    newDistplot = T,
    sc_dge_byfactor = list(),
    sc_dge_byfactor_msg = NULL,
    sc_dge_byfactor_novars_msg = NULL,
    allowPrecomputedClusters = F,
    datasetsRowsSelected = NULL,
    datasetNames = NULL,
    data_type = NULL,
    resType = NULL,
    showDataSummaryLinks = T,
    showQCLinks = F,
    showCorrelationLinks = T,
    correlation_showLabels = F,
    distance_showLabels = F,
    showClusteringLinks = F,
    enrichRLib = NULL,
    enrichRLib_sc = NULL,
    enrichRLib_manuallyselectcells = NULL,
    sc_showOverviewLinks = T,
    sc_showClusteringLinks = F,
    sc_showGeneExpressionLinks = F,
    sc_showDGELinks = F,
    showTabs = T,
    totalFigureHeight = NULL,
    totalFigureHeight_sc = NULL,
    resValues = NULL,
    sc_dge_factor_lfcThreshold = 0.25,
    sc_dge_factor_test = "wilcox",
    sc_dge_factor_minPct = 0.1,
    sc_dge_factor_minDiffPct = -1
  )

  count_content <- reactive({
    if(is.null(d$end_count)) return()
    fread(d$end_count$datapath, header = T) %>%
      mutate_all(funs(replace_na(.,0))) %>% 
      column_to_rownames(names(.)[1])
  })
  
  meta_content <- reactive({
    if(is.null(d$end_meta)) return()
    fread(d$end_meta$datapath, header = T) %>%
      column_to_rownames(names(.)[1])
  })
  
  # Dataset names in row at top of app
  # output$datasetNames <- renderText({
  output$datasetNames <- renderUI({
    req(SubmitData$data_type, d$datasetNames, input$tab_structure)
    txt <- paste(d$datasetNames, collapse = ", ")
    dataType <- "Bulk RNA-seq experiment"
    if(length(d$datasetNames) > 1) dataType <- "Bulk RNA-seq experiments"
    if(SubmitData$data_type == "Single-cell") {
      dataType <- "Single-cell RNA-seq experiment"
      if(length(d$datasetNames) > 1) dataType <- "Single-cell RNA-seq experiments"
    }
    if(SubmitData$data_type == "RASL") {
      dataType <- "RASL-seq experiment"
      if(length(d$datasetNames) > 1) dataType <- "RASL-seq experiments"
    }
    h4(strong(paste0(dataType, ":")), txt)
  })
  
  # Show selected resolution for single-cell in top bar
  output$sc_res <- renderUI({
    req(SubmitData$data_type, d$resType, d$res)
    if(SubmitData$data_type != "Single-cell") return()
    myRes <- "pre-computed"
    if(d$resType == "seurat_res") {
      myRes <- gsub("RNA_snn_res\\.", replacement = "", x = d$res)
    }
    h4(strong("Resolution: "), myRes)
  })
  
  # Button to download Seurat and scClustViz data
  output$download_sc <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type != "Single-cell" || is.null(d$inD) || is.null(d$seurat_sc)) return()
    actionButton("download_sc", label = "Download", icon = icon("download"))
  })
  
  output$downloadSeurat_msg <- renderUI({
    req(d$resType)
    txt <- "Download data in Seurat format as an RDS file which can be loaded into R."
    if(d$resType == "seurat_res") {
      txt <- paste0(
        "Download a ",
        strong("Seurat"), 
        " object containing data at the current resolution or a ",
        strong("List"),
        " of Seurat objects from all resolutions. Data will be saved as an",
        " RDS file which can be directly loaded in R."
      )
    }
    HTML(txt)
  })
  
  output$downloadSCV_msg <- renderUI({
    req(d$resType)
    txt <- paste0(
      "Download data in sCVdata format (scClustViz package) as an RDS file",
      " which can be loaded into R."
    )
    if(d$resType == "seurat_res") {
      txt <- paste0(
        "Download an ",
        strong("sCVdata"),
        " object (scClustViz package) containing data at the current resolution or a ",
        strong("List"),
        " of sCVdata objects from all resolutions. Data will be saved as an",
        " RDS file which can be directly loaded in R."
      )
    }
    HTML(txt)
  })
  
  output$downloadSeurat_format <- renderUI({
    req(d$resType)
    if(d$resType != "seurat_res") return()
    awesomeRadio(
      inputId = "downloadSeurat_format", 
      label = NULL,
      choices = c(
        "Seurat (current resolution)" = "single",
        "List (all resolutions)" = "multi"
      ),
      selected = "single",
      status = "success",
      width = "240px"
    )
  })
  
  output$downloadSCV_format <- renderUI({
    req(d$resType)
    if(d$resType != "seurat_res") return()
    awesomeRadio(
      inputId = "downloadSCV_format", 
      label = NULL,
      choices = c(
        "sCVdata (current resolution)" = "single",
        "List (all resolutions)" = "multi"
      ),
      selected = "single",
      status = "success",
      width = "240px"
    )
  })
  
  observeEvent(input$download_sc, {
    showModal(modalDialog(
      title = "Download single-cell data",
      size = "l",
      easyClose = T,
      fluidPage(
        fluidRow(h4(strong("Seurat"))),
        fluidRow(uiOutput("downloadSeurat_msg")),
        br(),
        fluidRow(
          div(
            style = "display: inline-block;",
            uiOutput("downloadSeurat_format"),
          ),
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = "downloadSeurat", 
              label = "Download Seurat",
              style = "vertical-align: bottom;"
            )
          )
        ),
        br(), br(),
        fluidRow(h4(strong("scClustViz"))),
        fluidRow(uiOutput("downloadSCV_msg")),
        br(),
        fluidRow(
          div(
            style = "display: inline-block;",
            uiOutput("downloadSCV_format"),
          ),
          div(
            style = "display: inline-block;",
            downloadButton(
              outputId = "downloadSCV", 
              label = "Download scClustViz",
              style = "vertical-align: bottom;"
            )
          )
        )
      )
    ))
  })
  
  output$downloadSeurat <- downloadHandler(
    filename = function() {
      if(!is.null(input$downloadSeurat_format) && 
         input$downloadSeurat_format == "multi" && d$resType == "seurat_res") {
        return("SeuratList.rds")
      }
      return("Seurat.rds")
    },
    content = function(file) {
      if(!is.null(input$downloadSeurat_format) &&
         input$downloadSeurat_format == "multi" &&  d$resType == "seurat_res") {
        saveRDS(seurat_only(), file = file)
      } else {
        saveRDS(d$inD, file = file)
      }
    }
  )
  
  output$downloadSCV <- downloadHandler(
    filename = function() {
      if(!is.null(input$downloadSCV_format) && 
         input$downloadSCV_format == "multi" && d$resType == "seurat_res") {
        return("sCVdataList.rds")
      }
      return("sCVdata.rds")
    },
    content = function(file) {
      if(!is.null(input$downloadSCV_format) && 
         input$downloadSCV_format == "multi" && d$resType == "seurat_res") {
        saveRDS(d$SCV, file = file)
      } else {
        saveRDS(d$SCV[[d$res]], file = file)
      }
    }
  )
  
  
  # Reset app after user clicks "Reset" and confirms
  observeEvent(input$reset, {
    showModal(modalDialog(
      easyClose = T,
      footer = NULL,
      fluidPage(
        p("Do you really want to reset the app? All loaded data and analyses will be erased."),
        p("Click anywhere outside this box to cancel or click ", strong("OK"), "to reset."),
        actionButton("reset_ok", label = "OK")
      )
    ))
  })
  
  observeEvent(input$reset_ok, {
    removeModal()
    shinyjs::refresh()
  })
  
  # Go back to Load Data page when user clicks "Update"
  observeEvent(input$update, {
    updateNavbarPage(inputId = "tab_structure", selected = "val1")
  })
  
  # Update button
  output$update <- renderUI({
    req(input$tab_structure)
    if(input$tab_structure == "val1") return()
    actionButton("update", label = "Update", icon = icon("edit"))
  })
  
  # Reset button
  output$reset <- renderUI({
    actionButton("reset", label = "Reset", icon = icon("redo"))
  })
  
  # Top bar with dataset names and reset button
  output$topBar <- renderUI({
    req(input$tab_structure, SubmitData)
    if(input$tab_structure == "val1" & is.null(SubmitData$selectedData)) return()
    
    div(
      style = "background-color: #eee; vertical-align: middle; margin-bottom: 20px;",
      fluidPage(
        # br(),
        div(
          style = "margin-top: 10px; margin-bottom: 10px;",
          fluidRow(
            column(
              width = 6,
              div(
                # style = "font-size: 20px; vertical-align: center;",
                uiOutput("datasetNames")
              )
            ),
            column(
              width = 3,
              uiOutput("sc_res")
            ),
            column(
              width = 3,
              align = "right",
              fluidRow(
                div(style = "display: inline-block;", uiOutput("download_sc")),
                div(style = "display: inline-block;", uiOutput("update")),
                div(style = "display: inline-block; margin-right: 15px;", uiOutput("reset"))
              )
            )
          )
        )
        # br()
      )
    )
  })
  
  # Bulk toggle for sidebar elements to only appear in certain tabs
  observe({
    # Reactive values for # genes to use in GSE/heatmap
    if(is.null(input$n_genes)) {
      d$n_genes <- 100
    } else {
      d$n_genes <- input$n_genes
    }
    
    if(!IsInteger(input$n_genes_heatmap) || input$n_genes_heatmap < 1) {
      if(!is.null(dgeout3()) && nrow(dgeout3()) > 1) {
        maxGenes <- nrow(dgeout3())
        d$n_genes_heatmap <- min(nrow(dgeout3()), 100)
      } else {
        d$n_genes_heatmap <- 100
      }
    } else {
      d$n_genes_heatmap <- input$n_genes_heatmap
    }
  })
  
  # scRNA-Seq toggle for sidebar elements to only appear in certain tabs
  observe({
    # Reactive values for # genes to use in GSE/heatmap
    if(is.null(input$n_genes_sc)) {
      d$n_genes_sc <- 100
    } else {
      d$n_genes_sc <- input$n_genes_sc
    }
    
    if(!IsInteger(input$n_genes_heatmap_sc) || input$n_genes_heatmap_sc < 1) {
      if(!is.null(sc_mytab_heatmap()) && nrow(sc_mytab_heatmap()) > 1) {
        maxGenes <- nrow(sc_mytab_heatmap())
        d$n_genes_heatmap_sc <- min(nrow(sc_mytab_heatmap()), 100)
      } else {
        d$n_genes_heatmap_sc <- 100
      }
    } else {
      d$n_genes_heatmap_sc <- input$n_genes_heatmap_sc
    }
    
    if(is.null(input$n_genes_manuallyselectcells_sc)) {
      d$n_genes_manuallyselectcells_sc <- 100
    } else {
      d$n_genes_manuallyselectcells_sc <- input$n_genes_manuallyselectcells_sc
    }
  })
  
  ############################################################  
  ############################################################
  ### UI OUTPUT FOR DGE ANALYSIS: BULK & scRNA-SEQ 
  ############################################################
  ############################################################
  
  # bulk rna-seq uiOutput for DGE Analysis tab
  output$showRunDGE <- renderUI({
    req(input$dge, SubmitData$data_type)
    if(SubmitData$data_type != "Bulk") return()
    myLab <- "Run DGE"
    if(input$dge == "runDGE") myLab <- strong(myLab)
    actionLink(
      inputId = "showRunDGE", 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })
  
  output$showOverview <- renderUI({
    req(SubmitData$data_type, dgeover(), input$dge)
    if(SubmitData$data_type == "Single-cell") return()
    if(is.null(dgeover())) return()
    if(d$dge_warning_msg != "") return()
    myLab <- "Overview"
    if(input$dge == "overview") myLab <- strong(myLab)
    actionLink(
      inputId = "showOverview", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$showVolcano <- renderUI({
    req(SubmitData$data_type, dgeout1(), input$dge)
    if(SubmitData$data_type == "Single-cell") return()
    if(is.null(dgeout1())) return()
    if(d$dge_warning_msg != "") return()
    myLab <- "Volcano plot"
    if(input$dge == "volcano") myLab <- strong(myLab)
    actionLink(
      inputId = "showVolcano", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$showGSE <- renderUI({
    req(SubmitData$data_type, input$dge)
    if(SubmitData$data_type == "Single-cell") return()
    myLab <- "Gene set enrichment"
    if(input$dge == "gse") myLab <- strong(myLab)
    actionLink(
      inputId = "showGSE", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$showHeatmap <- renderUI({
    req(SubmitData$data_type, input$dge)
    if(SubmitData$data_type == "Single-cell") return()
    myLab <- "Heatmap"
    if(input$dge == "heatmap") myLab <- strong(myLab)
    actionLink(
      inputId = "showHeatmap", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$showDGEClustering <- renderUI({
    req(input$dge)
    myLab <- "Clustering"
    if(input$dge == "clustering") myLab <- strong(myLab)
    actionLink(
      inputId = "showDGEClustering", 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })
  
  output$clusteringLinks <- renderUI({
    req(input$clustalg, clustout(), d$showClusteringLinks, input$clustering_tabsetPanel)
    if(length(clustout()) != 2 & length(clustout()) != 8) return()
    if(length(clustout()) == 8) {
      lab2 <- "WGCNA - dendrogram"
      lab3 <- "WGCNA - TOM"
      if(input$clustering_tabsetPanel == "clustPlotW02") lab2 <- strong(lab2)
      if(input$clustering_tabsetPanel == "clustPlotW03") lab3 <- strong(lab3)
      
      fluidPage(
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "showClustPlotW02", 
            label = lab2,
            style = "color: black;"
          )
        ),
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "showClustPlotW03", 
            label = lab3,
            style = "color: black;"
          )
        )
      )
    } else {
      fluidPage(
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "showClustPlotW01", 
            label = strong("K-medoids"),
            style = "color: black;"
          )
        )
      )
    }
  })
  
  observeEvent(input$showRunDGE, {
    d$showClusteringLinks <- F
    updateTabsetPanel(inputId = "dge", selected = "runDGE")
  })
  
  observeEvent(input$showOverview, {
    d$showClusteringLinks <- F
    updateTabsetPanel(inputId = "dge", selected = "overview")
  })
  
  observeEvent(input$showVolcano, {
    d$showClusteringLinks <- F
    updateTabsetPanel(inputId = "dge", selected = "volcano")
  })
  
  observeEvent(input$showGSE, {
    d$showClusteringLinks <- F
    updateTabsetPanel(inputId = "dge", selected = "gse")
  })
  
  observeEvent(input$showHeatmap, {
    d$showClusteringLinks <- F
    updateTabsetPanel(inputId = "dge", selected = "heatmap")
  })
  
  observeEvent(input$showDGEClustering, {
    d$showClusteringLinks <- ifelse(
      test = (input$dge == "clustering"), 
      yes = ifelse(
        test = d$showClusteringLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    updateTabsetPanel(inputId = "dge", selected = "clustering")
  })
  
  observeEvent(input$showClustPlotW01, {
    updateTabsetPanel(inputId = "clustering_tabsetPanel", selected = "clustPlotW01")
  })
  
  observeEvent(input$showClustPlotW02, {
    updateTabsetPanel(inputId = "clustering_tabsetPanel", selected = "clustPlotW02")
  })
  
  observeEvent(input$showClustPlotW03, {
    updateTabsetPanel(inputId = "clustering_tabsetPanel", selected = "clustPlotW03")
  })

  output$dge_bulk <- renderUI({
    # if(is.null(SubmitData$goqc) || SubmitData$goqc == 0) return()
    req(SubmitData$goqc, SubmitData$data_type)
    sidebarLayout(
      sidebarPanel(
        width = 2,
        style = "width: 190px;",
        fluidRow(div(style = "margin-left: 10px;", uiOutput("showRunDGE"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("showOverview"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("showVolcano"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("showHeatmap"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("showGSE"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("showDGEClustering"))),
        fluidRow(div(style = "margin-left: 4px;", uiOutput("clusteringLinks")))
      ),
      mainPanel(
        width = 10,
        tabsetPanel(
          id = "dge",
          type = "hidden",
          selected = ifelse(
            test = (SubmitData$data_type == "Bulk"), 
            yes = "runDGE", 
            no = "heatmap"
          ),
          tabPanelBody(
            value = "runDGE",
            fluidPage(
              fluidRow(
                div(style = "display: inline-block;", h3("DGE options", style = "margin-top: 0px;")),
                div(
                  style = "display: inline-block; margin-top: 0px; float: right;",
                  actionLink("overview_help_button", label = "Help")
                )
              ),
              fluidRow(
                column(
                  width = 4,
                  style = "background-color: #eee; width: 380px; height: 160px;",
                  fluidPage(
                    fluidRow(uiOutput("dgeexpsetup")),
                    fluidRow(
                      div(style = "display: inline-block; width: 160px;", uiOutput("dgemethod")),
                      div(style = "display: inline-block; width: 160px; margin-left: 10px;", uiOutput("dgeexpedgernorm"))
                    )
                  )
                ),
                column(
                  width = 2,
                  style = "background-color: #eee; width: 160px; height: 160px; margin-left: 20px;",
                  fluidPage(
                    fluidRow(div(style = "width: 120px;", uiOutput("dgepadjcutoff"))),
                    fluidRow(div(style = "width: 120px;", uiOutput("dgefcmin")))
                  )
                ),
                column(
                  width = 3,
                  style = "background-color: #eee; width: 280px; height: 160px; margin-left: 20px;",
                  fluidPage(
                    fluidRow(div(style = "width: 160px;", uiOutput("batchFactor"))),
                    fluidRow(uiOutput("housekeepingSelect"))
                  )
                )
              ),
              br(),
              conditionalPanel(
                condition = "input.dgeexpsetup == 'exp1'",
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp1a")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp1b")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp1c"))
                )
              ),
              conditionalPanel(
                condition = "input.dgeexpsetup == 'exp2'", 
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp2a")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp2c")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp2c5"))
                ),
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp2b")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp2d")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp2d5"))
                )
              ),
              conditionalPanel(
                condition = "input.dgeexpsetup == 'exp3'",
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp3a")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp3c")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp3b")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp3d"))
                )
              ),
              conditionalPanel(
                condition = "input.dgeexpsetup == 'exp4'",
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp4a")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp4c")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp4b")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp4d"))
                )
              ),
              conditionalPanel(
                condition = "input.dgeexpsetup == 'exp5'",
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp5a")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp5b")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp5c")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp5d"))
                )
              ),
              conditionalPanel(
                condition = "input.dgeexpsetup == 'exp6'",
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp6c")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp6d")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp6e")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp6f"))
                ),            
                fluidRow(
                  style = "background-color: #eee; width: 860px;",
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp6a")),
                  div(style = "display: inline-block; width: 200px; margin-left: 10px;", uiOutput("dgeexp6b"))
                )
              ),
              br(),
              fluidRow(
                style = "background-color: #eee; width: 860px;",
                div(style = "display: inline-block; margin-left: 10px; margin-top: 5px; margin-bottom: 10px;", uiOutput("dgeexpformhead")),
                div(
                  style = "display: inline-block; margin-left: 20px; margin-top: 5px; margin-bottom: 10px; vertical-align: top;",
                  conditionalPanel(condition = "input.dgeexpsetup == 'exp1'", uiOutput("dgeexpform1")),
                  conditionalPanel(condition = "input.dgeexpsetup == 'exp2'", uiOutput("dgeexpform2")),
                  conditionalPanel(condition = "input.dgeexpsetup == 'exp3'", uiOutput("dgeexpform3")),
                  conditionalPanel(condition = "input.dgeexpsetup == 'exp4'", uiOutput("dgeexpform4")),
                  conditionalPanel(condition = "input.dgeexpsetup == 'exp5'", uiOutput("dgeexpform5")),
                  conditionalPanel(
                    condition = "input.dgeexpsetup == 'exp6'",
                    uiOutput("dgeexpform6a"),
                    uiOutput("dgeexpform6b"),
                    uiOutput("dgeexpform6c"),
                    uiOutput("dgeexpform6d")
                  )
                ),
                div(style = "display: inline-block; margin-left: 40px; margin-top: 5px; margin-bottom: 10px; vertical-align: top;", uiOutput("godge"))
              ),
              br(),
              fluidRow(uiOutput("dge_expdesign_warning_ui")),
              fluidRow(uiOutput("dge_warning_msg_ui")),
              br()
            )
          ),
          tabPanelBody(
            value = "overview",
            fluidPage(
              fluidRow(uiOutput("dgeoverviewHeader")),
              fluidRow(div(style = "width: 600px;", uiOutput("dgeoverview"))),
              br(),
              fluidRow(uiOutput("dgeplot2Header")),
              fluidRow(plotlyOutput("dgeplot2", width = "600px", height = "400px")),
              fluidRow(
                div(style = "display:inline-block", uiOutput("dldgeoverpdf")),
                div(style = "display:inline-block", uiOutput("dldgeoverpng"))
              ),
              br()
            )
          ),
          tabPanelBody(
            value = "volcano",
            fluidPage(
              fluidRow(
                div(style = "display: inline-block; width: 280px;", uiOutput("dgemaincontrasts")),
                div(style = "display: inline-block; margin-left: 20px; vertical-align: top;", uiOutput("vistype")),
                actionLink(
                  inputId = "volcanoplots_help_button", 
                  label = "Help",
                  style = "display: inline-block; float: right;"
                )
              ),
              fluidRow(plotlyOutput("volplots")),
              br(),
              fluidRow(
                div(style = "display:inline-block;", uiOutput("dldgemavolpdf")),
                div(style = "display:inline-block; margin-left: 10px;", uiOutput("dldgemavolpng")),
              ),
              br(),
              fluidRow(uiOutput("mytableHeader")),
              fluidRow(DT::dataTableOutput("mytable")),
              br(),
              fluidRow(uiOutput("downloadall")),
              br()
            )
          ),
          tabPanelBody(
            value = "gse",
            fluidPage(
              shinyjs::useShinyjs(),
              div(
                style = "border: 1px solid #eee; padding: 20px;",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-top: 0px; margin-left: 10px;", 
                    h3("Gene set enrichment options", style = "margin-top: 0px;")
                  ),
                  div(
                    style = "display: inline-block; float: right; margin-right: 10px;",
                    actionLink("gse_help_button", label = "Help"),
                  )
                ),
                fluidRow(
                  div(
                    style = "display: inline-block; margin-left: 10px; vertical-align: top;",
                    uiOutput("dge_list")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("dgemaincontrasts_gse")
                  ),
                  div(style = "display: inline-block;", uiOutput("geneOptions_ui")),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("gene_list_filtering")
                  ),
                  div(style = "display: inline-block;", uiOutput("gene_list"))
                ),
                fluidRow(div(style = "margin-left: 10px;", uiOutput("enrichRLib"))),
              ),
              br(),
              div(uiOutput("dgeInfoSubmitWarning_ui")),
              hr(),
              uiOutput("tbl_func_enrHeader"),
              DT::dataTableOutput("func_enr_tbl"),
              br(),
              br(),
              textOutput("no_list"),
              br()
            )
          ),
          tabPanelBody(
            value = "heatmap",
            fluidPage(
              div(
                style = "border: 1px solid #eee; padding: 20px;",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-left: 10px;",
                    h3("Heatmap options", style = "margin-top: 0px;")
                  ),
                  div(
                    style = "display: inline-block; float: right; margin-right: 10px;",
                    actionLink("heatmap_help_button", label = "Help"),
                  )
                ),
                fluidRow(
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("dge_list_heatmap")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("dgemaincontrasts_heatmap")
                  ),
                  div(style = "display: inline-block;", uiOutput("geneOptions_heatmap_ui")),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("gene_list_filtering_heatmap")
                  ),
                  div(style = "display: inline-block;", uiOutput("gene_list_heatmap")),
                  div(style = "display: inline-block;", uiOutput("filterFactor_heatmap")),
                  div(style = "display: inline-block;", uiOutput("filterFactorSelection_heatmap"))
                ),
                fluidRow(
                  # style = "background-color: #eee;",
                  div(
                    style = "display: inline-block; vertical-align: top; width: 180px; margin-left: 10px;",
                    uiOutput("heatColors")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top; width: 300px; margin-left: 30px;",
                    fluidPage(
                      fluidRow(uiOutput("heatfactorlabel")),
                      fluidRow(uiOutput("contrast_groups_only"))
                    )
                  ),
                  div(
                    style = "display: inline-block; width: 100px; margin-top: 30px; margin-left: 30px;",
                    fluidPage(
                      fluidRow(uiOutput("row_type")),
                      fluidRow(uiOutput("scale_genes"))
                    )
                  ),
                  div(
                    style = "display: inline-block; width: 100px; vertical-align: top; margin-top: 30px; margin-left: 30px;",
                    uiOutput("col_type")
                  )
                )
              ),
              br(),
              div(uiOutput("dgeInfoSubmitWarning_heatmap_ui")),
              uiOutput("headheat"),
              withSpinner(
                uiOutput("heatplot1_ui"),
                type = 8, size = 1, color = "black", proxy.height = "300px"
              ),
              br(),
              div(
                div(style = "display: inline-block;", uiOutput("dlqcheatplot1pdf")),
                div(style = "display: inline-block;", uiOutput("dlqcheatplot1png")),
              ),
              br(),
              br(),
              br(),
              uiOutput("heatfactor"),
              plotlyOutput("heatplot2", width = 800, height = 600),
              br(),
              div(
                div(style = "display: inline-block;", uiOutput("dlqcheatplot2pdf")),
                div(style = "display: inline-block;", uiOutput("dlqcheatplot2png"))
              ),
              br()
            )
          ),
          tabPanelBody(
            value = "clustering",
            fluidPage(
              div(
                style = "border: 1px solid #eee; padding: 20px;",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-left: 10px;",
                    h3("Clustering options", style = "margin-top: 0px;")
                  ),
                  div(
                    style = "display: inline-block; float: right; margin-right: 10px;",
                    actionLink(inputId = "postdge_clustering_help_button", label = "Help"),
                  )
                ),
                fluidRow(
                  div(style = "display: inline-block; width: 140px; margin-left: 10px;", uiOutput("clustalg")),
                  div(
                    style = "display: inline-block; width: 100px; margin-left: 40px; vertical-align: top;", 
                    uiOutput("clustvarnumber")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("min_module_size")
                  ),
                  div(style = "display: inline-block; margin-left: 40px; margin-top: 25px; vertical-align: top;", uiOutput("goclust"))
                ),
                div(style = "color: red; margin-left: 10px;", textOutput("wgcna_warning_msg"))
              ),
              br()
            ),
            tabsetPanel(
              id = "clustering_tabsetPanel",
              type = "hidden",
              tabPanelBody(
                value = "clustPlotW01",
                fluidPage(
                  uiOutput("headclustplotW01"),
                  uiOutput("clustplotW01_ui"),
                  div(
                    div(style = "display:inline-block", uiOutput("downloadclustplotW01pdf")),
                    div(style = "display:inline-block", uiOutput("downloadclustplotW01png")),
                    div(style = "display:inline-block", uiOutput("downloadclustmodK"))
                  ),
                  br()
                )
              ),
              tabPanelBody(
                value = "clustPlotW02",
                fluidPage(
                  uiOutput("headclustplotW02"),
                  plotOutput("clustplotW02", width = 900, height = 400),
                  br(),
                  div(
                    div(
                      style = "display:inline-block",
                      uiOutput("downloadclustplotW02pdf")
                    ),
                    div(
                      style = "display:inline-block",
                      uiOutput("downloadclustplotW02png")
                    ),
                    # Moving this here temporarily
                    div(
                      style = "display:inline-block",
                      uiOutput("downloadclustmod")
                    )
                  ),
                  br()
                )
              ),
              tabPanelBody(
                value = "clustPlotW03",
                fluidPage(
                  uiOutput("headclustplotW03"),
                  plotOutput("clustplotW03", width = 700, height = 700),
                  br(),
                  div(
                    div(style = "display:inline-block", uiOutput("downloadclustplotW03pdf")),
                    div(style = "display:inline-block", uiOutput("downloadclustplotW03png"))
                  ),
                  br()
                )
              )
            )
          )
        )
      )
    )
  })

  # Nested sidebar actionLinks for single-cell
  output$sc_showSelectResolution <- renderUI({
    req(input$sc_dge_tabsetPanel, d$resType)
    if(d$resType != "seurat_res") return()
    myLab <- "Select resolution"
    if(input$sc_dge_tabsetPanel == "sc_dge_selectResolution") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showSelectResolution", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_showOverview <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Overview"
    if(input$sc_dge_tabsetPanel == "sc_dge_overview") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showOverview", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_overviewLinks <- renderUI({
    req(d$sc_showOverviewLinks, input$sc_dge_tabsetPanel, input$sc_dge_overview_tabsetPanel)
    if(input$sc_dge_tabsetPanel != "sc_dge_overview") return()
    clusterLab <- "Clusters"
    metadataLab <- "Metadata"
    if(input$sc_dge_overview_tabsetPanel == "sc_dge_overview_clusters") clusterLab <- strong(clusterLab)
    if(input$sc_dge_overview_tabsetPanel == "sc_dge_overview_metadata") metadataLab <- strong(metadataLab)
    
    fluidPage(
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "sc_dge_overview_showClusters", 
          label = clusterLab, 
          style = "color: black;"
        )
      ),
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "sc_dge_overview_showMetadata", 
          label = metadataLab, 
          style = "color: black;"
        )
      )
    ) 
  })
  
  output$sc_showGeneExpression <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Gene expression"
    if(input$sc_dge_tabsetPanel == "sc_dge_geneExpression") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showGeneExpression", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_geneExpressionLinks <- renderUI({
    req(d$sc_showGeneExpressionLinks, input$sc_dge_tabsetPanel, input$sc_geneExpression_tabsetPanel)
    if(input$sc_dge_tabsetPanel != "sc_dge_geneExpression") return()
    clusterLab <- "By cluster"
    cellLab <- "Cell distribution"
    if(input$sc_geneExpression_tabsetPanel == "sc_dge_geneExpression_cluster") clusterLab <- strong(clusterLab)
    if(input$sc_geneExpression_tabsetPanel == "sc_dge_geneExpression_cell") cellLab <- strong(cellLab)
    
    fluidPage(
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "sc_dge_geneExpression_showCluster", 
          label = clusterLab, 
          style = "color: black;"
        )
      ),
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "sc_dge_geneExpression_showCell", 
          label = cellLab, 
          style = "color: black;"
        )
      )
    ) 
  })
  
  output$sc_showDGE <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "DGE"
    if(input$sc_dge_tabsetPanel == "sc_dge_DGE") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showDGE", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_DGELinks <- renderUI({
    req(d$sc_showDGELinks, input$sc_dge_tabsetPanel, input$sc_dge_DGE_tabsetPanel)
    if(input$sc_dge_tabsetPanel != "sc_dge_DGE") return()
    clusterLab <- "DGE by cluster"
    customLab <- "Custom DGE"
    if(input$sc_dge_DGE_tabsetPanel == "sc_dge_DGE_cluster") clusterLab <- strong(clusterLab)
    if(input$sc_dge_DGE_tabsetPanel == "sc_dge_DGE_custom") customLab <- strong(customLab)
    
    fluidPage(
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "sc_dge_DGE_showCluster", 
          label = clusterLab, 
          style = "color: black;"
        )
      ),
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "sc_dge_DGE_showCustom", 
          label = customLab, 
          style = "color: black;"
        )
      )
    ) 
  })
  
  output$sc_showVolcano <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Volcano plot"
    if(input$sc_dge_tabsetPanel == "sc_dge_volcano") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showVolcano", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_showGSE <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Gene set enrichment"
    if(input$sc_dge_tabsetPanel == "sc_dge_gse") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showGSE", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_showHeatmap <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Heatmap"
    if(input$sc_dge_tabsetPanel == "sc_dge_heatmap") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showHeatmap", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_showManuallySelectCells <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Manually select cells"
    if(input$sc_dge_tabsetPanel == "sc_dge_manuallySelectCells") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showManuallySelectCells", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_showClustering <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Clustering"
    if(input$sc_dge_tabsetPanel == "sc_dge_clustering") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showClustering", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_clusteringLinks <- renderUI({
    req(input$clustalg_sc, clustout_sc(), d$sc_showClusteringLinks, input$sc_dge_tabsetPanel, 
        input$sc_clustering_tabsetPanel)
    if(input$sc_dge_tabsetPanel != "sc_dge_clustering") return()
    if(length(clustout_sc()) != 2 & length(clustout_sc()) != 8) return()
    if(length(clustout_sc()) == 8) {
      lab2 <- "WGCNA - dendrogram"
      lab3 <- "WGCNA - TOM"
      if(input$sc_clustering_tabsetPanel == "sc_clustPlotW02") lab2 <- strong(lab2)
      if(input$sc_clustering_tabsetPanel == "sc_clustPlotW03") lab3 <- strong(lab3)
      
      fluidPage(
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "sc_showClustPlotW02", 
            label = lab2,
            style = "color: black;"
          )
        ),
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "sc_showClustPlotW03", 
            label = lab3,
            style = "color: black;"
          )
        )
      )
    } else {
      fluidPage(
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "sc_showClustPlotW01", 
            label = strong("K-medoids"),
            style = "color: black;"
          )
        )
      )
    }
  })
  
  output$sc_showCombineClusters <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Merge clusters"
    if(input$sc_dge_tabsetPanel == "sc_dge_combineClusters") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showCombineClusters", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  output$sc_showGroupCells <- renderUI({
    req(input$sc_dge_tabsetPanel)
    myLab <- "Group cells by gene expression"
    if(input$sc_dge_tabsetPanel == "sc_dge_GroupCells") myLab <- strong(myLab)
    actionLink(
      inputId = "sc_showGroupCells", 
      label = myLab, 
      style = "font-size: 16px; color: black;"
    )
  })
  
  observeEvent(input$sc_showSelectResolution, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_selectResolution")
    # Reset d$inD and d$res to prevent errors after running DGE on 
    # Manually select cells page. 
    if(!is.null(d$res) && grepl(pattern = "^Comp:", x = d$res)) {
      if(d$resType == "seurat_res") {
        if(!is.null(input$overall_res)) {
          resTemp <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
          if(!(resTemp %in% names(seurat_only()))) resTemp <- names(seurat_only())[1]
          d$inD <- seurat_only()[[resTemp]]
          resTemp <- gsub(":.*$", replacement = "", x = input$overall_res)
          if(!(resTemp %in% names(d$SCV))) resTemp <- names(d$SCV)[!(grepl("^Comp:", names(d$SCV)))][1]
          d$res <- resTemp
        } else {
          resTemp <- names(seurat_only())[1] 
          d$inD <- seurat_only()[[resTemp]]
          resTemp <- names(d$SCV)[!(grepl("^Comp:", names(d$SCV)))][1]
          d$res <- resTemp
        }
      } else {
        d$inD <- seurat_only()
        d$res <- "Clust"
      }
    }
  })
  
  observeEvent(input$sc_showOverview, {
    # Reset d$inD and d$res to prevent errors after running DGE on 
    # Manually select cells page. 
    if(!is.null(d$res) && grepl(pattern = "^Comp:", x = d$res)) {
      if(d$resType == "seurat_res") {
        if(!is.null(input$overall_res)) {
          resTemp <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
          if(!(resTemp %in% names(seurat_only()))) resTemp <- names(seurat_only())[1]
          d$inD <- seurat_only()[[resTemp]]
          resTemp <- gsub(":.*$", replacement = "", x = input$overall_res)
          if(!(resTemp %in% names(d$SCV))) resTemp <- names(d$SCV)[!(grepl("^Comp:", names(d$SCV)))][1]
          d$res <- resTemp
        } else {
          resTemp <- names(seurat_only())[1] 
          d$inD <- seurat_only()[[resTemp]]
          resTemp <- names(d$SCV)[!(grepl("^Comp:", names(d$SCV)))][1]
          d$res <- resTemp
        }
      } else {
        d$inD <- seurat_only()
        d$res <- "Clust"
      }
    }
    
    d$sc_showOverviewLinks <- ifelse(
      test = (input$sc_dge_tabsetPanel == "sc_dge_overview"), 
      yes = ifelse(
        test = d$sc_showOverviewLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_overview")
  })
  
  observeEvent(input$sc_dge_overview_showClusters, {
    updateTabsetPanel(inputId = "sc_dge_overview_tabsetPanel", selected = "sc_dge_overview_clusters")
  })
  
  observeEvent(input$sc_dge_overview_showMetadata, {
    updateTabsetPanel(inputId = "sc_dge_overview_tabsetPanel", selected = "sc_dge_overview_metadata")
  })
  
  observeEvent(input$sc_showGeneExpression, {
    d$sc_showGeneExpressionLinks <- ifelse(
      test = (input$sc_dge_tabsetPanel == "sc_dge_geneExpression"), 
      yes = ifelse(
        test = d$sc_showGeneExpressionLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_geneExpression")
  })
  
  observeEvent(input$sc_dge_geneExpression_showCluster, {
    updateTabsetPanel(inputId = "sc_geneExpression_tabsetPanel", selected = "sc_dge_geneExpression_cluster")
  })
  
  observeEvent(input$sc_dge_geneExpression_showCell, {
    updateTabsetPanel(inputId = "sc_geneExpression_tabsetPanel", selected = "sc_dge_geneExpression_cell")
  })
  
  observeEvent(input$sc_showDGE, {
    
    # Reset d$inD and d$res to prevent errors after running DGE on 
    # Manually select cells page. 
    if(!is.null(d$res) && grepl(pattern = "^Comp:", x = d$res)) {
      if(d$resType == "seurat_res") {
        if(!is.null(input$overall_res)) {
          resTemp <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
          if(!(resTemp %in% names(seurat_only()))) resTemp <- names(seurat_only())[1]
          d$inD <- seurat_only()[[resTemp]]
          resTemp <- gsub(":.*$", replacement = "", x = input$overall_res)
          if(!(resTemp %in% names(d$SCV))) resTemp <- names(d$SCV)[!(grepl("^Comp:", names(d$SCV)))][1]
          d$res <- resTemp
        } else {
          resTemp <- names(seurat_only())[1] 
          d$inD <- seurat_only()[[resTemp]]
          resTemp <- names(d$SCV)[!(grepl("^Comp:", names(d$SCV)))][1]
          d$res <- resTemp
        }
      } else {
        d$inD <- seurat_only()
        d$res <- "Clust"
      }
    }
    
    d$sc_showDGELinks <- ifelse(
      test = (input$sc_dge_tabsetPanel == "sc_dge_DGE"), 
      yes = ifelse(
        test = d$sc_showDGELinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_DGE")
  })
  
  observeEvent(input$sc_dge_DGE_showCluster, {
    updateTabsetPanel(inputId = "sc_dge_DGE_tabsetPanel", selected = "sc_dge_DGE_cluster")
  })
  
  observeEvent(input$sc_dge_DGE_showCustom, {
    updateTabsetPanel(inputId = "sc_dge_DGE_tabsetPanel", selected = "sc_dge_DGE_custom")
  })
  
  observeEvent(input$sc_showVolcano, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_volcano")
  })
  
  observeEvent(input$sc_showGSE, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_gse")
  })
  
  observeEvent(input$sc_showHeatmap, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_heatmap")
  })
  
  observeEvent(input$sc_showManuallySelectCells, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_manuallySelectCells")
    if(!is.null(clustList()) && length(clustList()) > 0 && 
      !is.null(input$comps_manuallyselectcells_sc) && any(grepl("^Comp:", x = clustList()))) {
      d$res <- input$comps_manuallyselectcells_sc
    }
  })
  
  observeEvent(input$sc_showClustering, {
    d$sc_showClusteringLinks <- ifelse(
      test = (input$sc_dge_tabsetPanel == "sc_dge_clustering"), 
      yes = ifelse(
        test = d$sc_showClusteringLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_clustering")
  })
  
  observeEvent(input$sc_showClustPlotW01, {
    updateTabsetPanel(inputId = "sc_clustering_tabsetPanel", selected = "sc_clustPlotW01")
  })
  
  observeEvent(input$sc_showClustPlotW02, {
    updateTabsetPanel(inputId = "sc_clustering_tabsetPanel", selected = "sc_clustPlotW02")
  })
  
  observeEvent(input$sc_showClustPlotW03, {
    updateTabsetPanel(inputId = "sc_clustering_tabsetPanel", selected = "sc_clustPlotW03")
  })
  
  observeEvent(input$sc_showCombineClusters, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_combineClusters")
  })
  
  observeEvent(input$sc_showGroupCells, {
    updateTabsetPanel(inputId = "sc_dge_tabsetPanel", selected = "sc_dge_GroupCells")
  })
  
  # sc rna-seq uiOutput for DGE Analysis tab
  output$dge_sc <- renderUI({
    if(is.null(SubmitData$goqc) || SubmitData$goqc == 0) return()
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        width = 2,
        style = "width: 190px;",
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showSelectResolution"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showOverview"))),
        fluidRow(div(style = "margin-left: 4px;", uiOutput("sc_overviewLinks"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showGeneExpression"))),
        fluidRow(div(style = "margin-left: 4px;", uiOutput("sc_geneExpressionLinks"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showDGE"))),
        fluidRow(div(style = "margin-left: 4px;", uiOutput("sc_DGELinks"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showVolcano"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showHeatmap"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showGSE"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showManuallySelectCells"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showClustering"))),
        fluidRow(div(style = "margin-left: 4px;", uiOutput("sc_clusteringLinks"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showCombineClusters"))),
        fluidRow(div(style = "margin-left: 10px;", uiOutput("sc_showGroupCells")))
      ),
      mainPanel = mainPanel(
        width = 10,
        tabsetPanel(
          id = "sc_dge_tabsetPanel",
          type = "hidden",
          selected = ifelse(
            test = (!is.null(d$resType) && d$resType == "seurat_res"),
            yes = "sc_dge_selectResolution",
            no = "sc_dge_overview"
          ),
          tabPanelBody(
            value = "sc_dge_selectResolution",
            fluidPage(
              fluidRow(
                div(
                  style = "display: inline-block; margin-left: 10px;",
                  h3("Identify optimal resolution for clustering", style = "margin-top: 0px;")
                ),
                div(
                  style = "display: inline-block; float: right;",
                  actionLink("sc_resolution_help_button", label = "Help")
                )
              ),
              fluidRow(hr()),
              fluidRow(
                div(
                  style = "display: inline-block; width: 240px; margin-left: 10px;",
                  uiOutput("overall_res")
                ),
                div(
                  style = "display: inline-block; width: 100px; margin-left: 20px;",
                  uiOutput("red")
                ),
                div(
                  style = "display: inline-block; width: 200px; margin-left: 20px;",
                  uiOutput("select_res_factor")
                )
              ),
              fluidRow(
                div(
                  style = "margin-left: 10px;",
                  plotOutput("res_umap", height = 500, width = 500)
                )
              ),
              fluidRow(hr()),
              fluidRow(
                div(
                  style = "margin-left: 10px;", 
                  plotOutput("clustree_png", height = 600, width = 750)
                )
              ),
              br()
            )
          ),
          tabPanelBody(
            value = "sc_dge_overview",
            div(
              style = "display: inline-block; float: right;",
              actionLink(inputId = "sc_overview_help_button", label = "Help")
            ),
            tabsetPanel(
              id = "sc_dge_overview_tabsetPanel",
              type = "hidden",
              tabPanelBody(
                value = "sc_dge_overview_clusters",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-left: 10px;",
                    h3("Cluster exploration", style = "margin-top: 0px;")
                  )
                ),
                # hr(),
                fluidPage(
                  fluidRow(
                    column(
                      width = 6,
                      style = "width: 480px; height: 620px; border: 1px solid #eee; padding-left: 20px; padding-bottom: 20px;",
                      fluidRow(h4("Cluster separation")),
                      fluidRow(
                        div(
                          style = "display: inline-block; width: 340px;",
                          uiOutput("deType")
                        ),
                        div(
                          style = "display: inline-block; width: 100px; margin-left: 10px; vertical-align: top;",
                          uiOutput("FDRthresh1")
                        )
                      ),
                      fluidRow(
                        plotOutput(
                          outputId = "clustSep", 
                          height = "450px", 
                          # width = "500px",
                          width = "400px",
                          dblclick = "clustSep_dblclick",
                          brush = brushOpts(id = "clustSep_brush", resetOnNew = T)
                        )
                      ),
                      fluidRow(
                        downloadButton("clustSepSave", label = "Save as PNG")
                      )
                    ),
                    column(
                      width = 5,
                      offset = 1,
                      style = "width: 480px; height: 620px; border: 1px solid #eee; margin-left: 20px; padding-left: 20px; padding-bottom: 20px;",
                      fluidRow(h4("Silhouette plot")),
                      # fluidRow(plotOutput("sil", height = "600px", width = "450px")),
                      fluidRow(plotOutput("sil", height = "508px", width = "400px")),
                      br(),
                      fluidRow(downloadButton("silSave", label = "Save as PNG"))
                    )
                  ),
                  br()
                )
              ),
              tabPanelBody(
                value = "sc_dge_overview_metadata",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-left: 10px;",
                    h3("Metadata exploration", style = "margin-top: 0px;")
                  )
                ),
                fluidPage(
                  fluidRow(
                    column(
                      width = 6,
                      style = "width: 480px; height: 620px; border: 1px solid #eee; padding-left: 20px; padding-bottom: 20px;",
                      fluidRow(h4("Metadata relationships")),
                      fluidRow(
                        div(
                          style = "display: inline-block; width: 200px;", 
                          uiOutput("mdScatterX")
                        ),
                        div(
                          style = "display: inline-block; width: 200px; margin-left: 10px;", 
                          uiOutput("mdScatterY")
                        )
                      ),
                      fluidRow(
                        div(
                          style = "display: inline-block;", uiOutput("scatterLogX")
                        ),
                        div(
                          style = "display: inline-block;", uiOutput("scatterLogY")
                        )
                      ),
                      fluidRow(plotOutput("mdScatter", height = "400px", width = "400px")),
                      br(),
                      fluidRow(
                        downloadButton("mdScatterSave", label = "Save as PNG", style = "margin-left: 10px;")
                      )
                    ),
                    column(
                      width = 5,
                      offset = 1,
                      style = "width: 480px; height: 620px; border: 1px solid #eee; margin-left: 20px; padding-left: 20px; padding-bottom: 20px;",
                      fluidRow(h4("Metadata by cluster")),
                      fluidRow(
                        div(
                          style = "display: inline-block; width: 200px;", 
                          uiOutput("mdFactorData")
                        ),
                        div(
                          style = "display: inline-block; vertical-align: top;",
                          uiOutput("mdFactorOptsF")
                        ),
                        div(
                          style = "display: inline-block; vertical-align: top; margin-top: 30px;",
                          uiOutput("mdFactorOptsN")
                        )
                      ),
                      fluidRow(plotOutput("mdFactor", height = "400px", width = "400px")),
                      br(),
                      fluidRow(
                        downloadButton("mdFactorSave", label = "Save as PNG", style = "margin-left: 10px;")
                      )
                    )
                  ),
                  br()
                )
              )
            )
          ),
          tabPanelBody(
            value = "sc_dge_geneExpression",
            div(
              style = "display: inline-block; float: right;",
              actionLink("sc_geneexpression_help_button", label = "Help"),
            ),
            tabsetPanel(
              id = "sc_geneExpression_tabsetPanel",
              type = "hidden",
              tabPanelBody(
                value = "sc_dge_geneExpression_cluster",
                fluidPage(
                  fluidRow(h3("Gene expression by cluster", style = "margin-top: 0px;")),
                  fluidRow(
                    div(style = "display: inline-block; width: 150px;", uiOutput("cgSelect")),
                    div(
                      style = "display: inline-block; vertical-align: middle; margin-left: 20px; width: 100px;",
                      uiOutput("geneExpIncludeJitter"),
                    ),
                    div(
                      style = "display: inline-block; vertical-align: middle; margin-left: 20px; width: 100px;",
                      uiOutput("geneExpIncludeDR"),
                    )
                  ),
                  fluidRow(plotOutput("geneTest", height = "600px", width = "600px")),
                  fluidRow(uiOutput("geneTestSaveButton")),
                  br()
                )
              ),
              tabPanelBody(
                value = "sc_dge_geneExpression_cell",
                fluidPage(
                  fluidRow(h3("Cell distribution of genes of interest", style = "margin-top: 0px;")),
                  fluidRow(
                    div(style = "display: inline-block; width: 150px;", uiOutput("GOI1select")),
                    div(style = "display: inline-block; width: 150px; margin-left: 20px;", uiOutput("GOI_EmbType")),
                    div(style = "display: inline-block; width: 150px; margin-left: 20px;", uiOutput("GOI_EmbDimX")),
                    div(style = "display: inline-block; width: 150px; margin-left: 20px;", uiOutput("GOI_EmbDimY"))
                  ),
                  fluidRow(
                    div(style = "display: inline-block; width: 200px;", uiOutput("plotClust1")),
                    div(
                      style = "display: inline-block; vertical-align: top; margin-top: 30px; margin-left: 20px;", 
                      uiOutput("plotLabel1")
                    )
                  ),
                  fluidRow(plotOutput("goiPlot1", height = "600px", width = "600px")),
                  fluidRow(uiOutput("goiPlot1SaveButton")),
                  br()
                )
              )
            )
          ),
          tabPanelBody(
            value = "sc_dge_DGE",
            tabsetPanel(
              id = "sc_dge_DGE_tabsetPanel",
              type = "hidden",
              tabPanelBody(
                value = "sc_dge_DGE_cluster",
                fluidPage(
                  fluidRow(
                    h3("Differential gene expression by cluster", style = "display: inline-block; margin-top: 0px;"),
                    div(
                      style = "display: inline-block; vertical-align: top; float: right;",
                      actionLink("sc_dgeanalysis_help_button", label = "Help")
                    )
                  ),
                  fluidRow(
                    div(style = "display: inline-block; width: 100px; vertical-align: top;", uiOutput("DEclustNum")),
                    div(
                      style = "display: inline-block; width: 100px; margin-left: 20px;",
                      uiOutput("FDRthresh2")
                    ),
                    div(
                      style = "display: inline-block; width: 150px; margin-left: 20px; vertical-align: top;",
                      uiOutput("heatDEtype")
                    ),
                    div(
                      style = "display: inline-block; width: 200px; margin-left: 20px; vertical-align: top;",
                      uiOutput("DEgeneSlider")
                    )
                  ),
                  fluidRow(h4("Dotplot")),
                  fluidRow(uiOutput("dotplotWarning")),
                  fluidRow(plotOutput("dotplot", height = "500px", width = "1000px")),
                  fluidRow(
                    div(
                      style = "display: inline-block;", 
                      downloadButton("CGSsave0", label = "Download cluster gene stats")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 20px;",
                      downloadButton("deGeneSave", label = "Download DE results")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 20px;",
                      uiOutput("heatmapSave")
                    )
                  ),
                  fluidRow(hr()),
                  fluidRow(h3("DGE table", style = "margin-top: 0px;")),
                  fluidRow(h4("Filters")),
                  fluidRow(
                    div(
                      style = "display: inline-block; width: 180px;", 
                      uiOutput("sc_dgeByCluster_table_lfc")
                    ),
                    div(
                      style = "display: inline-block; width: 150px; margin-left: 20px;",
                      uiOutput("sc_dgeByCluster_table_padj")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 20px;", 
                      uiOutput("sc_dgeByCluster_table_load")
                    )
                  ),
                  fluidRow(uiOutput("dge_cluster_sc")),
                  br()
                )
              ),
              tabPanelBody(
                value = "sc_dge_DGE_custom",
                fluidPage(
                  fluidRow(
                    h3("Differential gene expression options", style = "display: inline-block; margin-top: 0px;"),
                    div(
                      style = "display: inline-block; float: right;", 
                      actionLink("sc_dge_factor_help_button", label = "Help")
                    )
                  ),
                  fluidRow(textOutput("sc_dge_byfactor_novars_msg"),style="color:red"),
                  fluidRow(
                    div(style = "display: inline-block; width: 150px;", uiOutput("sc_dge_factor_fact")),
                    div(
                      style = "display: inline-block; width: 150px; margin-left: 20px;", 
                      uiOutput("sc_dge_factor_group1")
                    ),
                    div(
                      style = "display: inline-block; width: 150px; margin-left: 20px;", 
                      uiOutput("sc_dge_factor_group2")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 20px; vertical-align: top; margin-top: 23px;", 
                      uiOutput("sc_dge_factor_submit")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 20px; vertical-align: top; margin-top: 23px;", 
                      uiOutput("sc_dge_factor_clear")
                    )
                  ),
                  fluidRow(uiOutput("sc_dge_factor_viewOpts")),
                  uiOutput("sc_dge_factor_opts"),
                  fluidRow(textOutput("sc_dge_byfactor_msg"),style="color:red"),
                  fluidRow(hr()),
                  fluidRow(uiOutput("sc_dge_factor_results_header")),
                  fluidRow(
                    div(style = "display: inline-block; width: 200px;", uiOutput("sc_dge_factor_tablefact")),
                    div(style = "display: inline-block; width: 400px; margin-left: 20px;", uiOutput("sc_dge_factor_comp"))
                  ),
                  fluidRow(uiOutput("sc_dge_factor_results")),
                  br()
                )
              )
            )
          ),
          tabPanelBody(
            value = "sc_dge_volcano",
            fluidPage(
              fluidRow(
                h3("Volcano plot", style = "display: inline-block; margin-top: 0px;"),
                div(
                  style = "display: inline-block; float: right;", 
                  actionLink("sc_volcanoplots_help_button", label = "Help")
                )
              ),
              fluidRow(
                div(
                  style = "display: inline-block; width: 220px; vertical-align: top;", 
                  uiOutput("scatterInput")
                ),
                div(
                  style = "display: inline-block; width: 180px; margin-left: 20px; vertical-align: top;", 
                  uiOutput("diffLabelType")
                ),
                div(
                  style = "display: inline-block; width: 140px; margin-left: 20px;",
                  fluidPage(
                    fluidRow(div(style = "width: 120px;", uiOutput("setScatterA"))),
                    fluidRow(div(style = "width: 140px;", downloadButton("CGSsaveA", label = "Download stats")))
                  )
                ),
                div(
                  style = "display: inline-block; width: 140px; margin-left: 20px;",
                  fluidPage(
                    fluidRow(div(style = "width: 120px;", uiOutput("setScatterB"))),
                    fluidRow(div(style = "width: 140px;", downloadButton("CGSsaveB", label = "Download stats")))
                  )
                ),
                div(style = "display: inline-block; width: 80px; margin-left: 20px; vertical-align: top;", uiOutput("diffLabelSelect"))
              ),
              br(),
              fluidRow(plotOutput("setScatter", width = "600px", height = "500px", click = "scatterClick")),
              fluidRow(
                div(style = "display: inline-block;", downloadButton("setScatterSave", label = "Save as PNG")),
                div(
                  style = "display: inline-block; margin-left: 20px;", 
                  downloadButton("setComparisonSave", label = "Download DE results")
                )
              ),
              fluidRow(hr()),
              fluidRow(h3("Differential gene expression table", style = "margin-top: 0px;")),
              fluidRow(h4("Filters")),
              fluidRow(
                div(
                  style = "display: inline-block; width: 180px;", 
                  uiOutput("vol_cluster_sc_lfc")
                ),
                div(
                  style = "display: inline-block; width: 150px; margin-left: 20px;",
                  uiOutput("vol_cluster_sc_padj")
                ),
                div(
                  style = "display: inline-block; margin-left: 20px;", 
                  uiOutput("vol_cluster_sc_load")
                )
              ),
              fluidRow(uiOutput("vol_cluster_sc")),
              br()
            )
          ),
          tabPanelBody(
            value = "sc_dge_gse",
            fluidPage(
              shinyjs::useShinyjs(),
              div(
                style = "border: 1px solid #eee; padding: 20px;",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-top: 0px; margin-left: 10px;", 
                    h3("Gene set enrichment options", style = "margin-top: 0px;")
                  ),
                  div(
                    style = "display: inline-block; float: right; margin-right: 10px;",
                    actionLink("sc_gse_help_button", label = "Help")
                  )
                ),
                fluidRow(
                  div(
                    style = "display: inline-block; vertical-align: top; margin-left: 10px;",
                    uiOutput("dge_list_sc")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("gse_customdge_fact")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("gse_customdge_contrast")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("sc_contrasts")
                  ),
                  div(style = "display: inline-block;", uiOutput("geneOptions_sc_ui")),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("sc_gene_list_filtering")
                  ),
                  div(style = "display: inline-block;", uiOutput("gene_list_sc"))
                ),
                fluidRow(div(style = "margin-left: 10px;", uiOutput("enrichRLib_sc"))),
              ),
              br(),
              div(
                div(style = "display: inline-block;", uiOutput("dge_info_sc")),
                div(style = "display: inline-block;", uiOutput("submit_fe_sc"))
              ),
              div(textOutput("gse_warning_sc"), style = "color: red;"),
              hr(),
              uiOutput("sc_tbl_func_enrHeader"),
              DT::dataTableOutput("sc_func_enr_tbl"),
              br(),
              br(),
              textOutput("sc_no_list"),
              br()
            )
          ),
          tabPanelBody(
            value = "sc_dge_heatmap",
            fluidPage(
              div(
                style = "border: 1px solid #eee; padding: 20px;",
                fluidRow(
                  h3("Heatmap options", style = "display: inline-block; margin-left: 10px; margin-top: 0px;"),
                  div(
                    style = "display: inline-block; float: right; margin-right: 10px;",
                    actionLink("sc_heatmap_help_button", label = "Help"),
                  )
                ),
                fluidRow(
                  div(
                    style = "display: inline-block; vertical-align: top; margin-left: 10px;",
                    uiOutput("dge_list_heatmap_sc")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("heatmap_customdge_fact")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("heatmap_customdge_contrast")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("sc_contrasts_heatmap")
                  ),
                  div(style = "display: inline-block;", uiOutput("geneOptions_heatmap_sc_ui")),
                  div(
                    style = "display: inline-block; vertical-align: top;", 
                    uiOutput("sc_gene_list_filtering_heatmap")
                  ),
                  div(style = "display: inline-block;", uiOutput("gene_list_heatmap_sc"))
                ),
                fluidRow(
                  div(
                    style = "display: inline-block; vertical-align: top; width: 180px; margin-left: 10px;",
                    uiOutput("heatColors_sc")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top; width: 300px; margin-left: 30px;",
                    fluidPage(
                      fluidRow(uiOutput("heatfactorlabel_sc")),
                      fluidRow(uiOutput("contrast_groups_only_sc"))
                    )
                  ),
                  div(
                    style = "display: inline-block; width: 100px; margin-top: 30px; margin-left: 30px;",
                    fluidPage(
                      fluidRow(uiOutput("row_type_sc")),
                      fluidRow(uiOutput("scale_genes_sc"))
                    )
                  ),
                  div(
                    style = "display: inline-block; width: 100px; vertical-align: top; margin-top: 30px; margin-left: 30px;",
                    uiOutput("col_type_sc")
                  )
                )
              ),
              br(),
              div(
                div(style = "display: inline-block;", uiOutput("dge_info_heatmap_sc")),
                div(style = "display: inline-block;", uiOutput("submit_list_sc"))
              ),
              div(textOutput("heatmap_warning_sc"), style = "color: red;"),
              hr(),
              uiOutput("headheat_sc"),
              withSpinner(
                uiOutput("heatplot1_sc_ui"),
                type = 8, size = 1, color = "black", proxy.height = "300px"
              ),
              br(),
              div(
                div(style = "display: inline-block;", uiOutput("dlqcheatplot1pdf_sc")),
                div(style = "display: inline-block;", uiOutput("dlqcheatplot1png_sc")),
              ),
              br(),
              br(),
              br(),
              uiOutput("heatfactor_sc"),
              plotlyOutput("heatplot2_sc", width = 800, height = 600),
              br(),
              div(
                div(style = "display:inline-block", uiOutput("dlqcheatplot2pdf_sc")),
                div(style = "display:inline-block", uiOutput("dlqcheatplot2png_sc"))
              ),
              br()
            )
          ),
          tabPanelBody(
            value = "sc_dge_manuallySelectCells",
            fluidPage(
              fluidRow(
                h3(
                  "Manually select cells for DGE", 
                  style = "display: inline-block; margin-top: 0px;"
                ),
                div(
                  style = "display: inline-block; float: right;",
                  actionLink("sc_manuallyselectcells_help_button", label = "Help")
                )
              ),
              fluidRow(
                div(
                  style = "display: inline-block; width: 100px;",
                  uiOutput("SelDE_EmbType")
                ),
                div(
                  style = "display: inline-block; width: 100px; margin-left: 20px;",
                  uiOutput("SelDE_EmbDimX")
                ),
                div(
                  style = "display: inline-block; width: 100px; margin-left: 20px;",
                  uiOutput("SelDE_EmbDimY")
                ),
                div(
                  style = "display: inline-block; width: 220px; margin-left: 20px; vertical-align: top;",
                  uiOutput("tsneSelDEcol")
                ),
                div(
                  style = "display: inline-block; width: 60px; margin-left: 20px; vertical-align: top; margin-top: 23px;",
                  uiOutput("plusFilt")
                ),
                div(
                  style = "display: inline-block; width: 100px; margin-left: 20px; vertical-align: top; margin-top: 23px;",
                  uiOutput("minusFilt")
                ),
                div(
                  style = "display: inline-block; width: 150px; margin-left: 20px; vertical-align: top; margin-top: 23px;",
                  uiOutput("MDfiltsRemoveAll")
                )
              ),
              fluidRow(uiOutput("MDfilts")),
              fluidRow(
                column(
                  width = 6, 
                  style = "width: 460px; border: 1px solid #eee; padding: 10px;",
                  fluidRow(
                    fluidPage(div(htmlOutput("textSetA"), style = "font-size: 16px;"))
                  ),
                  br(),
                  fluidRow(
                    div(
                      style = "display: inline-block; width: 155px;",
                      fluidPage(uiOutput("addCellsA"))
                    ),
                    div(
                      style = "display: inline-block; margin-left: 5px;",# width: 170px;",
                      uiOutput("removeCellsA")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 5px;",# width: 80px;",
                      uiOutput("clearA")
                    )
                  )
                ),
                column(
                  width = 5,
                  offset = 1,
                  style = "width: 460px; margin-left: 30px; border: 1px solid #eee; padding: 10px;",
                  fluidRow(
                    fluidPage(div(htmlOutput("textSetB"), style = "font-size: 16px;"))
                  ),
                  br(),
                  fluidRow(
                    div(
                      style = "display: inline-block; width: 155px;",
                      fluidPage(uiOutput("addCellsB"))
                    ),
                    div(
                      style = "display: inline-block; margin-left: 5px;",# width: 170px;",
                      uiOutput("removeCellsB")
                    ),
                    div(
                      style = "display: inline-block; margin-left: 5px;",# width: 80px;",
                      uiOutput("clearB")
                    )
                  )
                )
              ),
              fluidRow(textOutput("textOverlap"), style = "color: red;"),
              br(),
              fluidRow(
                div(
                  style = "display: inline-block; width: 240px;",
                  uiOutput("DEsetName")
                ),
                div(
                  style = "display: inline-block; width: 300px; margin-left: 20px;",
                  uiOutput("calcDE")
                )
              ),
              fluidRow(textOutput("calcText"), style = "color: red;"),
              fluidRow(textOutput("cellsHovered")),
              fluidRow(plotlyOutput("tsneSelDE", width = "500px", height = "500px")),
              fluidRow(uiOutput("manuallyselectcells_downloadDGE_ui")),
              br(),
              fluidRow(uiOutput("manuallyselectcells_downloadGSE_ui")),
              br()
            )
          ),
          tabPanelBody(
            value = "sc_dge_clustering",
            fluidPage(
              div(
                style = "border: 1px solid #eee; padding: 20px;",
                fluidRow(
                  div(
                    style = "display: inline-block; margin-left: 10px;",
                    h3("Clustering options", style = "margin-top: 0px;")
                  ),
                  div(
                    style = "display: inline-block; float: right; margin-right: 10px;",
                    actionLink(inputId = "sc_postdge_clustering_help_button", label = "Help"),
                  )
                ),
                fluidRow(
                  div(style = "display: inline-block; width: 140px; margin-left: 10px;", uiOutput("clustalg_sc")),
                  div(
                    style = "display: inline-block; width: 100px; margin-left: 40px; vertical-align: top;", 
                    uiOutput("clustvarnumber_sc")
                  ),
                  div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput("min_module_size_sc")
                  ),
                  div(style = "display: inline-block; margin-left: 40px; margin-top: 25px; vertical-align: top;", uiOutput("goclust_sc"))
                ),
                div(style = "color: red; margin-left: 10px;", textOutput("wgcna_warning_msg_sc"))
              ),
              br()
            ),
            tabsetPanel(
              id = "sc_clustering_tabsetPanel",
              type = "hidden",
              tabPanelBody(
                value = "sc_clustPlotW01",
                fluidPage(
                  uiOutput("headclustplotW01_sc"),
                  uiOutput("clustplotW01_sc_ui"),
                  div(
                    div(style = "display:inline-block", uiOutput("downloadclustplotW01pdf_sc")),
                    div(style = "display:inline-block", uiOutput("downloadclustplotW01png_sc")),
                    div(style = "display:inline-block", uiOutput("downloadclustmodK_sc"))
                  ),
                  br()
                )
              ),
              tabPanelBody(
                value = "sc_clustPlotW02",
                fluidPage(
                  uiOutput("headclustplotW02_sc"),
                  plotOutput("clustplotW02_sc", width = 900, height = 400),
                  br(),
                  div(
                    div(
                      style = "display:inline-block",
                      uiOutput("downloadclustplotW02pdf_sc")
                    ),
                    div(
                      style = "display:inline-block",
                      uiOutput("downloadclustplotW02png_sc")
                    ),
                    # Moving this here temporarily
                    div(
                      style = "display:inline-block",
                      uiOutput("downloadclustmod_sc")
                    )
                  ),
                  br()
                )
              ),
              tabPanelBody(
                value = "sc_clustPlotW03",
                fluidPage(
                  uiOutput("headclustplotW03_sc"),
                  plotOutput("clustplotW03_sc", width = 700, height = 700),
                  br(),
                  div(
                    div(style = "display:inline-block", uiOutput("downloadclustplotW03pdf_sc")),
                    div(style = "display:inline-block", uiOutput("downloadclustplotW03png_sc"))
                  ),
                  br()
                )
              )
            )
          ),
          tabPanelBody(
            value = "sc_dge_combineClusters",
            fluidPage(
              div(
                style = "border: 1px solid #eee; padding: 10px;",
                div(
                  h3("Merge clusters", style = "display: inline-block; margin-top: 0px;"),
                  div(
                    style = "display: inline-block; float: right;",
                    actionLink("sc_combineClusters_help_button", label = "Help")
                  )
                ),
                div(
                  div(
                    style = "display: inline-block; width: 200px; vertical-align: top;",
                    uiOutput("clusters_comb")
                  ),
                  div(
                    style = "display: inline-block; width: 100px; margin-left: 20px; vertical-align: top;",
                    uiOutput("username")
                  ),
                  div(
                    style = "display: inline-block; width: 200px; margin-left: 20px;",
                    uiOutput("note")
                  ),
                  div(
                    style = "display: inline-block; width: 200px; margin-left: 20px; vertical-align: top;",
                    uiOutput("table_name")
                  )
                ),
                uiOutput("comb_buttons_ui")
              ),
              fluidRow(plotlyOutput("comb_mdFactor_plotly", height = "600px", width = "600px")),
              br()
            )
          ),
          tabPanelBody(
            value = "sc_dge_GroupCells",
            fluidPage(
              div(
                style = "border: 1px solid #eee; padding: 10px;",
                div(
                  h3("Group cells by gene expression", style = "display: inline-block; margin-top: 0px;"),
                  div(
                    style = "display: inline-block; float: right;",
                    actionLink("sc_groupCells_help_button", label = "Help")
                  )
                ),
                div(
                  div(style = "display: inline-block; width: 200px; vertical-align: top;", uiOutput("gene")),
                  div(style = "display: inline-block; vertical-align: top;", uiOutput("metric")),
                  div(style = "display: inline-block;", uiOutput("varRange_ui"))
                ),
                div(
                  div(
                    style = "display: inline-block; width: 100px; vertical-align: top;",
                    uiOutput("username2")
                  ),
                  div(
                    style = "display: inline-block; width: 200px; margin-left: 20px;",
                    uiOutput("note2")
                  ),
                  div(
                    style = "display: inline-block; width: 200px; margin-left: 20px; vertical-align: top;",
                    uiOutput("table_name2")
                  ),
                  div(
                    style = "display: inline-block; margin-left: 20px; vertical-align: top; margin-top: 23px;", 
                    uiOutput("submit_mysql2")
                  )
                )
              ),
              br(),
              fluidRow(
                div(
                  style = "display: inline-block; width: 540px; height: 440px;", 
                  plotlyOutput("selectClustersByGOI", width = "540px", height = "440px")
                ),
                div(
                  style = "display: inline-block; width: 400px; margin-left: 30px; vertical-align: top;",
                  fluidPage(
                    fluidRow(div(style = "width: 150px;", uiOutput("setName"))),
                    fluidRow(uiOutput("addSet")),
                    br(),
                    fluidRow(uiOutput("addCells"))
                  )
                )
              ),
              br(),
              fluidRow(
                div(style = "display: inline-block;", uiOutput("download_new_dataset_GOI")),
                div(style = "display: inline-block; margin-left: 10px;", uiOutput("download_new_metadata_GOI"))
              ),
              br()
            )
          )
        )
      )
    )
  })
  
  output$dge <- renderUI({
    req(SubmitData$data_type, SubmitData$goqc)
    if(SubmitData$data_type == "Single-cell") {
      if(is.null(d$resType)) return()
      tagList(
        uiOutput("dge_sc")
      )
    } else {
      tagList(
        uiOutput("dge_bulk")
      )
    }
  })
  
  ###################################################################
  ###################################################################
  ### SECTION 01 - DATA LOAD (DATA)
  ###################################################################
  ###################################################################
  
  # LoadData module. SubmitData is a reactiveValues-type object.
  SubmitData <- callModule(module = LoadData, id = "load_data", maxSamples = 10000)
  
  observeEvent(SubmitData$data_type, {
    d$cts_out <- d$selectedData <- d$datasetNames <- NULL
  }, ignoreNULL = F, ignoreInit = T)
  
  observeEvent(SubmitData$select_datasets_table_rows_selected, {
    d$cts_out <- d$selectedData <- d$datasetNames <- NULL
  }, ignoreNULL = F, ignoreInit = T)
  
  observeEvent(SubmitData$selectedData, {
    d$cts_out <- d$selectedData <- d$datasetNames <- NULL
  }, ignoreNULL = F, ignoreInit = T)
  
  observeEvent(SubmitData$goqc, {
    if(is.null(SubmitData$goqc)) return()
    if(SubmitData$goqc == 0) return()
    d$resType <- SubmitData$resType
    d$resValues <- SubmitData$resValues
    d$datasetNames <- SubmitData$datasetNames
    d$drt <- SubmitData$drt
    d$pcaVarianceFraction <- SubmitData$pcaVarianceFraction
  }, priority = 1)
  
  ddsout <- eventReactive(SubmitData$goqc, {
    if(is.null(SubmitData$goqc)) return()
    if(SubmitData$goqc == 0) return()
    SubmitData$ddsout
  })
    
  ddstran <- eventReactive(SubmitData$goqc, {
    if(is.null(SubmitData$goqc)) return()
    if(SubmitData$goqc == 0) return()
    SubmitData$ddstran
  })
  
  seurat_only <- eventReactive(SubmitData$goqc, {
    if(is.null(SubmitData$goqc) || SubmitData$goqc == 0) return()
    SubmitData$seurat_only
  })
  
  seurat_sc <- eventReactive(SubmitData$goqc, {
    if(is.null(SubmitData$goqc) || SubmitData$goqc == 0) return()
    SubmitData$seurat_sc
  })
  
  
  ###################################################################
  ###################################################################
  ### SECTION 02 - QUALITY CONTROL (QC)
  ###################################################################
  ###################################################################
  
  output$showDataSummaryLinks <- renderUI({
    req(input$dataSummaryQCTabset)
    myLab <- "Data summary"
    if(input$dataSummaryQCTabset == "dataSummary") myLab <- strong(myLab)
    actionLink(
      inputId = "showDataSummaryLinks", 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })
  # observeEvents with ifelse to expand/collapse subtabs
  observeEvent(input$showDataSummaryLinks, {
    d$showDataSummaryLinks <- ifelse(
      test = (input$dataSummaryQCTabset == "dataSummary"), 
      yes = ifelse(
        test = d$showDataSummaryLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    d$showQCLinks <- F
    updateTabsetPanel(inputId = "dataSummaryQCTabset", selected = "dataSummary")
  })
  
  output$showQCLinks <- renderUI({
    req(input$dataSummaryQCTabset)
    myLab <- "Quality control"
    if(input$dataSummaryQCTabset == "qc") myLab <- strong(myLab)
    actionLink(
      inputId = "showQCLinks", 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })
  
  observeEvent(input$showQCLinks, {
    d$showQCLinks <- ifelse(
      test = (input$dataSummaryQCTabset == "qc"), 
      yes = ifelse(
        test = d$showQCLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    d$showDataSummaryLinks <- F
    updateTabsetPanel(inputId = "dataSummaryQCTabset", selected = "qc")
  })
  
  # renderUI with subtab links which can be expanded/collapsed
  output$dataSummaryLinks <- renderUI({
    req(d$showDataSummaryLinks, input$dataSummaryQCTabset, input$dataSummary_tabsetPanel)
    if(input$dataSummaryQCTabset != "dataSummary") return()
    
    countDataLab <- "Count data"
    metadataLab <- "Metadata"
    if(input$dataSummary_tabsetPanel == "dataSummary_subsetCountData") countDataLab <- strong(countDataLab)
    if(input$dataSummary_tabsetPanel == "dataSummary_metadata") metadataLab <- strong(metadataLab)

    fluidPage(
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "showSubsetCountData", 
          label = countDataLab,
          style = "color: black;" 
        )
      ),
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "showMetadata", 
          label = metadataLab,
          style = "color: black;" 
        )
      )
    )
  })
  
  output$qcLinks <- renderUI({
    req(d$showQCLinks, input$dataSummaryQCTabset, input$qc_tabsetPanel)
    if(input$dataSummaryQCTabset != "qc") return()
    
    boxAndWhiskerLab <- "Box and whisker"
    histogramLab <- "Histogram"
    totalReadsLab <- "Total reads"
    
    if(input$qc_tabsetPanel == "qc_boxAndWhisker") boxAndWhiskerLab <- strong(boxAndWhiskerLab)
    if(input$qc_tabsetPanel == "qc_histogram") histogramLab <- strong(histogramLab)
    if(input$qc_tabsetPanel == "qc_totalReads") totalReadsLab <- strong(totalReadsLab)
    
    fluidPage(
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "showBoxAndWhisker", 
          label = boxAndWhiskerLab,
          style = "color: black;"
        )
      ),
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "showHistogram", 
          label = histogramLab,
          style = "color: black;"
        )
      ),
      fluidRow(
        style = "margin-left: 6px;", 
        actionLink(
          inputId = "showTotalReads", 
          label = totalReadsLab,
          style = "color: black;"
        )
      )
    )
  })
  
  # observeEvents to update dataSummaryQCTabset tabsetPanel
  observeEvent(input$showSubsetCountData, {
    updateTabsetPanel(inputId = "dataSummary_tabsetPanel", selected = "dataSummary_subsetCountData")
  })
  
  observeEvent(input$showMetadata, {
    updateTabsetPanel(inputId = "dataSummary_tabsetPanel", selected = "dataSummary_metadata")
  })
  
  observeEvent(input$showBoxAndWhisker, {
    updateTabsetPanel(inputId = "qc_tabsetPanel", selected = "qc_boxAndWhisker")
  })
  
  observeEvent(input$showHistogram, {
    updateTabsetPanel(inputId = "qc_tabsetPanel", selected = "qc_histogram")
  })
  
  observeEvent(input$showTotalReads, {
    updateTabsetPanel(inputId = "qc_tabsetPanel", selected = "qc_totalReads")
  })
  
  # Info boxes for pre-/post-filtration gene counts and sample/cell counts
  output$pre <- renderUI({ 
    req(d$cts_out, ddsout())
    fluidPage(
      p(formatC(nrow(ddsout()[[4]]), big.mark = ","), 
        style = "font-size: 35px; margin-bottom: 0%;"),
      p("genes (pre-filtration)", 
        style = "white-space: nowrap; font-size: 20px;")
    )
  })
  
  output$post <- renderUI({ 
    req(d$cts_out, ddsout())
    fluidPage(
      p(formatC(nrow(ddsout()[[3]]), big.mark = ","), 
        style = "font-size: 35px; margin-bottom: 0%;"),
      p("genes (post-filtration)", 
        style = "white-space: nowrap; font-size: 20px;")
    )
  })
  
  output$samplecount <- renderUI({ 
    req(d$cts_out, ddsout(), SubmitData$data_type)
    txt <- "samples"
    if(SubmitData$data_type == "Single-cell") txt <- "cells"
    fluidPage(
      p(formatC(ncol(ddsout()[[3]]), big.mark = ","), 
        style = "font-size: 35px; margin-bottom: 0%;"),
      p(txt, style = "white-space: nowrap; font-size: 20px;")
    )
  })
  
  # QC - header (2) - file summary (count data)
  output$filesummarycts <- renderUI({
    req(d$cts_out)
    h4("Subset count data (first 5 rows and last 5 rows)")
  })
  
  # QC - download button for raw count data
  output$download_raw_counts <- renderUI({
    req(d$cts_out)
    downloadButton("download_raw_counts_button", label = "Download raw counts")
  })
  
  # QC - download file for raw count data
  output$download_raw_counts_button <- downloadHandler(
    filename = function() {
      paste0("Raw_Counts.csv")
    },
    content = function(file) {
      write.csv(
        d$cts_out,
        file = file,
        quote = F
      )
    }
  )
  
  # QC - download button for normalized count data
  output$download_normalized_data <- renderUI({
    downloadButton("download_normalized_data_button", label = "Download normalized data")
  })
  
  # QC - download file for normalized count data
  output$download_normalized_data_button <- downloadHandler(
    filename = function() {
      paste0("Normalized_Data.csv")
    },
    content = function(file) {
      if(SubmitData$data_type == "Bulk") {
        dat <- assay(ddstran()[[1]])
      } else if(SubmitData$data_type == "RASL") {
        dat <- ddsout()[[4]]
      } else {
        if(d$resType == "seurat_res") {
          dat <- seurat_only()[[1]]@assays$RNA@data
        } else {
          dat <- seurat_only()@assays$RNA@data
        }
      }
      write.csv(
        dat,
        file = file,
        quote = F
      )
    }
  )
  
  # QC - add data relevant to QC to object d
  observeEvent(SubmitData$goqc, {
    if(is.null(SubmitData$goqc) || SubmitData$goqc == 0) return()
    d$corplotClick <- d$heatmap1Click <- d$setNames <- d$goqc_warning <- 
      d$wgcna_warning <- d$wgcna_warning_sc <- d$wgcna_minModuleSizeTried <-
      d$wgcna_minModuleSizeTried_sc <- d$sc_dge_byfactor_msg <- NULL
    d$newHeat <- d$newCorplot <- d$newDistplot <- T
    d$sc_dge_byfactor <- d$cellSets <- d$addCellsOverlap <- d$addCellsClicks <- list()
    # Reset event_data("plotly_selected")
    runjs("Shiny.setInputValue('plotly_selected-A', null);") 

    d$newCluster <- T
    if(SubmitData$data_type == "Single-cell") {
      d$newExpSetup <- F
      d$a <- d$b <- NULL
      withProgress(message = "Processing...", value = 0, {
        incProgress(1/2)
        d$cts_out <- data.frame(ddsout()[[4]])
        incProgress(1/2)
      })
      
      if(d$resType == "seurat_res" && length(seurat_only()) == 0) {
        return()
      }
      
      sCVdL <- seurat_sc()
      d$seurat_sc <- sCVdL
      d$SCV <- sCVdL
      d$meta <- meta()
      d$clustList <- clustList()
      d$numClust <- numClust()
      if(d$resType == "seurat_res") {
        d$inD <- seurat_only()[[grep("^res", names(seurat_only()))[1]]]
        d$res <- names(d$SCV)[1]
        d$meta$seurat_clusters <- d$meta[[names(d$SCV)[1]]]
        d$meta$silWidth <- d$meta[[paste0(names(d$SCV)[1]), "_silWidth"]]
        # determine if any metadata factors can be used for DGE
        mdCols <- colnames(seurat_only()[[1]]@meta.data)
        mdCols <- mdCols[!(apply(seurat_only()[[1]]@meta.data, 2, function(j) is.numeric(j) | length(unique(j)) == 1))]
        if(length(mdCols) > 0) {
          mdCols <- grep("^nCount_RNA$|^nFeature_RNA$", x = mdCols, value = T, invert = T)
        }
        if(length(mdCols) > 0) {
          mdCols <- mdCols[sapply(mdCols, function(mdCol) sum(as.vector(table(seurat_only()[[1]]@meta.data[, mdCol])) >= 3) >= 2)]
        }
      } else {
        d$inD <- seurat_only()
        d$res <- "Clust"
        # Determine if any metadata factors can be used for DGE
        mdCols <- colnames(seurat_only()@meta.data)
        mdCols <- mdCols[!(apply(seurat_only()@meta.data, 2, function(j) is.numeric(j) | length(unique(j)) == 1))]
        if(length(mdCols) > 0) {
          mdCols <- grep("^nCount_RNA$|^nFeature_RNA$", x = mdCols, value = T, invert = T)
        }
        if(length(mdCols) > 0) {
          mdCols <- mdCols[sapply(mdCols, function(mdCol) sum(as.vector(table(seurat_only()@meta.data[, mdCol])) >= 3) >= 2)]
        }
      }
      if(length(mdCols) == 0) {
        d$sc_dge_byfactor_novars_msg <- "No valid columns in metadata with 3 or more cells from at least 2 groups."
      } else {
        d$sc_dge_byfactor_novars_msg <- NULL
      }
      d$MD <- getMD(d$inD)
      d$inD_orig <- d$inD
      symbolMap <<- NULL

    } else if(SubmitData$data_type == "Bulk") {
      d$cts_out <- data.frame(ddsout()[[4]])
      callModule(profiler, 
                 id = "profiler_bulk", 
                 seuratData = list(counts = ddsout()[[3]], metadata = ddsout()[[2]]))
      callModule(dimred, 
                 id = "dda_dimred_bulk", 
                 dat = list(data = ddstran()[[1]], 
                            metadata = ddsout()[[2]]))
      d$newExpSetup <- T
    } else if (SubmitData$data_type == "RASL") {
      d$cts_out <- data.frame(ddsout()[[4]])
      callModule(dimred, 
                 id = "dda_dimred_rasl", 
                 dat = list(data = ddsout()[[4]], 
                            metadata = ddsout()[[2]]))
      d$newExpSetup <- T
    }
    # Insert additional tabs and switch to File summary.
    if(d$showTabs) {
      showTab(inputId = "tab_structure", target = "summary", select = T)
      showTab(inputId = "tab_structure", target = "val2")
      showTab(inputId = "tab_structure", target = "val3")
      if(SubmitData$data_type != "RASL") {
        showTab(inputId = "tab_structure", target = "val4")
      }
      
    } else {
      updateTabsetPanel(inputId = "tab_structure", selected = "summary")
    }
    d$showTabs <- F
    d$showDataSummaryLinks <- T
    updateTabsetPanel(inputId = "dataSummaryQCTabset", selected = "dataSummary")
    updateTabsetPanel(inputId = "dataSummary_tabsetPanel", selected = "dataSummary_subsetCountData")
  }, priority = 0)
  
  # When d$cts_out is NULL (when new data are loaded), remove Summary tab
  observeEvent(d$cts_out, {
    if(is.null(d$cts_out)) {
      hideTab(inputId = "tab_structure", target = "val4")
      hideTab(inputId = "tab_structure", target = "val3")
      hideTab(inputId = "tab_structure", target = "val2")
      hideTab(inputId = "tab_structure", target = "summary")
      d$showTabs <- T
    } 
  }, ignoreNULL = F, ignoreInit = T)
  
  # DimRed tab (bulk or sc)
  output$dda_dimred <- renderUI({
    req(SubmitData$data_type, ddsout())
    if(SubmitData$data_type == "Single-cell") {
      tagList(
        dimredUI("dda_dimred_sc")
      )
    } else if(SubmitData$data_type == "Bulk") {
      if(is.null(ddstran())) return()
      tagList(
        dimredUI("dda_dimred_bulk")
      )
    } else if(SubmitData$data_type == "RASL") {
      if(is.null(ddsout())) return()
      tagList(
        dimredUI("dda_dimred_rasl")
      )
    }
  })
  
  # Profiler UI
  output$profiler <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type == "Single-cell") {
      if(is.null(d$inD)) return()
      # tagList(
        profilerUI("profiler_sc")
      # )
    } else if(SubmitData$data_type == "Bulk") {
      if(is.null(ddsout())) return()
      # tagList(
        profilerUI("profiler_bulk")
      # )
    } else {
      return()
    }
  })
  
  # WHAT IS THIS? IS IT NEEDED?
  # QC - update items in object d?
  observeEvent(input$goscdge, {
    d$a <- d$b <- NULL
    compNames <- NULL
    output$calcText <- renderText("")
  }, ignoreInit = T)
  
  # QC - count data sub table DT output
  output$fileoutputcts <- renderUI({
    req(d$cts_out)
    # if(is.null(d$cts_out)) return()
    DT::dataTableOutput("outputcts")
  })
  
  # QC - count data sub table DT
  output$outputcts <- DT::renderDataTable(server = FALSE, {
    req(d$cts_out)
    colsToShow <- 1:min(10, ncol(d$cts_out))
    add <- data.frame(d$cts_out[1,colsToShow])
    
    add <- add %>%
      tidyr::gather(sample, value) %>%
      mutate(value = "...") %>%
      tidyr::spread(sample, value)
    cts <- rbind(d$cts_out[c(1:5),colsToShow], add, d$cts_out[(nrow(d$cts_out) - 4):nrow(d$cts_out), colsToShow])
    rownames(cts)[6] <- ""
    cts <- cts %>% 
      rownames_to_column(var = 'Gene')
    
    DT::datatable(
      cts,
      rownames = F,
      options = list(
        lengthChange = F,
        bFilter = F,
        bInfo = F, 
        bPaginate = F,
        bOrder = F,
        ordering = F,
        autoWidth = FALSE,
        columnDefs = list(list(className = 'dt-left', width = '40%', targets = "_all")),
        dom = 'frtip',
        scrollX = TRUE
      ),
      class = "display"
    ) %>%
      formatStyle(columns = "Gene", "white-space"="nowrap")
  })

  # QC - header for metadata summary
  output$filesummarycoldata <- renderUI({
    req(d$cts_out)
    if(is.null(d$cts_out)) return() 
    h4("Metadata")
  })
  
  # QC - metadata DT
  output$fileoutputcoldata <- renderUI({
    req(d$cts_out)
    if(is.null(d$cts_out)) return()
    DT::dataTableOutput("outcol")
  })
  
  # QC - DT elements for metadata output
  output$outcol <- DT::renderDataTable(server = FALSE, {
    req(ddsout(), SubmitData$data_type)
    sampleType <- "cell"
    if(SubmitData$data_type %in% c("Bulk", "RASL")) sampleType <- "sample"
    DT::datatable({
      coldata <- ddsout()[[2]] %>%
        rownames_to_column(var = sampleType)
      coldata
    }, rownames = F,
    extensions = 'Buttons',
    options = list(
      autoWidth = F,
      columnDefs = list(list(className = 'dt-left', width = '40%', targets = "_all")),
      dom = 'Bfrtip',
      scrollX = F,
      buttons = list(list(extend = 'csv', 
                          filename = "Metadata",
                          text = 'Download metadata'))
    ), class = "display")
  })
  
  # QC - header for boxplot count distributions
  output$countbox <- renderUI({
    req(d$cts_out)
    h4("Count data distributions - box and whisker")
  })
  
  # QC - select factors to use for box plots and legend
  output$countsummarylabel <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | SubmitData$data_type == "Single-cell") return()
    choices = colnames(ddsout()[[2]])
    choices = choices[!(choices %in% c("nCount_RNA","nFeature_RNA"))]
    selectInput(inputId = "countsummarylabel", 
                label = "Choose factor(s) for plots", 
                choices = choices, 
                selected = choices[1],
                multiple = T)
  })
  
  # QC - option to rotate labels on boxplot
  output$boxplot_rotate <- renderUI({
    req(d$cts_out)
    if(is.null(d$cts_out)) return()
    prettyCheckbox("boxplot_rotate", label = "Rotate labels", value = T, 
                   status = "default", icon = icon("check"))
  })
  
  # QC - visualize - boxplot
  output$boxplot <- renderPlotly({
    req(SubmitData$data_type)
    if(is.null(d$cts_out) || is.null(input$boxplot_rotate)) return()
    if(SubmitData$data_type %in% c("Bulk", "RASL") && is.null(input$countsummarylabel)) return()

    if(SubmitData$data_type == "RASL") {
      tmp <- ddsout()[[4]]
      lab <- "log<sub>2</sub>(fold-changes)"
    } else {
      tmp <- ddstran()[[1]]
      tmp <- assay(tmp)
      lab <- ddstran()[[2]]
    }
    
    tickAngle <- 0
    if(input$boxplot_rotate) tickAngle <- -45
    
    if(SubmitData$data_type == "Single-cell") {
      if("orig.ident" %in% colnames(ddsout()[[2]])) {
        sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
        xAxisLabel <- "orig.ident"
      } else {
        sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
        xAxisLabel <- colnames(ddsout()[[2]])[1]
      }
    } else {
      sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                            1, function(i) paste(i, collapse = "_"))
      xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
    }
    names(sampleLabels) <- colnames(tmp)
    isolate({
      box <- as.data.frame(tmp)
      box <- tidyr::gather(box)
      box <- box %>%
        mutate(sample = sampleLabels[box$key]) %>%
        dplyr::group_by(key, sample) %>%
        summarise(value = mean(value, na.rm = T))
      box$sample <- factor(box$sample, 
                           levels = mixedsort(unique(box$sample)))
      len <- length(unique(box$sample))
      
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
      
      graphics.off()
      pdf(NULL)
      plot_ly(
        box,
        type = "box",
        y = ~value,
        color = ~sample,
        colors = hue_pal()(len),
        hoverlabel = label
      ) %>%
        layout(
          margin = list(b = 90),
          xaxis = list(title = xAxisLabel, tickangle = tickAngle),
          yaxis = list(title = paste("Average ", lab, sep = "")),
          font = list(family = "Noto Sans JP"),
          showlegend = T
        ) 
    })
  })
  
  # QC - download button for boxplot (PDF)
  output$dlqcboxplotpdf <- renderUI({
    req(SubmitData$data_type, d$cts_out)
    if(is.null(d$cts_out) | is.null(input$boxplot_rotate)) return()
    if(SubmitData$data_type %in% c("Bulk", "RASL") && is.null(input$countsummarylabel)) return()
    downloadButton("dlqcboxplotpdfimg", "Download static plot (PDF)")
  })
  
  # QC - download file for boxplot (PDF)
  output$dlqcboxplotpdfimg <- downloadHandler(
    filename =  function() {
      paste("QC-boxplot.pdf")
    },
    content = function(file) {
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
      }
      tickAngle <- 0
      if(input$boxplot_rotate) tickAngle <- 45
      if(SubmitData$data_type == "Single-cell") {
        if("orig.ident" %in% colnames(ddsout()[[2]])) {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
          xAxisLabel <- "orig.ident"
        } else {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
          xAxisLabel <- colnames(ddsout()[[2]])[1]
        }
      } else {
        sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                              1, function(i) paste(i, collapse = "_"))
        xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
      }
      names(sampleLabels) <- colnames(tmp)
      box <- as.data.frame(tmp)
      box <- tidyr::gather(box)
      box <- box %>%
        mutate(sample = sampleLabels[box$key]) %>%
        dplyr::group_by(key, sample) %>%
        summarise(value = mean(value, na.rm = T))
      box$sample <- factor(box$sample, 
                           levels = mixedsort(unique(box$sample)))
      len <- length(unique(box$sample))
      cols <- hue_pal()(len)
      p <- ggplot(box, aes(x = sample, y = value, fill = sample)) +
        stat_boxplot(geom = "errorbar") +
        geom_boxplot() + theme_classic() +
        xlab(xAxisLabel) + ylab(expression(average~log[2]~ (counts + 1))) + 
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(angle = tickAngle, vjust = .4, size = 14),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_manual(values = cols)
      CairoPDF(file, width = 12, height = 6.5)
      print(p)
      dev.off()
    }
  )
  
  # QC - download button for boxplot (PNG)
  output$dlqcboxplotpng <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(input$boxplot_rotate)) return()
    if(SubmitData$data_type %in% c("Bulk", "RASL") && is.null(input$countsummarylabel)) return()
    downloadButton("dlqcboxplotpngimg", "Download static plot (PNG)")
  })
  
  # QC - download file for boxplot (PNG)
  output$dlqcboxplotpngimg <- downloadHandler(
    filename =  function() {
      paste("QC-boxplot.png")
    },
    content = function(file) {
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
      }
      tickAngle <- 0
      if(input$boxplot_rotate) tickAngle <- 45
      if(SubmitData$data_type == "Single-cell") {
        if("orig.ident" %in% colnames(ddsout()[[2]])) {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
          xAxisLabel <- "orig.ident"
        } else {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
          xAxisLabel <- colnames(ddsout()[[2]])[1]
        }
      } else {
        sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                              1, function(i) paste(i, collapse = "_"))
        xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
      }
      names(sampleLabels) <- colnames(tmp)
      box <- as.data.frame(tmp)
      box <- tidyr::gather(box)
      box <- box %>%
        mutate(sample = sampleLabels[box$key]) %>%
        dplyr::group_by(key, sample) %>%
        summarise(value = mean(value, na.rm = T))
      box$sample <- factor(box$sample, 
                           levels = mixedsort(unique(box$sample)))
      len <- length(unique(box$sample))
      cols <- hue_pal()(len)
      p <- ggplot(box, aes(x = sample, y = value, fill = sample)) +
        stat_boxplot(geom = "errorbar") +
        geom_boxplot() + theme_classic() +
        xlab(xAxisLabel) + ylab(expression(average~log[2]~ (counts + 1))) + 
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(angle = tickAngle, vjust = .4, size = 14),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_manual(values = cols)
      png(file, width = 1200, height = 650)
      print(p)
      dev.off()
    }
  )
  
  # QC - header for histogram of count data
  output$counthist <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out)) return()
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    h4("Count data distributions - histogram")
  })
  
  # QC - histogram of count data
  output$hist <- renderPlotly({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type %in% c("Bulk", "RASL") && is.null(input$countsummarylabel)) return()

    if(SubmitData$data_type == "RASL") {
      tmp <- ddsout()[[4]]
      lab <- "log<sub>2</sub>(fold-changes)"
    } else {
      tmp <- ddstran()[[1]]
      tmp <- assay(tmp)
      lab <- ddstran()[[2]]
    }
    
    if(SubmitData$data_type == "Single-cell") {
      if("orig.ident" %in% colnames(ddsout()[[2]])) {
        sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
        xAxisLabel <- "orig.ident"
      } else {
        sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
        xAxisLabel <- colnames(ddsout()[[2]])[1]
      }
    } else {
      sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                            1, function(i) paste(i, collapse = "_"))
      xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
    }
    names(sampleLabels) <- colnames(tmp)
    isolate({
      hist <- as.data.frame(tmp)
      hist <- tidyr::gather(hist)
      hist <- hist %>%
        mutate(sample = sampleLabels[hist$key]) %>%
        dplyr::group_by(key, sample) %>%
        summarise(value = mean(value, na.rm = T))
      hist$sample <- factor(hist$sample, 
                            levels = mixedsort(unique(hist$sample)))
      len <- length(unique(hist$sample))
      
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
      
      graphics.off()
      pdf(NULL)
      plot_ly(
        hist,
        type = "histogram",
        x = ~value,
        color = ~sample,
        colors = hue_pal()(len), 
        hoverlabel = label
      ) %>%
        layout(
          xaxis = list(title = lab),
          yaxis = list(title = paste("Average Frequency")),
          font = list(family = "Noto Sans JP"),
          showlegend = T,
          bargap = 0.2
        )
    })
  })
  
  # QC - download button for histogram (PDF)
  output$dlqchistpdf <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type)) return()
    
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    downloadButton("dlqchistpdfimg", "Download static plot (PDF)")
  })
  
  # QC - download file for histogram (PDF)
  output$dlqchistpdfimg <- downloadHandler(
    filename =  function() {
      paste("QC-histogram.pdf")
    },
    content = function(file) {
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
      }
      if(SubmitData$data_type == "Single-cell") {
        if("orig.ident" %in% colnames(ddsout()[[2]])) {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
          xAxisLabel <- "orig.ident"
        } else {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
          xAxisLabel <- colnames(ddsout()[[2]])[1]
        }
      } else {
        sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                              1, function(i) paste(i, collapse = "_"))
        xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
      }
      names(sampleLabels) <- colnames(tmp)
      hist <- as.data.frame(tmp)
      hist <- tidyr::gather(hist)
      hist <- hist %>%
        mutate(sample = sampleLabels[hist$key]) %>%
        dplyr::group_by(key, sample) %>%
        summarise(value = mean(value, na.rm = T))
      hist$sample <- factor(hist$sample, 
                            levels = mixedsort(unique(hist$sample)))
      len <- length(unique(hist$sample))
      cols <- hue_pal()(len)
      p <- ggplot(hist, aes(x = value, fill = sample)) +
        geom_histogram(position = "dodge") + theme_classic() +
        xlab(expression(average~log[2]~ (counts + 1))) + ylab("Average frequency") + 
        theme(axis.text.x = element_text(angle = 0, vjust = .4)) +
        scale_fill_manual(values = cols) +
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0))
      CairoPDF(file, width = 12, height = 6.5)
      print(p)
      dev.off()
    }
  )
  
  # QC - download button for histogram (PNG)
  output$dlqchistpng <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type)) return()
    
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    downloadButton("dlqchistpngimg", "Download static plot (PNG)")  
  })
  
  # QC - download file for histogram (PNG)
  output$dlqchistpngimg <- downloadHandler(
    filename =  function() {
      paste("QC-histogram.png")
    },
    content = function(file) {
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
      }
      hist <- as.data.frame(tmp)
      if(SubmitData$data_type == "Single-cell") {
        if("orig.ident" %in% colnames(ddsout()[[2]])) {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
          xAxisLabel <- "orig.ident"
        } else {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
          xAxisLabel <- colnames(ddsout()[[2]])[1]
        }
      } else {
        sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                              1, function(i) paste(i, collapse = "_"))
        xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
      }
      names(sampleLabels) <- colnames(tmp)
      hist <- tidyr::gather(hist)
      hist <- hist %>%
        mutate(sample = sampleLabels[hist$key]) %>%
        dplyr::group_by(key, sample) %>%
        summarise(value = mean(value, na.rm = T))
      hist$sample <- factor(hist$sample, 
                            levels = mixedsort(unique(hist$sample)))
      len <- length(unique(hist$sample))
      cols <- hue_pal()(len)
      p <- ggplot(hist, aes(x = value, fill = sample)) +
        geom_histogram(position = "dodge") + theme_classic() +
        xlab(expression(average~log[2]~ (counts + 1))) + ylab("Average frequency") + 
        theme(axis.text.x = element_text(angle = 0, vjust = .4)) +
        scale_fill_manual(values = cols) + 
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0))
      png(file, width = 1200, height = 650)
      print(p)
      dev.off()
    }
  )
  
  # QC - header for total reads barplot
  output$counttotal <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type)) return()
    
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    h4("Total reads")
  })
  
  # QC - option to rotate labels on barplot
  output$barplot_rotate <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out)) return()
    prettyCheckbox("barplot_rotate", label = "Rotate labels", value = T, status = "default", icon = icon("check"))
  })
  
  # QC - total reads barplot
  output$barplot <- renderPlotly({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type) | is.null(input$barplot_rotate)) return()
    
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    
    if(SubmitData$data_type == "RASL") {
      tmp <- ddsout()[[4]]
    } else {
      tmp <- ddstran()[[1]]
      tmp <- assay(tmp)
    }
    tickAngle <- 0
    if(input$barplot_rotate) tickAngle <- -45
    if(SubmitData$data_type == "Single-cell") {
      if("orig.ident" %in% colnames(ddsout()[[2]])) {
        sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
        xAxisLabel <- "orig.ident"
      } else {
        sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
        xAxisLabel <- colnames(ddsout()[[2]])[1]
      }
    } else {
      sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                            1, function(i) paste(i, collapse = "_"))
      xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
    }
    names(sampleLabels) <- colnames(tmp)
    isolate({
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
        bar <- as.data.frame(tmp)
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
        dds <- ddsout()[[1]]
        bar <- as.data.frame(assay(dds))
      }
      bar <- colSums(bar)
      bar <- as.data.frame(t(bar))
      bar <- tidyr::gather(bar)
      bar <- bar %>%
        mutate(sample = sampleLabels[bar$key]) %>%
        dplyr::group_by(sample) %>%
        summarise(value = mean(value, na.rm = T))
      
      bar$sample <- factor(bar$sample, 
                           levels = mixedsort(bar$sample))
      len <- length(unique(bar$sample))
      
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
      
      graphics.off()
      pdf(NULL)
      plot_ly(
        bar,
        type = "bar",
        y = ~value,
        x = ~sample,
        color = ~sample,
        colors = hue_pal()(len),
        hoverlabel = label
      ) %>%
        layout(
          margin = list(b = 90),
          xaxis = list(title = xAxisLabel, tickangle = tickAngle),
          yaxis = list(title = "Counts"),
          font = list(family = "Noto Sans JP"),
          showlegend = T
        )
    })
  })
  
  # QC - download button for barplot (PDF)
  output$dlqcbarplotpdf <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type) | is.null(input$barplot_rotate)) return()
    
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    downloadButton("dlqcbarplotpdfimg", "Download static plot (PDF)")
  })
  
  # QC - download file for barplot (PDF)
  output$dlqcbarplotpdfimg <- downloadHandler(
    filename =  function() {
      paste("QC-barplot.pdf")
    },
    content = function(file) {
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
        bar <- as.data.frame(tmp)
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
        dds <- ddsout()[[1]]
        bar <- as.data.frame(assay(dds))
      }
      bar <- colSums(bar)
      bar <- as.data.frame(t(bar))
      tickAngle <- 0
      if(input$barplot_rotate) tickAngle <- 45
      if(SubmitData$data_type == "Single-cell") {
        if("orig.ident" %in% colnames(ddsout()[[2]])) {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
          xAxisLabel <- "orig.ident"
        } else {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
          xAxisLabel <- colnames(ddsout()[[2]])[1]
        }
      } else {
        sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                              1, function(i) paste(i, collapse = "_"))
        xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
      }
      names(sampleLabels) <- colnames(tmp)
      bar <- tidyr::gather(bar)
      bar <- bar %>%
        mutate(sample = sampleLabels[bar$key]) %>%
        dplyr::group_by(sample) %>%
        summarise(value = mean(value, na.rm = T))
      bar$sample <- factor(bar$sample, 
                           levels = mixedsort(bar$sample))
      len <- length(unique(bar$sample))
      cols <- hue_pal()(len)
      p <- ggplot(bar, aes(x = sample, y = value, fill = sample)) +
        geom_bar(stat = "identity") + theme_classic() +
        xlab(xAxisLabel) + ylab("Counts")  + 
        theme(axis.text.x = element_text(angle = tickAngle, vjust = .4)) +
        scale_fill_manual(values = cols) +
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_discrete(expand = c(0,0))
      CairoPDF(file, width = 12, height = 6.5)
      print(p)
      dev.off()
    }
  )
  
  # QC - download button for barplot (PNG)
  output$dlqcbarplotpng <- renderUI({
    req(d$cts_out, SubmitData$data_type)
    if(is.null(d$cts_out) | is.null(SubmitData$data_type) | is.null(input$barplot_rotate)) return()
    
    if(SubmitData$data_type %in% c("Bulk", "RASL") & is.null(input$countsummarylabel)) return()
    downloadButton("dlqcbarplotpngimg", "Download static plot (PNG)")    
  })
  
  # QC - download file for barplot (PNG)
  output$dlqcbarplotpngimg <- downloadHandler(
    filename =  function() {
      paste("QC-barplot.png")
    },
    content = function(file) {
      if(SubmitData$data_type == "RASL") {
        tmp <- ddsout()[[4]]
        bar <- as.data.frame(tmp)
      } else {
        tmp <- ddstran()[[1]]
        tmp <- assay(tmp)
        dds <- ddsout()[[1]]
        bar <- as.data.frame(assay(dds))
      }
      bar <- colSums(bar)
      bar <- as.data.frame(t(bar))
      tickAngle <- 0
      if(input$barplot_rotate) tickAngle <- 45
      if(SubmitData$data_type == "Single-cell") {
        if("orig.ident" %in% colnames(ddsout()[[2]])) {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), "orig.ident"])
          xAxisLabel <- "orig.ident"
        } else {
          sampleLabels <- as.character(ddsout()[[2]][colnames(tmp), 1])
          xAxisLabel <- colnames(ddsout()[[2]])[1]
        }
      } else {
        sampleLabels <- apply(ddsout()[[2]][colnames(tmp), input$countsummarylabel, drop = F], 
                              1, function(i) paste(i, collapse = "_"))
        xAxisLabel <- paste(input$countsummarylabel, collapse = " + ")
      }
      names(sampleLabels) <- colnames(tmp)
      bar <- tidyr::gather(bar)
      bar <- bar %>%
        mutate(sample = sampleLabels[bar$key]) %>%
        dplyr::group_by(sample) %>%
        summarise(value = mean(value, na.rm = T))
      bar$sample <- factor(bar$sample, 
                           levels = mixedsort(bar$sample))
      len <- length(unique(bar$sample))
      cols <- hue_pal()(len)
      p <- ggplot(bar, aes(x = sample, y = value, fill = sample)) +
        geom_bar(stat = "identity") + theme_classic() +
        xlab(xAxisLabel) + ylab("Counts")  + 
        theme(axis.text.x = element_text(angle = tickAngle, vjust = .4)) +
        scale_fill_manual(values = cols) +
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_discrete(expand = c(0,0))
      png(file, width = 1200, height = 650)
      print(p)
      dev.off()
    }
  )
  
  ###################################################################
  ###################################################################
  ### SECTION 03 - Discovery-Driven Analysis (DDA)
  ###################################################################
  ###################################################################

  output$showCorrelationLinks <- renderUI({
    req(input$dda_tabsetPanel)
    myLab <- "Correlation"
    if(input$dda_tabsetPanel == "dda_correlation") myLab <- strong(myLab)
    actionLink(
      inputId = "showCorrelationLinks", 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })
  
  # observeEvents with ifelse to expand/collapse subtabs
  observeEvent(input$showCorrelationLinks, {
    d$showCorrelationLinks <- ifelse(
      test = (input$dda_tabsetPanel == "dda_correlation"), 
      yes = ifelse(
        test = d$showCorrelationLinks,
        yes = F,
        no = T
      ), 
      no = T
    )
    d$showClustering <- F
    updateTabsetPanel(inputId = "dda_tabsetPanel", selected = "dda_correlation")
  })

  output$correlationLinks <- renderUI({
    req(d$showCorrelationLinks, SubmitData$data_type, input$dda_correlation_tabsetPanel)
    distanceLab <- "Distance"
    correlationLab <- "Correlation"
    if(input$dda_correlation_tabsetPanel == "dda_correlation_correlation") correlationLab <- strong(correlationLab)
    if(input$dda_correlation_tabsetPanel == "dda_correlation_distance") distanceLab <- strong(distanceLab)
    
    if(SubmitData$data_type %in% c("Bulk", "RASL")) {
      fluidPage(
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "showCorrelation", 
            label = correlationLab, 
            style = "color: black;"
          )
        ),
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "showDistance", 
            label = distanceLab, 
            style = "color: black;"
          )
        )
      )
    } else {
      fluidPage(
        fluidRow(
          style = "margin-left: 6px;", 
          actionLink(
            inputId = "showCorrelation", 
            label = correlationLab, 
            style = "color: black;"
          )
        )
      )
    }
  })
  
  output$showClustering <- renderUI({
    req(input$dda_tabsetPanel)
    myLab <- "Clustering"
    if(input$dda_tabsetPanel == "dda_clustering") myLab <- strong(myLab)
    actionLink(
      inputId = "showClustering", 
      style = "font-size: 16px; color: black;", 
      label = myLab
    )
  })
  # observeEvents to update dda_tabsetPanel tabsetPanel
  observeEvent(input$showCorrelation, {
    updateTabsetPanel(inputId = "dda_correlation_tabsetPanel", selected = "dda_correlation_correlation")
  })
  
  observeEvent(input$showDistance, {
    updateTabsetPanel(inputId = "dda_correlation_tabsetPanel", selected = "dda_correlation_distance")
  })
  
  observeEvent(input$showClustering, {
    d$showCorrelationLinks <- F
    updateTabsetPanel(inputId = "dda_tabsetPanel", selected = "dda_clustering")
  })
  
  # DDA - help for correlation heatmap
  observeEvent(input$cor_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Correlation")),
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
      includeMarkdown("markdown/help/correlation_help.md")
    ))
  })

  # DDA - reset plotly click in correlation plot
  observeEvent(input$build_cor_matrix, {
    d$corplotClick <- NULL
    d$newCorplot <- F
  }, priority = 10)
  
  # DDA - correlation matrix
  corout <- eventReactive(input$build_cor_matrix, {
    nSamplesDownsample <- min(500, input$cor_nSamples)
    if(SubmitData$data_type == "Bulk") {
      sampRows <- 1:nrow(ddstran()[[1]])
      sampCols <- 1:ncol(ddstran()[[1]])
      if(input$cor_nGenes < nrow(ddstran()[[1]])) {
        sampRows <- names(sort(rowMeans(assay(ddstran()[[1]])), decreasing = T))[1:input$cor_nGenes]
      }
      if(nSamplesDownsample < ncol(ddstran()[[1]])) {
        set.seed(input$cor_seed)
        sampCols <- sample(1:ncol(ddstran()[[1]]), size = nSamplesDownsample)
      }
      mat <- assay(ddstran()[[1]][sampRows, sampCols])
    } else if(SubmitData$data_type == "RASL") {
      sampRows <- 1:nrow(ddsout()[[4]])
      sampCols <- 1:ncol(ddsout()[[4]])
      if(input$cor_nGenes < nrow(ddsout()[[4]])) {
        sampRows <- names(sort(rowMeans(ddsout()[[4]]), decreasing = T))[1:input$cor_nGenes]
      }
      if(nSamplesDownsample < ncol(ddsout()[[4]])) {
        set.seed(input$cor_seed)
        sampCols <- sample(1:ncol(ddsout()[[4]]), size = nSamplesDownsample)
      }
      mat <- ddsout()[[4]][sampRows, sampCols]
      
    } else {
      if(d$resType == "seurat_res") {
        sampRows <- 1:nrow(seurat_only()[[1]])
        sampCols <- 1:ncol(seurat_only()[[1]])
        if(input$cor_nGenes < nrow(seurat_only()[[1]])) {
          sampRows <- names(sort(rowMeans(seurat_only()[[1]]@assays$RNA@data), decreasing = T))[1:input$cor_nGenes]
        }
        if(nSamplesDownsample < ncol(seurat_only()[[1]])) {
          set.seed(10)
          sampCols <- sort(sample(1:ncol(seurat_only()[[1]]), size = nSamplesDownsample), decreasing = F)
        }
        mat <- as.matrix(seurat_only()[[1]]@assays$RNA@data[sampRows, sampCols])
      } else {
        sampRows <- 1:nrow(seurat_only())
        sampCols <- 1:ncol(seurat_only())
        if(input$cor_nGenes < nrow(seurat_only())) {
          sampRows <- names(sort(rowMeans(seurat_only()@assays$RNA@data), decreasing = T))[1:input$cor_nGenes]
        }
        if(nSamplesDownsample < ncol(seurat_only())) {
          set.seed(10)
          sampCols <- sort(sample(1:ncol(seurat_only()), size = nSamplesDownsample), decreasing = F)
        }
        mat <- as.matrix(seurat_only()@assays$RNA@data[sampRows, sampCols])
      }
    }
    return(list(cor(mat), mat))
  })
  
  # DDA - distance matrix
  distout <- eventReactive(input$build_dist_matrix, {
    if(is.null(SubmitData$data_type) || SubmitData$data_type != "Bulk") return()
    nSamplesDownsample <- min(500, input$dist_nSamples)
    sampRows <- 1:nrow(ddstran()[[1]])
    sampCols <- 1:ncol(ddstran()[[1]])
    if(input$dist_nGenes < nrow(ddstran()[[1]])) {
      sampRows <- names(sort(rowMeans(assay(ddstran()[[1]])), decreasing = T))[1:input$dist_nGenes]
    }
    if(nSamplesDownsample < ncol(ddstran()[[1]])) {
      set.seed(input$dist_seed)
      sampCols <- sample(1:ncol(ddstran()[[1]]), size = nSamplesDownsample)
    }
    mat <- assay(ddstran()[[1]][sampRows, sampCols])

    withProgress(message = "Creating distance matrix...", value = 0, {
      incProgress(1/3)
      sdm <- as.matrix(stats::dist(t(mat)))
      incProgress(1/3)
      tooltips <- paste0(
        "<b>Dist:</b> ", round(sdm, 3)
      )
      tooltips <- matrix(
        tooltips,
        ncol = ncol(sdm),
        byrow = FALSE
      )
      dimnames(tooltips) <- dimnames(sdm)
      
      for(x in rownames(tooltips)) {
        for(y in colnames(tooltips)) {
          tooltips[x, y] <- paste0(
            tooltips[x, y],
            "<br />",
            "<b>X:</b> ",
            x,
            "<br />",
            "<b>Y:</b> ",
            y
          )
        }
      }
      incProgress(1/3)
    })

    return(list(sdm = sdm, tooltips = tooltips))
  })
  
  # DDA - help for distance heatmap
  observeEvent(input$dist_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Distance matrix")),
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
      includeMarkdown("markdown/help/distance_help.md")
    ))
  })
  
  observeEvent(input$build_dist_matrix, {
    d$newDistplot <- F
  })
  

  # DDA - specify number of samples/cells for correlation heatmap
  output$cor_nSamples <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type == "Bulk") {
      lab <- "Samples"
      maxSamples <- ncol(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      lab <- "Samples"
      maxSamples <- ncol(ddsout()[[4]])
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    tempVal <- min(500, maxSamples)
    numericInput("cor_nSamples", label = h3(lab), value = tempVal, min = 2, max = tempVal, step = 1)
  })

  observeEvent(input$cor_nSamples, {
    if(is.finite(input$cor_nSamples) && input$cor_nSamples >= 2) return()
    if(SubmitData$data_type == "Bulk") {
      lab <- "Samples"
      maxSamples <- ncol(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      lab <- "Samples"
      maxSamples <- ncol(ddsout()[[4]])
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    tempVal <- min(500, maxSamples)
    updateNumericInput(session = session, inputId = "cor_nSamples", value = tempVal)
  })
  
  # DDA - action button to use all cells or samples for
  # correlation plot
  output$cor_all_samples <- renderUI({
    actionButton("cor_all_samples", label = "All")
  })

  # DDA - select the number of samples you want to visualize in
  # correlation matrix for downsampling
  observeEvent(input$cor_all_samples, {
    if(is.null(input$cor_all_samples)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxSamples <- ncol(assay(ddstran()[[1]]))
      lab <- "Samples"
    } else if(SubmitData$data_type == "RASL") {
      maxSamples <- ncol(ddsout()[[4]])
      lab <- "Samples"
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "cor_nSamples", label = lab, value = maxSamples,
                       min = 2, max = maxSamples, step = 1)
  })

  # # DDA - action button to submit default # cells or samples
  output$cor_default_samples <- renderUI({
    actionLink("cor_default_samples", label = "Max")
  })

  # # DDA - random seed for downsampling cells/samples
  output$cor_seed <- renderUI({
    fluidRow(
      div(
        style = "display: inline-block; width: 100px;", 
        numericInput("cor_seed", label = NULL, value = 10, step = 1)
      ),
      div(style = "display: inline-block; margin-left: 10px;", h5("Seed"))
    )
  })

  # DDA - THIS LOOKS DUPLICATED AS IN LINE 3624, WE
  # SHOULD BE ABLE TO COMBINE THESE INTO ONE OBSERVEEVENT
  observeEvent(input$cor_default_samples, {
    if(is.null(input$cor_default_samples)) return()
    # if(is.null(ddstran())) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxSamples <- ncol(assay(ddstran()[[1]]))
      lab <- "Samples"
    } else if(SubmitData$data_type == "RASL") {
      maxSamples <- ncol(ddsout()[[4]])
      lab <- "Samples"
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "cor_nSamples", value = min(500, maxSamples))
  })

  # Correlation options header
  output$correlation_options_header <- renderUI({
    h3("Correlation options", style = "margin-top: 0px;")
  })
  
  # DDA - number of genes for correlation heatmap
  output$cor_nGenes <- renderUI({
    req(SubmitData$data_type)

    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    tempVal <- min(2000, maxGenes)
    numericInput(
      inputId = "cor_nGenes", 
      label = h3("Genes"), 
      value = tempVal, 
      min = 2, 
      max = maxGenes, 
      step = 1
    )
  })
  
  # Reset default if user inputs invalid values
  observeEvent(input$cor_nGenes, {
    if(is.finite(input$cor_nGenes) && input$cor_nGenes >= 2) return()

    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    tempVal <- min(2000, maxGenes)
    updateNumericInput(session = session, inputId = "cor_nGenes", value = tempVal)
  })

  # DDA - action button to use all genes
  output$cor_all_genes <- renderUI({
    actionLink("cor_all_genes", label = "All")
  })

  # DDA - THIS LOOKS DUPLICATED AS IN LINE 3684, WE
  # SHOULD BE ABLE TO COMBINE THESE INTO ONE OBSERVEEVENT
  observeEvent(input$cor_all_genes, {
    if(is.null(input$cor_all_genes)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "cor_nGenes", value = maxGenes)
  })

  # DDA - button to use default # genes
  output$cor_default_genes <- renderUI({
    actionLink("cor_default_genes", label = "Default")
  })

  # DDA - THIS LOOKS DUPLICATED AS IN LINE 3712, WE
  # SHOULD BE ABLE TO COMBINE THESE INTO ONE OBSERVEEVENT
  observeEvent(input$cor_default_genes, {
    if(is.null(input$cor_default_genes)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "cor_nGenes", value = min(2000, maxGenes))
  })

  # DDA - select input - correlation heatmap colors
  output$plotColors <- renderUI({
    selectInput(
      inputId = "plotColors",
      label = "Color palette",
      choices = c(
        "Default (Blue-white-red)" = "Default",
        "Viridis (Blue-green-yellow)" = "Viridis",
        "Green-yellow-red" = "Green-yellow-red"
      ),
      selected = "Default"
    )
  })

  # DDA - action button for correlation heatmap
  output$build_cor_matrix <- renderUI({
    actionButton("build_cor_matrix", "Build heatmap", icon = icon("space-shuttle"))
  })
  
  # DDA - specify number of samples/cells for distance heatmap
  output$dist_nSamples <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type == "Bulk") {
      lab <- "Samples"
      maxSamples <- ncol(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      lab <- "Samples"
      maxSamples <- ncol(ddsout()[[4]])
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    tempVal <- min(500, maxSamples)
    numericInput("dist_nSamples", label = h3(lab), value = tempVal, min = 2, max = tempVal, step = 1)
  })
  
  observeEvent(input$dist_nSamples, {
    if(is.finite(input$dist_nSamples) && input$dist_nSamples >= 2) return()
    if(SubmitData$data_type == "Bulk") {
      lab <- "Samples"
      maxSamples <- ncol(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      lab <- "Samples"
      maxSamples <- ncol(ddsout()[[4]])
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    tempVal <- min(500, maxSamples)
    updateNumericInput(session = session, inputId = "dist_nSamples", value = tempVal)
  })
  
  # DDA - action button to use all cells or samples for
  # distance plot
  output$dist_all_samples <- renderUI({
    actionButton("dist_all_samples", label = "All")
  })
  
  # DDA - select the number of samples you want to visualize in
  # distance matrix for downsampling
  observeEvent(input$dist_all_samples, {
    if(is.null(input$dist_all_samples)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxSamples <- ncol(assay(ddstran()[[1]]))
      lab <- "Samples"
    } else if(SubmitData$data_type == "RASL") {
      maxSamples <- ncol(ddsout()[[4]])
      lab <- "Samples"
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "dist_nSamples", label = lab, value = maxSamples,
                       min = 2, max = maxSamples, step = 1)
  })
  
  # # DDA - action button to submit default # cells or samples
  output$dist_default_samples <- renderUI({
    actionLink("dist_default_samples", label = "Max")
  })
  
  # # DDA - random seed for downsampling cells/samples
  output$dist_seed <- renderUI({
    fluidRow(
      div(
        style = "display: inline-block; width: 100px;", 
        numericInput("dist_seed", label = NULL, value = 10, step = 1)
      ),
      div(style = "display: inline-block; margin-left: 10px;", h5("Seed"))
    )
  })
  
  # DDA - THIS LOOKS DUPLICATED AS IN LINE 3624, WE
  # SHOULD BE ABLE TO COMBINE THESE INTO ONE OBSERVEEVENT
  observeEvent(input$dist_default_samples, {
    if(is.null(input$dist_default_samples)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxSamples <- ncol(assay(ddstran()[[1]]))
      lab <- "Samples"
    } else if(SubmitData$data_type == "RASL") {
      maxSamples <- ncol(ddsout()[[4]])
      lab <- "Samples"
    } else {
      lab <- "Cells"
      if(d$resType == "seurat_res") {
        maxSamples <- ncol(seurat_only()[[1]])
      } else {
        maxSamples <- ncol(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "dist_nSamples", value = min(500, maxSamples))
  })
  
  # Distance options header
  output$distance_options_header <- renderUI({
    h3("Distance matrix options", style = "margin-top: 0px;")
  })
  
  # DDA - number of genes for distance heatmap
  output$dist_nGenes <- renderUI({
    req(SubmitData$data_type)
    
    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    tempVal <- min(2000, maxGenes)
    numericInput(
      inputId = "dist_nGenes", 
      label = h3("Genes"), 
      value = tempVal, 
      min = 2, 
      max = maxGenes, 
      step = 1
    )
  })
  
  # Reset default if user inputs invalid values
  observeEvent(input$dist_nGenes, {
    if(is.finite(input$dist_nGenes) && input$dist_nGenes >= 2) return()
    
    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    tempVal <- min(2000, maxGenes)
    updateNumericInput(session = session, inputId = "dist_nGenes", value = tempVal)
  })
  
  # DDA - action button to use all genes
  output$dist_all_genes <- renderUI({
    actionLink("dist_all_genes", label = "All")
  })
  
  # DDA - THIS LOOKS DUPLICATED AS IN LINE 3684, WE
  # SHOULD BE ABLE TO COMBINE THESE INTO ONE OBSERVEEVENT
  observeEvent(input$dist_all_genes, {
    if(is.null(input$dist_all_genes)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "dist_nGenes", value = maxGenes)
  })
  
  # DDA - button to use default # genes
  output$dist_default_genes <- renderUI({
    actionLink("dist_default_genes", label = "Default")
  })
  
  # DDA - THIS LOOKS DUPLICATED AS IN LINE 3712, WE
  # SHOULD BE ABLE TO COMBINE THESE INTO ONE OBSERVEEVENT
  observeEvent(input$dist_default_genes, {
    if(is.null(input$dist_default_genes)) return()
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type == "Bulk") {
      maxGenes <- nrow(assay(ddstran()[[1]]))
    } else if(SubmitData$data_type == "RASL") {
      maxGenes <- nrow(ddsout()[[4]])
    } else {
      if(d$resType == "seurat_res") {
        maxGenes <- nrow(seurat_only()[[1]])
      } else {
        maxGenes <- nrow(seurat_only())
      }
    }
    updateNumericInput(session, inputId = "dist_nGenes", value = min(2000, maxGenes))
  })
  
  # DDA - select input - distance heatmap colors
  output$dist_plotColors <- renderUI({
    selectInput(
      inputId = "dist_plotColors",
      label = "Color palette",
      choices = c(
        "Default (Blue-white-red)" = "Default",
        "Viridis (Blue-green-yellow)" = "Viridis",
        "Green-yellow-red" = "Green-yellow-red"
      ),
      selected = "Default"
    )
  })
  
  # DDA - color schemes for distance matrix
  getColors_dist <-  eventReactive(input$build_dist_matrix, {
    if(input$dist_plotColors == "Viridis") {
      return(viridis(100))
    } else if(input$dist_plotColors == "Green-yellow-red") {
      return(rev(brewer.pal(n=11, name = "RdYlGn")))
    } else {
      return(colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100))
    }
  })
  
  # DDA - action button for distance heatmap
  output$build_dist_matrix <- renderUI({
    actionButton("build_dist_matrix", "Build heatmap", icon = icon("space-shuttle"))
  })

  # DDA - display data type (bulk or sc)
  output$dataTypeCor <- renderText({
    req(SubmitData$data_type)
    dataType = "Single-cell RNA-seq"
    if(SubmitData$data_type == "Bulk") dataType <- "Bulk RNA-seq"
    if(SubmitData$data_type == "RASL") dataType <- "RASL-seq"
    return(dataType)
  })

  # THIS IS REDUNDANT TO ABOVE IN LINE 3783
  # DDA - display RNA-seq data type (single-cell or bulk)
  output$dataTypeDimRed <- renderText({
    req(SubmitData$data_type)
    dataType = "Single-cell RNA-seq"
    if(SubmitData$data_type == "Bulk") dataType <- "Bulk RNA-seq"
    if(SubmitData$data_type == "RASL") dataType <- "RASL-seq"
    return(dataType)
  })

  # THE NEXT 3 INPUTS ARE REDUNDANT; WRAP INTO MODULE
  # DDA - display data type (bulk or sc)
  output$dataTypePCA <- renderText({
    req(SubmitData$data_type)
    dataType = "Single-cell RNA-seq"
    if(SubmitData$data_type == "Bulk") dataType <- "Bulk RNA-seq"
    if(SubmitData$data_type == "RASL") dataType <- "RASL-seq"
    return(dataType)
  })

  # DDA - display data type (bulk or sc)
  output$dataTypeTSNE <- renderText({
    req(SubmitData$data_type)
    dataType = "Single-cell RNA-seq"
    if(SubmitData$data_type == "Bulk") dataType <- "Bulk RNA-seq"
    if(SubmitData$data_type == "RASL") dataType <- "RASL-seq"
    return(dataType)
  })

  # DDA - display data type (bulk or sc)
  output$dataTypeUMAP <- renderText({
    req(SubmitData$data_type)
    dataType = "Single-cell RNA-seq"
    if(SubmitData$data_type == "Bulk") dataType <- "Bulk RNA-seq"
    if(SubmitData$data_type == "RASL") dataType <- "RASL-seq"
    return(dataType)
  })

  # DDA - checkbox to show/hide row and column labels
  output$showLabels <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type == "Single-cell") return()
    if(is.null(corout())) return()
    prettyCheckbox(
      inputId = "showLabels",
      label = "Show colum/row labels (may not align correctly for large heatmaps)",
      value = F,
      width = "100%",
      status = "default", 
      icon = icon("check")
    )
  })
  observeEvent(input$showLabels, {
    d$correlation_showLabels = input$showLabels
  })

  # DDA - correlation heatmap downsample message
  output$cor_downsample_msg <- renderUI({
    req(SubmitData$data_type, input$cor_nGenes, input$cor_nSamples)
    if(SubmitData$data_type %in% c("Bulk", "RASL")) {
      sampleType <- "samples"
      if(is.null(ddstran())) return()
      nGenes <- nrow(ddstran()[[1]])
      nSamples <- ncol(ddstran()[[1]])
    } else {
      sampleType <- "cells"
      if(is.null(seurat_only())) return()
      if(d$resType == "seurat_res") {
        nGenes <- nrow(seurat_only()[[1]])
        nSamples <- ncol(seurat_only()[[1]])
      } else {
        nGenes <- nrow(seurat_only())
        nSamples <- ncol(seurat_only())
      }
    }
    nSamplesDownsample <- min(500, input$cor_nSamples)
    if(nGenes == input$cor_nGenes & nSamples == nSamplesDownsample) return()
    if(input$cor_nGenes < nGenes & nSamplesDownsample < nSamples) {
      txt1 <- paste0("Downsampling genes from ", nGenes, " to ", input$cor_nGenes, ".")
      txt2 <- paste0("Downsampling ", sampleType, " from ", nSamples, " to ", nSamplesDownsample, ".")
      HTML(paste(txt1, txt2, sep = "<br/>"))
    } else if(input$cor_nGenes < nGenes) {
      txt <- paste0("Downsampling genes from ", nGenes, " to ", input$cor_nGenes, ".")
      HTML(txt)
    } else {
      txt <- paste0("Downsampling ", sampleType, " from ", nSamples, " to ", nSamplesDownsample, ".")
      HTML(txt)
    }
  })

  # DDA - color schemes for correlation matrix
  getColors <-  eventReactive(input$build_cor_matrix, {
    if (input$plotColors == "Viridis") {
      return(viridis(100))
    } else if (input$plotColors == "Green-yellow-red") {
      return(rev(brewer.pal(n=11, name = "RdYlGn")))
    } else {
      return(colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100))
    }
  })

  # DDA - header for correlation heatmap
  output$headcor <- renderUI({
    req(input$build_cor_matrix, !d$newCorplot, SubmitData$data_type)
    h3("Correlation heatmap")
  })

  # DDA - correlation matrix
  output$corplot1 <- renderPlotly({
    req(input$build_cor_matrix, !d$newCorplot, SubmitData$data_type)
    if(is.null(input$build_cor_matrix) || input$build_cor_matrix == 0 ||
       is.null(corout())) return()
    if(is.null(d$correlation_showLabels)) return()
    withProgress(message = "Creating correlation heatmap...", value = 0, {
      incProgress(1/3)

      # Add sample/cell names to tooltip
      tooltips <- paste0(
        "<b>R:</b> ", round(corout()[[1]], 3)
      )
      tooltips <- matrix(
        tooltips,
        ncol = ncol(corout()[[1]]),
        byrow = FALSE
      )
      dimnames(tooltips) <- dimnames(corout()[[1]])

      if(nrow(tooltips) <= 100) {
        tooltips2 <- outer(rownames(tooltips), rownames(tooltips), function(x, y) {
          paste(x, y, sep = "<br /><b>Y:</b> ")
        })
        tooltips[] <- paste0(tooltips, "<br /><b>X:</b> ", tooltips2)
      }
      
      xLab <- yLab <- NULL
      myMargins <- list(r = 450, b = 150)
      axisTitleText <- ifelse(SubmitData$data_type %in% c("Bulk", "RASL"), yes = "Samples", no = "Cells")
      axisTitle <- list(text = axisTitleText, standoff = 10, font = list(size = 20))
      xLayout <- list(title = axisTitle)
      yLayout <- list(title = axisTitle)
      if(SubmitData$data_type %in% c("Bulk", "RASL")) {
        if(d$correlation_showLabels) {
          xLab <- colnames(corout()[[1]])
          yLab <- rownames(corout()[[1]])
          maxChar <- max(nchar(yLab), na.rm = T)
          myMargins <- list(l = maxChar*5, r = 200)
          xLayout <- list(title = axisTitle, 
                          tickangle = 30, 
                          tickfont = list(size = 10))
          yLayout <- list(title = axisTitle, 
                          tickfont = list(size = 10))
        }
      }
      incProgress(1/3)

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
        x = xLab,
        y = yLab,
        z = corout()[[1]],
        colors = getColors(),
        type = "heatmap",
        text = tooltips,
        hoverinfo = "text",
        hoverlabel = label,
        source = "corplot",
      ) %>%
        layout(
          xaxis = xLayout,
          yaxis = yLayout,
          margin = myMargins,
          font = list(family = "Noto Sans JP")
          ) %>%
        colorbar(limits = c(-1,1))
      
    })
  })

  # DDA - download button for corplot1 (pdf)
  output$dlqccorplot1pdf <- renderUI({
    req(corout(), input$build_cor_matrix, !d$newCorplot)

    if(is.null(input$build_cor_matrix)) return()
    if(input$build_cor_matrix == 0) return()
    if(is.null(corout())) return()
    downloadButton("dlqccorplot1pdfimg", "Download static plot (PDF)")
  })

  # DDA - download file for corplot1 (pdf)
  output$dlqccorplot1pdfimg <- downloadHandler(
    filename =  function() {
      paste("cor-matrix.pdf")
    },
    content = function(file) {
      cor.mat <- corout()[[1]]
      cor.mat <- data.frame(cor.mat) %>%
        rownames_to_column("samp") %>%
        gather(sample, value, -samp)
      
      xLab <- colnames(corout()[[1]])
      yLab <- rownames(corout()[[1]])

      if(input$plotColors == "Viridis") {
        cols <- viridis(100)
      } else if(input$plotColors == "Green-yellow-red") {
        cols <- rev(brewer.pal(n=11, name = "RdYlGn"))
      } else {
        cols <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)
      }

      showLabels <- F
      if(SubmitData$data_type %in% c("Bulk", "RASL")) {
        if(d$correlation_showLabels) {
          showLabels <- T
        }
      }
      if(showLabels) {
        p <- ggplot(cor.mat, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples") +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text.x = element_text(angle = 90, vjust = .4)) +
          scale_fill_gradientn(colours = cols, limits = c(-1, 1))
      } else {
        p <- ggplot(cor.mat, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples") +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()) +
          scale_fill_gradientn(colours = cols, limits = c(-1, 1))
      }

      pdf(file, width = 7, height = 6.5)
      print(p)
      dev.off()
    }
  )

  # DDA - download button for corplot1 (png)
  output$dlqccorplot1png <- renderUI({
    req(corout(), input$build_cor_matrix, !d$newCorplot)
    
    
    if(is.null(input$build_cor_matrix)) return()
    if(input$build_cor_matrix == 0) return()
    if(is.null(corout())) return()
    downloadButton("dlqccorplot1pngimg", "Download static plot (PNG)")
  })

  # DDA - download file for corplot1 (png)
  output$dlqccorplot1pngimg <- downloadHandler(
    filename =  function() {
      paste("cor-matrix.png")
    },
    content = function(file) {
      cor.mat <- corout()[[1]]
      cor.mat <- data.frame(cor.mat) %>%
        rownames_to_column("samp") %>%
        gather(sample, value, -samp)

      if (input$plotColors == "Viridis") {
        cols <- viridis(100)
      } else if (input$plotColors == "Green-yellow-red") {
        cols <- rev(brewer.pal(n=11, name = "RdYlGn"))
      } else {
        cols <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)
      }

      showLabels <- F
      if(SubmitData$data_type %in% c("Bulk", "RASL")) {
        if(d$correlation_showLabels) {
          showLabels <- T
        }
      }
      if(showLabels) {
        p <- ggplot(cor.mat, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples") +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text.x = element_text(angle = 90, vjust = .4)) +
          scale_fill_gradientn(colours = cols, limits = c(-1, 1))
      } else {
        p <- ggplot(cor.mat, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples") +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()) +
          scale_fill_gradientn(colours = cols, limits = c(-1, 1))
      }

      png(file, width = 800, height = 750)
      print(p)
      dev.off()
    }
  )

  # DDA - header for correlation scatterplot
  output$corplot2Header <- renderUI({
    req(d$corplotClick)
    h3("Correlation scatterplot")
  })

  # DDA - store correlation plot point click in object d
  observe({
    d$corplotClick <- event_data("plotly_click", source = "corplot")
  })

  # DDA - correlation sample vs sample (cell vs cell) scatterplot
  output$corplot2 <- renderPlotly({
    req(corout(), d$corplotClick)
    
    
    if(is.null(corout())) return()
    if(is.null(d$corplotClick)) return()
    withProgress(message = "Rendering count plot...", value = 0, {
      incProgress()
      s.cor <- d$corplotClick
      validate(
        need(
          s.cor != "",
          message = "Click on one of the heatmap cells to view this plot!"
        )
      )
      cts.tran <- corout()[[2]]
      cts.tran <- as.data.frame(cts.tran)
      z <- s.cor[["z"]]
      # x and y are actually reversed, so need to flip them
      x <- s.cor[["y"]]
      y <- s.cor[["x"]]
      if(is.integer(x)) x <- colnames(cts.tran)[x]
      if(is.integer(y)) y <- colnames(cts.tran)[y]
      myTitle <- paste0(x, " vs. ", y)
      if(sum(nchar(c(x,y))) > 50) myTitle <- paste0(x, " vs.\n", y)

      # Add gene expression to tooltip
      geneLabs <- paste0("<b>Gene:</b> ", rownames(cts.tran))
      xLab <- paste0("<b>", x, ":</b> ")
      xLabs <- paste0(xLab, formatC(cts.tran[, x]))
      yLab <- paste0("<b>", y, ":</b> ")
      yLabs <- paste0(yLab, formatC(cts.tran[, y]))

      tooltips <- paste(geneLabs, xLabs, yLabs, sep = "<br />")

      m <- list(
        l = 50,
        r = 50,
        b = 100,
        t = 100,
        pad = 4
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
        data = cts.tran,
        type = "scatter",
        mode = "markers",
        x = cts.tran[, x],
        y = cts.tran[, y],
        text = tooltips,
        marker = list(size = 2),
        hoverinfo = "text",
        hoverlabel = label
      ) %>%
        layout(
          title = myTitle,
          font = list(size = 12, family = "Noto Sans JP"),
          xaxis = list(title = x),
          yaxis = list(title = y),
          margin = m
        )
    })
  })

  # DDA - download button corplot2 (pdf)
  output$dlqccorplot2pdf <- renderUI({
    
    if(is.null(d$corplotClick)) {
      return()
    } else {
      downloadButton(
        "dlqccorplot2pdfimg",
        "Download static plot (PDF)"
      )
    }
  })

  # DDA - download file corplot2 (pdf)
  output$dlqccorplot2pdfimg <- downloadHandler(
    filename =  function() {
      paste("cor-scatterplot.pdf")
    },
    content = function(file) {
      cts.tran <- corout()[[2]]
      cts.tran <- as.data.frame(cts.tran)
      s.cor <- event_data("plotly_click", source = "corplot")
      z = s.cor[["z"]]
      x <- s.cor[["x"]]
      y <- s.cor[["y"]]
      if(is.integer(x)) x <- colnames(cts.tran)[x]
      if(is.integer(y)) y <- colnames(cts.tran)[y]
      myTitle <- paste0(y, " vs. ", x)
      if(sum(nchar(c(x,y))) > 50) myTitle <- paste0(x, " vs.\n", y)

      pdf(file, width = 8, height = 6.5)
      p <- ggplot(cts.tran, aes(x = cts.tran[, x], y = cts.tran[, y])) +
        geom_point(color = "dodgerblue4") + theme_classic() +
        xlab(x) + ylab(y) + ggtitle(myTitle) +
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(angle = 90, vjust = .4))
      print(p)
      dev.off()
    }
  )

  # DDA - download button corplot2 (png)
  output$dlqccorplot2png <- renderUI({
    if(is.null(d$corplotClick)) return()
    downloadButton(
      "dlqccorplot2pngimg",
      "Download static plot (PNG)"
    )
  })

  # DDA - download file corplot2 (png)
  output$dlqccorplot2pngimg <- downloadHandler(
    filename =  function() {
      paste("cor-scatterplot.png")
    },
    content = function(file) {
      cts.tran <- corout()[[2]]
      cts.tran <- as.data.frame(cts.tran)
      s.cor <- event_data("plotly_click", source = "corplot")
      z = s.cor[["z"]]
      x <- s.cor[["x"]]
      y <- s.cor[["y"]]
      if(is.integer(x)) x <- colnames(cts.tran)[x]
      if(is.integer(y)) y <- colnames(cts.tran)[y]
      myTitle <- paste0(y, " vs. ", x)
      if(sum(nchar(c(x,y))) > 50) myTitle <- paste0(x, " vs.\n", y)

      png(file, width = 800, height = 750)
      p <- ggplot(cts.tran, aes(x = cts.tran[, x], y = cts.tran[, y])) +
        geom_point(color = "dodgerblue4") + theme_classic() +
        xlab(x) + ylab(y) + ggtitle(myTitle) +
        theme(text = element_text(family = "noto-sans-jp"),
              axis.text.x = element_text(angle = 90, vjust = .4))
      print(p)
      dev.off()
    }
  )

  # DDA - header correlation analysis
  output$headcor2 <- renderUI({
    req(SubmitData$data_type, distout())
    if(SubmitData$data_type == "Single-cell") return()
    h3("Sample distance matrix")
  })

  # DDA - checkbox to show/hide row and column labels
  output$dist_labels <- renderUI({
    req(SubmitData$data_type, distout())
    if(!(SubmitData$data_type %in% c("Bulk", "RASL"))) return()
    prettyCheckbox(
      inputId = "dist_labels",
      label = "Show colum/row labels (may not align correctly for large heatmaps)",
      value = F, status = "default", icon = icon("check")
    )
  })
  observeEvent(input$dist_labels, {
    d$distance_showLabels <- input$dist_labels
  })
  
  # DDA - distance heatmap downsample message
  output$dist_downsample_msg <- renderUI({
    req(SubmitData$data_type, input$dist_nGenes, input$dist_nSamples)
    if(SubmitData$data_type %in% c("Bulk", "RASL")) {
      sampleType <- "samples"
      if(is.null(ddstran())) return()
      nGenes <- nrow(ddstran()[[1]])
      nSamples <- ncol(ddstran()[[1]])
    } else {
      sampleType <- "cells"
      if(is.null(seurat_only())) return()
      if(d$resType == "seurat_res") {
        nGenes <- nrow(seurat_only()[[1]])
        nSamples <- ncol(seurat_only()[[1]])
      } else {
        nGenes <- nrow(seurat_only())
        nSamples <- ncol(seurat_only())
      }
    }
    nSamplesDownsample <- min(500, input$dist_nSamples)
    if(nGenes == input$dist_nGenes && nSamples == nSamplesDownsample) return()
    if(input$dist_nGenes < nGenes && nSamplesDownsample < nSamples) {
      txt1 <- paste0("Downsampling genes from ", nGenes, " to ", input$dist_nGenes, ".")
      txt2 <- paste0("Downsampling ", sampleType, " from ", nSamples, " to ", nSamplesDownsample, ".")
      HTML(paste(txt1, txt2, sep = "<br/>"))
    } else if(input$dist_nGenes < nGenes) {
      txt <- paste0("Downsampling genes from ", nGenes, " to ", input$dist_nGenes, ".")
      HTML(txt)
    } else {
      txt <- paste0("Downsampling ", sampleType, " from ", nSamples, " to ", nSamplesDownsample, ".")
      HTML(txt)
    }
  })

  # DDA - distance matrix (bulk only)
  output$corplot3 <- renderPlotly({
    req(input$build_dist_matrix, !d$newDistplot, getColors_dist(), SubmitData$data_type, distout())
    if(is.null(input$build_dist_matrix) || input$build_dist_matrix == 0 ||
       is.null(d$distance_showLabels) || SubmitData$data_type != "Bulk") return()

    sdm <- distout()[["sdm"]]
    tooltips <- distout()[["tooltips"]]
    
    xLab <- yLab <- NULL
    myMargins <- list(r = 450, b = 150)
    axisTitleText <- "Samples"
    axisTitle <- list(text = axisTitleText, standoff = 10, font = list(size = 20))
    xLayout <- list(title = axisTitle)
    yLayout <- list(title = axisTitle)
    if(d$distance_showLabels) {
      xLab <- colnames(sdm)
      yLab <- rownames(sdm)
      maxChar <- max(nchar(yLab), na.rm = T)
      myMargins <- list(l = maxChar*5, r = 200)
      xLayout <- list(title = axisTitle, 
                      tickangle = 30, 
                      tickfont = list(size = 10))
      yLayout <- list(title = axisTitle, 
                      tickfont = list(size = 10))
    }    
    cols <- getColors_dist()
    
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
      x = xLab,
      y = yLab,
      z = sdm,
      colors = cols,
      type = "heatmap",
      text = tooltips,
      hoverinfo = "text", 
      hoverlabel = label
    ) %>%
      layout(
        xaxis = xLayout,
        yaxis = yLayout,
        margin = myMargins,
        font = list(family = "Noto Sans JP")
      )
  })

  # DDA - download button - sample distance matrix (PDF)
  output$dlqcorplot3pdf <- renderUI({
    req(SubmitData$data_type, distout(), getColors_dist())
    if(!(SubmitData$data_type %in% c("Bulk", "RASL"))) return()
    downloadButton("dlqcorplot3pdfimg", "Download plot (PDF)")
  })

  # DDA - download file - sample distance matrix (PDF)
  output$dlqcorplot3pdfimg <- downloadHandler(
    filename =  function() {
      return("sample-dists.pdf")
    },
    content = function(file) {
      sdm <- distout()[["sdm"]]
      sdm <- data.frame(sdm) %>%
        rownames_to_column("samp") %>%
        gather(sample, value, -samp)
      cols <- getColors_dist()

      if(d$distance_showLabels) {
        p <- ggplot(sdm, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples")  +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text.x = element_text(angle = 90, vjust = .4),
                text = element_text(size = 8)) +
          scale_fill_gradientn(colours = cols)
      } else {
        pdf(file, width = 7, height = 6.5)
        p <- ggplot(sdm, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples")  +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()) +
          scale_fill_gradientn(colours = cols)
      }
      pdf(file, width = 8, height = 6.5)
      print(p)
      dev.off()
    }
  )

  # DDA - download button - sample distance matrix (PNG)
  output$dlqcorplot3png <- renderUI({
    req(SubmitData$data_type, distout(), getColors_dist())
    if(!(SubmitData$data_type %in% c("Bulk", "RASL"))) return() 
    downloadButton("dlqcorplot3pngimg", "Download plot (PNG)")
  })

  # DDA - download file - sample distance matrix (PNG)
  output$dlqcorplot3pngimg <- downloadHandler(
    filename =  function() {
      return("sample-dists.png")
    },
    content = function(file) {
      sdm <- distout()[["sdm"]]
      sdm <- data.frame(sdm) %>%
        rownames_to_column("samp") %>%
        gather(sample, value, -samp)
      cols <- getColors_dist()
      
      if(d$distance_showLabels) {
        p <- ggplot(sdm, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples")  +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text.x = element_text(angle = 90, vjust = .4, size = 10),
                axis.text.y = element_text(size = 12)) +
          scale_fill_gradientn(colours = cols)
      } else {
        p <- ggplot(sdm, aes(x = sample, y = samp)) +
          geom_tile(aes(fill = value)) + theme_classic() +
          xlab("Samples") + ylab("Samples")  +
          theme(text = element_text(family = "noto-sans-jp"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()) +
          scale_fill_gradientn(colours = cols)
      }
      png(file, width = 800, height = 750)
      print(p)
      dev.off()
    }
  )
  
  ###################################################################
  ###################################################################
  ### SECTION 04 - DIFFERENTIAL GENE EXPRESSION ANALYSIS (DGE)
  ###################################################################
  ###################################################################

  ### SECTION 03 - BULK DGE
    
  ######### HELP BUTTONS - ALL FOR BULK RNA-SEQ######################
    
  # BULK-DGE-OVER - overview help button
  observeEvent(input$overview_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Differential gene expression")),
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
      includeMarkdown("markdown/help/runDGE_help.md")
    ))
  })

  # BULK-DGE-VOL - volcano plot help button
  observeEvent(input$volcanoplots_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Volcano plot")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "m",
      easyClose = T,
      includeMarkdown("markdown/help/volcano_help.md")
    ))
  })
  
  # BULK-DGE-GSE - GSE help button
  observeEvent(input$gse_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Gene set enrichment")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html#gene-set-enrichment", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/gse_help.md")
    ))
  })

  # BULK-DGE-CLUS - clustering help button
  observeEvent(input$postdge_clustering_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Clustering")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "m",
      easyClose = T,
      includeMarkdown("markdown/help/clustering_help.md")
    ))
  })
  
  ######### OVERVIEW TAB CODE: BULK RNA-SEQ######################
    
  # BULK-DGE-OVER - select algorithm for dge
  output$dgeexpsetup <- renderUI({
    selectInput(
      inputId = "dgeexpsetup",
      label = "Experimental design",
      choices = c(
        "Two-group comparisons" = "exp1",
        "Multiple factor comparisons (factorial)" = "exp2",
        "Classical interaction design" = "exp3",
        "Additive models - paired or blocking" = "exp4",
        "Main effects" = "exp5",
        "Main effects with grouping factors" = "exp6"
      ),
      selected = "exp1",
      multiple = FALSE
    )
  })
  
  # BULK-DGE-OVER - padj value cut-off
  output$dgepadjcutoff <- renderUI({
    textInput(
      inputId = "dgepadjcutoff",
      label = withMathJax("Adj. \\(p\\)-val cutoff"),
      value = 0.05
    )
  })
  
  # BULK-DGE-OVER - fold-change cut-off
  output$dgefcmin <- renderUI({
    textInput(
      inputId = "dgefcmin",
      label = "Min. fold change",
      value = 1
    )
  })
  
  # BULK-DGE-OVER - choose analytical methods: DESeq2, edgeR, limma
  output$dgemethod <- renderUI({
    req(input$dgeexpsetup)
    myChoices <- c("DESeq2" = "deseq", "edgeR" = "edger", "limma-voom" = "limma")
    if(input$dgeexpsetup %in% c("exp5", "exp6")) {
      myChoices <- c("DESeq2" = "deseq", "edgeR" = "edger")
    }
    selectInput(
      inputId = "dgemethod",
      label = "DGE method",
      choices = myChoices,
      selected = "deseq"
    )
  })
  
  # Returns T/F whether all required inputs are filled for the chosen 
  # experimental design. 
  CheckExpDesign <- reactive({
    req(input$dgeexpsetup)
    nInputs <- switch(
      input$dgeexpsetup,
      exp1 = 3,
      exp2 = 6,
      exp3 = 2,
      exp4 = 2,
      exp5 = 4,
      exp6 = 2
    )
    
    requiredInputs <- switch(
      input$dgeexpsetup,
      exp1 = c(input$dgeexp1a, input$dgeexp1b, input$dgeexp1c),
      exp2 = c(input$dgeexp2a, input$dgeexp2b,
               input$dgeexp2c, input$dgeexp2c5, input$dgeexp2d, 
               input$dgeexp2d5),
      exp3 = c(input$dgeexp3a, input$dgeexp3b),
      exp4 = c(input$dgeexp4a, input$dgeexp4b),
      exp5 = c(input$dgeexp5a, input$dgeexp5b, input$dgeexp5c, input$dgeexp5d),
      exp6 = c(input$dgeexp6a, input$dgeexp6c)
    )
    if(length(requiredInputs) != nInputs) return(F)
    if(any(requiredInputs == "")) return(F)
    return(T)
  })
  
  # Warning message if data is insufficient for experimental design
  output$dge_expdesign_warning <- renderText("Data insufficient for this design.")
  output$dge_expdesign_warning_ui <- renderUI({
    if(is.null(CheckExpDesign())) return()
    if(CheckExpDesign()) return()
    span(textOutput("dge_expdesign_warning"), style = "color:red;")
  })
  
  # BULK-DGE-OVER - action button to submit DGE analysis
  output$godge <- renderUI({
    if(is.null(CheckExpDesign()) || !CheckExpDesign() ||
       is.null(IsFullRank()) || !IsFullRank() ||
       is.null(EnoughSamples()) || !EnoughSamples()) return()
    actionButton(inputId = "godge", label = "Submit", icon = icon("space-shuttle"))
  })
  
  # BULK-DGE-DE - dge warning
  output$dge_warning_msg <- renderText(d$dge_warning_msg)
  output$dge_warning_msg_ui <- renderUI({
    if(d$dge_warning_msg == "") return()
    span(textOutput("dge_warning_msg"), style = "color:red;")
  })
  
  # BULK-DGE-OVER - batch correction factor
  output$batchFactor <- renderUI({
    req(ddsout())
    myChoices <- c("None", colnames(ddsout()[[2]]))
    mySelected <- "None"
    if("Dataset" %in% myChoices) mySelected <- "Dataset"
    selectInput(
      inputId = "batchFactor", 
      label = "Batch correction factor", 
      choices = myChoices, 
      selected = mySelected,
    )
  })
  
  # BULK-DGE-OVER - whether to include housekeeping genes effect in DE model
  output$housekeepingSelect <- renderUI({
    req(ddsout())
    prettyCheckbox(
      inputId = "housekeepingSelect",
      label = "Control for housekeeping genes",
      value = FALSE,
      status = "default",
      icon = icon("check")
    )
  })
  
  # BULK-DGE-OVER - exp. setup 1 - two group comparisons - factor choice
  output$dgeexp1a <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp1") return()
    myChoices <- colnames(ddsout()[[2]])
    myChoices <- myChoices[myChoices != input$batchFactor]
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp1a",
      label = "Factor",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 1 - two group comparisons - choose comp a
  output$dgeexp1b <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp1") return()
    if(is.null(input$dgeexp1a) || !(input$dgeexp1a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp1a]])))
    }
    selectInput(
      inputId = "dgeexp1b", 
      label = "Group 1",
      choices = myChoices 
    )
  })
  
  # BULK-DGE-OVER - exp. setup 1 - two group comparisons - choose comp b
  output$dgeexp1c <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp1") return()
    if(is.null(input$dgeexp1a) || !(input$dgeexp1a %in% colnames(ddsout()[[2]])) ||
       is.null(input$dgeexp1b)) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp1a]])))
      myChoices <- myChoices[myChoices != input$dgeexp1b]
    }
    if(length(myChoices) == 0) {
      myChoices <- ""
      mySelected <- NULL
    } else {
      myChoices <- c("Rest", myChoices)
      mySelected <- "Rest"
    }
    selectInput(
      inputId = "dgeexp1c", 
      label = "Group 2",
      choices = myChoices,
      selected = mySelected
    )
  })
  
  # BULK-DGE-OVER - exp. setup 2 - mult. group comparisons - factor choice A
  output$dgeexp2a <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp2") return()
    myChoices <- colnames(ddsout()[[2]])
    myChoices <- myChoices[myChoices != input$batchFactor]
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp2a",
      label = "Factor A",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 2 - mult. group comparisons - group comb. choices
  output$dgeexp2c <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp2") return()
    if(is.null(input$dgeexp2a) || !(input$dgeexp2a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp2a]])))
    }
    selectInput(
      inputId = "dgeexp2c", 
      label = "Factor A group 1",
      choices = myChoices 
    )
  })
  
  # BULK-DGE-OVER - exp. setup 2 - mult. group comparisons - group comb. choices
  output$dgeexp2c5 <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp2") return()
    if(is.null(input$dgeexp2a) || !(input$dgeexp2a %in% colnames(ddsout()[[2]])) ||
       is.null(input$dgeexp2c)) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp2a]])))
      myChoices <- myChoices[myChoices != input$dgeexp2c]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp2c5",
      label = "Factor A group 2",
      choices = myChoices,
    )
  })
  
  # BULK-DGE-OVER - exp. setup 2 - mult. group comparisons - factor choice B
  output$dgeexp2b <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp2") return()
    myChoices <- colnames(ddsout()[[2]])
    if(is.null(input$dgeexp2a)) {
      myChoices <- ""
    } else {
      myChoices <- myChoices[myChoices != input$dgeexp2a]
      myChoices <- myChoices[myChoices != input$batchFactor]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp2b",
      label = "Factor B",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 2 - mult. group comparisons - group comb. choices
  output$dgeexp2d <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp2") return()
    if(is.null(input$dgeexp2b) || !(input$dgeexp2b %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp2b]])))
      myChoices <- myChoices[myChoices != input$batchFactor]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp2d", 
      label = "Factor B group 1",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 2 - mult. group comparisons - group comb. choices
  output$dgeexp2d5 <- renderUI({
    req(ddsout())
    if(is.null(input$dgeexp2b) || !(input$dgeexp2b %in% colnames(ddsout()[[2]])) ||
       is.null(input$dgeexp2d)) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp2b]])))
      myChoices <- myChoices[myChoices != input$dgeexp2d]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp2d5",
      label = "Factor B group 2",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 3 - interaction - factor choice A
  output$dgeexp3a <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp3") return()
    myChoices <- colnames(ddsout()[[2]])
    myChoices <- myChoices[myChoices != input$batchFactor]
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp3a",
      label = "Factor A",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 3 - interaction - factor choice B
  output$dgeexp3b <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp3") return()
    if(is.null(input$dgeexp3a) || !(input$dgeexp3a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- colnames(ddsout()[[2]])
      myChoices <- myChoices[myChoices != input$dgeexp3a]
      myChoices <- myChoices[myChoices != input$batchFactor]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp3b",
      label = "Factor B",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 3 - interaction - reference level for factor A
  output$dgeexp3c <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp3") return()
    if(is.null(input$dgeexp3a) || !(input$dgeexp3a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp3a]])))
    }
    selectInput(
      inputId = "dgeexp3c",
      label = "Factor A reference",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 3 - interaction - reference level for factor B
  output$dgeexp3d <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp3") return()
    if(is.null(input$dgeexp3b) || !(input$dgeexp3b %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp3b]])))
    }
    selectInput(
      inputId = "dgeexp3d",
      label = "Factor B reference",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 4 - additive model - blocking factor
  output$dgeexp4a <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp4") return()
    myChoices <- colnames(ddsout()[[2]])
    myChoices <- myChoices[myChoices != input$batchFactor]
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp4a",
      label = "Blocking factor",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 4 - additive model - treatment factor
  output$dgeexp4b <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp4") return()
    if(is.null(input$dgeexp4a) || !(input$dgeexp4a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- colnames(ddsout()[[2]])
      myChoices <- myChoices[myChoices != input$dgeexp4a]
      myChoices <- myChoices[myChoices != input$batchFactor]
    }
    if(length(myChoices) == 0) myChoices <- ""
    
    selectInput(
      inputId = "dgeexp4b",
      label = "Treatment factor",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 4 - additive model - reference level for blocking factor
  output$dgeexp4c <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp4") return()
    if(is.null(input$dgeexp4a) || !(input$dgeexp4a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp4a]])))
    }
    selectInput(
      inputId = "dgeexp4c",
      label = "Blocking factor reference",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 4 - additive model - reference level for treatment factor
  output$dgeexp4d <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp4") return()
    if(is.null(input$dgeexp4b) || !(input$dgeexp4b %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp4b]])))
    }
    selectInput(
      inputId = "dgeexp4d",
      label = "Treatment factor reference",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 5 - main effect (ME) - choose ME
  output$dgeexp5a <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp5") return()
    myChoices <- colnames(ddsout()[[2]])
    myChoices <- myChoices[myChoices != input$batchFactor]
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp5a",
      label = "Main effect factor",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 5 - ME - choose ME reference
  output$dgeexp5b <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp5") return()
    if(is.null(input$dgeexp5a) || !(input$dgeexp5a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp5a]])))
    }
    
    selectInput(
      inputId = "dgeexp5b",
      label = "Main effect reference",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 5 - ME - choose contrasts
  output$dgeexp5c <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp5") return()
    if(is.null(input$dgeexp5a) || !(input$dgeexp5a %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp5a]])))
    }
    selectInput(
      inputId = "dgeexp5c", 
      label = "Group 1",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 5 - ME - choose contrasts
  output$dgeexp5d <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp5") return()
    if(is.null(input$dgeexp5a) || !(input$dgeexp5a %in% colnames(ddsout()[[2]])) ||
       is.null(input$dgeexp5c)) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp5a]])))
      myChoices <- myChoices[myChoices != input$dgeexp5c]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp5d", 
      label = "Group 2",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 6 - ME + group fact - choose ME
  output$dgeexp6a <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(input$dgeexp6c) || !(input$dgeexp6c %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- colnames(ddsout()[[2]])
      myChoices <- myChoices[myChoices != input$batchFactor]
      myChoices <- myChoices[myChoices != input$dgeexp6c]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp6a",
      label = "Grouping factor",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 6 - ME + group fact - choose group fact
  output$dgeexp6b <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(input$dgeexp6c) || !(input$dgeexp6c %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      metaDF <- ddsout()[[2]]
      if(is.null(input$dgeexp6d) || sum(metaDF[[input$dgeexp6c]] == input$dgeexp6d) < 3) {
        myChoices <- ""
      } else {
        indices <- which(metaDF[[input$dgeexp6c]] == input$dgeexp6d)
        myChoices <- mixedsort(unique(as.character(metaDF[[input$dgeexp6a]][indices])))
      }
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp6b",
      label = "Grouping factor level",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 6 - ME + group fact - choose ME reference
  output$dgeexp6c <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp6") return()
    myChoices <- colnames(ddsout()[[2]])
    myChoices <- myChoices[myChoices != input$batchFactor]
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp6c",
      label = "Main effect factor",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 6 - ME + group fact - choose group fact level
  output$dgeexp6d <- renderUI({
    req(ddsout(), input$dgeexpsetup)
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(input$dgeexp6c) || !(input$dgeexp6c %in% colnames(ddsout()[[2]]))) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp6c]])))
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp6d",
      label = "Main effect reference",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 5 - ME - choose contrasts
  output$dgeexp6e <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(input$dgeexp6c) || !(input$dgeexp6c %in% colnames(ddsout()[[2]])) ||
       is.null(input$dgeexp6d)) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp6c]])))
      myChoices <- myChoices[myChoices != input$batchFactor]
      myChoices <- myChoices[myChoices != input$dgeexp6d]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp6e", 
      label = "Group 1",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup 5 - ME - choose contrasts
  output$dgeexp6f <- renderUI({
    req(ddsout(), input$dgeexpsetup, input$batchFactor)
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(input$dgeexp6c) || !(input$dgeexp6c %in% colnames(ddsout()[[2]])) ||
       is.null(input$dgeexp6e)) {
      myChoices <- ""
    } else {
      myChoices <- mixedsort(unique(as.character(ddsout()[[2]][[input$dgeexp6c]])))
      myChoices <- myChoices[myChoices != input$batchFactor]
      myChoices <- myChoices[myChoices != input$dgeexp6e]
    }
    if(length(myChoices) == 0) myChoices <- ""
    selectInput(
      inputId = "dgeexp6f", 
      label = "Group 2",
      choices = myChoices
    )
  })
  
  # BULK-DGE-OVER - exp. setup - formula - header
  output$dgeexpformhead <- renderUI({
    
    if (is.null(input$dgeexpsetup)) return()
    if (input$dgeexpsetup != "exp1" & ncol(ddsout()[[2]]) < 2) return()
    if (is.null(CheckExpDesign())) return()
    if (!CheckExpDesign()) return()
    h5(strong("Linear model:"), style = "margin-top: 0px;")
  })
  
  # BULK-DGE-OVER - exp. setup - formula - formula 1
  output$dgeexpform1 <- renderUI({
    req(input$dgeexpsetup, input$batchFactor)
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp1") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    batchString <- ""
    if(input$housekeepingSelect) batchString <- "Batch + "
    if(input$batchFactor != "None") batchString <- paste0(batchString, input$batchFactor, " + ")
    
    code(paste0(" ~ ", batchString, input$dgeexp1a))
  })
  
  # BULK-DGE-OVER - exp. setup - formula - formula 2
  output$dgeexpform2 <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp2") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    batchString <- "" 
    tmp <- colnames(ddsout()[[2]])
    if(input$housekeepingSelect) batchString <- "Batch + "
    if(input$batchFactor != "None") {
      batchString <- paste0(batchString, input$batchFactor, " + ")
      tmp <- tmp[tmp != input$batchFactor]
    }
    
    if(length(tmp) < 2) return()
    code(paste0(" ~ ", batchString, input$dgeexp2a, "_", input$dgeexp2b))
  })
  
  # BULK-DGE-OVER - exp. setup 3 - interaction - formula layout
  output$dgeexpform3 <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp3") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    batchString <- "" 
    tmp <- colnames(ddsout()[[2]])
    if(input$housekeepingSelect) batchString <- "Batch + "
    if(input$batchFactor != "None") {
      batchString <- paste0(batchString, input$batchFactor, " + ")
      tmp <- tmp[tmp != input$batchFactor]
    }
    
    if(length(tmp) < 2) return()
    code(
      paste0(
        " ~ ",
        batchString,
        input$dgeexp3a,
        " + ",
        input$dgeexp3b,
        " + ",
        input$dgeexp3a,
        ":",
        input$dgeexp3b
      )
    )
  })
  
  # BULK-DGE-OVER - exp. setup 4 - added effects - formula layout
  output$dgeexpform4 <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp4") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    batchString <- "" 
    tmp <- colnames(ddsout()[[2]])
    if(input$housekeepingSelect) batchString <- "Batch + "
    if(input$batchFactor != "None") {
      batchString <- paste0(batchString, input$batchFactor, " + ")
      tmp <- tmp[tmp != input$batchFactor]
    }
    
    if(length(tmp) < 2) return()
    code(paste0(" ~ ", batchString, input$dgeexp4a, " + ", input$dgeexp4b))
  })
  
  # BULK-DGE-OVER - exp. setup 5 - ME - formula layout
  output$dgeexpform5 <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp5") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    batchString <- "" 
    if(input$housekeepingSelect) batchString <- "Batch + "
    if(input$batchFactor != "None") {
      batchString <- paste0(batchString, input$batchFactor, " + ")
    }
    code(paste0(" ~ ", batchString, input$dgeexp5a, " (as main effect)"))
  })
  
  # BULK-DGE-OVER - exp. setup 6A - ME + group fact. - formula layout
  output$dgeexpform6a <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    batchString <- "" 
    if(input$housekeepingSelect) batchString <- "Batch + "
    if(input$batchFactor != "None") {
      batchString <- paste0(batchString, input$batchFactor, " + ")
    }
    code(paste0(" ~ ", batchString, input$dgeexp6c, " (as main effect)"))
  })
  
  # BULK-DGE-OVER - exp. setup 6B - ME + group fact. - formula layout
  output$dgeexpform6b <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    code(paste0("* Limited by: "))
  })
  
  # BULK-DGE-OVER - exp. setup 6C - ME + group fact. - formula layout
  output$dgeexpform6c <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    code(paste("...factor: ", input$dgeexp6a))
  })
  
  # BULK-DGE-OVER - exp. setup 6D - ME + group fact. - formula layout
  output$dgeexpform6d <- renderUI({
    req(input$dgeexpsetup, input$batchFactor, ddsout())
    
    if(is.null(input$dgeexpsetup)) return()
    if(input$dgeexpsetup != "exp6") return()
    if(is.null(CheckExpDesign())) return()
    if(is.null(input$housekeepingSelect)) return()
    if(!CheckExpDesign()) return()    
    
    code(paste("...factor level: ", input$dgeexp6b))
  })
  
  # BULK-DGE-OVER - edgeR normalization option
  output$dgeexpedgernorm <- renderUI({
    req(input$dgemethod)
    if(input$dgemethod != "edger") return()
    selectInput(
      inputId = "dgeexpedgernorm",
      label = "Normalization type",
      choices = c(
        "TMM" = "TMM",
        "RLE" = "RLE",
        "upperquartile" = "upperquartile",
        "none" = "none"
      )
    )
  })
  
  # BULK-DGE-OVER - warning messages for data inputs for dge analysis
  observe({
    if(!is.null(EnoughSamples())) {
      if(EnoughSamples()) {
        d$newExpSetup <- F
      } else {
        d$dge_warning_msg <- "Not enough samples from factor levels."
        return()
      }
    }
    
    if(!is.null(IsFullRank())) {
      if(IsFullRank()) {
        d$newExpSetup <- F
      } else {
        d$dge_warning_msg <- "Experimental design is not full rank."
        return()
      }
    }
    
    d$dge_warning_msg <- ""
    
  })
  
  # BULK-DGE-OVER - store dge main contrasts in object d
  observeEvent(input$dgemaincontrasts, {
    if(is.null(d$dgeContrast) || input$dgemaincontrasts != d$dgeContrast) {
      d$dgeContrast <- input$dgemaincontrasts
    }
  })
  
  # BULK-DGE-OVER - store dge main contrasts gse in object d
  observeEvent(input$dgemaincontrasts_gse, {
    if(is.null(d$dgeContrast) || input$dgemaincontrasts_gse != d$dgeContrast) {
      d$dgeContrast <- input$dgemaincontrasts_gse
    }
  })
  
  # BULK-DGE-OVER - choose main contrasts for dge (sidepar panel)
  output$dgemaincontrasts <- renderUI({
    req(input$godge, dgeout1())
    myChoices <- colnames(dgeout1()[[2]])
    mySelected <- myChoices[1]
    if(!is.null(d$dgeContrast) && d$dgeContrast %in% myChoices) {
      mySelected <- d$dgeContrast
    }
    
    selectInput(
      inputId = "dgemaincontrasts",
      label = "Contrast",
      choices = myChoices,
      selected = mySelected,
      width = "400px"
    )
  })
  
  # BULK-DGE-OVER - user input for model matrix
  mod.matrix <- eventReactive(input$godge, {
    mod.matrix <- input$mod.matrix
    mod.matrix <- as.matrix(read.csv(
      mod.matrix$datapath,
      header = TRUE,
      row.names = 1
    ))
    return(list(mod.matrix))
  })
  
  # Returns T/F whether current selection for bulk DGE is full rank
  IsFullRank <- reactive({
    req(SubmitData$data_type, ddsout()[[2]], input$dgeexpsetup, input$batchFactor)
    if(SubmitData$data_type != "Bulk") return()
    myFactors <- switch(
      input$dgeexpsetup,
      exp1 = input$dgeexp1a,
      exp2 = c(input$dgeexp2a, input$dgeexp2b),
      exp3 = c(input$dgeexp3a, input$dgeexp3b),
      exp4 = c(input$dgeexp4a, input$dgeexp4b),
      exp5 = input$dgeexp5a,
      exp6 = c(input$dgeexp6a, input$dgeexp6c)
    )
    if(input$batchFactor != "None") myFactors <- c(myFactors, input$batchFactor)
    
    if(!all(myFactors %in% colnames(ddsout()[[2]]))) return()
    
    fullRankCheck <- CheckFullRank(ddsout()[[2]][, myFactors, drop = F])
    if(is.null(fullRankCheck)) return(T)
    return(F)
  })
  
  # Returns T/F whether current selection for bulk DGE has enough samples
  # in each combination of factor level in experimental design
  EnoughSamples <- reactive({
    req(SubmitData$data_type, ddsout()[[2]], CheckExpDesign(), input$dgeexpsetup, input$batchFactor)
    if(SubmitData$data_type != "Bulk") return()
    metaDF <- ddsout()[[2]]

    if(input$dgeexpsetup == "exp1") {
      if(input$dgeexp1c == "Rest") {
        metaDF[[input$dgeexp1a]] <- as.character(metaDF[[input$dgeexp1a]])
        metaDF[[input$dgeexp1a]][metaDF[[input$dgeexp1a]] != input$dgeexp1b] <- "Rest"
        metaDF[[input$dgeexp1a]] <- factor(metaDF[[input$dgeexp1a]], levels = c(input$dgeexp1b, "Rest"))
      }
      checkList <- list(c(input$dgeexp1b, input$dgeexp1c))
      names(checkList) <- input$dgeexp1a
    } else if(input$dgeexpsetup == "exp2") {
      checkList <- list(c(input$dgeexp2c, input$dgeexp2c5), c(input$dgeexp2d, input$dgeexp2d5))
      names(checkList) <- c(input$dgeexp2a, input$dgeexp2b)
    } else if(input$dgeexpsetup == "exp3") {
      checkList <- list(as.character(unique(metaDF[, input$dgeexp3a])), as.character(unique(metaDF[, input$dgeexp3b])))
      names(checkList) <- c(input$dgeexp3a, input$dgeexp3b)
    } else if(input$dgeexpsetup == "exp4") {
      checkList <- list(as.character(unique(metaDF[, input$dgeexp4a])), as.character(unique(metaDF[, input$dgeexp4b])))
      names(checkList) <- c(input$dgeexp4a, input$dgeexp4b)
    } else if(input$dgeexpsetup == "exp5") {
      checkList <- list(unique(c(input$dgeexp5b, input$dgeexp5c, input$dgeexp5d)))
      names(checkList) <- input$dgeexp5a
    } else if(input$dgeexpsetup == "exp6") {
      checkList <- list(input$dgeexp6b, as.character(unique(metaDF[(metaDF[, input$dgeexp6a] == input$dgeexp6b), input$dgeexp6c])))
      names(checkList) <- c(input$dgeexp6a, input$dgeexp6c)
    }
    if(!all(names(checkList) %in% colnames(metaDF))) return()
    for(myFact in names(checkList)) {
      if(!all(checkList[[myFact]] %in% metaDF[[myFact]])) return()
    }
    
    includedRownames <- rownames(metaDF)[metaDF[[names(checkList)[1]]] %in% checkList[[names(checkList)[1]]]]
    if(length(checkList) > 1) {
      for(fact in names(checkList)[-1]) {
        includedRownames <- intersect(includedRownames, rownames(metaDF)[metaDF[[fact]] %in% checkList[[fact]]])
      }
    }
    metaDF <- metaDF[includedRownames, , drop = FALSE]
    
    if(input$batchFactor != "None") {
      if(length(unique(metaDF[[input$batchFactor]])) == 1) return(F)
      checkList[[input$batchFactor]] <- unique(metaDF[[input$batchFactor]])
    }
    CheckMultiFactorLevels(df = metaDF, factLevelsList = checkList, minCount = 2)
  })
  
  # BULK-DGE-OVER - reactive expression for storing dge analysis
  dgeout1 <- eventReactive(input$godge, {
    req(input$batchFactor)
    cts <- ddsout()[[3]]
    coldata <- ddsout()[[2]]
    batchVar <- w1 <- NULL
    if(input$batchFactor != "None") batchVar <- input$batchFactor
    
    if(!is.null(input$housekeepingSelect)) {
      if(input$housekeepingSelect) w1 <- CalcW1(cts)
    }
    if (input$dgemethod == "limma") {
      if (input$dgeexpsetup == "exp1") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp1b, "_VS_", input$dgeexp1c, sep = "")
          if(input$dgeexp1c == "Rest") {
            coldata[[input$dgeexp1a]] <- as.character(coldata[[input$dgeexp1a]])
            coldata[[input$dgeexp1a]][coldata[[input$dgeexp1a]] != input$dgeexp1b] <- "Rest"
            coldata[[input$dgeexp1a]] <- factor(coldata[[input$dgeexp1a]], levels = c(input$dgeexp1b, "Rest"))
          }
          de.genes <- limma.exp1(
            fact = input$dgeexp1a,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp2") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp2c, "_", input$dgeexp2d, "_VS_", 
                       input$dgeexp2c5, "_", input$dgeexp2d5, sep = "")
          de.genes <- limma.exp2(
            fact1 = input$dgeexp2a,
            fact2 = input$dgeexp2b,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp3") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          de.genes <- limma.exp3(
            fact1 = input$dgeexp3a,
            fact2 = input$dgeexp3b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp3c,
            fact2.rlvl = input$dgeexp3d,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp4") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          de.genes <- limma.exp4(
            fact1 = input$dgeexp4a,
            fact2 = input$dgeexp4b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp4c,
            fact2.rlvl = input$dgeexp4d,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } 
    } else if (input$dgemethod == "edger") {
      if (input$dgeexpsetup == "exp1") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp1b, "_VS_", input$dgeexp1c, sep = "")
          if(input$dgeexp1c == "Rest") {
            coldata[[input$dgeexp1a]] <- as.character(coldata[[input$dgeexp1a]])
            coldata[[input$dgeexp1a]][coldata[[input$dgeexp1a]] != input$dgeexp1b] <- "Rest"
            coldata[[input$dgeexp1a]] <- factor(coldata[[input$dgeexp1a]], levels = c(input$dgeexp1b, "Rest"))
          }
          de.genes <- edger.exp1(
            fact = input$dgeexp1a,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            norm = input$dgeexpedgernorm,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp2") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp2c, "_", input$dgeexp2d,
                       "_VS_", input$dgeexp2c5, "_", input$dgeexp2d5, sep = "")
          de.genes <- edger.exp2(
            fact1 = input$dgeexp2a,
            fact2 = input$dgeexp2b,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            norm = input$dgeexpedgernorm,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp3") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp3(
            fact1 = input$dgeexp3a,
            fact2 = input$dgeexp3b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp3c,
            fact2.rlvl = input$dgeexp3d,
            norm = input$dgeexpedgernorm,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp4") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp4(
            fact1 = input$dgeexp4a,
            fact2 = input$dgeexp4b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp4c,
            fact2.rlvl = input$dgeexp4d,
            norm = input$dgeexpedgernorm,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp5") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp5(
            fact = input$dgeexp5a,
            fact.levl = input$dgeexp5b,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp5c,
            norm = input$dgeexpedgernorm,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp6") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp6(
            me.fact = input$dgeexp6c,
            me.levl = input$dgeexp6d,
            gp.fact = input$dgeexp6a,
            gp.levl = input$dgeexp6b,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp6e,
            norm = input$dgeexpedgernorm,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } 
    } else if (input$dgemethod == "deseq") {
      if (input$dgeexpsetup == "exp1") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp1b, "_VS_", input$dgeexp1c, sep = "")
          if(input$dgeexp1c == "Rest") {
            coldata[[input$dgeexp1a]] <- as.character(coldata[[input$dgeexp1a]])
            coldata[[input$dgeexp1a]][coldata[[input$dgeexp1a]] != input$dgeexp1b] <- "Rest"
            coldata[[input$dgeexp1a]] <- factor(coldata[[input$dgeexp1a]], levels = c(input$dgeexp1b, "Rest"))
          }
          de.genes <- deseq.exp1(
            fact = input$dgeexp1a,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp2") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp2c, "_", input$dgeexp2d,
                       "_VS_", input$dgeexp2c5, "_", input$dgeexp2d5, sep = "")
          de.genes <- deseq.exp2(
            fact1 = input$dgeexp2a,
            fact2 = input$dgeexp2b,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp3") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          de.genes <- deseq.exp3(
            fact1 = input$dgeexp3a,
            fact2 = input$dgeexp3b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp3c,
            fact2.rlvl = input$dgeexp3d,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp4") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          de.genes <- deseq.exp4(
            fact1 = input$dgeexp4a,
            fact2 = input$dgeexp4b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp4c,
            fact2.rlvl = input$dgeexp4d,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp5") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp5c, "_VS_", input$dgeexp5d, sep = "")

          de.genes <- deseq.exp5(
            fact = input$dgeexp5a,
            fact.levl = input$dgeexp5b,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp6") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          con <- paste(input$dgeexp6e, "_VS_", input$dgeexp6f, sep = "")
          de.genes <- deseq.exp6(
            me.fact = input$dgeexp6c,
            me.levl = input$dgeexp6d,
            gp.fact = input$dgeexp6a,
            gp.levl = input$dgeexp6b,
            cts = cts,
            coldata = coldata,
            perm.h = con,
            batchVar = batchVar,
            w1 = w1
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      }
    }
    return(list(fit.cont, fit.names))
  })
  
  # ISOLATE DOES THE IDENTICAL THING - DO WE EVEN NEED THIS ADD'L CODE?
  # BULK-DGE-OVER - prevent app from trying to load DGE 
  # results immediately after experimental 
  # setup is changed
  observeEvent(input$dgeexpsetup, {
    d$newExpSetup <- T
  })
  observeEvent(input$fdr, {
    d$newExpSetup <- T
    d$fdr <- as.numeric(input$fdr)
  })
  observeEvent(input$dgepadjcutoff, {
    d$newExpSetup <- T
  })
  observeEvent(input$dgefcmin, {
    d$newExpSetup <- T
  })
  observeEvent(SubmitData$goqc, {
    if(is.null(SubmitData$goqc)) return()
    if(SubmitData$goqc == 0) return()
    d$newExpSetup <- T
  })
  
  # BULK-DGE-OVER - reactive expression for storing dge contrast table
  dgeout2 <- reactive({
    req(ddsout(), dgeout1(), input$dgeexpsetup, d$dgeContrast, input$dgeexpsetup)
    expset <- input$dgeexpsetup
    contTable <- getContTable(
      de.genes = dgeout1()[[1]],
      coef = d$dgeContrast,
      cts = ddsout()[[3]],
      expset = expset,
      design = dgeout1()[[2]],
      fact = input$dgeexp1a,
      fact5 = input$dgeexp5a,
      fact6 = input$dgeexp6c
    )
    if(is.null(contTable)) return()
    return(list(contTable))
  })
  
  # BULK-DGE-OVER - reactive expression for filtering data based on cut-off
  dgeout3 <- reactive({
    req(dgeout2(), input$dgepadjcutoff, input$dgefcmin)
    tmp <- dgeout2()[[1]]
    if(any(!(c("log2FoldChange", "padj") %in% colnames(tmp)))) return()
    tmp <- tmp[!(is.na(tmp$log2FoldChange) | is.na(tmp$padj)), ]
    padj <- as.numeric(input$dgepadjcutoff)
    lfc <- as.numeric(input$dgefcmin)
    tmp <- tmp[abs(tmp$log2FoldChange) >= lfc, ]
    tmp <- tmp[tmp$padj <= padj, ]
    rownames(tmp) <- tmp$id
    return(tmp)
  })
  
  # BULK-DGE-OVER - reactive expression to generate overview data
  dgeover0 <- eventReactive(input$godge, {
    if(d$dge_warning_msg != "") return()
    expset <- input$dgeexpsetup

    perm <- dgeout1()[[2]]
    perm <- colnames(perm)
    cont.ls <- list()
    for(i in perm) {
      cont.ls[[i]] <- contTable <- getContTable(
        de.genes = dgeout1()[[1]],
        coef = i,
        cts = ddsout()[[3]],
        expset = expset,
        design = dgeout1()[[2]],
        fact = input$dgeexp1a,
        fact5 = input$dgeexp5a,
        fact6 = input$dgeexp6c
      )
    }
    return(list(cont.ls))
  })
  
  # BULK-DGE-OVER - reactive expression for sign dge
  dgeover <- eventReactive(input$godge, {
    req(dgeover0(), input$dgepadjcutoff, input$dgefcmin)
    p <- as.numeric(input$dgepadjcutoff)
    lf <- as.numeric(input$dgefcmin)
    
    test <- dgeOverTbl(cont.ls = dgeover0()[[1]],
                       lf = lf,
                       p = p)
    return(list(test))
  })
  
  # BULK-DGE-OVER - dge overview table front-end
  observeEvent(input$godge, {
    d$dge_warning_msg <- ""
    d$gse_warning_msg <- ""
    d$heatmap_warning_msg <- ""
    d$newGSE <- T
    d$newHeat <- T
    d$newExpSetup <- F
    
    output$dgeoverview <- renderUI({
      req(input$godge)
      if (d$dge_warning_msg != "") return()
      DT::dataTableOutput("dgeoverviewtable")
    })
    
    output$dgeoverviewHeader <- renderUI({
      req(input$godge)
      if(d$dge_warning_msg != "") return()
      h3("DGE regulation table", style = "margin-top: 0px;")
    })
    d$dgeContrast <- colnames(dgeout1()[[2]])[1]
    updateTabsetPanel(inputId = "dge", selected = "overview")
  })
  
  # BULK-DGE-OVER - store experimental design info in object d after DGE analysis is run
  observeEvent(input$godge, {
    if(input$dgeexpsetup == "exp1") {
      d$bulkExpFactorInfo <- list(c(input$dgeexp1b, input$dgeexp1c))
      names(d$bulkExpFactorInfo) <- input$dgeexp1a
    } else if(input$dgeexpsetup == "exp2") {
      d$bulkExpFactorInfo <- list(c(input$dgeexp2c, input$dgeexp2c5),
                                  c(input$dgeexp2d, input$dgeexp2d5))
      names(d$bulkExpFactorInfo) <- c(input$dgeexp2a, input$dgeexp2b)
    } else if(input$dgeexpsetup == "exp3") {
      d$bulkExpFactorInfo <- NULL
    } else if(input$dgeexpsetup == "exp4") {
      d$bulkExpFactorInfo <- NULL
    } else if(input$dgeexpsetup == "exp5") {
      d$bulkExpFactorInfo <- list(c(input$dgeexp5c, input$dgeexp5d))
      names(d$bulkExpFactorInfo) <- input$dgeexp5a
    } else if(input$dgeexpsetup == "exp6") {
      d$bulkExpFactorInfo <- list(input$dgeexp6b,
                                  c(input$dgeexp6e, input$dgeexp6f))
      names(d$bulkExpFactorInfo) <- c(input$dgeexp6a, input$dgeexp6c)
    }
  })
  
  # BULK-DGE-OVER - dge overview table backend
  output$dgeoverviewtable <- DT::renderDataTable(server = FALSE, {
    req(dgeover())
    DT::datatable({
      cont <- dgeover()[[1]]
      colnames(cont) <- c("Comparison", "Regulation", "IDs")
      group1 = gsub("^(.*)_VS_.*", replacement = "\\1", x = cont$Comparison)
      cont$Regulation = paste0(cont$Regulation, " in ", group1)
      cont
    }, rownames = F,
    extensions = 'Buttons',
    options = list(
      autoWidth = FALSE,
      columnDefs = list(list(className = "dt-left", width = "100px", targets = "_all")),
      dom = 'Bt',
      buttons = list(list(extend = 'csv', 
                          filename = "Bulk RNA-Seq DGE analysis",
                          text = 'Download DGE analysis'))
    ), class = "display")
  })
  
  # BULK-DGE-OVER - download button (dge overview table)
  output$dldgeoverviewtbl <- renderUI({
    req(input$godge, req(dgeover()[[1]]))
    if (d$dge_warning_msg != "") return()
    downloadButton("dldgeoverviewtbl2", "Download table")
  })
  
  # BULK-DGE-OVER - Download dge file
  output$dldgeoverviewtbl2 <- downloadHandler(
    filename = function() {
      paste0("dge-overview.csv")
    },
    content = function(file) {
      cont <- dgeover()[[1]]
      colnames(cont) <- c("Comparison", "Regulation", "IDs")
      group1 = gsub("^(.*)_VS_.*", replacement = "\\1", x = cont$Comparison)
      cont$Regulation = paste0(cont$Regulation, " in ", group1)
      write.csv(
        cont,
        file,
        row.names = FALSE
      )
    }
  )
  
  # BULK-DGE-OVER - dge overview plot
  observeEvent(input$godge, {
    output$dgeplot2 <- renderPlotly({
      if (d$dge_warning_msg != "") return()
      input$godge
      comp <- dgeover()[[1]]
      
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
      tickAngle <- 0
      bottomMargin <- 10
      maxChars <- max(nchar(unique(as.character(comp$contrast))))
      if(maxChars > 20 & length(unique(comp$contrast)) > 1) {
        tickAngle <- -45
        bottomMargin <- 150
      }
      plot_ly(
        comp,
        type = "bar",
        y = ~ value,
        x = ~ contrast,
        colors = hue_pal()(2),
        color = ~ variable,
        hoverlabel = label
      ) %>%
        layout(
          margin = list(b = bottomMargin),
          xaxis = list(title = "", tickangle = tickAngle),
          yaxis = list(title = "Number of genes"),
          font = list(family = "Noto Sans JP")
        )
    }) 
    
    # BULK-DGE-OVER - header for DGE regulation barplot
    output$dgeplot2Header <- renderUI({
      req(input$godge)
      if(d$dge_warning_msg != "") return()
      h3("DGE regulation barplot")
    })
    
  })
  
  # BULK-DGE-OVER - download button - dge overview (PDF)
  output$dldgeoverpdf <- renderUI({
    if(d$dge_warning_msg != "") return()
    req(input$godge) 
    downloadButton("dldgeoverpdfimg", "Download static plot (PDF)")
  })
  
  # BULK-DGE-OVER - download file - dge overview (PDF)
  output$dldgeoverpdfimg <- downloadHandler(
    filename =  function() {
      paste("dge-overview.pdf")
    },
    content = function(file) {
      pdf(file,
          width = 8,
          height = 6.5,
          onefile = FALSE) # open the pdf device
      dgeOverPlot(comp = dgeover()[[1]])
      dev.off()
    }
  )
  
  # BULK-DGE-OVER - Show download button - DGE Overview (PNG)
  output$dldgeoverpng <- renderUI({
    if(d$dge_warning_msg != "") return()
    req(input$godge)
    downloadButton("dldgeoverpngimg", "Download static plot (PNG)")
  })
  
  # BULK-DGE-OVER - download plot - dge overview (PNG)
  output$dldgeoverpngimg <- downloadHandler(
    filename =  function() {
      paste("dge-overview.png")
    },
    content = function(file) {
      png(file, width = 800, height = 650)
      dgeOverPlot(comp = dgeover()[[1]])
      dev.off()
    }
  )
  
  ######### VOLCANO PLOTS TAB CODE: BULK RNA-SEQ#####################
  
  # BULK-DGE-VOL - select plot type: ma or volcano
  output$vistype <- renderUI({
    radioGroupButtons(
      inputId = "vistype",
      label = div(style = "margin-top: 0px;", "Plot type"),
      choices = c("Volcano plot" = "volplot", "MA plot" = "maplot"),
      selected = "volplot"
    )
  })
  
  # WHAT IS THIS DOING? PLEASE VERIFY
  # BULK-DGE-VOL - Shared data
  share <- reactive({
    SharedData$new(dgeout3())
  })
  
  myCounter2 <- reactiveVal(0)
  # BULK-DGE-VOL - volcano plots, with warning if no dge
  output$volplots <- renderPlotly({
    req(input$godge, dgeout1(), !d$newExpSetup, d$dgeContrast, input$vistype)
    if(!(d$dgeContrast %in% colnames(dgeout1()[[2]]))) return()

    visType <- input$vistype
    check <- dgeout3()
    check <- nrow(check)
    validate(need(
      expr = check > 1,
      message = paste(
        "Note:",
        "It seems that you have no differentially",
        "expressed IDs. There are zero differentially",
        "expressed points to plot. You may still download",
        "static images of the entire dataset as either",
        "PDF or PNG files."
      )
    ))
    s <- input$mytable_rows_selected
    
    group1 <- strsplit(d$dgeContrast, split = "_VS_")[[1]][1]

    datobj <- dgeout3() %>%
      mutate(fold = ifelse(log2FoldChange > 0, paste0("up in ", group1),
                           ifelse(log2FoldChange < 0, paste0("down in ", group1), NA)))
    tooltips <- paste0(
      "<b>Gene:</b> ",
      datobj$id,
      "<br />",
      "<b>LFC:</b> ",
      round(datobj$log2FoldChange, 3),
      "<br />",
      "<b>PADJ:</b> ",
      round(datobj$padj, 3),
      "<br />",
      "<b>BM:</b> ",
      round(datobj$baseMean, 3),
      "<br />",
      "<b>", d$dgeContrast, ":</b> ",
      datobj$fold,
      "<br />"
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
    
    if (visType == "maplot") {
      if (!length(s)) {
        data <- dgeout3() %>%
          mutate(color = ifelse(log2FoldChange < 0, "#00FF00",
                                ifelse(log2FoldChange > 0, "#FF0000", NA)))
        myColors <- NULL
        if(sum(data$log2FoldChange < 0) > 0) myColors = c(myColors, "#00FF00")
        if(sum(data$log2FoldChange > 0) > 0) myColors = c(myColors, "#FF0000")
        
        isolate({
          p <- data %>%
            plot_ly(
              x = ~ log10(baseMean),
              y = ~ log2FoldChange,
              type = "scatter",
              mode = "markers",
              color = ~ color,
              colors = myColors,
              text = tooltips,
              hoverlabel = label,
              hoverinfo = "text",
              name = "Unfiltered"
            ) %>%
            layout(
              showlegend = FALSE,
              xaxis = list(title = "log<sub>10</sub>(baseMean)", tickfont = list(size = 14), titlefont = list(size = 16)),
              yaxis = list(title = "log<sub>2</sub>(fold change)", tickfont = list(size = 14), titlefont = list(size = 16)),
              font = list(family = "Noto Sans JP")
            ) %>%
            highlight(
              "plotly_selected",
              color = I("royalblue1"),
              selected = attrs_selected(name = "Filtered")
            ) %>% 
            add_annotations(
              x = max(log10(data$baseMean)),
              xanchor = "right",
              y = max(data$log2FoldChange)*1.01,
              yanchor = "bottom",
              text = paste0("Up in ", group1),
              showarrow = F,
              font = list(color = "#FF0000", size = 18)
            ) %>% 
            add_annotations(
              x = max(log10(data$baseMean)),
              xanchor = "right",
              y = min(data$log2FoldChange)*1.01,
              yanchor = "top",
              text = paste0("Down in ", group1),
              showarrow = F,
              font = list(color = "#00FF00", size = 18)
            )
        })
      } else if (length(s)) {
        data <- dgeout3() %>%
          mutate(color = ifelse(log2FoldChange < 0, "#00FF00",
                                ifelse(log2FoldChange > 0, "#FF0000", NA)))
        myColors <- NULL
        if(sum(data$log2FoldChange < 0) > 0) myColors = c(myColors, "#00FF00")
        if(sum(data$log2FoldChange > 0) > 0) myColors = c(myColors, "#FF0000")
        
        isolate({
          pp <- data %>%
            plot_ly() %>%
            add_trace(
              x = ~ log10(baseMean),
              y = ~ log2FoldChange,
              type = "scatter",
              mode = "markers",
              color = ~ color,
              colors = myColors,
              text = tooltips,
              hoverlabel = label,
              hoverinfo = "text",
              name = "Unfiltered"
            ) %>%
            layout(
              showlegend = FALSE,
              xaxis = list(title = "log<sub>10</sub>(baseMean)", tickfont = list(size = 14), titlefont = list(size = 16)),
              yaxis = list(title = "log<sub>2</sub>(fold change)", tickfont = list(size = 14), titlefont = list(size = 16)),
              font = list(family = "Noto Sans JP")
            ) %>% 
            add_annotations(
              x = max(log10(data$baseMean)),
              xanchor = "right",
              y = max(data$log2FoldChange)*1.01,
              yanchor = "bottom",
              text = paste0("Up in ", group1),
              showarrow = F,
              font = list(color = "#FF0000", family = "Noto Sans JP", size = 18)
            ) %>% 
            add_annotations(
              x = max(log10(data$baseMean)),
              xanchor = "right",
              y = min(data$log2FoldChange)*1.01,
              yanchor = "top",
              text = paste0("Down in ", group1),
              showarrow = F,
              font = list(color = "#00FF00", family = "Noto Sans JP", size = 18)
            )
          # selected data
          pp <- add_trace(
            pp,
            data = data[s, , drop = FALSE],
            x = ~ log10(baseMean),
            y = ~ log2FoldChange,
            type = "scatter",
            mode = "markers",
            size = I(8),
            name = "Filtered",
            marker = list(
              color = "yellow",
              size = 8,
              opacity = 1,
              line = list(
                color = "black",
                width = 2
              )
            )
          )
          return(pp)
        })
      }
    } else if (visType == "volplot") {
      if (!length(s)) {
        data <- dgeout3() %>%
          mutate(color = ifelse(log2FoldChange < 0, "#00FF00",
                                ifelse(log2FoldChange > 0, "#FF0000", NA)))
        myColors <- NULL
        if(sum(data$log2FoldChange < 0) > 0) myColors = c(myColors, "#00FF00")
        if(sum(data$log2FoldChange > 0) > 0) myColors = c(myColors, "#FF0000")

        isolate({
          p <- data %>%
            plot_ly(
              x = ~ log2FoldChange,
              y = ~ -log10(pvalue),
              type = "scatter",
              mode = "markers",
              color = ~ color,
              colors = myColors,
              text = tooltips,
              hoverlabel = label,
              hoverinfo = "text",
              name = "Unfiltered"
            ) %>%
            layout(
              showlegend = FALSE,
              xaxis = list(title = "log<sub>2</sub>(fold change)", tickfont = list(size = 14), titlefont = list(size = 16)),
              yaxis = list(title = "-log<sub>10</sub>(p value)", tickfont = list(size = 14), titlefont = list(size = 16)),
              font = list(family = "Noto Sans JP")
            ) %>%
            highlight(
              "plotly_selected",
              color = I("royalblue1"),
              selected = attrs_selected(name = "Filtered")
            ) %>%
            add_annotations(
              x = 0.25,
              xanchor = "left",
              y = max(-log10(data$pvalue[data$pvalue > 0]))*1.01,
              yanchor = "bottom",
              text = paste0("Up in ", group1),
              showarrow = F,
              font = list(color = "#FF0000", size = 18)
            ) %>%
            add_annotations(
              x = -0.25,
              xanchor = "right",
              y = max(-log10(data$pvalue[data$pvalue > 0]))*1.01,
              yanchor = "bottom",
              text = paste0("Down in ", group1),
              showarrow = F,
              font = list(color = "#00FF00", size = 18)
            )

        })
      } else if (length(s)) {
        data <- dgeout3() %>%
          mutate(color = ifelse(log2FoldChange < 0, "#00FF00",
                                ifelse(log2FoldChange > 0, "#FF0000", NA)))
        myColors <- NULL
        if(sum(data$log2FoldChange < 0) > 0) myColors = c(myColors, "#00FF00")
        if(sum(data$log2FoldChange > 0) > 0) myColors = c(myColors, "#FF0000")

        isolate({
          pp <- data %>%
            plot_ly() %>%
            add_trace(
              x = ~ log2FoldChange,
              y = ~ -log10(pvalue),
              type = "scatter",
              mode = "markers",
              color = ~ color,
              colors = myColors,
              text = tooltips,
              hoverlabel = label,
              hoverinfo = "text",
              name = "Unfiltered"
            ) %>%
            layout(
              showlegend = FALSE,
              xaxis = list(title = "log<sub>2</sub>(fold change)", tickfont = list(size = 14), titlefont = list(size = 16)),
              yaxis = list(title = "-log<sub>10</sub>(p value)", tickfont = list(size = 14), titlefont = list(size = 16)),
              font = list(family = "Noto Sans JP")
            ) %>%
            add_annotations(
              x = 0.25,
              xanchor = "left",
              y = max(-log10(data$pvalue))*1.01,
              yanchor = "bottom",
              text = paste0("Up in ", group1),
              showarrow = F,
              font = list(color = "#FF0000", family = "Noto Sans JP", size = 18)
            ) %>%
            add_annotations(
              x = -0.25,
              xanchor = "right",
              y = max(-log10(data$pvalue))*1.01,
              yanchor = "bottom",
              text = paste0("Down in ", group1),
              showarrow = F,
              font = list(color = "#00FF00", family = "Noto Sans JP", size = 18)
            )
          
          # selected data
          pp <- add_trace(
            pp,
            data = data[s, , drop = FALSE],
            x = ~ log2FoldChange,
            y = ~ -log10(pvalue),
            type = "scatter",
            mode = "markers",
            size = I(8),
            name = "Filtered",
            marker = list(
              color = "yellow",
              size = 8,
              opacity = 1,
              line = list(
                color = "black",
                width = 2
              )
            )
          )
        })
      }
    }
  })
  
  # BULK-DGE-VOL - download button - ma or volano plot (PDF)
  output$dldgemavolpdf <- renderUI({
    req(input$godge, !d$newExpSetup, d$dgeContrast)
    if(!(d$dgeContrast %in% colnames(dgeout1()[[2]]))) return()
    
    validate(
      need(input$vistype != "", "")
    )
    if (input$vistype == "maplot") {
      # DGE - Show download button - MA plot (PDF)
      downloadButton("dldgemapdfimg", "Download MA plot (PDF)")
    } else if (input$vistype == "volplot") {
      # DGE - Show download button - Volcano plot (PDF)
      downloadButton("dldgevolpdfimg", "Download volcano plot (PDF)")
    }
  })
  
  # BULK-DGE-VOL - download button - ma or volano plot (PNG)
  output$dldgemavolpng <- renderUI({
    req(input$godge, !d$newExpSetup, d$dgeContrast)
    if(!(d$dgeContrast %in% colnames(dgeout1()[[2]]))) return()
    
    validate(
      need(input$vistype != "", "")
    )
    if (input$vistype == "maplot") {
      # DGE - Show download button - MA plot (PNG)
      downloadButton("dldgemapngimg", "Download MA plot (PNG)")
    } else if (input$vistype == "volplot") {
      # DGE - Show download button - Volcano plot (PNG)
      downloadButton("dldgevolpngimg", "Download volcano plot (PNG)")
    }
  })
  
  # BULK-DGE-VOL - download file - ma plot (PDF)
  output$dldgemapdfimg <- downloadHandler(
    filename =  function() {
      paste("dge-ma-plot.pdf")
    },
    content = function(file) {
      pdf(file, width = 8, height = 6.5, onefile = FALSE)
      dgeMAPlot(
        dgeout2 = dgeout2()[[1]],
        p = as.numeric(input$dgepadjcutoff),
        l = as.numeric(input$dgefcmin),
        cont = d$dgeContrast
      )
      dev.off()
    }
  )
  
  # BULK-DGE-VOL - download file - ma plot (PNG)
  output$dldgemapngimg <- downloadHandler(
    filename =  function() {
      paste("dge-ma-plot.png")
    },
    content = function(file) {
      png(file, width = 800, height = 650)
      dgeMAPlot(
        dgeout2 = dgeout2()[[1]],
        p = as.numeric(input$dgepadjcutoff),
        l = as.numeric(input$dgefcmin),
        cont = d$dgeContrast
      )
      dev.off()
    }
  )
  
  # BULK-DGE-VOL - download file - volcano plot (PDF)
  output$dldgevolpdfimg <- downloadHandler(
    filename =  function() {
      paste("dge-vol-plot.pdf")
    },
    content = function(file) {
      pdf(file, width = 8, height = 6.5, onefile = FALSE)
      dgeVolPlot(
        dgeout2 = dgeout2()[[1]],
        p = as.numeric(input$dgepadjcutoff),
        l = as.numeric(input$dgefcmin),
        cont = d$dgeContrast
      )
      dev.off()
    }
  )
  
  # BULK-DGE-VOL - download file - volcano plot (PNG)
  output$dldgevolpngimg <- downloadHandler(
    filename =  function() {
      paste("dge-vol-plot.png")
    },
    content = function(file) {
      png(file, width = 850, height = 600, res = 125)
      dgeVolPlot(
        dgeout2 = dgeout2()[[1]],
        p = as.numeric(input$dgepadjcutoff),
        l = as.numeric(input$dgefcmin),
        cont = d$dgeContrast
      )
      dev.off()
    }
  )
  
  # BULK-DGE-VOL - dge table header equation
  output$mytableHeader <- renderUI({
    req(input$godge, !d$newExpSetup, input$dgepadjcutoff, input$dgefcmin)
    txt <- paste('abs(\\(\\log_{2}\\)foldchange) \\(\\gt\\)',
                 input$dgefcmin,
                 '& adjusted p-value \\(\\leq\\)',
                 input$dgepadjcutoff)
    
    tagList(
      h3("DGE table"),
      withMathJax(p(txt))
    )
  })
  
  # BULK-DGE-VOL - filtered dge table
  output$mytable <- DT::renderDataTable(server = FALSE, {
    req(input$godge, !d$newExpSetup, d$dgeContrast)
    if(!(d$dgeContrast %in% colnames(dgeout1()[[2]]))) return()
    isolate({
      tmp <- dgeout3() %>% mutate_if(is.numeric, round, digits = 4)
      tmp <- tmp[complete.cases(tmp), ]
      m2 <- tmp[share()$selection(),]
      dt <- DT::datatable({
        tmp
      }, rownames = F,
      extensions = 'Buttons',
      options = list(
        autoWidth = FALSE,
        columnDefs = list(list(className = "dt-left", width = "100px", targets = "_all")),
        dom = 'Bfrtip',
        buttons = list(list(extend = 'csv', 
                            filename = "Bulk RNA-Seq DGE filtered analysis",
                            text = 'Download DGE filtered analysis')),
        language = list(zeroRecords = "No genes passing thresholds") 
      ), class = "display")
      if (NROW(m2) == 0) {
        dt
      } else {
        DT::formatStyle(
          dt,
          "id",
          target = "row",
          color = DT::styleEqual(m2$id, rep("white", length(m2$id))),
          backgroundColor = DT::styleEqual(m2$id,
                                           rep("darkgray", length(m2$id)))
        )
      }
    })
  })
  
  # BULK-DGE-VOL - download button for all data (not-filtered)
  output$downloadall <- renderUI({
    if(is.null(input$godge) || input$godge == 0 || d$newExpSetup) return()
    downloadButton("downallData", label = "Download all data")
  })
  
  # BULK-DGE-VOL - download all data file (not filtered)
  output$downallData <- downloadHandler(
    filename = function() {
      paste0(d$dgeContrast
             , "_all_results.csv")
    },
    content = function(file) {
      write.csv(dgeout2()[[1]],
                file, row.names = FALSE)
    }
  )
  
  ######### GSE & HEATMAP TAB CODE: BULK RNA-SEQ#####################
  
  # BULK-DGE-GSE - filtered dge dataset
  mytab <- reactive({
    req(input$godge, d$dgeContrast, dgeout3(), !d$newExpSetup)
    tmp <- dgeout3() %>% mutate_if(is.numeric, round, digits = 4)
    tmp <- tmp[complete.cases(tmp), ]
    m2 <- tmp[share()$selection(),]
  })
  
  # BULK-DGE-GSE - load list of enrichr libs as default
  output$enrichRLib <- renderUI({
    enrichRLibs <- scan("./data/enrichrlibrary.txt", what = "character", quiet = T)
    selectInput(
      inputId = "enrichRLib", 
      label = "EnrichR libraries", 
      choices = enrichRLibs, 
      selected = enrichRLibs[1:6], 
      multiple = T, 
      width = "600px"
    )
  })
  observe({
    if(is.null(input$enrichRLib)) {
      d$enrichRLib <- scan("./data/enrichrlibrary.txt", what = "character", quiet = T)[1:6]
    } else {
      d$enrichRLib <- input$enrichRLib
    }
  })
  
  # BULK-DGE-GSE - choose contrasts for gse
  output$dgemaincontrasts_gse <- renderUI({
    req(input$godge, input$dge_list, dgeout1())
    if(input$dge_list == "Custom") return()
    tmp <- dgeout1()[[2]]
    tmp <- colnames(tmp)
    selected = NULL
    if(!is.null(d$dgeContrast)) {
      if(d$dgeContrast %in% tmp) selected = d$dgeContrast
    }
    div(
      style = "margin-left: 20px;",
      selectInput(
        inputId = "dgemaincontrasts_gse",
        label = "Contrast",
        choices = tmp,
        selected = selected,
        width = "300px"
      )        
    )
  })
  
  # BULK-DGE-GSE - store heatmap contrasts in object d 
  observeEvent(input$dgemaincontrasts_heatmap, {
    if(is.null(d$dgeContrast) || input$dgemaincontrasts_heatmap != d$dgeContrast) {
      d$dgeContrast <- input$dgemaincontrasts_heatmap
    }
  })
  
  # BULK-DGE-GSE - select type of gene list to use
  output$dge_list <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type == "Bulk") {
      if(is.null(input$godge) || input$godge == 0) return()
      myChoices <- c("DGE filtered", "Custom")
      mySelected <- "DGE filtered"
    } else {
      myChoices <- c("All", "Custom")
      mySelected <- "All"
    }

    div(
      style = "width: 110px;",
      awesomeRadio(
        inputId = "dge_list", 
        label = "Input gene list",
        choices = myChoices,
        selected = mySelected,
        status = "success"
      )
    )
  })   
  
  # BULK-DGE-GSE - copy/paste box of genes for targeted geneset enrich
  output$gene_list <- renderUI({
    if(!is.null(input$dge_list) && input$dge_list != "Custom") return()
    myStyle <- "margin-left: 30px;"
    if(is.null(input$dge_list)) myStyle <- NULL
    
    div(
      style = myStyle,
      textAreaInput(
        inputId = "gene_list", 
        label = "Paste genes", 
        value = "", 
        width = "150px"
      )
    )
  })
  
  # BULK-DGE-GSE - select type of geneset enrichment (which genes to include)
  output$gene_list_filtering <- renderUI({
    if(is.null(d$dgeContrast) || is.null(input$godge) || input$godge == 0) return()
    if(is.null(input$dge_list) || input$dge_list != "DGE filtered") return()

    if(grepl(".*_VS_.*", d$dgeContrast)) {
      group1 <- gsub("^(.*)_VS_.*", replacement = "\\1", x = d$dgeContrast)
      group2 <- gsub(".*_VS_(.*)$", replacement = "\\1", x = d$dgeContrast)
      tempChoices <- c(1:3)
      names(tempChoices) <- c("All DE genes", 
                              paste0("Up in ", group1),
                              paste0("Up in ", group2))
    } else {
      tempChoices <- c(1:3)
      names(tempChoices) <- c("All DE genes", 
                              paste0("Up in ", input$dgemaincontrasts_gse),
                              paste0("Down in ", input$dgemaincontrasts_gse))
    }
    div(
      style = "margin-left: 20px;",
      awesomeRadio(
        inputId = "dge_filter", 
        label = "Gene list filtering", 
        choices = tempChoices,
        selected = 1, 
        status = "success"
      )
    )
  })
  
  # THESE ARE ALL REDUNDANT - WRAP INTO A SINGLE OBSERVEEVENT
  # BULK-DGE-GSE - store updated filtered data size in object d
  observe({
    if(is.null(gsetran1()) || is.null(gsetran1()[[1]])) {
      if(!is.null(input$dge_list) && input$dge_list == "DGE filtered") {
        d$dge_info <- "No DE genes."
      } else if(!is.null(input$dge_list) && input$dge_list == "Custom") {
        d$dge_info <- "No genes found in data."
      }
    } else {
      dims = dim(gsetran1()[[1]])
      d$dge_info <- paste0("Data size: ", dims[1], " genes")
    }
  })
  
  # THESE ARE ALL REDUNDANT - WRAP INTO A SINGLE OBSERVEEVENT
  # BULK-DGE-Heatmap - store updated filtered data size in object d
  observe({
    if(is.null(heattran1()) || is.null(heattran1()[[1]])) {
      if(!is.null(input$dge_list_heatmap) && input$dge_list_heatmap == "DGE filtered") {
        d$dge_info_heatmap <- "No DE genes."
      } else if(!is.null(input$dge_list_heatmap) && input$dge_list_heatmap == "Custom") {
        d$dge_info_heatmap <- "No genes found in data."
      } else if(!is.null(input$gene_list_heatmap) && input$gene_list_heatmap != "") {
        d$dge_info_heatmap <- "No genes found in data."
      } else {
        d$dge_info_heatmap <- NULL
      }
    } else {
      dims <- dim(heattran1()[[1]])
      geneTxt <- " genes, "
      if(dims[1] == 1) geneTxt <- " gene, "
      
      d$dge_info_heatmap <- paste0("Data size: ", dims[1], geneTxt, dims[2], " samples")
    }
  })
  
  # BULK-DGE-GSE - select top number of genes to include
  output$top_n_genes <- renderUI({
    req(input$godge, input$dge_list)
    if(input$dge_list == "Custom") return()
    div(
      style = "width: 150px;",
      textInput(
        inputId = "n_genes", 
        label = "Top n genes", 
        value = 100, 
        placeholder = "ALL"
      )
    )
  })
  
  # BULK-DGE-GSE - select all genes to use in geneset enrichment
  output$use_all_genes <- renderUI({
    req(input$godge, input$dge_list)
    if(input$dge_list == "Custom") return()
    div(
      prettyCheckbox(
        inputId = "fea_all_genes", 
        label = "All genes passing filters", 
        value = F, 
        status = "default", 
        icon = icon("check")
      )
    )
  })
  
  output$geneOptions_ui <- renderUI({
    if(is.null(input$godge) || input$godge == 0 || is.null(input$dge_list) || 
       input$dge_list == "Custom") return()
    div(
      style = "width: 180px; margin-left: 20px;", 
      fluidPage(
        fluidRow(uiOutput("top_n_genes")),
        fluidRow(uiOutput("use_all_genes"))
      )
    )
  })
  
  # BULK-DGE-GSE - set reactive values of all genes (if it changes) in object d
  observeEvent(input$fea_all_genes, {
    d$fea_all_genes <- input$fea_all_genes
  })
  
  # BULK-DGE-GSE - set reactive values of all genes (if it changes) in object d
  observeEvent(input$fea_all_genes_heatmap, {
    d$fea_all_genes_heatmap <- input$fea_all_genes_heatmap
  })
  
  # BULK-DGE-GSE - print info about DGE results (e.g., number of genes and
  # samples) before user starts gse processing
  output$dge_info <- renderUI({
    if(!is.null(input$dge_list) && input$dge_list == "Custom" &&
       (is.null(input$gene_list) || input$gene_list == "")) return()
    h4(d$dge_info, style = "margin-top: 0px; margin-right: 30px;")
  })
  
  # BULK-DGE-GSE - action button to run gse
  output$submit_fe <- renderUI({
    if(is.null(gsetran1()) || is.null(gsetran1()[[1]]) || nrow(gsetran1()[[1]]) == 0) return()
    actionButton("submit_fe", "Run GSE", icon = icon("space-shuttle"))
  })
  
  # Show/hide submit_fe button depending on whether data are ready
  observe({
    if(is.null(gsetran1())) {
      shinyjs::hide(id = "submit_fe")
    } else {
      shinyjs::show(id = "submit_fe")
    }
  })
  
  # BULK-DGE-GSE - set warning as empty if the gse data is empty
  observeEvent(input$submit_fe, {
    if(is.null(d$enrichRLib) || length(d$enrichRLib) == 0 ||
       d$enrichRLib[1] == "") {
      d$gse_warning_msg <- "Please select at least one EnrichR library."
    } else if(d$n_genes != "" & suppressWarnings(is.na(as.integer(na.omit(d$n_genes))))) {
      d$gse_warning_msg <- "Please provide an integer in 'Use top n genes'."
    } else if(!is.null(input$dge_list) && input$dge_list == "Custom") {
      if(length(unlist(strsplit(input$gene_list, split = "\\s+|,\\s?"))) == 0) {
        d$gse_warning_msg <- "No genes were submitted."
      }
    } else if(is.null(gsetran1())) {
      d$gse_warning_msg <- "No DE genes."
    } else {
      d$gse_warning_msg <- ""
      d$newGSE <- F
    }
    
  })
  
  # BULK-DGE-GSE - reactive expression for warnings for gse
  func_enrich <- eventReactive(input$submit_fe, {
      withProgress(message = "Computing functional enrichment gene lists...", value = 0, {
        incProgress(0.3)
        isolate({
          if(!is.null(input$dge_list) && input$dge_list == "DGE filtered") {
            dge <- dgeout3()
            if("rowname" %in% colnames(dge)) colnames(dge)[colnames(dge) == "rowname"] = "id"
            uni_genes <- unique(dge$id)
            if(input$dge_filter == 2) dge <- dge[dge$log2FoldChange > 0, ]
            if(input$dge_filter == 3) dge <- dge[dge$log2FoldChange < 0, ]
            dge <- dge[order(1/abs(dge$log2FoldChange), dge$pvalue, decreasing = F), ]
            # load library list
            libs <- d$enrichRLib
            if(d$n_genes != "" & !(input$fea_all_genes)) dge <- dge[1:min(as.integer(d$n_genes), nrow(dge)), ]
            # convert genes to a vector
            genes <- dge %>% pull(id)
            d$genes_debug <- genes
            # run enrichR w/selected libraries
            enr <- enrichr(genes, libs)
            enr <- enr[sapply(enr, function(x) nrow(x) > 0)]
            if (length(enr) == 0) {
              return(NULL)
            } else {
              enrCheck <- sapply(enr, function(i) {
                if(ncol(i) == 1) {
                  if(colnames(i)[1] == "X.html.") return(FALSE)
                }
                return(TRUE)
              })
              if(all(!enrCheck)) {
                d$gse_warning_msg <- "EnrichR is unavailable. Please try again later."
                return()
              }
              
              # return significant genes
              enr <- getSigTerms(enr, libs)
              enr <- enr[sapply(enr, function(x) !is.null(x))]
              if (length(enr) == 0) {
                return(NULL)
              } else {
                nam <- unlist(lapply(1:length(enr), function(x)
                  unique(enr[[x]]$libName)))
                names(enr) <- nam
                incProgress(0.6)
                # function for removing redundant go terms
                enr <- rbindlist(enr)
                enr <- enr[order(-score)]
                enr.go <- enr[grepl(pattern = "GO", x = libName) & grepl(pattern = "\\(GO", x = term), , drop = F]
                if (nrow(enr.go) > 0) {
                  enr.go[,GOID := tstrsplit(term, "\\(GO")[2]]
                  enr.go[,GOID := paste0('GO', GOID)]
                  enr.go[,GOID := gsub("\\)", '', GOID)]
                  enr.go[,term.short := tstrsplit(term, "\\(")[1]]
                  enr.go[,term.short := trimws(term.short, which="right")]
                  enr.go <- enr.go %>%
                    mutate(term = paste(term.short, " (", GOID, ")", sep = "")) %>%
                    dplyr::select(-c(term.short, GOID))
                  enr.rest <- enr[!grepl('GO', libName),]
                  enr <- rbind(enr.rest, enr.go)
                }
                enr <- enr[order(-score)]
                enr <- enr %>%
                  mutate(pval = round(pval, digits = 3),
                         adjPval = round(adjPval, digits = 3),
                         Z_score = round(Z_score, digits = 3),
                         score = round(score, digits = 3))
                names(enr) <- c("Library name", "Library rank", "Gene count",
                                "Term", "Overlap", "P-value", "Adjusted p-value",
                                "Old p-value", "Old adjusted p-value", "Z-score",
                                "Score", "Gene list")
                enr <- enr %>% 
                  arrange(`P-value`)
                enr <- enr[, -c("Old p-value", "Old adjusted p-value"), with = FALSE]
                return(enr)
              }
            }
          } else if(is.null(input$dge_list) || input$dge_list == "Custom") {
            genes <- unlist(strsplit(input$gene_list, split = "\\s+|,\\s?"))
            if(length(genes) == 0) return()
            genes <- toupper(genes)
            cts <- ddsout()[[3]]
            inc <- attr(cts, "dimnames")[[1]][attr(cts, "dimnames")[[1]] %in% genes]
            genes <- genes[genes %in% inc]
            if(length(genes) == 0) return()
            # load library list
            libs <- d$enrichRLib
            # run enrichR w/selected libraries
            enr <- enrichr(genes, libs)
            enr <- enr[sapply(enr, function(x) nrow(x) > 0)]
            if(length(enr) == 0) return() 
            enrCheck <- sapply(enr, function(i) {
              if(ncol(i) == 1) {
                if(colnames(i)[1] == "X.html.") return(FALSE)
              }
              return(TRUE)
            })
            if(all(!enrCheck)) {
              d$gse_warning_msg <- "EnrichR is unavailable. Please try again later."
              return()
            }
            # return significant genes
            enr <- getSigTerms(enr, libs)
            enr <- enr[sapply(enr, function(x) !is.null(x))]
            if(length(enr) == 0) return()
            nam <- unlist(lapply(1:length(enr), function(x) unique(enr[[x]]$libName)))
            names(enr) <- nam
            incProgress(0.6)
            # function for removing redundant go terms
            enr <- rbindlist(enr)
            enr <- enr[order(-score)]
            enr.go <- enr[grepl(pattern = "GO", x = libName) & grepl(pattern = "\\(GO", x = term), , drop = F]
            if(nrow(enr.go) > 0) {
              enr.go[,GOID := tstrsplit(term, "\\(GO")[2]]
              enr.go[,GOID := paste0('GO', GOID)]
              enr.go[,GOID := gsub("\\)", '', GOID)]
              enr.go[,term.short := tstrsplit(term, "\\(")[1]]
              enr.go[,term.short := trimws(term.short, which="right")]
              enr.go <- enr.go %>%
                mutate(term = paste(term.short, " (", GOID, ")", sep = "")) %>%
                dplyr::select(-c(term.short, GOID))
              enr.rest <- enr[!grepl('GO', libName),]
              enr <- rbind(enr.rest, enr.go)
            }
            enr <- enr[order(-score)]
            enr <- enr %>%
              mutate(pval = round(pval, digits = 3),
                     adjPval = round(adjPval, digits = 3),
                     Z_score = round(Z_score, digits = 3),
                     score = round(score, digits = 3))
            names(enr) <- c("Library name", "Library rank", "Gene count",
                            "Term", "Overlap", "P-value", "Adjusted p-value",
                            "Old p-value", "Old adjusted p-value", "Z-score",
                            "Score", "Gene list")
            enr <- enr %>% 
              arrange(`P-value`)
            enr <- enr[, -c("Old p-value", "Old adjusted p-value"), with = FALSE]
            return(enr)
          }
        })
      })
    # }
  })
  
  # BULK-DGE-GSE - header for gse table
  output$tbl_func_enrHeader <- renderUI({
    req(func_enrich(), !d$newGSE, input$submit_fe)
    if(!is.null(d$gse_warning_msg) && d$gse_warning_msg != "") return()
    h3("Gene set enrichment table", style = "margin-top: 0px;")
  })
  
  output$func_enr_tbl <- DT::renderDataTable(server = FALSE, {
    if(is.null(func_enrich())) return()
    DT::datatable({
      func_enrich()
    },
    rownames = FALSE,
    extensions = 'Buttons',
    options = list(
      autoWidth = FALSE,
      columnDefs = list(list(className = 'dt-left', width = '40%', targets = "_all")),
      dom= 'Bfrtip',
      buttons = list(list(extend = 'csv', 
                          filename = "Bulk_GSE",
                          text = 'Download GSE')),
      scrollX = TRUE
    ), class = "display")
  })
  
  # BULK-DGE-GSE - if no gse is found, print messages
  output$no_list <- renderText({
    req(input$godge, input$submit_fe)
    if(d$gse_warning_msg != "" | d$newExpSetup) return(NULL)
    if (input$dge_list == "DGE filtered" &&
        is.null(func_enrich())) {
      "There are no significant enrichments of common annotated biological features."
    } else if (input$dge_list == "Custom" &&
               is.null(func_enrich())) {
      "There are no significant enrichments of common annotated biological features for the submitted genes."
    }
  })
  
  # BULK-DGE-GSE - paste user-selected gene list here
  output$gene_targ <- renderUI({
    req(input$godge, input$dge_list)
    if (input$dge_list == "Custom") {
      textAreaInput("gene_targ", "Paste genes", value = "" , width = "100px")
    }
  })

  # BULK-DGE-GSE - Heatmap help button
  observeEvent(input$heatmap_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Heatmap")),
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
      includeMarkdown("markdown/help/heatmap_help.md")
    ))
  })
  
  # BULK-DGE-Heatmap - select type of gene list to use
  output$dge_list_heatmap <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type == "Bulk") {
      if(is.null(input$godge) || input$godge == 0) {
        myChoices <- "Custom"
        mySelected <- "Custom"
      } else {
        myChoices <- c("DGE filtered", "Custom")
        mySelected <- "DGE filtered"
      }
    } else {
      myChoices <- c("All", "Custom")
      mySelected <- "All"
    }
    div(
      style = "width: 110px; margin-left: 10px;",
      awesomeRadio(
        inputId = "dge_list_heatmap", 
        label = "Input gene list",
        choices = myChoices,
        selected = mySelected,
        status = "success"
      )
    )
  }) 
  
  # BULK-DGE-Heatmap - choose contrasts for heatmap
  output$dgemaincontrasts_heatmap <- renderUI({
    req(input$godge, input$dge_list_heatmap, dgeout1())
    if(input$dge_list_heatmap != "DGE filtered") return()
    tmp <- dgeout1()[[2]]
    tmp <- colnames(tmp)
    selected = NULL
    if(!is.null(d$dgeContrast)) {
      if(d$dgeContrast %in% tmp) selected = d$dgeContrast
    }
    div(
      style = "margin-left: 20px;",
      selectInput(
        inputId = "dgemaincontrasts_heatmap",
        label = "Contrast",
        choices = tmp,
        selected = selected,
        width = "300px"
      )        
    )
  })
  
  # BULK-DGE-Heatmap - copy/paste box of genes for heatmap
  output$gene_list_heatmap <- renderUI({
    if(!is.null(input$dge_list_heatmap) && input$dge_list_heatmap != "Custom") return()
    myStyle <- "margin-left: 30px;"
    if(is.null(input$dge_list_heatmap)) myStyle <- NULL
    div(
      style = myStyle,
      textAreaInput(
        inputId = "gene_list_heatmap", 
        label = "Paste genes", 
        value = "", 
        width = "150px"
      )
    )
  })
  
  # RASL - Choose a factor to use to filter log2 fold-changes (e.g., Treatment, Ratio)
  output$filterFactor_heatmap <- renderUI({
    req(SubmitData$data_type)
    if(SubmitData$data_type != "RASL") return()
    if(is.null(ddsout())) return()
    myChoices <- c("None", colnames(ddsout()[[2]]))
    mySelected <- "None"
    div(
      style = "width: 200px; margin-left: 30px; vertical-align: top;",
      selectInput(
        inputId = "filterFactor_heatmap",
        label = "Filter samples by factor",
        choices = myChoices,
        selected = mySelected
      )
    )
  })
  
  output$filterFactorSelection_heatmap <- renderUI({
    req(SubmitData$data_type, input$filterFactor_heatmap)
    if(SubmitData$data_type != "RASL") return()
    if(is.null(input$filterFactor_heatmap) || input$filterFactor_heatmap == "None") return()
    if(is.null(ddsout())) return()
    meta <- ddsout()[[2]]
    myChoices <- mixedsort(unique(meta[, input$filterFactor_heatmap]))
    mySelected <- myChoices[1]
    myLabel <- paste0(input$filterFactor_heatmap, " level")
    div(
      style = "width: 200px; margin-left: 30px; vertical-align: top;",
      selectInput(
        inputId = "filterFactorSelection_heatmap",
        label = myLabel,
        choices = myChoices,
        selected = mySelected,
        multiple = T
      )
    )
  })
  
  # BULK-DGE-Heatmap - select which genes to include
  output$gene_list_filtering_heatmap <- renderUI({
    if(is.null(d$dgeContrast) || is.null(input$godge) || input$godge == 0) return()
    if(is.null(input$dge_list_heatmap) || input$dge_list_heatmap != "DGE filtered") return()
    if(grepl(".*_VS_.*", d$dgeContrast)) {
      group1 <- gsub("^(.*)_VS_.*", replacement = "\\1", x = d$dgeContrast)
      group2 <- gsub(".*_VS_(.*)$", replacement = "\\1", x = d$dgeContrast)
      tempChoices <- c(1:3)
      names(tempChoices) <- c("All DE genes", 
                              paste0("Up in ", group1),
                              paste0("Up in ", group2))
    } else {
      tempChoices <- c(1:3)
      names(tempChoices) <- c("All DE genes", 
                              paste0("Up in ", input$dgemaincontrasts_heatmap),
                              paste0("Down in ", input$dgemaincontrasts_heatmap))
    }
    div(
      style = "margin-left: 20px;",
      awesomeRadio(
        inputId = "dge_filter_heatmap", 
        label = "Gene list filtering", 
        choices = tempChoices,
        selected = 1, 
        status = "success"
      )
    )
  })
  
  # BULK-DGE-Heatmap - select top number of genes to include
  output$top_n_genes_heatmap <- renderUI({
    req(input$godge, input$dge_list_heatmap)
    if(input$dge_list_heatmap == "Custom") return()
    if(is.null(dgeout3())) return()
    maxGenes <- nrow(dgeout3())
    myValue <- min(maxGenes, 100)
    div(
      style = "width: 150px;",
      numericInput(
        inputId = "n_genes_heatmap", 
        label = "Top n genes", 
        min = 1,
        max = maxGenes,
        value = myValue, 
        step = 1
      )
    )
  })
  
  observeEvent(input$n_genes_heatmap, {
    if(
      !IsInteger(input$n_genes_heatmap) ||
      input$n_genes_heatmap < 1 || 
      input$n_genes_heatmap > nrow(dgeout3())
    ) {
      if(!is.null(dgeout3()) && nrow(dgeout3()) > 1) {
        myValue <- min(nrow(dgeout3()), 100)
      } else {
        myValue <- 100
      }
      updateNumericInput(session = session, inputId = "n_genes_heatmap", value = myValue)
    }
  })
  
  # BULK-DGE-Heatmap - select all genes to use in heatmap
  output$use_all_genes_heatmap <- renderUI({
    req(input$godge, input$dge_list_heatmap)
    if(input$dge_list_heatmap == "Custom") return()
    div(
      prettyCheckbox(
        inputId = "fea_all_genes_heatmap", 
        label = "All genes passing filters", 
        value = F, 
        status = "default", 
        icon = icon("check")
      )
    )
  })
  
  output$geneOptions_heatmap_ui <- renderUI({
    if(is.null(input$godge) || input$godge == 0 || is.null(input$dge_list_heatmap) || 
       input$dge_list_heatmap == "Custom") return()
    div(
      style = "width: 180px; margin-left: 20px;", 
      fluidPage(
        fluidRow(uiOutput("top_n_genes_heatmap")),
        fluidRow(uiOutput("use_all_genes_heatmap"))
      )
    )
  })
  
  # BULK-DGE-Heatmap - set reactive values of all genes (if it changes) in object d
  observeEvent(input$fea_all_genes_heatmap, {
    d$fea_all_genes_heatmap <- input$fea_all_genes_heatmap
  })
  
  # BULK-DGE-Heatmap - print info about DGE results (e.g., number of genes and
  # samples) before user starts heatmap processing
  output$dge_info_heatmap <- renderUI({
    if(!is.null(input$dge_list_heatmap) && input$dge_list_heatmap == "Custom" &&
       (is.null(input$gene_list_heatmap) || input$gene_list_heatmap == "")) return()
    if(is.null(d$dge_info_heatmap) || d$dge_info_heatmap == "") return()
    h4(d$dge_info_heatmap, style = "margin-top: 0px; margin-right: 30px;")
  })
  
  # BULK-DGE-GSE - color palette for heatmap
  output$heatColors <- renderUI({
    selectInput(
      inputId = "heatColors",
      label = "Color palette",
      choices = c(
        'Default',
        'Red-white-blue',
        'Viridis',
        'Green-yellow-red'
      ),
      selected = 'Default'
    )
  })
  
  # BULK-DGE-Heatmap - cluster row selection
  output$row_type <- renderUI({
    prettyCheckbox("row_type", label = "Cluster genes", 
                  value = T, status = "default", icon = icon("check"))
  })
  
  # BULK-DGE-Heatmap - cluster column selection
  output$col_type <- renderUI({
    prettyCheckbox("col_type", label = "Cluster samples",
                  value = T, status = "default", icon = icon("check"))
  })
  
  # SC-DGE-Heatmap - scale genes option
  output$scale_genes <- renderUI({
    prettyCheckbox("scale_genes", label = "Scale genes",
                   value = T, status = "default", icon = icon("check"))
  })
  
  
  # BULK-DGE-GSE - gse warning
  output$gse_warning <- renderText(d$gse_warning_msg)
  
  output$dgeInfoSubmitWarning_ui <- renderUI({
    if((is.null(d$dge_info) || d$dge_info == "") && 
       (is.null(gsetran1()) || is.null(gsetran1()[[1]]) || nrow(gsetran1()[[1]]) == 0) &&
       (is.null(d$gse_warning_msg) || d$gse_warning_msg == "")) {
      return()
    }
    tagList(
      div(
        div(style = "display: inline-block;", uiOutput("dge_info")),
        div(style = "display: inline-block;", uiOutput("submit_fe"))
      ),
      div(textOutput("gse_warning"),style="color:red"),
      hr(),
    )
  })
  
  # BULK-DGE-GSE - heatmap warning
  output$heatmap_warning <- renderText(d$heatmap_warning_msg)
  
  output$dgeInfoSubmitWarning_heatmap_ui <- renderUI({
    tagList(
      div(
        div(style = "display: inline-block;", uiOutput("dge_info_heatmap")),
        div(style = "display: inline-block;", uiOutput("submit_list"))
      ),
      div(textOutput("heatmap_warning"), style = "color: red;"),
      hr()
    )
  })
  
  # BULK-DGE-GSE - action button for building the heatmap
  output$submit_list <- renderUI({
    req(heattran1())
    actionButton(
      inputId = "submit_list", 
      label = "Build heatmap", 
      icon = icon("space-shuttle")
    )
  })
  
  # BULK-DGE-Heatmap - genes for heatmap
  genelist_heatmap <- eventReactive(input$submit_list, { 
    req(input$gene_list_heatmap)
    input$gene_list_heatmap
  })
  
  # BULK-DGE-GSE - genes not found in user-selected list
  observeEvent(input$submit_fe, {
    d$gse_warning_msg <- ""
    d$newGSE <- F
    if(is.null(gsetran1())) {
      d$gse_warning_msg <- "No DE genes."
      return()
    }
    if(is.null(input$dge_list) || input$dge_list == "Custom") {
      if(length(unlist(strsplit(input$gene_list, split = "\\s+|,\\s?"))) == 0 |
         input$gene_list == "") {
        d$gse_warning_msg <- "No genes were submitted."
      }
    }
    output$miss_genesa <- renderUI({
      if(!is.null(input$dge_list) && input$dge_list == "DGE filtered") return()
      if(is.null(input$gene_list) || input$gene_list == "") {
        d$gse_warning_msg <- "No genes were submitted."
        return()
      }
      genes <- unlist(strsplit(input$gene_list, split = "\\s+|,\\s?"))
      genes <- toupper(genes)
      cts <- ddsout()[[3]]
      inc <- attr(cts, "dimnames")[[1]]
      genes <- genes[!genes %in% inc]
      if(length(genes) == 0) {
        d$gse_warning_msg <- ""
        return()
      }
      d$gse_warning_msg <- "One or more genes were not found."
      genes <- paste(genes, collapse = ', ')
      textAreaInput("miss_genes", "The following genes were not found in analysis:", 
                    value = genes, width = "100px")
    })
  })
  
  # BULK-DGE-GSE - header for heatmap: gene selection directions
  output$headheat <- renderUI({
    req(heattran2(), !d$newHeat, input$submit_list)
    h3("Interactive heatmap (click on cells)")
  })
  
  # BULK-DGE-GSE - use only samples in groups from contrast
  output$contrast_groups_only <- renderUI({
    req(input$godge, input$dgeexpsetup, input$dge_list_heatmap)
    if(input$dge_list_heatmap != "DGE filtered") return()
    if(!(input$dgeexpsetup %in% c("exp1","exp2","exp5","exp6"))) return()
    prettyCheckbox(
      inputId = "contrast_groups_only", 
      label = "Use samples from contrast only", 
      value = F, 
      status = "default", 
      icon = icon("check")
    )
  })
  
  
  # BULK-DGE-GSE - select factors to use for sample labels and legend for heatmap
  output$heatfactorlabel <- renderUI({
    req(ddsout(), SubmitData$data_type)
    myChoices <- colnames(ddsout()[[2]])
    mySelected <- myChoices[1]
    if(SubmitData$data_type == "RASL" && !is.null(input$filterFactor_heatmap) &&
       input$filterFactor_heatmap != "None" && !is.null(input$filterFactorSelection_heatmap) &&
       input$filterFactor_heatmap %in% myChoices && length(myChoices) >= 2) {
      mySelected <- myChoices[myChoices != input$filterFactor_heatmap][1]
      myChoices <- c(mySelected, myChoices[myChoices != mySelected])
    }
    selectInput(inputId = "heatfactorlabel", 
                label = "Choose factor(s) for labeling", 
                choices = myChoices, 
                selected = mySelected,
                multiple = T)
  })
  
  # BULK-DGE-GSE - color schemes for heatmap
  heatCols <-  reactive({
    # req(input$heatColors, !d$newExpSetup)
    req(input$heatColors)
    if(input$heatColors == "Red-white-blue") {
      cols <- colorRampPalette(c('red', 'white', 'blue'))(100)
    } else if (input$heatColors == "Viridis") {
      cols <- viridis(100)
    } else if (input$heatColors == "Green-yellow-red") {
      cols <- rev(brewer.pal(n=11, name = "RdYlGn"))
    } else {
      cols <- rev(brewer.pal(n = 11, name = "RdBu"))
    }
  })
  
  # Returns samples to be used in heatmap
  heatmapsamples <- reactive({
    req(ddsout(), ddstran())
    samples <- colnames(ddstran()[[1]])
    if(!is.null(d$bulkExpFactorInfo) && !is.null(input$contrast_groups_only) &&
       input$contrast_groups_only && !is.null(input$dge_list_heatmap) &&
       input$dge_list_heatmap == "DGE filtered" && !(input$dgeexpsetup == "exp1" && input$dgeexp1c == "Rest")) {
      keep <- rep.int(F, times = length(samples))
      names(keep) <- samples
      for(i in names(d$bulkExpFactorInfo)) {
        keep[as.character(ddsout()[[2]][, i]) %in% d$bulkExpFactorInfo[[i]]] <- T
      }
      samples <- names(keep)[keep]
    }
    return(samples)
  })
  
  # Returns genes to be used in heatmap
  heatmapgenes <- reactive({
    req(SubmitData$data_type)
    if(SubmitData$data_type != "Bulk" || is.null(ddstran())) return()

    if(is.null(input$dge_list_heatmap) || input$dge_list_heatmap == "Custom") {
      if(is.null(input$gene_list_heatmap) || input$gene_list_heatmap == "") return()
      genes <- toupper(unlist(strsplit(input$gene_list_heatmap, split = "\\s+|,\\s?")))
      if(length(genes) == 0) return()
      genes <- genes[genes %in% toupper(rownames(ddstran()[[1]]))]
      if(length(genes) == 0) return()
      return(genes)
    }
    if(is.null(dgeout3()) || nrow(dgeout3()) == 0 || is.null(input$dge_filter_heatmap)) return()
    dge <- dgeout3()
    if(input$dge_filter_heatmap == 2) dge <- dge[dge$log2FoldChange > 0, ]
    if(input$dge_filter_heatmap == 3) dge <- dge[dge$log2FoldChange < 0, ]
    dge <- dge[order(1/abs(dge$log2FoldChange), dge$pvalue, decreasing = F), ]
    if(d$n_genes_heatmap != "") {
      if(!is.null(input$fea_all_genes_heatmap) && !(input$fea_all_genes_heatmap)) {
        if(suppressWarnings(is.na(as.integer(na.omit(d$n_genes_heatmap))))) return()
        dge <- dge[1:min(as.integer(d$n_genes_heatmap), nrow(dge)), , drop = F]
      }
    }
    return(rownames(dge))
  })
  
  heattran1 <- reactive({
    if(is.null(ddstran()) || nrow(ddstran()[[1]]) == 0) return()
    if(is.null(heatmapgenes()) || length(heatmapgenes()) == 0) return()
    if(is.null(heatmapsamples()) || length(heatmapsamples()) == 0) return()
    if(!all(heatmapgenes() %in% rownames(ddstran()[[1]]))) return()
    if(!(all(heatmapsamples() %in% colnames(ddstran()[[1]])))) return()
    list(assay(ddstran()[[1]])[heatmapgenes(), heatmapsamples(), drop = F])
  })
  
  # BULK-DGE-GSE - pull genes for gse
  gsetran1 <- reactive({
    req(SubmitData$data_type)
    if(is.null(d$cts_out)) return()
    if(is.null(ddsout())) return()
    if(SubmitData$data_type != "RASL" && is.null(ddstran())) return()
    if(!is.null(input$dge_list) && input$dge_list == "Custom" &&
       (is.null(input$gene_list) || input$gene_list == "")) return()
    
    if(SubmitData$data_type == "RASL") {
      samples <- colnames(ddsout()[[4]])
    } else {
      samples <- colnames(assay(ddstran()[[1]]))
    }
    
    if(is.null(input$dge_list) || input$dge_list == "Custom") {
      if(is.null(input$gene_list) || input$gene_list == "") return()
      genes <- unlist(strsplit(input$gene_list, split = "\\s+|,\\s?"))
      genes <- toupper(genes)
      cts <- ddsout()[[3]]
      inc <- toupper(attr(cts, "dimnames")[[1]])[toupper(attr(cts, "dimnames")[[1]]) %in% genes]
      genes <- genes[genes %in% inc]
      if(length(genes) <= 1) return()
      
      dds.counts <- ddsout()[[3]]
      if(SubmitData$data_type == "RASL") {
        heat.counts <- ddsout()[[4]]
      } else {
        heat.counts <- assay(ddstran()[[1]])
      }
      topID <- order(rowMeans(dds.counts), decreasing = TRUE)
      heat.mat <- heat.counts[topID, samples]
      heat.mat <- heat.mat[toupper(rownames(heat.mat)) %in% genes,]
      if (length(genes) < 50 & length(genes) > 1) {
        heat.mat <- heat.mat[1:length(genes), ,drop = FALSE]
      }
      return(list(heat.mat))
      
    } else if(input$dge_list == "DGE filtered") {
      
      if(is.null(dgeout3()) || is.null(input$dge_filter) || 
         is.null(d$n_genes) || is.null(input$fea_all_genes)) return()
      
      dge <- dgeout3()
      dge_filter = input$dge_filter
      logfc_col = "log2FoldChange"
      pval_col = "pvalue"
      if(dge_filter == 2) dge <- dge[dge[, logfc_col] > 0, ]
      if(dge_filter == 3) dge <- dge[dge[, logfc_col] < 0, ]
      dge <- dge[order(1/abs(dge[, logfc_col]), dge[, pval_col], decreasing = F), ]
      if(d$n_genes != "" & !(input$fea_all_genes)) dge <- dge[1:min(as.integer(d$n_genes), nrow(dge)), ]
      uni_genes <- unique(dge$id)
      if (length(uni_genes) == 0) return()
      # convert genes to a vector
      genes <- dge %>% pull(id)
      if (length(genes) <= 1) return()
      dds.counts <- ddsout()[[3]]
      heat.counts <- ddstran()[[1]]
      heat.counts <- assay(heat.counts)
      topID <- order(rowMeans(dds.counts), decreasing = TRUE)
      heat.mat <- heat.counts[topID, samples]
      heat.mat <- heat.mat[rownames(heat.mat) %in% genes,]
      if (length(genes) < 50 & length(genes) > 1) {
        heat.mat <- heat.mat[1:length(genes), ,drop = FALSE]
      }
      return(list(heat.mat))
    }
  }) 
  
  # BULK-DGE-Heatmap - reactive expression for heatmap data
  heattran2 <- eventReactive(input$submit_list, {
    rescaled_mat <- heattran1()[[1]]
    rescaled_mat_unclipped <- rescaled_mat
    
    keyLabel <- "Normalized expression"
    if(!is.null(input$scale_genes) && input$scale_genes) {
      keyLabel <- "Row z-score"
      rescaled_mat <- t(scale(t(rescaled_mat)))
      rescaled_mat_unclipped <- rescaled_mat
      # Clip values like Seurat
      rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat > 2.5), values = 2.5)
      rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat < -2.5), values = -2.5)
      # re-scale matrix -2 to 2 as per CM's request
      rescaled_mat <- scales::rescale(rescaled_mat, to = c(-2, 2))
      rescaled_mat_unclipped <- scales::rescale(rescaled_mat_unclipped, to = c(-2, 2))
    } else {
      # Clip values like Seurat
      rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat > 6), values = 6)
      rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat < -2.5), values = -2.5)
    }
    
    # replace NAs with 0s
    rescaled_mat[is.na(rescaled_mat)] <-  0
    rescaled_mat_unclipped[is.na(rescaled_mat_unclipped)] <-  0
    # remove rows with all 0s
    rescaled_mat <- rescaled_mat[which(rowSums(rescaled_mat) != 0), , drop = F]
    rescaled_mat_unclipped <- rescaled_mat_unclipped[rownames(rescaled_mat), colnames(rescaled_mat), drop = F]
    
    samples <- colnames(rescaled_mat)
    metaDF <- ddsout()[[2]]
    if(!is.null(input$dge_list_heatmap) && input$dge_list_heatmap == "DGE filtered" && 
       input$dgeexpsetup == "exp1" && input$dgeexp1c == "Rest") {
      metaDF[[input$dgeexp1a]] <- as.character(metaDF[[input$dgeexp1a]])
      metaDF[[input$dgeexp1a]][metaDF[[input$dgeexp1a]] != input$dgeexp1b] <- "Rest"
      metaDF[[input$dgeexp1a]] <- factor(metaDF[[input$dgeexp1a]], levels = c(input$dgeexp1b, "Rest"))
    }
    
    sampleLabels <- apply(metaDF[samples, input$heatfactorlabel, drop = F], 
                          1, function(i) paste(i, collapse = "_"))
    sampleLabels <- factor(sampleLabels, levels = mixedsort(unique(sampleLabels), decreasing = T))
    
    len <- length(unique(sampleLabels))
    cols <- hue_pal()(len)
    pal <- unlist(sapply(1:len, function(x) {
      col <- cols[x]
      names(col) <- unique(sampleLabels)[x]
      return(col)
    }))
    
    customHoverMat <- matrix(paste0(
      matrix(paste0("value: ", as.matrix(round(rescaled_mat_unclipped, digits = 5)), "<br>")),
      matrix(paste0("gene: ", rownames(rescaled_mat_unclipped), "<br>"), nrow = nrow(rescaled_mat_unclipped), ncol = ncol(rescaled_mat_unclipped), byrow = F)
    ), nrow = nrow(rescaled_mat_unclipped), ncol = ncol(rescaled_mat_unclipped))
    
    rowDend <- colDend <- F
    if(input$row_type && nrow(rescaled_mat) >= 2) {
      rowDist <- stats::dist(rescaled_mat_unclipped)
      rowHClust <- hclust(rowDist)
      rowDend <- as.dendrogram(rowHClust)
      if(nrow(rescaled_mat_unclipped) <= 2000) {
        rowDend <- seriate_dendrogram(dend = rowDend, x = rowDist, method = "OLO")
      }
    } 
    
    legendLabels <- unique(sampleLabels)
    if(length(legendLabels) == 1) {
      legendHeight <- 80
    } else {
      legendHeight <- ifelse(test = (length(legendLabels) == 2), yes = 60, no = 40) * length(legendLabels)
    }
    
    heatmapHeight <- max(round(800 * (nrow(heattran1()[[1]]) / 30)), 267)
    if(nrow(heattran1()[[1]]) > 50) heatmapHeight <- 600
    
    # Keep dendrogram and side colorbar same height when heatmap height changes
    subplotHeights <- c(40) / heatmapHeight
    subplotHeights <- c(subplotHeights, 1 - sum(subplotHeights))
    
    if(input$col_type && nrow(rescaled_mat) >= 2) {
      # Keep dendrogram and side colorbar same height when heatmap height changes
      subplotHeights <- c(80, 40) / heatmapHeight
      subplotHeights <- c(subplotHeights, 1 - sum(subplotHeights))
      colDist <- stats::dist(t(rescaled_mat_unclipped))
      colHClust <- hclust(colDist)
      colDend <- as.dendrogram(colHClust)
      if(ncol(rescaled_mat_unclipped) <= 2000) {
        colDend <- seriate_dendrogram(dend = colDend, x = colDist, method = "OLO")
      }
    } else {
      rescaled_mat <- ReorderColsByFactor(x = rescaled_mat, fact = sampleLabels)
      rescaled_mat_unclipped <- rescaled_mat_unclipped[rownames(rescaled_mat), colnames(rescaled_mat), drop = F]
    }
    if(nrow(rescaled_mat) > 50) {
      ticks <- c(FALSE, FALSE)
    } else {
      ticks <- c(FALSE, TRUE)
    }
    
    cols <- heatCols()
    
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
    
    colorbarLabel <- paste(input$heatfactorlabel, collapse = "_")
    df <- data.frame("a" = sampleLabels)
    names(df) <- colorbarLabel
    
    # Keep color scale key same height when heatmap height changes
    colorscaleKeyHeight <- 200 / heatmapHeight
    
    # Color scale and legend y position
    legendYPos <- subplotHeights[length(subplotHeights)]
    
    bottomMargin <- max(d$totalFigureHeight - heatmapHeight, 0)
    if(!is.dendrogram(rowDend) || !is.dendrogram(colDend)) bottomMargin <- 0
    
    return(list(
      rescaled_mat = rescaled_mat,
      rescaled_mat_unclipped = rescaled_mat_unclipped,
      cols = cols,
      ticks = ticks,
      customHoverMat = customHoverMat,
      df = df,
      pal = pal,
      colorscaleKeyHeight = colorscaleKeyHeight,
      legendYPos = legendYPos,
      keyLabel = keyLabel,
      rowDend = rowDend,
      colDend = colDend,
      heatmapHeight = heatmapHeight,
      subplotHeights = subplotHeights,
      colorbarLabel = colorbarLabel,
      legendHeight = legendHeight,
      label = label,
      bottomMargin = bottomMargin,
      sampleLabels = sampleLabels
    ))
  })
  
  # BULK-DGE-GSE - warning messages for heatmap
  observeEvent(input$submit_list, {
    # Reset event_data("plotly_selected")
    runjs("Shiny.setInputValue('plotly_selected-A', null);") 
    d$heatmap1Click <- NULL
    d$newHeat <- F
    if(d$n_genes_heatmap != "" & suppressWarnings(is.na(as.integer(na.omit(d$n_genes_heatmap))))) {
      d$heatmap_warning_msg <- "Please provide an integer in 'Use top n genes'."
    } else if(is.null(input$dge_list_heatmap) || input$dge_list_heatmap == "Custom") {
      if(length(unlist(strsplit(genelist_heatmap(), split = "\\s+|,\\s?"))) == 0) {
        d$heatmap_warning_msg <- "No genes were submitted."
      }
    } else if(is.null(heattran1())) {
      d$heatmap_warning_msg <- "No DE genes."
    }
    if(!is.null(heattran1())) d$heatmap_warning_msg <- ""
    
    if(is.null(heattran1())) return()
    
    samples <- colnames(heattran1()[[1]])
    sampleLabels <- apply(ddsout()[[2]][samples, input$heatfactorlabel, drop = F], 
                          1, function(i) paste(i, collapse = "_"))
    
    legendLabels <- unique(sampleLabels)
    
    # Legend height
    if(length(legendLabels) == 1) {
      legendHeight <- 80
    } else {
      legendHeight <- ifelse(test = (length(legendLabels) == 2), yes = 60, no = 40) * length(legendLabels)
    }
    
    heatmapHeight <- max(round(800 * (nrow(heattran1()[[1]]) / 30)), 267)
    if(nrow(heattran1()[[1]]) > 50) heatmapHeight <- 600
    
    if(input$col_type) {
      if(legendHeight > heatmapHeight) {
        totalFigureHeight <- legendHeight + 120
      } else {
        totalFigureHeight <- heatmapHeight
      }
    } else {
      if(legendHeight > heatmapHeight) {
        totalFigureHeight <- legendHeight + 40
      } else {
        totalFigureHeight <- heatmapHeight
      }
    }
    d$totalFigureHeight <- totalFigureHeight
  })
  
  # BULK-DGE-GSE - heatmap
  output$heatplot1 <- renderPlotly({
    if(is.null(heattran2())) return()
    isolate({
      p <- heatmaply(
        x = heattran2()$rescaled_mat,
        plot_method = "plotly",
        colors = heattran2()$cols,
        dend_hoverinfo = F,
        custom_hovertext = heattran2()$customHoverMat,
        showticklabels = heattran2()$ticks,
        col_side_colors = heattran2()$df,
        col_side_palette = heattran2()$pal,
        colorbar_len = heattran2()$colorscaleKeyHeight,
        key.title = heattran2()$keyLabel,
        colorbar_xpos = -0.3,
        colorbar_ypos = heattran2()$legendYPos,
        colorbar_yanchor = "top",
        fontsize_row = 14,
        fontsize_col = 14,
        Rowv = heattran2()$rowDend,
        Colv = heattran2()$colDend,
        revC = T,
        width = 1000,
        height = heattran2()$heatmapHeight,
        subplot_heights = heattran2()$subplotHeights
      ) %>%
        colorbar(
          tickfont = list(size = 14, family = "Noto Sans JP"),
          titlefont = list(size = 16, family = "Noto Sans JP"),
          which = 2
        ) %>%
        colorbar(
          title = heattran2()$colorbarLabel,
          len = heattran2()$legendHeight,
          lenmode = "pixels",
          y = heattran2()$legendYPos,
          yanchor = "top",
          tickfont = list(size = 14, family = "Noto Sans JP"),
          titlefont = list(size = 16, family = "Noto Sans JP"),
          which = 1
        ) %>%
        layout(
          font = list(family = "Noto Sans JP", size = 14),
          hoverlabel = heattran2()$label,
          legend = list(font = list(size = 14, family = "Noto Sans JP")),
          margin = list(b = heattran2()$bottomMargin)
        ) 
      
      if(is.dendrogram(heattran2()$colDend) && is.dendrogram(heattran2()$rowDend)) {
        p$x$data[[3]]$hoverinfo <- "none"
      } else if(!is.dendrogram(heattran2()$colDend)) {
        p$x$data[[1]]$text <- "none"
      } else {
        p$x$data[[2]]$text <- "none"
        p$x$layout$yaxis$showticklabels <- F
      }
      
      isolate(return(p))
    })
  })
  
  output$heatplot1_ui <- renderUI({
    plotlyOutput("heatplot1", width = 1000, height = d$totalFigureHeight)
  })
  
  
  # BULK-DGE-GSE - download button - heatmap (PDF)
  output$dlqcheatplot1pdf <- renderUI({
    req(input$submit_list, !d$newHeat)
    if(is.null(heattran2())) return()
    downloadButton("dlqcheatplot1pdfimg", "Download static plot (PDF)")
  })
  
  # BULK-DGE-GSE - download file - heatmap (PDF)
  output$dlqcheatplot1pdfimg <- downloadHandler(
    filename =  function() {
      paste("qc-heatmap.pdf")
    },
    content = function(file) {
      pdf(file = file, width = 7.5, height = 7.5, onefile = F) 
      qcHeatMap(
        heat = heattran2()$rescaled_mat,
        color = heattran2()$cols,
        rows = heattran2()$rowDend,
        col = heattran2()$colDend,
        nam = heattran2()$sampleLabels,
        keyLabel = heattran2()$keyLabel,
        colorLabel = heattran2()$colorbarLabel,
        pal = heattran2()$pal,
        type = "pdf"
      )
      dev.off()
    }
  )
  
  # BULK-DGE-GSE - download button - heatmap (PNG)
  output$dlqcheatplot1png <- renderUI({
    req(input$submit_list, !d$newHeat)
    if(is.null(heattran2())) return()
    downloadButton("dlqcheatplot1pngimg", "Download static plot (PNG)")
  })
  
  
  # BULK-DGE-GSE - download file - heatmap (PNG)
  output$dlqcheatplot1pngimg <- downloadHandler(
    filename =  function() {
      paste("qc-heatmap.png")
    },
    content = function(file) {
      png(file, width = 1200, height = 850)
      qcHeatMap(
        heat = heattran2()$rescaled_mat,
        color = heattran2()$cols,
        rows = heattran2()$rowDend,
        col = heattran2()$colDend,
        nam = heattran2()$sampleLabels,
        keyLabel = heattran2()$keyLabel,
        colorLabel = heattran2()$colorbarLabel,
        pal = heattran2()$pal,
        type = "png"
      )
      dev.off()
    }
  )
  
  # BULK-DGE-GSE - choose factor for heatmap
  output$heatfactor <- renderUI({
    req(input$submit_list, !d$newHeat, d$heatmap1Click)
    s <- d$heatmap1Click
    validate(need(s != "",
                  message = ""))
    if(is.null(s)) return()
    tmp <- ddsout()[[2]]
    selectInput(inputId = "heatfactor",
                label = "Choose factor",
                choices = colnames(tmp))
  })
  
  # BULK-DGE-GSE - store heatmap plotly click in object d
  observe({
    d$heatmap1Click <- event_data("plotly_click", source = "A")
  })

  # BULK-DGE-GSE - count plot from heatmap
  output$heatplot2 <- renderPlotly({
    req(SubmitData$data_type)
    if(is.null(input$heatfactor) || is.null(d$newHeat) || d$newHeat ||
       is.null(input$submit_list) || input$submit_list == 0 ||
       is.null(input$scale_genes)) return()
    s <- d$heatmap1Click
    if(is.null(s)) return()
    
    if(is.dendrogram(heattran2()$rowDend)) {
      gene <- rev(labels(heattran2()$rowDend))[s$y]
    } else {
      gene <- rev(rownames(heattran2()$rescaled_mat))[s$y]
    }
    
    rc.data <- counts(ddsout()[[1]])
    
    test <- getGenes(
      rc.data = rc.data,
      id = gene,
      coldata = ddsout()[[2]],
      type = "Bulk"
    )
    
    m <- list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
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
    
    fact <- test[, input$heatfactor]
    y <- test[, "counts"]
    
    x_ord <- factor(fact, levels = mixedsort(unique(fact), decreasing = F))
    pal <- hue_pal()(length(levels(x_ord)))
    names(pal) <- levels(x_ord)
    
    tooltips <- paste0("<b>Sample:</b> ",
                       test$Sample,
                       "<br />",
                       "<b>Counts:</b> ",
                       round(test$counts, 3))
    
    plot_ly(
      data = test,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      x = x_ord,
      y = y,
      color = x_ord,
      colors = pal,
      text = tooltips,
      marker = list(size = 7),
      hoverinfo = "text",
      hoverlabel = label
    ) %>%
      layout(font = list(size = 12, family = "Noto Sans JP"),
             title = paste0(gene, " counts"),
             xaxis = list(title = input$heatfactor),
             yaxis = list(title = "Normalized counts"),
             margin = m
      )
  })
  
  # BULK-DGE-GSE - download button - heat counts (PDF)
  output$dlqcheatplot2pdf <- renderUI({
    req(input$submit_list, !d$newHeat, d$heatmap1Click)
    downloadButton("dlqcheatplot2pdfimg", "Download static plot (PDF)")
  })
  
  # BULK-DGE-GSE - download file - heat counts (PDF)
  output$dlqcheatplot2pdfimg <- downloadHandler(
    filename =  function() {
      paste("qc-heat-counts.pdf")
    },
    content = function(file) {
      s <- d$heatmap1Click
      if(is.null(s)) return()
      
      if(is.dendrogram(heattran2()$rowDend)) {
        gene <- rev(labels(heattran2()$rowDend))[s$y]
      } else {
        gene <- rev(rownames(heattran2()$rescaled_mat))[s$y]
      }
      
      rc.data <- counts(ddsout()[[1]])
      
      test <- getGenes(
        rc.data = rc.data,
        id = gene,
        coldata = ddsout()[[2]],
        type = "Bulk"
      )

      fact <- test[, input$heatfactor]
      fact <- factor(fact, levels = mixedsort(unique(fact), decreasing = F))
      y <- test[, "counts"]

      pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
      qcHeatCount(
        data = test,
        fact = fact,
        var = y,
        xaxis = input$heatfactor,
        title = gene,
        type = "bulk"
      )
      dev.off()
    }
  )
  
  # BULK-DGE-GSE - download button - heat counts (PNG)
  output$dlqcheatplot2png <- renderUI({
    req(input$submit_list, !d$newHeat, d$heatmap1Click)
    downloadButton("dlqcheatplot2pngimg", "Download static plot (PNG)")
  })
  
  # BULK-DGE-GSE - download plot - heat counts (PNG)
  output$dlqcheatplot2pngimg <- downloadHandler(
    filename =  function() {
      paste("qc-heat-counts.png")
    },
    content = function(file) {
      s <- d$heatmap1Click
      if(is.null(s)) return()
      
      if(is.dendrogram(heattran2()$rowDend)) {
        gene <- rev(labels(heattran2()$rowDend))[s$y]
      } else {
        gene <- rev(rownames(heattran2()$rescaled_mat))[s$y]
      }
      
      rc.data <- counts(ddsout()[[1]])
      
      test <- getGenes(
        rc.data = rc.data,
        id = gene,
        coldata = ddsout()[[2]],
        type = "Bulk"
      )

      fact <- test[, input$heatfactor]
      fact <- factor(fact, levels = mixedsort(unique(fact), decreasing = F))
      y <- test[, "counts"]
      
      png(file, width = 800, height = 750)
      qcHeatCount(
        data = test,
        fact = fact,
        var = y,
        xaxis = input$heatfactor,
        title = gene,
        type = "bulk"
      )
      dev.off()
    }
  )
  
  ######### CLUSTERING TAB CODE: BULK RNA-SEQ#####################
  
  # BULK-DGE-CLUS - actionbutton to submit clustering
  output$goclust <- renderUI({
    actionButton(
      inputId = "goclust",
      label = "Launch clustering analysis"
    )
  })
  
  # BULK-DGE-CLUS - wgcna warning
  output$wgcna_warning_msg <- renderText(d$wgcna_warning)
  
  # BULK-DGE-CLUS - set cluster in d
  observeEvent(input$goclust, {
    d$newCluster <- F
  })
  
  # BULK-DGE-CLUS - header for clustering
  output$headclust <- renderUI({
    h4("Clustering analysis")
  })
  
  # BULK-DGE-CLUS - choose number of variable genes
  observeEvent(SubmitData$data_type, {
    output$clustvarnumber <- renderUI({
      textInput(
        inputId = "clustvarnumber",
        label = div(style = "width: 200px;", "Top variable genes"),
        value = 500
      )
    })
  })
  
  # BULK-DGE-CLUS - choose clustering algorithm
  observeEvent(SubmitData$data_type, {
    output$clustalg <- renderUI({
      selectInput(
        inputId = "clustalg",
        label = div(style = "width: 200px;", "Clustering algorithm"),
        choices = c(
          "WGCNA" = "wgcna",
          "K-Medoids" = "kmed"
        ),
        selected = "WGCNA"
      )
    })
  })
  
  # BULK-DGE-CLUS - set min module size
  output$min_module_size <- renderUI({
    req(SubmitData$data_type, input$clustalg)
    if(SubmitData$data_type == "Single-cell") return()
    if(input$clustalg != "wgcna") return()
    div(
      style = "width: 100px; margin-left: 40px;",
      numericInput(
        inputId = "min_module_size",
        label = div(style = "width: 180px;", "Min. module size"),
        value = 30,
        min = 8,
        step = 1
      )
    )
  })
  
  # BULK-DGE-CLUS - update min module size for wgcna
  observeEvent(input$min_module_size, {
    if(is.null(input$min_module_size)) return()
    if(is.numeric(input$min_module_size)) {
      if(input$min_module_size >= 8) return()
    } 
    updateNumericInput(
      inputId = "min_module_size", 
      value = 30
    )
  })  
  
  # BULK-DGE-CLUS - reactive - perform clustering get 
  # variable counts (WIP - very mess atm...)
  clustout <- eventReactive(input$goclust, {
    clust_t0 <- as.numeric(Sys.time())
    num <- input$clustvarnumber
    if(SubmitData$data_type == "Bulk") {
      cts <- ddsout()[[1]]
      cts <- assay(cts)
      tran <- ddstran()[[1]]
      tran <- assay(tran)
    } else {
      cts <- ddsout()[[4]]
      tran <- ddsout()[[4]]
    }
    topID <- order(rowVars(cts), decreasing = TRUE)
    cts.var <- tran[topID, ]
    dds_mat <- cts.var[1:num, ]
    if (input$clustalg == "wgcna") {
      withProgress(message = "Running WGCNA...", value = 0, {
        incProgress(1/3)
        enableWGCNAThreads(max(1, parallel::detectCores()-2))
        dds_mat <- as.matrix(dds_mat)
        gene.names <- sort(rownames(dds_mat))
        datExpr <- t(dds_mat)
        # Create an object called "datTraits" that contains your
        # trait data
        if(SubmitData$data_type == "Bulk") {
          datTraits <- colData(ddsout()[[1]])
        } else {
          datTraits <- ddsout()[[2]]
        }
        nLevels <- apply(datTraits, 2, function(col) length(unique(col)))
        datTraits <- datTraits[, nLevels > 1 & nLevels <= 12, drop = F]
        A <- adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
        k <- as.numeric(apply(A,2,sum))-1 # standardized connectivity
        Z.k <- scale(k)
        thresholdZ.k <- -2.5 # often -2.5
        outlierColor <- ifelse(Z.k<thresholdZ.k,"red","black")
        sampleTree <- flashClust(stats::as.dist(1-A), method = "average")
        
        traitColors <- datTraits
        for(i in 1:ncol(traitColors)) {
          if(i %% 2 == 0) {
            traitColors[, i] <- labels2colors(traitColors[, i], 
                                              colorSeq = brewer.pal(12, "Set3"))
          } else {
            traitColors[, i] <- labels2colors(traitColors[, i],
                                              colorSeq = brewer.pal(12, "Paired"))
          }
        }
        
        dimnames(traitColors)[[2]] <- paste(names(datTraits))
        datColors <- data.frame(outlier = outlierColor, traitColors)
        incProgress(1/3)
        # TOM analysis - (computationally expensive)
        softPower <- 18
        adjacency <- adjacency(datExpr, power = softPower, type = "signed") #specify network type
        TOM <- TOMsimilarity(adjacency, TOMType = "signed") # specify network type
        dissTOM <- 1 - TOM
        geneTree <- flashClust(stats::as.dist(dissTOM), method="average")
        # This sets the minimum number of genes to cluster into a module
        minModuleSize <- input$min_module_size
        dynamicMods <- cutreeDynamic(
          dendro = geneTree,
          distM = dissTOM,
          deepSplit = 2,
          pamRespectsDendro = FALSE,
          minClusterSize = minModuleSize
        )
        if(length(unique(dynamicMods)) == 1) {
          if(minModuleSize == 8) return()
          minModuleSizes <- (minModuleSize-1):8
          for(i in minModuleSizes) {
            dynamicMods <- cutreeDynamic(
              dendro = geneTree,
              distM = dissTOM,
              deepSplit = 2,
              pamRespectsDendro = FALSE,
              minClusterSize = i
            )
            if(length(unique(dynamicMods)) > 1) {
              d$wgcna_minModuleSizeTried <- i
              break
            }
          }
          if(length(unique(dynamicMods)) == 1) {
            d$wgcna_minModuleSizeTried <- 8
            return()
          }
        } else {
          d$wgcna_minModuleSizeTried <- NULL
        }
        dynamicColors <- labels2colors(dynamicMods, colorSeq = brewer.pal(12, "Set3"))
        MEList <- moduleEigengenes(
          datExpr,
          colors = dynamicColors,
          softPower = 18
        )
        MEs <- MEList$eigengenes
        MEDiss <- 1 - cor(MEs)
        METree <- flashClust(stats::as.dist(MEDiss), method = "average")
        # set a threhold for merging modules. In this example we are
        # not merging so MEDissThres=0.0
        MEDissThres <- 0.0
        merge <- mergeCloseModules(
          datExpr,
          dynamicColors,
          cutHeight = MEDissThres,
          verbose = 3
        )
        mergedColors <- merge$colors
        mergedMEs <- merge$newMEs
        # Set the diagonal of the dissimilarity to NA
        diag(dissTOM) = NA;
        # Export modules to data frame
        module_colors= setdiff(unique(dynamicColors), "grey")
        modlist <- list()
        for (color in module_colors){
          module = gene.names[which(dynamicColors == color)]
          modlist[[color]] <- list(
            gene = module
          )
        }
        moddf <- data.frame(unlist(modlist))
        moddf$module <- gsub("\\..*", "", row.names(moddf))
        colnames(moddf)[1] <- "gene"
        rownames(moddf) <- seq_len(nrow(moddf))
        moddf$gene <- as.character(moddf$gene)
        moddf$module <- as.factor(moddf$module)
        sampleDF <- data.frame(
          sample = sampleTree$labels,
          outlier = datColors$outlier
        )
        disableWGCNAThreads()
        incProgress(1/3)
        return(
          list(
            sampleTree,
            datColors,
            geneTree,
            dynamicColors,
            mergedColors,
            dissTOM,
            moddf,
            sampleDF
          )
        )
      })
    } else if (input$clustalg == "kmed") {
      withProgress(message = "Running k-medoids...", value = 0, {
        incProgress(1/2)
        num <- as.matrix(dds_mat)
        mrwdist <- distNumeric(num, num, method = "mrw")
        # Detect cores of machine
        nclust <- parallel::detectCores() - 1
        message("Using ", nclust, " threads...")
        result <- fastkmed(mrwdist, ncluster = nclust, iterate = 50)
        # a simple and fast k-medoids function for bootstrap evaluation
        parkboot <- function(x, nclust) {
          res <- fastkmed(x, nclust, iterate = 50)
          return(res$cluster)
        }
        fastkmedboot <- clustboot(mrwdist, nclust = nclust, parkboot, nboot = 50)
        # consensus matrix
        wardorder <- function(x, nclust) {
          res <- fastcluster::hclust(x, method = "ward.D2")
          member <- cutree(res, nclust)
          return(member)
        }
        consensusfastkmed <- consensusmatrix(fastkmedboot, nclust = nclust, wardorder)
        # data frame generation
        output <- data.frame(
          gene_id = rownames(num),
          cluster = result$cluster
        )
        rownames(output) <- seq_len(nrow(output))
        output$gene_id <- as.character(output$gene_id)
        output$cluster <- as.factor(output$cluster)
        return(
          list(
            consensusfastkmed,
            output
          )
        )
        incProgress(2/2)
      })
    } 
  })
  
  # BULK-DGE-CLUS - create d$legendList and output DT datatables for WGCNA legend
  observeEvent(input$goclust, {
    if(input$clustalg != "wgcna") {
      updateTabsetPanel(inputId = "clustering_tabsetPanel", selected = "clustPlotW01")
      d$showClusteringLinks <- T
      return()
    }
    if(is.null(clustout())) {
      d$wgcna_warning <- "No gene modules found after setting 'Min. module size' to 8."
      return()
    }
    if(!is.null(d$wgcna_minModuleSizeTried)) {
      d$wgcna_warning <- paste0("Gene modules found after lowering 'Min. module size' to ",
                                d$wgcna_minModuleSizeTried, ".")
    } else {
      d$wgcna_warning <- NULL
    }
    d$showClusteringLinks <- T
    updateTabsetPanel(inputId = "clustering_tabsetPanel", selected = "clustPlotW02")
  })
  
  # BULK-DGE-CLUS - header for sample dendrogram or kmed
  output$headclustplotW01 <- renderUI({
    req(input$goclust, clustout(), !d$newCluster)
    isolate({
      if (input$goclust == 0) {
        return()
      } else if (is.null(clustout())) {
        return()
      } else if (input$clustalg == "wgcna") {
        return()
      } else if (input$clustalg == "kmed") {
        if (length(clustout()) == 2) {
          h4(strong("K-medoids - consensus matrix heatmap"))
        }
      } 
    })
  })
  
  # BULK-DGE-CLUS - dendrogram wgcna sample dendrogram plot
  output$clustplotW01 <- renderPlot({
    req(input$goclust, clustout())
    if(input$goclust == 0) return()
    if(length(clustout()) != 2) return()
    
    # K-medoids
    clustheatmap(clustout()[[1]], title = "K-medoids consensus matrix heatmap")
  })
  
  # Use renderUI to avoid blank area for sample dendrogram
  output$clustplotW01_ui <- renderUI({
    req(input$goclust, clustout())
    if(input$goclust == 0) return()
    if(length(clustout()) != 2) return()
    fluidRow(column(12, align = "left", plotOutput("clustplotW01", width = 900, height = 500)))
  })
  
  
  # BULK-DGE-CLUS - sample dendrogram download button (PNG)
  output$downloadclustplotW01png <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "kmed") {
        if (length(clustout()) == 2) {
          downloadButton(
            "downloadclustplotK01pngimg",
            "Download plot (PNG)"
          )
        }
      } 
    })
  })
  
  # BULK-DGE-CLUS - sample dendrogram download button (PDF)
  output$downloadclustplotW01pdf <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "kmed") {
        if (length(clustout()) == 2) {
          downloadButton(
            "downloadclustplotK01pdfimg",
            "Download plot (PDF)"
          )
        }
      } 
    })
  })
  
  # BULK-DGE-CLUS - download file for consensus matrix
  output$downloadclustplotK01pngimg <- downloadHandler(
    filename = function() {
      paste("kmed-consensus-matrix.png")
    },
    content = function(file) {
      req(clustout())
      png(file, width = 800, height = 600)
      sample <- clustout()[[1]]
      output <- clustheatmap(sample, "K-medoids consensus matrix heatmap")
      print(output)
      dev.off()
    }
  )
  
  # BULK-DGE-CLUS - header for gene dendrogram
  output$headclustplotW02 <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          h4(strong("WGCNA - gene dendrogram"))
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - gene dendrogram for wgcna
  output$clustplotW02 <- renderPlot({
    validate(need(input$goclust, ""))
    req(input$goclust)
    isolate({
      if (input$goclust == 0) {
        return()
      } else if (input$clustalg == "wgcna") {
        validate(
          need(input$goclust != "", "")
        )
        if (length(clustout()) == 8) {
          withProgress(message = "Creating gene dendrogram...", value = 0, {
            incProgress(1/2)
            geneTree <- clustout()[[3]]
            dynamicColors <- clustout()[[4]]
            datColors <- data.frame("Dynamic tree cut" = dynamicColors, check.names = F)
            incProgress(1/2)
            GGDend(geneTree, 
                   colorDF = datColors, 
                   leafLabels = F, 
                   touchZero = F)
          })
        }
        else {
          return()
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - gene dendrogram download button (PNG)
  output$downloadclustplotW02png <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        downloadButton(
          "downloadclustplotW02pngimg",
          "Download plot (PNG)"
        )
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - header to specify where to download gene modules
  output$headclustmoddown <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          h4(strong("WGCNA - download modules"))
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - gene dendrogram download file (PNG)
  output$downloadclustplotW02pngimg <- downloadHandler(
    filename = function() {
      paste("wgcna-gene-dendrogram.png")
    },
    content = function(file) {
      geneTree <- clustout()[[3]]
      dynamicColors <- clustout()[[4]]
      datColors <- data.frame("Dynamic tree cut" = dynamicColors, check.names = F)
      png(file, width = 800, height = 400)
      g <- GGDend(geneTree, 
                  colorDF = datColors, 
                  leafLabels = F, 
                  touchZero = F, 
                  plotTitle = "Gene dendrogram and module colors",
                  plotTitleSize = 30)
      print(g)
      dev.off()
    }
  )
  
  # BULK-DGE-CLUS - gene dendrogram download button (PDF)
  output$downloadclustplotW02pdf <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        downloadButton(
          "downloadclustplotW02pdfimg",
          "Download plot (PDF)"
        )
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - gene dendrogram download file (PDF)
  output$downloadclustplotW02pdfimg <- downloadHandler(
    filename = function() {
      paste("wgcna-gene-dendrogram.pdf")
    },
    content = function(file) {
      geneTree <- clustout()[[3]]
      dynamicColors <- clustout()[[4]]
      datColors <- data.frame("Dynamic tree cut" = dynamicColors, check.names = F)
      pdf(file, width = 8, height = 6.5)
      g <- GGDend(geneTree, 
                  colorDF = datColors, 
                  leafLabels = F, 
                  touchZero = F, 
                  plotTitle = "Gene dendrogram and module colors",
                  plotTitleSize = 30)
      print(g)
      dev.off()
      
    }
  )
  
  # BULK-DGE-CLUS - header for sample dendrogram
  output$headclustplotW03 <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          h4(strong("WGCNA - topological overlap matrix"))
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - TOM plot
  output$clustplotW03 <- renderPlot({
    validate(need(input$goclust, ""))
    req(input$goclust, clustout())
    isolate({
      if(input$goclust == 0) return()
      if(input$clustalg != "wgcna") return()
      if(length(clustout()) != 8) return()
      withProgress(message = "Creating TOM plot...", value = 0, {
        incProgress(1/2)
        geneTree <- clustout()[[3]]
        dynamicColors <- clustout()[[4]]
        dissTOM <- clustout()[[6]]
        incProgress(1/2)
        GGTom(geneTree, colors = dynamicColors, distMat = dissTOM, plotTitle = NULL)
      })
    })
  })
  
  # BULK-DGE-CLUS - gene module download file
  output$downloadclustmod2 <- downloadHandler(
    filename = function() {
      paste("wgcna-gene-modules.csv")
    },
    content = function(file) {
      moddf <- clustout()[[7]]
      colorPalette <- c("teal","lightyellow","lavender",
                        "salmon","blue","orange","green",
                        "lightpink","gray","violet","mint","yellow")
      names(colorPalette) <- brewer.pal(12, "Set3")
      moddf$module <- colorPalette[moddf$module]
      write.csv(moddf, file, row.names = FALSE)
    }
  )
  
  # BULK-DGE-CLUS - download button for TOM plot (PNG)
  output$downloadclustplotW03png <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          downloadButton(
            "downloadclustplotW03pngimg",
            "Download plot (PNG)"
          )
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - download button for TOM plot (PDF)
  output$downloadclustplotW03pdf <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          downloadButton(
            "downloadclustplotW03pdfimg",
            "Download plot (PDF)"
          )
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - TOM file download (PNG)
  output$downloadclustplotW03pngimg <- downloadHandler(
    filename = function() {
      paste("wgcna-tom-plot.png")
    },
    content = function(file) {
      geneTree <- clustout()[[3]]
      dynamicColors <- clustout()[[4]]
      dissTOM <- clustout()[[6]]
      png(file, width = 600, height = 400)
      g <- GGTom(geneTree, colors = dynamicColors, 
                 distMat = dissTOM)
      print(g)
      dev.off()
    }
  )
  
  # BULK-DGE-CLUS - download button gene module
  output$downloadclustmod <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          downloadButton(
            "downloadclustmod2",
            "Download gene modules (CSV)"
          )
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - TOM file download (PDF)
  output$downloadclustplotW03pdfimg <- downloadHandler(
    filename = function() {
      paste("wgcna-tom-plot.pdf")
    },
    content = function(file) {
      geneTree <- clustout()[[3]]
      dynamicColors <- clustout()[[4]]
      dissTOM <- clustout()[[6]]
      pdf(file, width = 8, height = 8)
      g <- GGTom(geneTree, colors = dynamicColors, 
                 distMat = dissTOM)
      print(g)
      dev.off()
      
    }
  )
  
  # BULK-DGE-CLUS - download button sample module 
  output$downloadclustsample <- renderUI({
    req(clustout())
    isolate({
      if (input$clustalg == "wgcna") {
        if (length(clustout()) == 8) {
          downloadButton(
            "downloadclustsample2",
            "Download sample modules (CSV)"
          )
        }
      } else {
        return()
      }
    })
  })
  
  # BULK-DGE-CLUS - sample module file download
  output$downloadclustsample2 <- downloadHandler(
    filename = function() {
      paste("wgcna-sample-modules.csv")
    },
    content = function(file) {
      sampledf <- clustout()[[8]]
      if(SubmitData$data_type == "Bulk") {
        metaDF <- as.data.frame(colData(ddsout()[[1]]))
      } else {
        metaDF <- ddsout()[[2]]
      }
      metaDF <- metaDF[sampledf$sample, ]
      sampledf <- cbind(sampledf, metaDF)
      sampledf <- sampledf[, colnames(sampledf) != "outlier"]
      write.csv(sampledf, file, row.names = FALSE)
    }
  )
  
  # BULK-DGE-CLUS - download file for consensus matrix
  output$downloadclustplotK01pdfimg <- downloadHandler(
    filename = function() {
      paste("kmed-consensus-matrix.pdf")
    },
    content = function(file) {
      pdf(file, width = 8, height = 6)
      consensusfastkmed <- clustout()[[1]]
      p <- clustheatmap(
        consensusfastkmed,
        "K-medoids consensus matrix heatmap"
      )
      print(p)
      dev.off()
    }
  )
  
  # BULK-DGE-CLUS - download button for gene module
  output$downloadclustmodK <- renderUI({
    req(clustout(), !d$newCluster)
    if(length(clustout()) != 2) return()
    downloadButton("downloadclustmodK2", label = "Download clusters (CSV)")
  })
  
  # BULK-DGE-CLUS - gene module file download
  output$downloadclustmodK2 <- downloadHandler(
    filename = function() {
      paste("kmed-gene-clusters.csv")
    },
    content = function(file) {
      clustdf <- clustout()[[2]]
      write.csv(clustdf, file, row.names = FALSE, col.names = TRUE)
    }
  )

  ### SECTION 03 - SCRNA-SEQ DGE
  
  ######### HELP BUTTONS - ALL FOR SCRNA-SEQ######################
  
  # SC-DGE-OR - help for overview tab
  observeEvent(input$sc_resolution_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Resolution selection")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html#clustree", 
            style = "color: black !important;"
          )
        )
      ),
      size = "m",
      easyClose = T,
      includeMarkdown("markdown/help/resolution_selection_help.md")
    ))
  })
  
  # SC-DGE-OVER - help for overview tab
  observeEvent(input$sc_overview_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Overview")),
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
      includeMarkdown("markdown/help/sc_overview_help.md")
    ))
  })
  
  # SC-DGE-GE - help for gene expression tab
  observeEvent(input$sc_geneexpression_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Gene expression")),
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
      includeMarkdown("markdown/help/sc_geneExpression_help.md")
    ))
  })
  
  # SC-DGE-DGE - help for dge tab
  observeEvent(input$sc_dgeanalysis_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Differential gene expression (DGE) by cluster")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html#dge", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/sc_dgeByCluster_help.md")
    ))
  })
  
  # SC-DGE-CUST - help for custom dge tab
  observeEvent(input$sc_dge_factor_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Custom differential gene expression (DGE)")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html#dge", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/sc_dgeCustom_help.md")
    ))
  })

  # SC-DGE-VOL - help for volcano tab
  observeEvent(input$sc_volcanoplots_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Volcano plot")),
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
      includeMarkdown("markdown/help/sc_volcano_help.md")
    ))
  })
  
  # SC-DGE-GSE - help for gse & heatmap tab
  observeEvent(input$sc_gse_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Gene set enrichment (GSE)")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html#gene-set-enrichment-1", 
            style = "color: black !important;"
          )
        )
      ),
      size = "l",
      easyClose = T,
      includeMarkdown("markdown/help/sc_gse_help.md")
    ))
  })
  
  # SC-DGE-MSC - help for manually select cells tab
  observeEvent(input$sc_manuallyselectcells_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Manually select cells for differential gene expression (DGE)")),
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
      includeMarkdown("markdown/help/sc_manuallyselectcells_help.md")
    ))
  })
  
  # SC-DGE-CLUS - help for clustering
  observeEvent(input$sc_postdge_clustering_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Clustering")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "m",
      easyClose = T,
      includeMarkdown("markdown/help/clustering_help.md")
    ))
  })
  
  # SC-DGE-CCC - help for combining clusters or cells
  observeEvent(input$sc_combineClusters_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Merge clusters")),
        div(
          style = "display: inline-block; float: right;",
          a(
            "Tutorial", target = "_blank", href = "ncats-SEQUIN-tutorial.html", 
            style = "color: black !important;"
          )
        )
      ),
      size = "m",
      easyClose = T,
      includeMarkdown("markdown/help/sc_mergeClusters_help.md")
    ))
  })

  # SC-DGE-CCC - help for combining clusters or cells
  observeEvent(input$sc_groupCells_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Group cells by gene expression")),
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
      includeMarkdown("markdown/help/sc_groupCells_help.md")
    ))
  })
  
  ######### SELECT OPTIMAL RESOLUTION: SCRNA-SEQ ############
  
  res_spec <- reactive({
    req(d$resType)
    if (d$resType == "seurat_res") return()
    meta <- ddsout()[[1]]
    meta <- data.frame(colData(meta))
    meta <- length(unique(meta$seurat_clusters))
    d$res <- "Clust"
    d$MD <- getMD(d$inD)
    return(meta)
  })
  
  # SC-DGE-OR - select resolution for umap plot
  output$overall_res <- renderUI({
    req(d$resType)
    if(d$resType != "seurat_res") return()
    choices <- clustList()
    choices <- choices[!(grepl("^Comp:", x = choices))]
    selectInput(
      inputId = "overall_res",
      label = "Resolution", 
      choices = choices
    )
  })
  
  # SC-DGE-OR - specify resolution for seurat obj
  # in object d
  observeEvent(input$overall_res,{
    d$sc_dge_byfactor <- NULL
    if(d$resType == "seurat_res") {
      resTemp <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
      if(!(resTemp %in% names(seurat_only()))) resTemp <- names(seurat_only())[1]
      d$inD <- seurat_only()[[resTemp]]
      resTemp <- gsub(":.*$", replacement = "", x = input$overall_res)
      if(!(resTemp %in% names(d$SCV))) resTemp <- names(d$SCV)[1]
      d$res <- resTemp
      d$meta$seurat_clusters <- d$meta[[resTemp]]
      d$meta$silWidth <- d$meta[[paste0(resTemp, "_silWidth")]]
    } 
    d$MD <- getMD(d$inD)
  })
  
  # Update sc profiler and dimred modules when d$inD changes (i.e., when resolution changes)
  observeEvent(d$inD, {
    callModule(profiler, id = "profiler_sc", seuratData = d$inD)
    callModule(dimred, 
               id = "dda_dimred_sc", 
               dat = d$inD,
               dimredMethods = c("pca", "umap", "tsne"))
  })
  
  # SC-DGE-OR - return clustree plot for multiple seurat resolutions
  output$clustree_png <- renderPlot({
    req(seurat_only(), d$resValues, d$resType)
    if(d$resType != "seurat_res" || length(d$resValues) < 3 ||
       length(seurat_only()) < 3) return()
    
    seur <- seurat_only()
    resNames <- names(seur)
    only_m <- lapply(resNames, function(resName) {
      excludeCols <- grep(
        pattern = "RNA_snn_res.", x = colnames(seur[[resName]]@meta.data), value = T
      )
      excludeCols <- c("seurat_clusters", excludeCols[excludeCols != paste0("RNA_snn_", resName)])
      m <- seur[[resName]]@meta.data[, !(colnames(seur[[resName]]@meta.data) %in% excludeCols)]
    })
    
    only_m <- only_m %>% reduce(left_join)
    
    seur[[1]]@meta.data <- only_m
    col_select <- names(seur[[1]]@meta.data)[grepl("res.", names(seur[[1]]@meta.data))]
    cols <- unique(gsub("[0-9].*", "", col_select))
    p <- clustree(seur[[1]], prefix = cols, scale_node_text = TRUE) +
      theme(text = element_text(family = "noto-sans-jp", size = 18),
            legend.text = element_text(family = "noto-sans-jp", size = 18))
    print(p)
  })
  
  output$red <- renderUI({
    req(seurat_only())
    selectInput(
      inputId = "red", 
      label = "Method", 
      choices = c("PCA" = "pca", "tSNE" = "tsne", "uMAP" = "umap"),
      selected = "PCA"
    )
  })
  
  output$select_res_factor <- renderUI({
    req(input$overall_res, d$resType, d$res, seurat_only())
    if(d$resType != "seurat_res") return()
    if(grepl("^Comp:", x = d$res)) return()
    
    res <- gsub("RNA_snn_", "", d$res)
    myChoices <- colnames(seurat_only()[[res]]@meta.data)
    myChoices <- myChoices[!(myChoices %in% c("nCount_RNA","nFeature_RNA"))]
    myChoices <- myChoices[!(apply(seurat_only()[[res]]@meta.data[, myChoices, drop = F], 2, is.numeric))]
    if("seurat_clusters" %in% myChoices) {
      myChoices <- c("seurat_clusters", myChoices[myChoices != "seurat_clusters"])
    }
    mySelected <- myChoices[1]
    selectInput(
      inputId = "select_res_factor", 
      label = "Grouping factor",
      choices = myChoices, 
      selected = mySelected)
  })
  
  # SC-DGE-OR - return umap plot by resolution
  output$res_umap <- renderPlot({
    req(input$overall_res, input$select_res_factor, d$resType, d$res, seurat_only())
    if(d$resType != "seurat_res") return()
    if(grepl("^Comp:", x = d$res)) return()

    seur <- seurat_only()
    res <- gsub("RNA_snn_", "", d$res)
    seur <- seur[[res]]
    if(is.null(input$select_res_factor) || !(input$select_res_factor %in% colnames(seur@meta.data))) return()
    xLab <- "PC 1"
    yLab <- "PC 2"
    if(input$red == "umap") {
      xLab <- "UMAP 1"
      yLab <- "UMAP 2"
    }
    if(input$red == "tsne") {
      xLab <- "tSNE 1"
      yLab <- "tSNE 2"
    }
    if(input$red == "pca") {
      xLab <- "PC 1"
      yLab <- "PC 2"
    }
    red <- input$red
    DimPlot(seur, reduction = red, group.by = input$select_res_factor) + 
      labs(title = "") +
      xlab(xLab) + ylab(yLab) +
      theme(
        plot.title = element_text(hjust = 0.5, family = "noto-sans-jp",  size = 18),
        text = element_text(family = "noto-sans-jp",  size = 18),
        legend.text = element_text(family = "noto-sans-jp", size = 18),
        aspect.ratio = 1
      )
  })
  
  ######### OVERVIEW TAB CODE: SCRNA-SEQ#####################
  
  # SC-DGE-OVER - calc the length of clusters in seurat obj
  numClust <- reactive({
    # for specified seurat resolutions
    if (d$resType == "seurat_res") {
      sapply(1:length(d$SCV), function(x)
        length(levels(Clusters(d$SCV[[x]]))))
      # based on input metadata file
    } else {
      length(unique(d$SCV[[1]]@Clusters))
    }
  })

  # SC-DGE-OVER - calc the total number of clusters for ea
  # resolution
  clustList <- reactive({
    temp <- names(d$SCV)
    temp <- unlist(sapply(1:length(temp), function(x)
      paste0(temp[x], ": ", numClust()[x]," clusters")))
    return(temp)
  })

  # SC-DGE-OVER - drop-down menu for the type of boxplots per cluster
  output$deType <- renderUI({
    temp_types <- list(
      "# positive DEGs per cluster to nearest cluster" = "DEneighb",
      "# positive DEGs per cluster to all other clusters" = "DEmarker",
      "Silhouette widths" = "silWidth"
    )
    selectInput(
      "deType", label = "Cluster separation metric",
      choices = temp_types,
      width = "100%")
  })


  # SC-DGE-OVER - update lasso selection tool for inter-cluster
  # DE boxplots
  observeEvent(input$clustSep_dblclick, {
    brush <- input$clustSep_brush
    if (!is.null(brush)) {
      d$x <- c(brush$xmin, brush$xmax)
      d$y <- c(brush$ymin, brush$ymax)
    } else {
      d$x <- NULL
      d$y <- NULL
    }
  })

  # SC-DGE-OVER - header for cluster separation boxplot
  output$clustSepHeader <- renderUI({
    req(d$SCV)
    h4("Cluster separation boxplot")
  })

  # SC-DGE-OVER - create DE boxplots by seurat resolution
  output$clustSep <- renderPlot({
    req(d$seurat_sc, d$res, input$deType, input$FDRthresh1)
    if(is.null(d$seurat_sc[[d$res]])) return()
    
    print(
      plot_clustSep(
        sCVdL = d$seurat_sc,
        DEtype = input$deType,
        FDRthresh = input$FDRthresh1,
        res = d$res,
        Xlim = d$x,
        Ylim = d$y,
        size = 1.25
      )
    )
  })

  # SC-DGE-OVER - specify FDR threshold for boxplots
  output$FDRthresh1 <- renderUI({
    req(input$deType)
    if(!(input$deType %in% c("DEneighb","DEmarker"))) return()
    numericInput(
      inputId = "FDRthresh1", 
      label = "FDR",
      value = 0.05, 
      min = 0,
      max = 1,
      step = 0.001
    )
  })
  
  observeEvent(input$FDRthresh1, {
    if(is.finite(input$FDRthresh1) && input$FDRthresh1 >= 0 && input$FDRthresh1 <= 1) return()
    updateNumericInput(session = session, inputId = "FDRthresh1", value = 0.05)
  })

  # SC-DGE-OVER - create silhouette plots by seurat resolution
  output$sil <- renderPlot({
    req(d$res, d$seurat_sc[[d$res]])
    if(grepl("^Comp", x = d$res)) return(NULL)
    plot_sil(d$seurat_sc[[d$res]], size = 1.25)
  })

  # SC-DGE-OVER - download handler for saving boxplots by png
  output$clustSepSave <- downloadHandler(
    filename = "Boxplots.png",
    content = function(file) {
      if(d$resType == "seurat_res") {
        grDevices::png(file, height = 600, units = "px", res = 150)
        print(plot_clustSep(sCVdL = d$seurat_sc,
                            DEtype = input$deType,
                            FDRthresh = input$FDRthresh1,
                            res = names(d$seurat_sc)[grepl(gsub("\\:.*", "", input$overall_res), names(d$seurat_sc))],
                            Xlim = d$x,
                            Ylim = d$y,
                            size = 0.75))
        grDevices::dev.off()
      } else {
        grDevices::png(file, height = 600, units = "px", res = 150)
        print(plot_clustSep(sCVdL = d$seurat_sc,
                            DEtype = input$deType,
                            FDRthresh = input$FDRthresh1,
                            res = "Clust",
                            Xlim = d$x,
                            Ylim = d$y,
                            size = 0.75))
        grDevices::dev.off()
      }
    }
  )

  # SC-DGE-OVER - download handler for silhouette plot by png
  output$silSave <- downloadHandler(
    filename = "Silhouette_plots.png",
    content=function(file) {
      if(length(d$SCV) > 1) {
        grDevices::png(file, height = 600, units = "px", res = 150)
        sc <- d$seurat_sc[grepl(gsub("\\:.*", "", input$overall_res), names(d$seurat_sc))]
        print(plot_sil(sc[[1]], size = 0.75))
        grDevices::dev.off()
      } else {
        grDevices::png(file, height = 600, units = "px", res = 150)
        print(plot_sil(d$seurat_sc[[1]], size = 0.75))
        grDevices::dev.off()
      }
    }
  )
  
  # SC-DGE-OVER - create metadata from seurat obj 
  # & generate silhouette data
  meta <- eventReactive(SubmitData$goqc, {
    req(SubmitData$data_type)
    if(is.null(SubmitData$goqc)) return()
    if(SubmitData$goqc == 0) return()
    if (SubmitData$data_type %in% c("Bulk", "RASL")) return()
    withProgress(message = "Building metadata file...", value = 0, {
      incProgress(1/3)
      # for several seurat resolution types
      if(d$resType == "seurat_res") {
        
        seur_obj <- lapply(seurat_only(), function(x) {
          y <- x@meta.data
          y$Sample <- rownames(y)
          y$seurat_clusters <- NULL
          return(y)
        })
        nam <- paste(names(d$seurat_sc), "_silWidth", sep = "")
        sil <- lapply(1:length(d$seurat_sc), function(x) {
          y <- d$seurat_sc[[x]]@Silhouette
          y <- data.frame("sil" = y[, 3])
          names(y) <- nam[x]
          return(y)
        })
        seur_obj <- lapply(1:length(seur_obj), function(x)
          cbind(seur_obj[[x]], sil[[x]])
        )
        cols <- lapply(1:length(seur_obj), function(x) {
          resColumn <- paste0("RNA_snn_", names(seurat_only())[x])
          silColumn <- paste0(resColumn, "_silWidth")
          names(seur_obj[[x]])[names(seur_obj[[x]]) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", resColumn, silColumn, "Sample", "sample")]
        })
          
        seur_obj <- lapply(1:length(seur_obj), function(x) {
          y <- seur_obj[[x]][, names(seur_obj[[x]]) %in% cols[[x]]]
          resColumn <- paste0("RNA_snn_", names(seurat_only())[x])
          silColumn <- paste0(resColumn, "_silWidth")
          grp <- names(y)[names(y) %in% c(resColumn, silColumn, "orig.ident", "Sample", "sample")]
          # grp <- grp[!grepl("silWidth", grp)]
          add <- y %>%
            dplyr::group_by_at(vars(one_of(grp))) %>%
            dplyr::summarise(`Number of Cells` = n())
          names(add)[ncol(add)] <- paste(resColumn, "_number of cells", sep = "")
          y <- y %>%
            full_join(., add)
        })
        seur_obj <- seur_obj %>% reduce(full_join)
        if("nCount_RNA" %in% colnames(seur_obj)) {
          seur_obj[, "nCount_RNA"] <- as.numeric(seur_obj[, "nCount_RNA"])
        }
        if("nFeature_RNA" %in% colnames(seur_obj)) {
          seur_obj[, "nFeature_RNA"] <- as.numeric(seur_obj[, "nFeature_RNA"])
        }
        if(any(grepl("silWidth", x = colnames(seur_obj)))) {
          for(i in grep("silWidth", x = colnames(seur_obj))) {
            seur_obj[, i] = as.numeric(seur_obj[, i])
          }
        }
        return(seur_obj)
        # for selected metadata file
      } else {
        seur_obj <- seurat_only()@meta.data
        seur_obj$Sample <- rownames(seur_obj) 
        cols <- names(seur_obj)[grepl("res|Cluster|orig.ident|seurat", names(seur_obj))]
        nam <- paste(names(d$seurat_sc), "_silWidth", sep = "")
        sil <- d$seurat_sc[[1]]@Silhouette
        sil <- data.frame("sil" = sil[, 3])
        names(sil) <- nam
        seur_obj <- cbind(seur_obj, sil)
        add <- seur_obj %>%
          dplyr::group_by_at(vars(one_of(cols))) %>%
          dplyr::summarise(`Number of Cells` = n())
        seur_obj <- seur_obj %>%
          full_join(., add)
        if("nCount_RNA" %in% colnames(seur_obj)) {
          seur_obj[, "nCount_RNA"] <- as.numeric(seur_obj[, "nCount_RNA"])
        }
        if("nFeature_RNA" %in% colnames(seur_obj)) {
          seur_obj[, "nFeature_RNA"] <- as.numeric(seur_obj[, "nFeature_RNA"])
        }
        if(any(grepl("silWidth", x = colnames(seur_obj)))) {
          for(i in grep("silWidth", x = colnames(seur_obj))) {
            seur_obj[, i] = as.numeric(seur_obj[, i])
          }
        }
        return(seur_obj)
      }
    })
  })

  # SC-DGE-OVER - create x input data for scatterplot
  # using metadata
  output$mdScatterX <- renderUI({
    req(d$meta, d$resValues)
    if(is.null(d$meta)) return()
    
    # resValues <- paste0("RNA_snn_res.", seq(from = 0.4, to = 2.8, by = 0.4))
    resValues <- paste0("RNA_snn_res.", d$resValues)
    silValues <- paste0(resValues, "_silWidth")
    namesExclude <- c("Sample", "sample", "Number of Cells", resValues, silValues)
    myChoices <- colnames(d$meta)[!(colnames(d$meta) %in% namesExclude)]
    myChoices <- grep(pattern = "_number.of.cells$", x = myChoices, value = T, invert = T)
    if(length(myChoices) == 0) return()
    mySelected <- myChoices[1]
    
    selectInput(
      inputId = "mdScatterX",
      label = "X axis", 
      choices = myChoices, 
      selected = mySelected
    )
  })

  # SC-DGE-OVER - create y input data for scatterplot
  # using metadata
  output$mdScatterY <- renderUI({
    req(d$meta, input$mdScatterX, d$resValues)
    if(is.null(d$meta)) return()
    
    # resValues <- paste0("RNA_snn_res.", seq(from = 0.4, to = 2.8, by = 0.4))
    resValues <- paste0("RNA_snn_res.", d$resValues)
    silValues <- paste0(resValues, "_silWidth")
    namesExclude <- c("Sample", "sample", "Number of Cells", resValues, 
                      silValues, input$mdScatterX)
    myChoices <- colnames(d$meta)[!(colnames(d$meta) %in% namesExclude)]
    myChoices <- grep(pattern = "_number.of.cells$", x = myChoices, value = T, invert = T)
    if(length(myChoices) == 0) return()
    mySelected <- myChoices[1]
    
    selectInput(
      inputId = "mdScatterY",
      label = "Y axis", 
      choices = myChoices, 
      selected = mySelected
    )
  })
  
  # SC-DGE-OVER - option to log transform x and y axes
  output$scatterLogX <- renderUI({
    req(d$meta, input$mdScatterX, input$mdScatterY)
    if(!(all(c(input$mdScatterX, input$mdScatterY) %in% colnames(d$meta)))) return()
    if(!is.numeric(d$meta[[input$mdScatterX]]) || any(d$meta[[input$mdScatterX]] <= 0)) return()
    div(
      style = "width: 80px; margin-right: 10px;",
      prettyCheckbox(
        inputId = "scatterLogX",
        label = "Log x axis", 
        value = F,
        status = "default",
        icon = icon("check")
      )
    )
  })
  
  output$scatterLogY <- renderUI({
    req(d$meta, input$mdScatterX, input$mdScatterY)
    if(!(all(c(input$mdScatterX, input$mdScatterY) %in% colnames(d$meta)))) return()
    if(!is.numeric(d$meta[[input$mdScatterY]]) || any(d$meta[[input$mdScatterY]] <= 0)) return()
    div(
      style = "width: 80px;",
      prettyCheckbox(
        inputId = "scatterLogY",
        label = "Log y axis", 
        value = F,
        status = "default",
        icon = icon("check")
      )
    )
  })
  
  # SC-DGE-OVER - create input data for barplot
  # using metadata
  output$mdFactorData <- renderUI({
    req(d$meta, d$resValues)
    # resValues <- paste0("RNA_snn_res.", seq(from = 0.4, to = 2.8, by = 0.4))
    resValues <- paste0("RNA_snn_res.", d$resValues)
    silValues <- paste0(resValues, "_silWidth")
    namesExclude <- c("Sample", "sample", "Cluster", "seurat_clusters", 
                      "Number of Cells", resValues, silValues)
    myChoices <- colnames(d$meta)[!(colnames(d$meta) %in% namesExclude)]
    myChoices <- grep(pattern = "_number.of.cells$", x = myChoices, value = T, invert = T)
    if(length(myChoices) == 0) return()
    mySelected <- myChoices[1]

    selectInput(
      inputId = "mdFactorData",
      label = "Metadata",
      choices = myChoices,
      selected = mySelected
    )
  })

  # SC-DGE-OVER - specify absolute/relative for metadata barplot
  output$mdFactorOptsF <- renderUI({
    req(d$meta, input$mdFactorData)
    if(!(input$mdFactorData %in% colnames(d$meta))) return()
    if(is.numeric(d$meta[[input$mdFactorData]])) return()
    div(
      style = "width: 200px; margin-left: 10px;",
      awesomeRadio(
        inputId = "mdFactorOptsF", 
        label = "Factor counts per cluster", 
        choices = c(Absolute = "absolute", Relative = "relative"),
        selected = "absolute",
        status = "default"
      )
    )
  })

  # SC-DGE-OVER - specify log scale y for metadata boxplot
  output$mdFactorOptsN <- renderUI({
    req(d$meta, input$mdFactorData)
    if(!(input$mdFactorData %in% colnames(d$meta))) return()
    if(!is.numeric(d$meta[[input$mdFactorData]])) return()
    if(grepl(pattern = "_silWidth$", x = input$mdFactorData)) return()
    if(input$mdFactorData == "Number of Cells") return()
    if(any(d$meta[[input$mdFactorData]] <= 0)) return()
    div(
      style = "width: 200px; margin-left: 10px;",
      prettyCheckbox(
        inputId = "mdFactorOptsN",
        label = "Log scale", 
        value = F,
        status = "default",
        icon = icon("check")
      )
    )
  })
  
  # SC-DGE-OVER - create scatterplot from metadata
  # based on selected x & y inputs
  output$mdScatter <- renderPlot({
    req(d$meta, input$mdScatterX, input$mdScatterY)
    if(!(all(c(input$mdScatterX, input$mdScatterY) %in% colnames(d$meta)))) return()
    
    logX <- logY <- F
    if(!is.null(input$scatterLogX)) logX <- input$scatterLogX
    if(!is.null(input$scatterLogY)) logY <- input$scatterLogY
    
    xRotate <- F
    if(!is.numeric(d$meta[[input$mdScatterX]]) &&
       max(nchar(as.character(unique(d$meta[[input$mdScatterX]])))) > 2) {
      xRotate <- T
    }

    plot_mdCompare_ggplot(
      df = d$meta, x = input$mdScatterX, y = input$mdScatterY, logX = logX, 
      logY = logY, xRotate = xRotate
    )
  })

  # SC-DGE-OVER - create barplot based on selected
  # metadata inputs
  output$mdFactor <- renderPlot({
    req(input$mdFactorData, d$meta)
    if(is.null(d$meta) || !all(c("seurat_clusters", input$mdFactorData) %in% colnames(d$meta))) return()
    
    logY <- F
    if(!is.null(input$mdFactorOptsN)) logY <- input$mdFactorOptsN
    
    plotType <- "absolute"
    if(!is.null(input$mdFactorOptsF)) plotType <- input$mdFactorOptsF
    
    if(is.numeric(d$meta[[input$mdFactorData]])) {
      plot_mdCompare_ggplot_box(
        df = d$meta, x = "seurat_clusters", y = input$mdFactorData, 
        logY = logY, color = T
      )
    } else {
      plot_mdCompare_ggplot_stackedbar(
        df = d$meta, x = "seurat_clusters", y = input$mdFactorData,
        plotType = plotType
      )
    }
  })

  # SC-DGE-OVER - download handler for scatterplot
  # using metadata
  output$mdScatterSave <- downloadHandler(
    filename=function() {
      paste0(gsub("^X|[_.]","",input$mdScatterX),"_vs_",
             gsub("^X|[_.]","",input$mdScatterY), ".png")
    },
    content = function(file) {
      logX <- logY <- F
      if(!is.null(input$scatterLogX)) logX <- input$scatterLogX
      if(!is.null(input$scatterLogY)) logY <- input$scatterLogY
      
      xRotate <- F
      if(!is.numeric(d$meta[[input$mdScatterX]]) &&
         max(nchar(as.character(unique(d$meta[[input$mdScatterX]])))) > 2) {
        xRotate <- T
      }
      
      p <- plot_mdCompare_ggplot(
        df = d$meta, x = input$mdScatterX, y = input$mdScatterY, 
        logX = logX, logY = logY, sizeFactor = 2, xRotate = xRotate
      )
      ggsave(file, plot = p, width = 5, height = 5)    
    }
  )

  # SC-DGE-OVER - create download handler for barplot
  # based on metadata
  output$mdFactorSave <- downloadHandler(
    filename=function() {
      if(input$mdFactorData == "seurat_clusters") {
        "Num_cells_PerCluster.png"
      } else {
        paste0(gsub("^X|[_.]", "", input$mdFactorData), "_PerCluster.png")
      }
    },
    content = function(file) {
      logY <- F
      if(!is.null(input$mdFactorOptsN)) logY <- input$mdFactorOptsN
      
      plotType <- "absolute"
      if(!is.null(input$mdFactorOptsF)) plotType <- input$mdFactorOptsF
      
      if(is.numeric(d$meta[[input$mdFactorData]])) {
        p <- plot_mdCompare_ggplot_box(
          df = d$meta, x = "seurat_clusters", y = input$mdFactorData, 
          logY = logY, color = T, sizeFactor = 2
        )
      } else {
        p <- plot_mdCompare_ggplot_stackedbar(
          df = d$meta, x = "seurat_clusters", y = input$mdFactorData,
          plotType = plotType, sizeFactor = 2
        )
      }      
      ggsave(file, plot = p, width = 5, height = 5)  
    }
  )

  ######### GENE EXPRESSION TAB CODE: SCRNA-SEQ#####################

  # SC-DGE-GE - select gene from seurat_sc
  output$cgSelect <- renderUI({
    req(d$seurat_sc, d$res)
    if (is.null(d$seurat_sc)) {
      # for multiple seurat resolutions
    } else if (length(d$seurat_sc) > 1 & !is.null(d$res)) {
      sc <- d$seurat_sc
      num <- sapply(1:length(sc), function(x)
        length(levels(Clusters(sc[[x]]))))
      temp <- names(sc)
      temp <- unlist(sapply(1:length(temp), function(x)
        paste0(temp[x], ": ", num[x]," clusters")))
      names(sc) <- temp
      if(sum(names(sc) %in% input$overall_res) == 0) return()
      sc <- sc[names(sc) %in% input$overall_res]
      temp_choices <- rownames(sc[[1]]@ClustGeneStats[[1]])
      if (is.null(names(temp_choices))) {
        temp_choices <- sort(temp_choices)
      } else {
        temp_choices <- temp_choices[order(names(temp_choices))]
      }
      # based on metadata
    } else {
      sc <- d$seurat_sc
      temp <- names(sc)
      num <- length(levels(Clusters(sc[[1]])))
      temp <- paste0(temp, ": ", num," clusters")
      names(sc) <- temp
      temp_choices <- rownames(sc[[1]]@ClustGeneStats[[1]])
      if (is.null(names(temp_choices))) {
        temp_choices <- sort(temp_choices)
      } else {
        temp_choices <- temp_choices[order(names(temp_choices))]
      }
    }
    selectizeInput("cgGene", choices = temp_choices, label = "Gene",
                   options = list(maxOptions = 1000000))
  })

  # SC-DGE-GE - header for gene expression by cluster
  output$geneexpressionbycluster_header <- renderUI({
    h3("Gene expression by cluster")
  })

  # SC-DGE-GE - runs hclust on output of scClustViz's DEdist function
  # in reactive so it is only run once per session.
  GEBoxplotHClust <- reactive({
    req(d$seurat_sc, d$resType)
    ind <- 1
    if(d$resType == "seurat_res") {
      if(is.null(d$res)) return()
      ind <- d$res
    }
    hclust(d = stats::as.dist(m = DEdist(d$seurat_sc[[ind]])), method = "single")
  })
  
  # SC-DGE-GE - returns detection rates for genes in all clusters. 
  # Returns a dataframe with clusters as column numbers and genes in rows.
  DR <- reactive({
    req(d$seurat_sc, d$resType)
    ind <- 1
    if(d$resType == "seurat_res") {
      if(is.null(d$res)) return()
      ind <- d$res
    }
    statsList <- ClustGeneStats(d$seurat_sc[[ind]])
    df <- as.data.frame(sapply(statsList, function(x) x$DR), row.names = rownames(statsList[[1]]))
    return(df)
  })
  
  # Include jitter option
  output$geneExpIncludeJitter <- renderUI({
    prettyCheckbox(
      inputId = "geneExpIncludeJitter",
      label = "Include jitter",
      value = T,
      status = "default",
      icon = icon("check")
    )
  })
  
  # Include detection rate option
  output$geneExpIncludeDR <- renderUI({
    prettyCheckbox(
      inputId = "geneExpIncludeDR",
      label = "Include detection rate",
      value = T,
      status = "default",
      icon = icon("check")
    )
  })
  
  # SC-DGE-GE - boxplots for gene expression comparison
  output$geneTest <- renderPlot({
    req(input$cgGene, GEBoxplotHClust(), d$seurat_sc, d$inD, d$resType)
    if(is.null(input$geneExpIncludeJitter) || is.null(input$geneExpIncludeDR)) return()
    if(!(input$cgGene %in% rownames(d$inD))) return()
    ind <- 1
    if(d$resType == "seurat_res") {
      if(is.null(d$res)) return()
      ind <- d$res
    }

    drData <- NULL
    if(input$geneExpIncludeDR) {
      if(is.null(DR())) return()
      drData <- data.frame(Cluster = colnames(DR()), DR = as.numeric(DR()[input$cgGene, ]))
    }

    plot_GEboxplot_ggplot(
      seuratData = d$inD, 
      scvData = d$seurat_sc[[ind]], 
      gene = input$cgGene,
      hcData = GEBoxplotHClust(),
      plotJitter = input$geneExpIncludeJitter,
      drData = drData,
      sizeFactor = 1
    )
  })
  
  # SC-DGE-GE -download button for gene expression boxplot
   output$geneTestSaveButton <- renderUI({
     downloadButton("geneTestSave", label = "Save as PNG")
  })

  # SC-DGE-GE - download handler for gene expression boxplots
  output$geneTestSave <- downloadHandler(
    filename=function() {
      paste0("Boxplot_", input$cgGene, ".png")
    },
    content=function(file) {
      ind <- 1
      if(d$resType == "seurat_res") {
        if(is.null(d$res)) return()
        ind <- d$res
      }
      
      drData <- NULL
      if(input$geneExpIncludeDR) {
        if(is.null(DR())) return()
        drData <- data.frame(Cluster = colnames(DR()), DR = as.numeric(DR()[input$cgGene, ]))
      }
      
      p <- plot_GEboxplot_ggplot(
        seuratData = d$inD, 
        scvData = d$seurat_sc[[ind]], 
        gene = input$cgGene,
        hcData = GEBoxplotHClust(),
        plotJitter = input$geneExpIncludeJitter,
        drData = drData,
        sizeFactor = 2
      )
      
      ggsave(file, plot = p, width = 5, height = 5)
    }
  )

  # SC-DGE-GE - header for distr of goi
  output$celldistributionHeader <- renderUI({
    h3("Cell distribution of genes of interest")
  })

  # SC-DGE-GE - goi embedding type selection
  output$GOI_EmbType <- renderUI({
    temp_embs <- hasEmb(d$inD)
    temp_embs <- temp_embs[
      sapply(temp_embs,function(X) ncol(getEmb(d$inD, X))) >= 2 &
        sapply(temp_embs,function(X) nrow(getEmb(d$inD, X))) == nrow(getMD(d$inD))
      ]
    temp_embs <- toupper(temp_embs)
    temp_embs <- gsub("TSNE", "tSNE", temp_embs)
    selectInput("GOI_EmbType", label = "Cell embedding",
                choices = temp_embs,
                selected = temp_embs[temp_embs %in% c("tSNE","UMAP")][1])
  })

  # SC-DGE-GE - goi x-axis
  output$GOI_EmbDimX <- renderUI({
    req(input$GOI_EmbType)
    selectInput("GOI_EmbDimX", label = "x-axis",
                choices = gsub("_", " ", colnames(getEmb(d$inD, input$GOI_EmbType))),
                selected = gsub("_", " ", colnames(getEmb(d$inD, input$GOI_EmbType))[1]))
  })

  # SC-DGE-GE - goi y-axis
  output$GOI_EmbDimY <- renderUI({
    req(input$GOI_EmbType)
    selectInput("GOI_EmbDimY", label = "y-axis",
                choices = gsub("_", " ", colnames(getEmb(d$inD, input$GOI_EmbType))),
                selected = gsub("_", " ", colnames(getEmb(d$inD, input$GOI_EmbType))[2]))
  })

  # SC-DGE-GE - goi select gene via name
  output$GOI1select <- renderUI({
    req(d$inD_orig)
    selectizeInput("goi1", label = "Gene", choices = rownames(d$inD_orig),
                   selected = rownames(d$inD_orig)[1], multiple = F,
                   options = list(maxOptions = 1000000))

  })

  # SC-DGE-GE - select type of dimensional reduction plot
  output$MD_EmbType <- renderUI({
    temp_embs <- hasEmb(seurat_only())
    temp_embs <- temp_embs[
      sapply(temp_embs,function(X) ncol(getEmb(seurat_only(), X))) >= 2 &
        sapply(temp_embs,function(X) nrow(getEmb(seurat_only(), X))) == nrow(getMD(seurat_only()))
      ]
    selectInput("MD_EmbType",label="Cell embedding",
                choices=toupper(temp_embs),
                selected=toupper(temp_embs)[toupper(temp_embs) %in% c("TSNE","UMAP")][1])
  })

  # SC-DGE-GE - dim red x-axis number
  output$MD_EmbDimX <- renderUI({
    selectInput("MD_EmbDimX", label = "x-axis",
                choices = gsub("_", " ", colnames(getEmb(seurat_only(), input$MD_EmbType))),
                selected = gsub("_", " ", colnames(getEmb(seurat_only(), input$MD_EmbType))[1]))
  })

  # SC-DGE-GE - dim red y-axis number
  output$MD_EmbDimY <- renderUI({
    selectInput("MD_EmbDimY", label = "y-axis",
                choices = gsub("_", " ", colnames(getEmb(seurat_only(), input$MD_EmbType))),
                selected = gsub("_", " ", colnames(getEmb(seurat_only(), input$MD_EmbType))[2]))
  })

  # SC-DGE-GE - specify clusters or goi overlay on plot
  output$plotClust1 <- renderUI({
    awesomeRadio("plotClust1",inline=F,label="Plot",selected="goi",
                 choices=list("Gene expression overlay"="goi", "Clusters"="clust"),
                 status = "success")
  })

  # SC-DGE-GE - overlay cluster labels on plot
  output$plotLabel1 <- renderUI({
    prettyCheckbox("plotLabel1", label = "Include cluster labels (style as above)", value = T, status = "default", icon = icon("check"))
  })

  # SC-DGE-GE - action button to plot gene exp by cluster or goi
  output$GOI1go <- renderUI({
    actionButton("GOI1go", "Search", icon = icon("search"))
  })

  # SC-DGE-GE - specify tSNE labels
  output$tsneLabels <- renderUI({
    if (length(d$res) > 0) {
      temp_choices <- list("Cluster annotations"="ClusterNames",
                           "Cluster annotations (label all)"="ClusterNamesAll",
                           "Cluster numbers"="Clusters")
      if (all(unique(attr(Clusters(d$SCV[[d$res]]),"ClusterNames")) == "")) {
        temp_choices <- temp_choices[3]
      } else if (grepl("^Comp:",d$res)) {
        temp_choices <- temp_choices[-1]
        names(temp_choices)[1] <- "Cluster annotations"
      }
      awesomeRadio("tsneLabels", "Labels:", inline = T, choices = temp_choices, status = "success")
    }
  })

  # SC-DGE-GE - tSNE plot for given cluster or goi + error msgs
  output$goiPlot1 <- renderPlot({
    req(d$inD, input$goi1, input$GOI_EmbType, input$plotClust1, 
        input$GOI_EmbDimX, input$GOI_EmbDimY)
    if(is.null(input$plotLabel1)) return()
    
    if(!(input$goi1 %in% rownames(d$inD))) return()
    
    xDim <- as.numeric(strsplit(input$GOI_EmbDimX, split = " ")[[1]][2])
    yDim <- as.numeric(strsplit(input$GOI_EmbDimY, split = " ")[[1]][2])

    plot_goi_ggplot(
      seuratData = d$inD, gene = input$goi1, method = tolower(input$GOI_EmbType),
      clusterLabels = input$plotLabel1, plotType = input$plotClust1, 
      sizeFactor = 1, dims = c(xDim, yDim)
    )
  })
  
  # SC-DGE-GE - download button for tSNE plot based on gene name
  output$goiPlot1SaveButton <- renderUI({
    req(input$plotClust1)
    if(is.null(input$plotClust1)) return()
    downloadButton("goiPlot1Save", label = "Save as PNG")
  })

  # NEED TO ADJUST TEXT SIZE IN DOWNLOAD - TOO SMALL
  # SC-DGE-GE - download tSNE file based on gene name
  output$goiPlot1Save <- downloadHandler(
    filename = function() {
      if(is.null(input$goi1) | input$plotClust1 == "clust") {
        paste0(input$GOI_EmbType, ".png")
      } else {
        paste0(input$GOI_EmbType, "_", input$goi1, ".png")
      }
    },
    content=function(file) {
      xDim <- as.numeric(strsplit(input$GOI_EmbDimX, split = " ")[[1]][2])
      yDim <- as.numeric(strsplit(input$GOI_EmbDimY, split = " ")[[1]][2])
      
      p <- plot_goi_ggplot(
        seuratData = d$inD, gene = input$goi1, method = tolower(input$GOI_EmbType),
        clusterLabels = input$plotLabel1, plotType = input$plotClust1, 
        sizeFactor = 2, dims = c(xDim, yDim)
      )
      
      ggsave(file, plot = p, width = 7, height = 7)
    }
  )

  ######### DGE ANALYSIS TAB CODE: SCRNA-SEQ#####################

  # SC-DGE-DGE - data input for dge analyis
  sc_mytab <- reactive({
    if(is.null(d$SCV) || is.null(d$res) || is.null(input$sc_contrasts) || 
       is.null(d$SCV[[d$res]]) || is.null(input$sc_gene_list_filtering)) return()
    if(!(input$sc_contrasts %in% names(d$SCV[[d$res]]@DEcombn))) return()
    df <- d$SCV[[d$res]]@DEcombn[names(d$SCV[[d$res]]@DEcombn) %in% input$sc_contrasts][[1]]
    if(nrow(df) == 0) return()
    df <- df[order(1/abs(df$logGER), df$pVal, decreasing = F), , drop = F]
    if(input$sc_gene_list_filtering == 2) df <- df[df$logGER > 0, , drop = F]
    if(input$sc_gene_list_filtering == 3) df <- df[df$logGER < 0, , drop = F]
    if(nrow(df) == 0) return()
    return(df)
  })

  # SC-DGE-DGE - select type of dotplot comparison
  output$heatDEtype <- renderUI({
    temp <- list(
      "DE vs rest" = "DEvsRest",
      
      # WHY IS THIS HASHED OUT??? PLEASE UPDATE.
      # "Marker genes" = "DEmarker",
      "DE vs neighbor" = "DEneighb"
    )
    awesomeRadio("dotplotDEtype",
                 "Dotplot genes",
                 choices = temp,
                 selected = "DEvsRest",
                 status = "success")
  })

  # SC-DGE-DGE - set false discovery rate
  output$FDRthresh2 <- renderUI({
    numericInput(
      inputId = "FDRthresh2", 
      label = "FDR",
      value = 0.05, 
      min = 0,
      max = 1,
      step = 0.001
    )
  })
  
  observeEvent(input$FDRthresh2, {
    if(is.finite(input$FDRthresh2) && input$FDRthresh2 >= 0 && input$FDRthresh2 <= 1) return()
    updateNumericInput(session = session, inputId = "FDRthresh2", value = 0.05)
  })
  

  # SC-DGE-DGE - create a slider to specify the # of genes to show in dotplot
  # based on min and max of dge
  output$DEgeneSlider <- renderUI({
    req(input$dotplotDEtype, input$DEclustNum, d$res, d$SCV)
    if(!(d$res %in% names(d$SCV))) return()
    if(!(input$DEclustNum %in% Clusters(d$SCV[[d$res]]))) return()
    if(is.null(DEgenes())) return()
    
    numericInput(
      inputId = "DEgeneCount", 
      label = "# genes per cluster", 
      value = 5,
      min = 1, 
      max = max(sapply(DEgenes(), length)),
      step = 1
    )
  })
  
  observeEvent(input$DEgeneCount, {
    if(is.finite(input$DEgeneCount) && !is.null(DEgenes()) &&
       input$DEgeneCount >= 1 && input$DEgeneCount <= max(sapply(DEgenes(), length))) return()
    updateNumericInput(session = session, inputId = "DEgeneCount", value = 5)
  })

  # SC-DGE-DGE - select cluster to analyze for dge
  output$DEclustNum <- renderUI({
    req(d$resType, d$SCV)
    ind <- 1
    if(d$resType == "seurat_res") {
      if(is.null(d$res) || !(d$res %in% names(d$SCV))) return()
      ind <- d$res
    }
    myChoices <- mixedsort(levels(Clusters(d$SCV[[ind]])))
    selectInput(inputId = "DEclustNum", label = "Cluster", choices = myChoices)
  })

  # SC-DGE-DGE - text file for gene summary stats for selected cluster
  output$CGSsave0 <- downloadHandler(
    filename = function() {
      paste0("Cluster_",
             input$DEclustNum,"_stats.txt")
    },
    content = function(file) {
      if (d$resType == "seurat_res") {
        res <- input$overall_res
        res <- gsub(".*_|\\:.*", "", res)
        res <- d$SCV[grepl(res, names(d$SCV))]
        outTable <- ClustGeneStats(res[[1]])[[input$DEclustNum]][, c("MGE","DR","MDGE")]
        write.table(outTable, file,
                    quote = F, sep = "\t",
                    row.names = T, col.names = NA)
      } else {
        outTable <- ClustGeneStats(d$SCV[[1]])[[input$DEclustNum]][, c("MGE","DR","MDGE")]
        write.table(outTable, file,
                    quote = F, sep = "\t",
                    row.names = T, col.names = NA)
      }
    }
  )

  # SC-DGE-DGE - download handler for text file for scClustViz analysis
  # for selected cluster
  output$deGeneSave <- downloadHandler(
    filename = function() {
      paste0(input$dotplotDEtype, "_", input$DEclustNum, ".txt")
    },
    content = function(file) {
      if (d$resType == "seurat_res") {
        res <- input$overall_res
        res <- gsub(".*_|\\:.*", "", res)
        res <- d$SCV[grepl(res, names(d$SCV))]
        outTable <- switch(
          EXPR = input$dotplotDEtype,
          DEvsRest = DEvsRest(res[[1]])[[input$DEclustNum]],
          DEmarker = DEmarker(res[[1]])[[input$DEclustNum]],
          DEneighb = DEneighb(res[[1]])[[input$DEclustNum]]
        )
        write.table(outTable, file,
                    quote = F, sep = "\t",
                    row.names = T, col.names = NA)
      } else {
        outTable <- switch(
          EXPR = input$dotplotDEtype,
          DEvsRest = DEvsRest(d$SCV[[1]])[[input$DEclustNum]],
          DEmarker = DEmarker(d$SCV[[1]])[[input$DEclustNum]],
          DEneighb = DEneighb(d$SCV[[1]])[[input$DEclustNum]]
        )
        write.table(outTable, file,
                    quote = F, sep = "\t",
                    row.names = T, col.names = NA)
      }
    }
  )

  # SC-DGE-DGE - reactive exp for de analysis for
  # selected cluster, de type and fdr
  DEgenes <- reactive({
    req(input$dotplotDEtype, input$FDRthresh2, d$resType, d$SCV, d$res)
    if(d$resType == "seurat_res") {
      if(is.null(d$res)) return()
      if(!(d$res %in% names(d$SCV))) return()
      dotplotDEgenes(sCVd = d$SCV[[d$res]],
                     DEtype = input$dotplotDEtype,
                     FDRthresh = input$FDRthresh2)
    } else {
      dotplotDEgenes(sCVd = d$SCV[[1]],
                     DEtype = input$dotplotDEtype,
                     FDRthresh = input$FDRthresh2)
    }
  })

  # SC-DGE-DGE - download handler for dge dotplot
  output$heatmapSave <- renderUI({
    req(d$seurat_sc, d$resType, DEgenes(), input$DEgeneCount, input$DEclustNum)
    if(is.null(DEgenes()) || length(DEgenes()[[as.character(input$DEclustNum)]]) == 0) return()
    if(d$resType == "seurat_res" && is.null(d$res)) return()
    downloadButton("heatmapSaveButton", label = "Save as PNG")
  })
  
  output$heatmapSaveButton <- downloadHandler(
    filename = function() {
      paste0("Dotplot_", input$dotplotDEtype, ".png")
    },
    content = function(file) {
      ind <- 1
      if(d$resType == "seurat_res") {
        if(is.null(d$res)) return()
        ind <- d$res
      }
      
      nGenes <- min(input$DEgeneCount, 69)
      
      p <- plot_deDotplot_ggplot(
        scvData = d$seurat_sc[[ind]],
        deGenes = DEgenes(),
        nGenes = nGenes,
        clust = input$DEclustNum,
        sizeFactor = 2
      )
      
      ggsave(file, plot = p, width = 13, height = 6)
    }
  )

  # SC-DGE-DGE - header for dge dotplot
  output$dotplotHeader <- renderUI({
    h4("Dotplot")
  })

  # SC-DGE-DGE - dotplot based on dge genes for selected cluster
  output$dotplotWarning <- renderUI({
    req(DEgenes(), input$DEclustNum)
    if(!is.null(DEgenes()) && length(DEgenes()[[as.character(input$DEclustNum)]]) > 0) return()
    HTML("No differentially expressed genes based on the chosen parameters.")
  })
  
  output$dotplot <- renderPlot({
    req(d$seurat_sc, d$resType, DEgenes(), input$DEgeneCount, input$DEclustNum)
    if(is.null(DEgenes()) || length(DEgenes()[[as.character(input$DEclustNum)]]) == 0) return()
    ind <- 1
    if(d$resType == "seurat_res") {
      if(is.null(d$res)) return()
      ind <- d$res
    }
    
    nGenes <- min(input$DEgeneCount, 69)

    plot_deDotplot_ggplot(
      scvData = d$seurat_sc[[ind]],
      deGenes = DEgenes(),
      nGenes = nGenes,
      clust = input$DEclustNum
    )
  })
  
  # SC-DGE-DGE - dge table header
  output$dge_cluster_sc_header <- renderUI({
    req(d$SCV)
    h3("DGE table")
  })

  # SC-DGE-DGE - header for dge table filters
  output$dge_cluster_sc_filterheader <- renderUI({
    req(d$SCV)
    h4("Filters")
  })

  # SC-DGE-DGE - abs log2FC threshold for table
  output$sc_dgeByCluster_table_lfc <- renderUI({
    numericInput("sc_dgeByCluster_table_lfc", label = "Abs. log2 fold-change",
                 value = 0.25, min = 0)
  })

  # SC-DGE-DGE - adj p-val threshold for table
  output$sc_dgeByCluster_table_padj <- renderUI({
    numericInput("sc_dgeByCluster_table_padj", label = "Adj. p-value",
                 value = 0.05, min = 0, max = 1)
  })

  # SC-DGE-DGE - action button to load dge table
  output$sc_dgeByCluster_table_load <- renderUI({
    req(input$DEclustNum, d$resType, d$res, d$SCV, input$dotplotDEtype,
        seurat_only())
    actionButton("sc_dgeByCluster_table_load", label = "Load table")
  })

  # SC-DGE-DGE - dge by cluster table
  output$dge_cluster_sc <- renderUI({
    req(sc_dge_cluster_table())
    DT::dataTableOutput("dge_cluster_sc_tbl")
  })

  # SC-DGE-DGE - reactive expression for dge cluster table
  sc_dge_cluster_table <- eventReactive(input$sc_dgeByCluster_table_load, {
    withProgress(message = "Creating DGE table...", value = 0, {
      incProgress(1/3)
      if (d$resType == "seurat_res") {
        res <- gsub(":.*$", replacement = "", x = input$overall_res)
        if(!(res %in% names(d$SCV))) return()
        clus_table <- ClustGeneStats(d$SCV[[res]])[[as.character(input$DEclustNum)]][, c("MGE","DR","MDGE")]
        if(is.null(clus_table)) return()
        out_table <- switch(
          EXPR = input$dotplotDEtype,
          DEvsRest = DEvsRest(d$SCV[[res]])[[as.character(input$DEclustNum)]],
          DEmarker = DEmarker(d$SCV[[res]])[[as.character(input$DEclustNum)]],
          DEneighb = DEneighb(d$SCV[[res]])[[as.character(input$DEclustNum)]]
        )
        if(is.null(out_table)) return()
        res <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
        if(!(res %in% names(seurat_only()))) return()
        seur <- seurat_only()[[res]]
        res <- gsub(":.*$", replacement = "", x = input$overall_res)
        Idents(seur) <- res
        avg <- data.frame(AverageExpression(seur))
        avg <- avg %>%
          rownames_to_column(., var = "Gene")
      } else {
        clus_table <- ClustGeneStats(d$SCV[[1]])[[input$DEclustNum]][, c("MGE","DR","MDGE")]
        if(is.null(clus_table)) return()
        out_table <- switch(
          EXPR = input$dotplotDEtype,
          DEvsRest = DEvsRest(d$SCV[[1]])[[as.character(input$DEclustNum)]],
          DEmarker = DEmarker(d$SCV[[1]])[[as.character(input$DEclustNum)]],
          DEneighb = DEneighb(d$SCV[[1]])[[as.character(input$DEclustNum)]]
        )
        if(is.null(out_table)) return()
        seur <- seurat_only()
        Idents(seur) <- "seurat_clusters"
        avg <- data.frame(AverageExpression(seur))
        avg <- avg %>%
          rownames_to_column(., var = "Gene")
      }
      incProgress(1/3)
      out_table <- out_table %>%
        rownames_to_column(., var = "Gene")
      clus_table <- clus_table %>%
        rownames_to_column(., var = "Gene")
      allDF <- out_table %>%
        left_join(., clus_table, by = "Gene") %>%
        left_join(., avg, by = "Gene") %>%
        mutate_if(., is.numeric, round, 3)
      names(allDF) <- gsub("logGER", "logFC", names(allDF))
      names(allDF) <- gsub("Wstat", "Wilcoxon statistic", names(allDF))
      names(allDF) <- gsub("pVal", "P-value", names(allDF))
      names(allDF) <- gsub("dDR", "Difference in detection rate", names(allDF))
      names(allDF) <- gsub("RNA\\.", "Avg exp cluster ", names(allDF))
      names(allDF) <- gsub("Gene", "GeneId", names(allDF))
      nCells <- ncol(seur)
      if (input$dotplotDEtype == "DEvsRest") {
        nam <- names(allDF)[grepl("Avg exp cluster ", names(allDF))]
        nam <- nam[!grepl(input$DEclustNum, nam)]
        allDF <- allDF %>%
          dplyr::select(-nam)

        rest <- seur@assays$RNA@data[, Idents(seur) != input$DEclustNum]
        nam <- rownames(rest)
        meanRest <- unlist(mclapply(1:nrow(rest), function(x)
          meanLogX(rest[x], ncell = ncol(seur), ex = exp(1)), mc.cores = 12))
        meanRest <- exp(1)^meanRest
        names(meanRest) <- nam
        allDF <- data.frame(GeneId = allDF$GeneId,
                            round(exp(1)^allDF$MGE, digits = 3),
                            Mean_Rest = round(meanRest[allDF$GeneId], digits = 3),
                            p_val_adj = allDF$FDR,
                            log2FC = round(log2(exp(1)^allDF$logFC), digits = 3))
        allDF <- allDF[as.numeric(allDF$p_val_adj) <= input$sc_dgeByCluster_table_padj &
                         abs(as.numeric(allDF$log2FC)) > input$sc_dgeByCluster_table_lfc, ]
        colnames(allDF)[2] <- paste0("Mean_",input$DEclustNum)
        allDF <- allDF[order(as.numeric(allDF$log2FC), decreasing = T), ]
      } else if (input$dotplotDEtype == "DEneighb") {
        nam <- names(allDF)[grepl("Avg exp cluster ", names(allDF))]
        sel <- paste(unlist(strsplit(gsub("logFC_", "", names(allDF)[2]), "-")), collapse = "|")
        nam <- nam[!grepl(sel, nam)]
        allDF <- allDF %>%
          dplyr::select(-nam)
        neighborName <- gsub("^.*\\|(.*)$", replacement = "\\1", x = sel)
        neigh <- seur@assays$RNA@data[, Idents(seur) == neighborName]
        nam <- rownames(neigh)
        meanNeighbor <- apply(seur@assays$RNA@data[, Idents(seur) == neighborName], 1,
                              function(x) meanLogX(x, ncell = nCells, ex = exp(1)))
        meanNeighbor <- exp(1)^meanNeighbor
        names(meanNeighbor) <- nam
        allDF <- data.frame(GeneId = allDF$GeneId,
                            round(exp(1)^allDF$MGE, digits = 3),
                            round(meanNeighbor[allDF$GeneId], digits = 3),
                            p_val_adj = allDF[, grep("FDR_", x = colnames(allDF))],
                            log2FC = round(log2(exp(1)^allDF[, grep("logFC_", x = colnames(allDF))]), digits = 3))
        allDF <- allDF[as.numeric(allDF$p_val_adj) <= input$sc_dgeByCluster_table_padj &
                         abs(as.numeric(allDF$log2FC)) > input$sc_dgeByCluster_table_lfc, ]
        colnames(allDF)[c(2,3)] <- c(paste0("Mean_", input$DEclustNum), paste0("Mean_", neighborName))
        allDF <- allDF[order(as.numeric(allDF$log2FC), decreasing = T), ]
      }
      return(allDF)
    })
  })

  # SC-DGE-DGE - dge table
  output$dge_cluster_sc_tbl <- DT::renderDataTable(server = FALSE, {
    req(sc_dge_cluster_table())
    df <- sc_dge_cluster_table()
    colnames(df)[colnames(df) == "GeneId"] = "Gene"
    colnames(df)[colnames(df) == "p_val_adj"] = "P-value (adj)"
    colnames(df)[colnames(df) == "log2FC"] = "Log2 fold-change"
    DT::datatable({
      df
    },
    rownames = F,
    extensions = 'Buttons',
    options = list(
      autoWidth = TRUE,
      columnDefs = list(list(className = "dt-left", width = "100px", targets = "_all")),
      dom = 'Bfrtip',
      buttons = list(list(extend = 'csv',
                          filename = paste(input$dotplotDEtype, "_", input$DEclustNum, sep = ""),
                          text = 'Download DGE by cluster')),
      language = list(zeroRecords = "No genes passing thresholds",
                      emptyTable = "")
    ), class = "display")
  })

  ######### Custom DGE tab: SCRNA-SEQ DATA ONLY ############

  # SC-DGE-CUST - warning message if there are factors w/o at
  # least two groups with 3 or more cells
  output$sc_dge_byfactor_novars_msg <- renderText({
    req(d$sc_dge_byfactor_novars_msg)
    return(d$sc_dge_byfactor_novars_msg)
  })

  # SC-DGE-CUST - select factor from metadata
  output$sc_dge_factor_fact <- renderUI({
    req(d$inD)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    
    myChoices <- colnames(d$inD@meta.data)
    myChoices <- myChoices[!(apply(d$inD@meta.data, 2, function(j) is.numeric(j) | length(unique(j)) == 1))]
    myChoices <- myChoices[!(myChoices %in% c("nCount_RNA", "nFeature_RNA"))]
    myChoices <- myChoices[sapply(myChoices, function(myChoice) sum(as.vector(table(d$inD@meta.data[, myChoice])) >= 3) >= 2)]
    
    if(length(myChoices) == 0) return()

    selectInput(
      inputId = "sc_dge_factor_fact",
      label = "Factor",
      choices = myChoices
    )
  })

  # SC-DGE-CUST - select factor group 1 for dge
  output$sc_dge_factor_group1 <- renderUI({
    req(d$inD, input$sc_dge_factor_fact)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(!(input$sc_dge_factor_fact %in% colnames(d$inD@meta.data))) return()

    myChoices <- mixedsort(unique(as.character(d$inD@meta.data[, input$sc_dge_factor_fact])))
    myChoices <- myChoices[sapply(myChoices, function(myChoice) sum(d$inD@meta.data[, input$sc_dge_factor_fact] == myChoice) >= 3)]

    if(length(myChoices) < 2) return()
    myChoices <- c("All (pairwise)" = "All", myChoices)
    selectInput("sc_dge_factor_group1", label = "Group 1", choices = myChoices)
  })

  # SC-DGE-CUST - select factor group 2 for dge
  output$sc_dge_factor_group2 <- renderUI({
    req(d$inD, input$sc_dge_factor_fact, input$sc_dge_factor_group1)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(!(input$sc_dge_factor_fact %in% colnames(d$inD@meta.data))) return()
    
    myChoices <- mixedsort(unique(d$inD@meta.data[, input$sc_dge_factor_fact]))
    
    myChoices <- isolate(as.character(myChoices[myChoices != input$sc_dge_factor_group1]))
    myChoices <- c("Rest", myChoices)
    selectInput("sc_dge_factor_group2", label = "Group 2", choices = myChoices, selected = "Rest")
  })

  # SC-DGE-CUST - view optional parameters
  output$sc_dge_factor_viewOpts <- renderUI({
    req(input$sc_dge_factor_fact, input$sc_dge_factor_group1,
        input$sc_dge_factor_group2)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    prettyCheckbox("sc_dge_factor_viewOpts", label = "View options", value = F, status = "default", icon = icon("check"))
  })

  output$sc_dge_factor_opts <- renderUI({
    req(input$sc_dge_factor_fact, input$sc_dge_factor_group1,
        input$sc_dge_factor_group2)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(is.null(input$sc_dge_factor_viewOpts) || !(input$sc_dge_factor_viewOpts)) return()
    tagList(
      div(
        style = "display: inline-block; width: 150px;", 
        numericInput(
          inputId = "sc_dge_factor_lfcThreshold", 
          label = "LFC threshold",
          min = 0, 
          step = 0.01,
          value = isolate(d$sc_dge_factor_lfcThreshold)
        )
      ),
      div(
        style = "display: inline-block; width: 150px; margin-left: 20px; vertical-align: top;", 
        selectInput(
          inputId = "sc_dge_factor_test", 
          label = "Test",
          choices = c(
            "wilcox (default)" = "wilcox",
            "bimod",
            "roc",
            "t",
            "negbinom",
            "poisson",
            "LR",
            "MAST",
            "DESeq2"
          ),
          selected = isolate(d$sc_dge_factor_test)
        )
      ),
      div(
        style = "display: inline-block; width: 150px; margin-left: 20px;", 
        numericInput(
          inputId = "sc_dge_factor_minPct", 
          label = "Min. pct",
          min = 0, 
          max = 1, 
          step = 0.1,
          value = isolate(d$sc_dge_factor_minPct)
        )
      ),
      div(
        style = "display: inline-block; width: 150px; margin-left: 20px;", 
        numericInput(
          inputId = "sc_dge_factor_minDiffPct", 
          label = "Min. diff pct",
          min = -1, 
          max = 1, 
          step = 0.1,
          value = isolate(d$sc_dge_factor_minDiffPct)
        )
      ),
      div(
        style = "display: inline-block; margin-left: 20px;",
        actionButton(inputId = "sc_dge_factor_resetDefault", label = "Reset defaults")
      )
    )
  })
  
  observe({
    if(is.null(input$sc_dge_factor_test)) return()
    d$sc_dge_factor_test <- input$sc_dge_factor_test
  })

  observe({
    if(is.null(input$sc_dge_factor_lfcThreshold)) return()
    d$sc_dge_factor_lfcThreshold <- input$sc_dge_factor_lfcThreshold
  })

  observe({
    if(is.null(input$sc_dge_factor_minPct)) return()
    d$sc_dge_factor_minPct <- input$sc_dge_factor_minPct
  })

  observe({
    if(is.null(input$sc_dge_factor_minDiffPct)) return()
    d$sc_dge_factor_minDiffPct <- input$sc_dge_factor_minDiffPct
  })

  # SC-DGE-CUST - update lfc thresh, factor, min pct and min diff pct
  observeEvent(input$sc_dge_factor_resetDefault, {
    updateNumericInput(session, inputId = "sc_dge_factor_lfcThreshold", value = 0.25)
    updateSelectInput(session, inputId = "sc_dge_factor_test", selected = "wilcox")
    updateNumericInput(session, inputId = "sc_dge_factor_minPct", value = 0.1)
    updateNumericInput(session, inputId = "sc_dge_factor_minDiffPct", value = -1)
  })

  # SC-DGE-CUST - action button to submit parameters for dge analysis
  output$sc_dge_factor_submit <- renderUI({
    req(input$sc_dge_factor_fact, input$sc_dge_factor_group1,
        input$sc_dge_factor_group2)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    actionButton("sc_dge_factor_submit", label = "Submit")
  })

  # SC-DGE-CUST - action button to clear results
  output$sc_dge_factor_clear <- renderUI({
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(length(d$sc_dge_byfactor) == 0) return()
    actionButton("sc_dge_factor_clear", label = "Clear results")
  })

  # SC-DGE-CUST - update object d factor for dge
  observeEvent(input$sc_dge_factor_clear, {
    d$sc_dge_byfactor <- NULL
  })

  # SC-DGE-CUST - warning if comparison has already been run
  output$sc_dge_byfactor_msg <- renderText({
    req(d$sc_dge_byfactor_msg)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    return(d$sc_dge_byfactor_msg)
  })

  # SC-DGE-CUST - based on inputs run custom dge analysis
  observeEvent(input$sc_dge_factor_submit, {
    withProgress(message = "Running DGE...", value = 0, {
      incProgress(1/3)
      
      d$sc_dge_byfactor_msg <- NULL      
      myResult <- SeuratDGE(
        seuratData = d$inD,
        fact = input$sc_dge_factor_fact,
        group1 = input$sc_dge_factor_group1,
        group2 = input$sc_dge_factor_group2,
        logfc.threshold = d$sc_dge_factor_lfcThreshold,
        test.use = d$sc_dge_factor_test,
        min.pct = d$sc_dge_factor_minPct,
        min.diff.pct = d$sc_dge_factor_minDiffPct
      )
      incProgress(1/3)
      myList <- d$sc_dge_byfactor[[input$sc_dge_factor_fact]]
      myList <- myList[!(names(myList) %in% names(myResult))]
      myList <- c(myResult, myList)
      d$sc_dge_byfactor[[input$sc_dge_factor_fact]] <- myList
      incProgress(1/3)
    })
  })

  # SC-DGE-CUST - header for custom dge results
  output$sc_dge_factor_results_header <- renderUI({
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(length(d$sc_dge_byfactor) == 0) return()
    h3("Differential gene expression results", style = "margin-top: 0px;")
  })

  # SC-DGE-CUST - select factor for custom dge table
  output$sc_dge_factor_tablefact <- renderUI({
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(length(d$sc_dge_byfactor) == 0) return()
    selectInput(
      inputId = "sc_dge_factor_tablefact",
      label = "Factor",
      choices = names(d$sc_dge_byfactor)
    )
  })

  # SC-DGE-CUST - select comparison for custom dge table
  output$sc_dge_factor_comp <- renderUI({
    req(input$sc_dge_factor_tablefact)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(length(d$sc_dge_byfactor) == 0) return()
    if(is.null(d$sc_dge_byfactor[[input$sc_dge_factor_tablefact]])) return()
    selectInput(
      inputId = "sc_dge_factor_comp",
      label = "Comparison",
      choices = names(d$sc_dge_byfactor[[input$sc_dge_factor_tablefact]])
    )
  })

  # SC-DGE-CUST - dge table output
  output$sc_dge_factor_results <- renderUI({
    req(input$sc_dge_factor_tablefact, input$sc_dge_factor_comp)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(length(d$sc_dge_byfactor) == 0) return()
    if(is.null(d$sc_dge_byfactor[[input$sc_dge_factor_tablefact]][[input$sc_dge_factor_comp]])) return()
    DT::dataTableOutput("sc_dge_factor_results_tbl")
  })

  # SC-DGE-CUST - dge table
  output$sc_dge_factor_results_tbl <- DT::renderDataTable(server = FALSE, {
    req(input$sc_dge_factor_tablefact, input$sc_dge_factor_comp)
    if(!is.null(d$sc_dge_byfactor_novars_msg)) return()
    if(length(d$sc_dge_byfactor) == 0) return()
    if(is.null(d$sc_dge_byfactor[[input$sc_dge_factor_tablefact]][[input$sc_dge_factor_comp]])) return()

    DT::datatable({
      d$sc_dge_byfactor[[input$sc_dge_factor_tablefact]][[input$sc_dge_factor_comp]] %>%
          mutate_if(., is.numeric, formatC)
    },
    rownames = F,
    extensions = 'Buttons',
    options = list(
      autoWidth = TRUE,
      columnDefs = list(list(className = "dt-left", width = "100px", targets = "_all")),
      dom = 'Bfrtip',
      buttons = list(list(extend = 'csv',
                          filename = paste0("DGE_", input$sc_dge_factor_tablefact,
                                            "_", input$sc_dge_factor_comp),
                          text = 'Download table')),
      language = list(zeroRecords = "No differentially expressed genes for this comparison.")
    ), class = "display")
  })

  ######### Volcano Plots TAB CODE: SCRNA-SEQ#####################
  
  # SC-DGE-VOL - Reactive function that returns  clusters from d$SCV[[d$res]]
  SCVClusters <- reactive({
    req(d$SCV, d$res)
    resTemp = d$res
    if(!(resTemp %in% names(d$SCV))) resTemp = names(d$SCV)[1]
    Clusters(d$SCV[[resTemp]])
  })

  # DGE - header for volcano plot
  output$setScatterHeader <- renderUI({
    h3("Volcano plot")
  })

  # SC-DGE-VOL - scatterplot comparing 2 clusters
  output$setScatter <- renderPlot({
    req(d$SCV, d$res, input$ssA, input$ssB)
    if(!(d$res %in% names(d$SCV)) || input$ssA == input$ssB ||
       !(input$ssA %in% levels(Clusters(d$SCV[[d$res]]))) ||
       !(input$ssB %in% levels(Clusters(d$SCV[[d$res]])))) return()

    plot_compareClusts2(
      sCVd = d$SCV[[d$res]],
      clA = input$ssA,
      clB = input$ssB,
      dataType = input$scatterInput,
      labType = input$diffLabelType,
      labTypeDiff = input$diffLabelChoice,
      labNum = input$diffCount,
      labGenes = rownames(d$SCV[[d$res]])
    )
  })

  # SC-DGE-VOL - volcano dge table
  output$vol_cluster_scHeader <- renderUI({
    req(d$SCV)
    h3("DGE table")
  })

  # SC-DGE-VOL - header for table filters
  output$vol_cluster_sc_filterheader <- renderUI({
    req(d$SCV)
    h4("Filters")
  })

  # SC-DGE-VOL - abs log2fc threshold for table
  output$vol_cluster_sc_lfc <- renderUI({
    numericInput("vol_cluster_sc_lfc", label = "Abs. log2 fold-change",
                 value = 0.25, min = 0)
  })

  # SC-DGE-VOL - adj p-val threshold for table
  output$vol_cluster_sc_padj <- renderUI({
    numericInput("vol_cluster_sc_padj", label = "Adj. p-value",
                 value = 0.05, min = 0, max = 1)
  })

  # SC-DGE-VOL - action button to load dge table
  output$vol_cluster_sc_load <- renderUI({
    req(d$SCV, input$ssA, input$ssB, d$res, DEcombn(d$SCV[[d$res]]),
        seurat_only(),
        ClustGeneStats(d$SCV[[d$res]])[[input$ssA]][,c("MGE","DR","MDGE")],
        ClustGeneStats(d$SCV[[d$res]])[[input$ssB]][,c("MGE","DR","MDGE")])
    if(!isTruthy(d$res)) return(NULL)
    if(!(d$res %in% names(d$SCV))) return(NULL)
    if(!(input$ssA %in% names(ClustGeneStats(d$SCV[[d$res]])))) return(NULL)
    if(!(input$ssB %in% names(ClustGeneStats(d$SCV[[d$res]])))) return(NULL)
    actionButton("vol_cluster_sc_load", label = "Load table")
  })

  # SC-DGE-VOL - reactive exp to return dge table
  vol_cluster_sc_table <- eventReactive(input$vol_cluster_sc_load, {
    a <- ClustGeneStats(d$SCV[[d$res]])[[input$ssA]][,c("MGE","DR","MDGE")]  %>%
      rownames_to_column(., var = "Gene")
    names(a)[2:4] <- c(paste("MGE_", input$ssA, sep = ""),
                       paste("DR_", input$ssA, sep = ""),
                       paste("MDGE_", input$ssA, sep = ""))
    b <- ClustGeneStats(d$SCV[[d$res]])[[input$ssB]][,c("MGE","DR","MDGE")] %>%
      rownames_to_column(., var = "Gene")
    names(b)[2:4] <- c(paste("MGE_", input$ssB, sep = ""),
                       paste("DR_", input$ssB, sep = ""),
                       paste("MDGE_", input$ssB, sep = ""))

    temp <- c(paste(input$ssA, input$ssB, sep = "-"), paste(input$ssB, input$ssA, sep = "-"))
    tempName <- temp[temp %in% names(DEcombn(d$SCV[[d$res]]))]
    if(!isTruthy(tempName)) return(NULL)
    if(!(tempName %in% names(DEcombn(d$SCV[[d$res]])))) return(NULL)

    out_table <- DEcombn(d$SCV[[d$res]])[[tempName]] %>%
      rownames_to_column(., var = "Gene")
    if (d$resType == "seurat_res") {
      seur <- seurat_only()[[1]]
      res <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
      if(!isTruthy(res)) return(NULL)
      if(sum(grepl(res, colnames(seur@meta.data)) == F) == 0 ) return(NULL)
      Idents(seur) <- res
      avg <- data.frame(AverageExpression(seur))
      avg <- avg %>%
        rownames_to_column(., var = "Gene")
    } else {
      seur <- seurat_only()
      Idents(seur) <- "Cluster"
      avg <- data.frame(AverageExpression(seur))
      avg <- avg %>%
        rownames_to_column(., var = "Gene")
    }

    all <- out_table %>%
      left_join(., a, by = "Gene") %>%
      left_join(., b, by = "Gene")  %>%
      left_join(., avg, by = "Gene") %>%
      mutate_if(., is.numeric, round, 3)
    names(all) <- gsub("logGER", "logFC", names(all))
    names(all) <- gsub("Wstat", "Wilcoxon statistic", names(all))
    names(all) <- gsub("pVal", "P-value", names(all))
    names(all) <- gsub("dDR", "Difference in detection rate", names(all))
    names(all) <- gsub("RNA\\.", "Avg exp cluster ", names(all))
    names(all) <- gsub("Gene", "GeneId", names(all))

    nam <- names(all)[grepl("Avg exp cluster ", names(all))]
    nam <- nam[!grepl(input$ssA, nam)]
    nam <- nam[!grepl(input$ssB, nam)]
    meanCol1 <- paste0("MGE_", input$ssA)
    meanCol2 <- paste0("MGE_", input$ssB)
    all <- all %>%
      dplyr::select(-nam)
    all <- all[, c("GeneId", meanCol1, meanCol2, "FDR","logFC")]
    colnames(all) <- c("GeneId", paste0("Mean_", input$ssA), paste0("Mean_", input$ssB),
                       "p_val_adj", "logFC")
    all <- data.frame(GeneId = all$GeneId,
                      round(exp(1)^all[, paste0("Mean_",input$ssA)], digits = 3),
                      round(exp(1)^all[, paste0("Mean_",input$ssB)], digits = 3),
                      p_val_adj = all$p_val_adj,
                      log2FC = round(log2(exp(1)^all$logFC), digits = 3))
    all <- all[as.numeric(all$p_val_adj) <= input$vol_cluster_sc_padj &
                 abs(as.numeric(all$log2FC)) > input$vol_cluster_sc_lfc, ]
    colnames(all)[c(2,3)] <- c(paste0("Mean_",input$ssA), paste0("Mean_",input$ssB))
    all <- all[order(as.numeric(all$log2FC), decreasing = T), ]
    return(all)
  })

  # SC-DGE-VOL - dge table output
  output$vol_cluster_sc <- renderUI({
    req(vol_cluster_sc_table())
    DT::dataTableOutput("vol_cluster_sc_tbl")
  })

  # SC-DGE-VOL - dge table
  output$vol_cluster_sc_tbl <- DT::renderDataTable(server = FALSE, {
    df <- vol_cluster_sc_table()
    colnames(df)[colnames(df) == "GeneId"] = "Gene"
    colnames(df)[colnames(df) == "p_val_adj"] = "P-value (adj)"
    colnames(df)[colnames(df) == "log2FC"] = "Log2 fold-change"
    DT::datatable({
      # vol_cluster_sc_table()
      df
    },
    rownames = F,
    extensions = 'Buttons',
    options = list(
      autoWidth = TRUE,
      columnDefs = list(list(className = "dt-left", width = "100px", targets = "_all")),
      dom = 'Bfrtip',
      buttons = list(list(extend = 'csv',
                          filename = paste(input$ssA, "_vs_", input$ssB, sep = ""),
                          text = 'Download DGE of selected clusters')),
      language = list(zeroRecords = "No genes passing thresholds")
    ), class = "display")
  })

  # SC-DGE-VOL - select cluster A
  output$setScatterA <- renderUI({
    req(SCVClusters(), d$res)
    if(length(d$res) == 0 || !(d$res %in% names(d$SCV))) return()
    temp_cl <- levels(SCVClusters())
    temp_sel <- temp_cl[1]
    selectInput(
      inputId = "ssA",
      label = "Cluster A",
      choices = temp_cl,
      selected = temp_sel
    )
  })

  # SC-DGE-VOL - select cluster B
  output$setScatterB <- renderUI({
    req(SCVClusters(), d$res, input$ssA)
    if(length(d$res) == 0 || !(d$res %in% names(d$SCV))) return()
    temp_cl <- levels(SCVClusters())
    temp_cl <- temp_cl[!(temp_cl == input$ssA)]
    temp_sel <- temp_cl[1]
    selectInput(
      inputId = "ssB",
      label = "Cluster B",
      choices = temp_cl,
      selected = temp_sel
    )
  })

  # SC-DGE-VOL - select type of volcano plot
  output$scatterInput <- renderUI({
    selectInput(
      inputId = "scatterInput", 
      label = "Plot type",
      choices = c(
        "Gene expression ratio" = "logGER",
        "Detection rate difference" = "dDR",
        "Gene expression difference" = "MGE"
      ),
      selected = "logGER"
    )
  })

  # SC-DGE-VOL - select type of volcano plot to display
  output$diffLabelType <- renderUI({
    selectInput(
      inputId = "diffLabelType", 
      label = "Label",
      choices = c("Largest fold-changes" = "diff", "Smallest p-values" = "de")
    )
  })

  # SC-DGE-VOL - upate cluster B based on selected cluster A
  observeEvent(input$ssA, {
    temp_cl <- levels(SCVClusters())
    temp_cl = temp_cl[!(temp_cl == input$ssA)]
    temp_sel = temp_cl[1]
    updateSelectInput(session, "ssB",
                      choices = temp_cl, selected = temp_sel)
  })

  # SC-DGE-VOL - download cluster gene stats A
  output$CGSsaveA <- downloadHandler(
    filename=function() { paste0("GeneStatsCluster_",input$ssA,".txt") },
    content=function(file) {
      outTable <- ClustGeneStats(d$SCV[[d$res]])[[input$ssA]][,c("MGE","DR","MDGE")]
      write.table(outTable,file,quote=F,sep="\t",row.names=T,col.names=NA)
    }
  )

  # SC-DGE-VOL - download cluster gene stats B
  output$CGSsaveB <- downloadHandler(
    filename=function() { paste0("GeneStatsCluster_",input$ssB,".txt") },
    content=function(file) {
      outTable <- ClustGeneStats(d$SCV[[d$res]])[[input$ssB]][,c("MGE","DR","MDGE")]
      write.table(outTable,file,quote=F,sep="\t",row.names=T,col.names=NA)
    }
  )

  # SC-DGE-VOL - select gene expression or detection rate
  output$diffLabelChoice <- renderUI({
    req(input$scatterInput, input$diffLabelType)
    if (input$scatterInput == "GERvDDR" & input$diffLabelType == "diff") {
      awesomeRadio("diffLabelChoice",label="Axis of difference:",
                   choices=c("Gene expression"="logGER",
                             "Detection rate"="dDR"),
                   status = "success")
    }
  })

  # SC-DGE-VOL - select number of most different genes to label
  output$diffLabelSelect <- renderUI({
    req(input$diffLabelType)
    numericInput(
      inputId = "diffCount",
      label = div(style = "width: 150px;", "# top genes to label"),
      value = 5,
      min = 1,
      max = 100,
      step = 1,
      width = "100%"
    )
  })

  # SC-DGE-VOL - download handler for volcano plot comp cluster
  # A and B
  output$setScatterSave <- downloadHandler(
    filename=function() {
      paste0(gsub("^X|[_.]", "", input$scatterInput), "_",
             gsub("^X|[_.]","",make.names(input$ssA)),
             "_vs_", gsub("^X|[_.]","",make.names(input$ssB)), ".png")
    },
    content=function(file) {
      if(length(d$res) > 0) {
        p <- plot_compareClusts2(
          sCVd = d$SCV[[d$res]],
          clA = input$ssA,
          clB = input$ssB,
          dataType = input$scatterInput,
          labType = input$diffLabelType,
          labTypeDiff = input$diffLabelChoice,
          labNum = input$diffCount,
          labGenes = GOI(),
          sizeFactor = 2
        )
        ggsave(file, plot = p, width = 6, height = 5)
      }
    }
  )

  # SC-DGE-VOL - download handler for comp cluster A and B data
  output$setComparisonSave <- downloadHandler(
    filename=function() {
      temp <- c(paste(input$ssA,input$ssB,sep="-"),paste(input$ssB,input$ssA,sep="-"))
      tempName <- temp[temp %in% names(DEcombn(d$SCV[[d$res]]))]
      return(paste0("DEcombn_",tempName,".txt"))
    },
    content=function(file) {
      temp <- c(paste(input$ssA,input$ssB,sep="-"),paste(input$ssB,input$ssA,sep="-"))
      tempName <- temp[temp %in% names(DEcombn(d$SCV[[d$res]]))]
      return(write.table(DEcombn(d$SCV[[d$res]])[[tempName]],
                         file,quote=F,sep="\t",row.names=T,col.names=NA))
    }
  )

  ######### GENE SET ENRICHMENT & HEATMAPS TAB CODE: SCRNA-SEQ#####################

  # SC-DGE-GSE - title for gse + heatmap
  output$gse_title_sc <- renderText({
    "GSE & heatmap"
  })

  # SC-DGE-GSE - select enrichr libraries
  output$enrichRLib_sc <- renderUI({
    enrichRLibs <- scan("./data/enrichrlibrary.txt", what = "character", quiet = T)
    selectInput(
      inputId = "enrichRLib_sc", 
      label = "EnrichR libraries",
      choices = enrichRLibs, 
      selected = enrichRLibs[1:6], 
      multiple = T,
      width = "600px"
    )
  })
  
  observe({
    if(is.null(input$enrichRLib_sc)) {
      d$enrichRLib_sc <- scan("./data/enrichrlibrary.txt", what = "character", quiet = T)[1:6]
    } else {
      d$enrichRLib_sc <- input$enrichRLib_sc
    }
  })
  
  # SC-DGE-GSE - contrasts for cluster comps
  output$sc_contrasts <- renderUI({
    req(d$SCV, input$dge_list_sc)
    if(input$dge_list_sc != "Cluster DGE") return()
    div(
      style = "margin-left: 6px;",
      selectInput(
        inputId = "sc_contrasts", 
        label = div(style = "width: 200px;", "Cluster contrast"),
        choices = names(d$SCV[[d$res]]@DEcombn), 
        multiple = F,
        width = "120px"
      )
    )
  })

  # SC-DGE-GSE - select input DGE filter list or user selected genes
  output$dge_list_sc <- renderUI({
    myChoices <- c("Cluster DGE", "Manual")
    if(!is.null(d$sc_dge_byfactor) && length(d$sc_dge_byfactor) > 0) {
      myChoices <- c("Cluster DGE", "Custom DGE", "Manual")
    }
    div(
      style = "width: 130px;",
      selectInput(
        inputId = "dge_list_sc",
        label = "Input gene list",
        choices = myChoices
      )
    )
  })
  
  # SC-DGE-CUST - select factor for custom dge table
  output$gse_customdge_fact <- renderUI({
    req(input$dge_list_sc)
    if(input$dge_list_sc != "Custom DGE" || !is.null(d$sc_dge_byfactor_novars_msg) ||
       length(d$sc_dge_byfactor) == 0) return()
    div(
      style = "width: 150px; margin-left: 10px;",
      selectInput(
        inputId = "gse_customdge_fact",
        label = "Factor",
        choices = names(d$sc_dge_byfactor)
      )
    )
  })
  
  # SC-DGE-CUST - select comparison for custom dge table
  output$gse_customdge_contrast <- renderUI({
    if(is.null(input$dge_list_sc) || input$dge_list_sc != "Custom DGE" ||
       !is.null(d$sc_dge_byfactor_novars_msg) || length(d$sc_dge_byfactor) == 0 ||
       is.null(input$gse_customdge_fact) || is.null(d$sc_dge_byfactor[[input$gse_customdge_fact]])) return()
    div(
      style = "width: 200px; margin-left: 10px;",
      selectInput(
        inputId = "gse_customdge_contrast",
        label = "Contrast",
        choices = names(d$sc_dge_byfactor[[input$gse_customdge_fact]])
      )
    )
  })

  # SC-DGE-GSE - copy/paste box for user selected genes
  output$gene_list_sc <- renderUI({
    req(input$dge_list_sc)
    if(input$dge_list_sc != "Manual") return()
    div(
      style = "margin-left: 30px;",
      textAreaInput(
        inputId = "gene_list_sc", 
        label = "Paste genes", 
        value = "" , 
        width = "150px"
      )
    )
  })

  # SC-DGE-GSE - select to complete gse on all, up in group 1 or group 2
  output$sc_gene_list_filtering <- renderUI({
    req(input$dge_list_sc)
    if(!(input$dge_list_sc %in% c("Cluster DGE", "Custom DGE")) || 
       (is.null(input$sc_contrasts) && is.null(input$gse_customdge_contrast))) return()
    tempChoices <- 1:3
    if(input$dge_list_sc == "Cluster DGE") {
      if(is.null(input$sc_contrasts)) return()
      myContrasts <- input$sc_contrasts
    } else {
      if(is.null(input$gse_customdge_contrast)) return()
      myContrasts <- input$gse_customdge_contrast
    }
    if(grepl(".*_VS_.*", myContrasts)) {
      group1 <- gsub("^(.*)_VS_.*", replacement = "\\1", x = myContrasts)
      group2 <- gsub(".*_VS_(.*)$", replacement = "\\1", x = myContrasts)
    } else {
      group1 <- gsub("^(.*)-.*", replacement = "\\1", myContrasts)
      group2 <- gsub("^.*-(.*)$", replacement = "\\1", myContrasts)
    }
    names(tempChoices) <- c("All DE genes", paste0("Up in ", group1), paste0("Up in ", group2))
    div(
      style = "margin-left: 10px;",
      selectInput(
        inputId = "sc_gene_list_filtering",
        label = "Gene list filtering",
        choices = tempChoices,
        width = "200px"
      )
    )
  })

  # SC-DGE-GSE - specify # of top genes to complete gse
  output$top_n_genes_sc <- renderUI({
    req(input$dge_list_sc)
    if(input$dge_list_sc == "Manual") return()
    div(
      style = "width: 150px;",
      textInput(
        inputId = "n_genes_sc", 
        label = "Top n genes", 
        value = 100, 
        placeholder = "ALL"
      )
    )
  })

  # SC-DGE-GSE - use all genes passing filter for gse
  output$use_all_genes_sc <- renderUI({
    req(input$dge_list_sc)
    if(input$dge_list_sc == "Manual") return()
    prettyCheckbox(
      inputId = "all_genes_sc", 
      label = "All genes passing filters", 
      value = F, 
      status = "default", 
      icon = icon("check")
    )
  })

  output$geneOptions_sc_ui <- renderUI({
    req(input$dge_list_sc)
    if(input$dge_list_sc == "Manual") return()
    div(
      style = "width: 150px; margin-left: 8px;", 
      fluidPage(
        fluidRow(uiOutput("top_n_genes_sc")),
        fluidRow(uiOutput("use_all_genes_sc"))
      )
    )
  })
  
  # SC-DGE-GSE - reactive exp for storing data in object d
  observeEvent(input$all_genes_sc, {
    d$fea_all_genes <- input$all_genes_sc
  })

  # SC-DGE-GSE - store updated gse data size in object d
  observe({
    if(is.null(gsegenes_sc())) {
      if(!is.null(input$dge_list_sc) && input$dge_list_sc == "Cluster DGE") {
        d$dge_info_sc <- "No DE genes."
      } else if(!is.null(input$dge_list_sc) && input$dge_list_sc == "Manual") {
        d$dge_info_sc <- "No genes found in data."
      }
    } else {
      nGenes <- length(gsegenes_sc())
      if(nGenes == 1) {
        d$dge_info_sc <- "1 gene found. At least 2 are required."
      } else {
        d$dge_info_sc <- paste0("Data size: ", nGenes, " genes")
      }
    }
  })
  
  # SC-DGE-GSE - output msgs after running gse
  output$dge_info_sc <- renderUI({
    if(!is.null(input$dge_list_sc) && input$dge_list_sc == "Manual" &&
       (is.null(input$gene_list_sc) || input$gene_list_sc == "")) return()
    h4(d$dge_info_sc, style = "margin-top: 0px; margin-right: 30px;")
  })

  # SC-DGE-GSE - action button to submit gse
  output$submit_fe_sc <- renderUI({
    actionButton("submit_fe_sc", "Run GSE", icon = icon("space-shuttle"))
  })
  
  # Show/hide submit_fe button depending on whether data are ready
  observe({
    if(is.null(gsegenes_sc()) || length(gsegenes_sc()) < 2) {
      shinyjs::hide(id = "submit_fe_sc")
    } else {
      shinyjs::show(id = "submit_fe_sc")
    }
  })
  

  # SC-DGE-GSE - specify colors for heatmap
  output$heatColors_sc <- renderUI({
    
    selectInput(
      inputId = "heatColors_sc",
      label = "Color palette",
      choices = c(
        'Default',
        'Red-white-blue',
        'Viridis',
        'Green-yellow-red'
      ),
      selected = 'Default'
    )
  })

  # SC-DGE-Heatmap - action button to build heatmap
  output$submit_list_sc <- renderUI({
    req(heattran1_sc())
    actionButton(
      inputId = "submit_list_sc", 
      label = "Build heatmap", 
      icon = icon("space-shuttle")
    )
  })

  # SC-DGE-Heatmap - cluster rows
  output$row_type_sc <- renderUI({
    prettyCheckbox("row_type_sc", label = "Cluster genes",
                  value = T, status = "default", icon = icon("check"))
  })

  # SC-DGE-GSE - cluster cols
  output$col_type_sc <- renderUI({
    prettyCheckbox("col_type_sc", label = "Cluster samples",
                  value = T, status = "default", icon = icon("check"))
  })
  
  # SC-DGE-Heatmap - scale genes option
  output$scale_genes_sc <- renderUI({
    prettyCheckbox("scale_genes_sc", label = "Scale genes",
                   value = T, status = "default", icon = icon("check"))
  })

  # SC-DGE-GSE - warning msgs for gse
  output$gse_warning_sc <- renderText(d$gse_warning_msg)

  # SC-DGE-GSE - missing genes (ie not in dataset) for user selected genes in gse
  observeEvent(input$submit_fe_sc, {
    d$newExpSetup <- F
    if(is.null(sc_mytab())) {
      d$gse_warning_msg <- "No DE genes."
    } else {
      d$gse_warning_msg <- ""
    }
    if(input$dge_list_sc == "Manual") {
      if(length(unlist(strsplit(input$gene_list_sc, split = "\\s+|,\\s?"))) == 0 |
         input$gene_list_sc == "") {
        d$gse_warning_msg <- "No genes were submitted."
      }
    }
  })

  # SC-DGE-GSE - reactive expression for gse
  sc_func_enrich <- eventReactive(input$submit_fe_sc, {
    d$gse_warning_msg <- ""
    withProgress(message = "Computing functional enrichment gene lists...", value = 0, {
      incProgress(1/4)
      isolate({
        enr <- enrichr(genes = gsegenes_sc(), databases = d$enrichRLib_sc)
        incProgress(1/4)
        enr <- enr[sapply(enr, function(x) nrow(x) > 0)]
        if(length(enr) == 0) return()
        enrCheck <- sapply(enr, function(i) {
          if(ncol(i) == 1) {
            if(colnames(i)[1] == "X.html.") return(FALSE)
          }
          return(TRUE)
        })
        if(all(!enrCheck)) {
          d$gse_warning_msg <- "EnrichR is unavailable. Please try again later."
          return()
        }
        # return significant genes
        enr <- getSigTerms(enr = enr, libs = d$enrichRLib_sc)
        incProgress(1/4)
        enr <- enr[sapply(enr, function(x) !is.null(x))]
        enr <- rbindlist(enr)
        if(length(enr) == 0) return()
        # function for removing redundant go terms
        enr <- enr[order(-score)]
        enr.go <- enr[grepl(pattern = "GO", x = libName) & grepl(pattern = "\\(GO", x = term), , drop = F]
        if(nrow(enr.go) > 0) {
          enr.go[,GOID := tstrsplit(term, "\\(GO")[2]]
          enr.go[,GOID := paste0('GO', GOID)]
          enr.go[,GOID := gsub("\\)", '', GOID)]
          enr.go[,term.short := tstrsplit(term, "\\(")[1]]
          enr.go[,term.short := trimws(term.short, which="right")]
          enr.go <- enr.go %>%
            mutate(term = paste(term.short, " (", GOID, ")", sep = "")) %>%
            dplyr::select(-c(term.short, GOID))
          enr.rest <- enr[!grepl('GO', libName),]
          enr <- rbind(enr.rest, enr.go)
          enr <- enr[order(-score)]
        }
        enr <- enr %>%
          mutate(pval = round(pval, digits = 3),
                 adjPval = round(adjPval, digits = 3),
                 Z_score = round(Z_score, digits = 3),
                 score = round(score, digits = 3))
        names(enr) <- c("Library name", "Library rank", "Gene count",
                        "Term", "Overlap", "P-value", "Adjusted p-value",
                        "Old p-value", "Old adjusted p-value", "Z-score",
                        "Score", "Gene list")
        enr <- enr %>%
          arrange(`P-value`)
        enr <- enr[, -c("Old p-value", "Old adjusted p-value"), with = FALSE]
        incProgress(1/4)
        return(enr)
      })
    })
  })

  # SC-DGE-GSE - header for gse
  output$sc_tbl_func_enrHeader <- renderUI({
    req(input$submit_fe_sc, sc_func_enrich())
    if(d$gse_warning_msg != "") return(NULL)
    if (is.null(sc_func_enrich())) {
      return(NULL)
    } else if (!is.null(sc_func_enrich())) {
      h3("Gene set enrichment table", style = "margin-top: 0px;")
    }
  })

  # SC-DGE-GSE - gse table output
  output$sc_func_enr_tbl <- DT::renderDataTable(server = FALSE, {
    if(is.null(sc_func_enrich())) return()
    DT::datatable({
      sc_func_enrich()
    }, rownames = F,
    extensions = 'Buttons',
    options = list(
      autoWidth = FALSE,
      columnDefs = list(list(className = 'dt-left', width = '40%', targets = "_all")),
      dom = 'Bfrtip',
      scrollX = TRUE,
      buttons = list(list(extend = 'csv',
                          filename = "scRNA-Seq GSE analysis",
                          text = 'Download GSE analysis'))
    ), class = "display")
  })

  # SC-DGE-GSE - warning text for specifying no sign gse
  output$sc_no_list <- renderText({
    if(d$gse_warning_msg != "") return(NULL)
    if(input$dge_list_sc == "Cluster DGE" && is.null(sc_func_enrich())) {
      "There are no significant enrichments of common annotated biological features."
    } else if(input$dge_list_sc == "Manual" && is.null(sc_func_enrich())) {
      "There are no significant enrichments of common annotated biological features for the submitted genes."
    }
  })

  # Single-cell heatmap help
  observeEvent(input$sc_heatmap_help_button, {
    showModal(modalDialog(
      title = div(
        div(style = "display: inline-block; ", strong("Heatmap")),
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
      includeMarkdown("markdown/help/sc_heatmap_help.md")
    ))
  })
  
  # SC-DGE-GSE - contrasts for cluster comps
  output$sc_contrasts_heatmap <- renderUI({
    req(d$SCV, input$dge_list_heatmap_sc)
    if(input$dge_list_heatmap_sc != "Cluster DGE") return()
    div(
      style = "margin-left: 6px;",
      selectInput(
        inputId = "sc_contrasts_heatmap", 
        label = div(style = "width: 200px;", "Cluster contrast"),
        choices = names(d$SCV[[d$res]]@DEcombn), 
        multiple = F,
        width = "120px"
      )
    )
  })
  
  # SC-DGE-Heatmap - select input DGE filter list or user selected genes
  output$dge_list_heatmap_sc <- renderUI({
    myChoices <- c("Cluster DGE", "Manual")
    if(!is.null(d$sc_dge_byfactor) && length(d$sc_dge_byfactor) > 0) {
      myChoices <- c("Cluster DGE", "Custom DGE", "Manual")
    }
    div(
      style = "width: 130px;",
      selectInput(
        inputId = "dge_list_heatmap_sc",
        label = "Input gene list",
        choices = myChoices
      )
    )
  })
  
  # SC-DGE-CUST - select factor for custom dge table
  output$heatmap_customdge_fact <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(input$dge_list_heatmap_sc != "Custom DGE" || !is.null(d$sc_dge_byfactor_novars_msg) ||
       length(d$sc_dge_byfactor) == 0) return()
    div(
      style = "width: 150px; margin-left: 10px;",
      selectInput(
        inputId = "heatmap_customdge_fact",
        label = "Factor",
        choices = names(d$sc_dge_byfactor)
      )
    )
  })
  
  # SC-DGE-CUST - select comparison for custom dge table
  output$heatmap_customdge_contrast <- renderUI({
    if(is.null(input$dge_list_heatmap_sc) || input$dge_list_heatmap_sc != "Custom DGE" ||
       !is.null(d$sc_dge_byfactor_novars_msg) || length(d$sc_dge_byfactor) == 0 ||
       is.null(input$heatmap_customdge_fact) || is.null(d$sc_dge_byfactor[[input$heatmap_customdge_fact]])) return()
    div(
      style = "width: 200px; margin-left: 10px;",
      selectInput(
        inputId = "heatmap_customdge_contrast",
        label = "Contrast",
        choices = names(d$sc_dge_byfactor[[input$heatmap_customdge_fact]])
      )
    )
  })
  
  # SC-DGE-Heatmap - copy/paste box for user selected genes
  output$gene_list_heatmap_sc <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(input$dge_list_heatmap_sc != "Manual") return()
    div(
      style = "margin-left: 30px;",
      textAreaInput(
        inputId = "gene_list_heatmap_sc", 
        label = "Paste genes", 
        value = "" , 
        width = "150px"
      )
    )
  })
  
  # SC-DGE-Heatmap - select to create heatmap from all, up in group 1 or group 2
  output$sc_gene_list_filtering_heatmap <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(!(input$dge_list_heatmap_sc %in% c("Cluster DGE", "Custom DGE")) || 
       (is.null(input$sc_contrasts_heatmap) && is.null(input$heatmap_customdge_contrast))) return()
    tempChoices <- 1:3
    if(input$dge_list_heatmap_sc == "Cluster DGE") {
      if(is.null(input$sc_contrasts_heatmap)) return()
      myContrasts <- input$sc_contrasts_heatmap
    } else {
      if(is.null(input$heatmap_customdge_contrast)) return()
      myContrasts <- input$heatmap_customdge_contrast
    }
    if(grepl(".*_VS_.*", myContrasts)) {
      group1 <- gsub("^(.*)_VS_.*", replacement = "\\1", x = myContrasts)
      group2 <- gsub(".*_VS_(.*)$", replacement = "\\1", x = myContrasts)
    } else {
      group1 <- gsub("^(.*)-.*", replacement = "\\1", myContrasts)
      group2 <- gsub("^.*-(.*)$", replacement = "\\1", myContrasts)
    }
    names(tempChoices) <- c("All DE genes", paste0("Up in ", group1), paste0("Up in ", group2))
    div(
      style = "margin-left: 10px;",
      selectInput(
        inputId = "sc_gene_list_filtering_heatmap",
        label = "Gene list filtering",
        choices = tempChoices,
        width = "200px"
      )
    )
  })
  
  # SC-DGE-Heatmap - specify # of top genes for heatmap
  output$top_n_genes_heatmap_sc <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(input$dge_list_heatmap_sc == "Manual") return()
    if(is.null(sc_mytab_heatmap())) return()
    maxGenes <- nrow(sc_mytab_heatmap())
    myValue <- min(maxGenes, 100)
    div(
      style = "width: 150px;",
      numericInput(
        inputId = "n_genes_heatmap_sc", 
        label = "Top n genes", 
        min = 1,
        max = maxGenes,
        value = myValue, 
        step = 1
      )
    )
  })
  
  observeEvent(input$n_genes_heatmap_sc, {
    if(
      !IsInteger(input$n_genes_heatmap_sc) ||
      input$n_genes_heatmap_sc < 1 || 
      input$n_genes_heatmap_sc > nrow(sc_mytab_heatmap())
    ) {
      if(!is.null(sc_mytab_heatmap()) && nrow(sc_mytab_heatmap()) > 1) {
        myValue <- min(nrow(sc_mytab_heatmap()), 100)
      } else {
        myValue <- 100
      }
      updateNumericInput(session = session, inputId = "n_genes_heatmap_sc", value = myValue)
    }
  })
  
  # SC-DGE-Heatmap - use all genes passing filter for heatmap
  output$use_all_genes_heatmap_sc <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(input$dge_list_heatmap_sc == "Manual") return()
    prettyCheckbox(
      inputId = "all_genes_heatmap_sc", 
      label = "All genes passing filters", 
      value = F, 
      status = "default", 
      icon = icon("check")
    )
  })
  
  output$geneOptions_heatmap_sc_ui <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(input$dge_list_heatmap_sc == "Manual") return()
    div(
      style = "width: 150px; margin-left: 8px;", 
      fluidPage(
        fluidRow(uiOutput("top_n_genes_heatmap_sc")),
        fluidRow(uiOutput("use_all_genes_heatmap_sc"))
      )
    )
  })
  
  # SC-DGE-GSE - use only cells in groups from contrast
  output$contrast_groups_only_sc <- renderUI({
    req(input$dge_list_heatmap_sc)
    if(!(input$dge_list_heatmap_sc %in% c("Cluster DGE", "Custom DGE"))) return()
    prettyCheckbox(
      inputId = "contrast_groups_only_sc", 
      label = "Use cells from contrast only", 
      value = F, 
      status = "default", 
      icon = icon("check")
    )
  })
  
  # SC-DGE-GSE - reactive exp for storing data in object d
  observeEvent(input$all_genes_heatmap_sc, {
    d$fea_all_genes <- input$all_genes_heatmap_sc
  })
  
  # SC-DGE-GSE - store updated gse data size in object d
  observe({
    if(is.null(heattran1_sc()[[1]])) {
      if(!is.null(input$dge_list_heatmap_sc) && input$dge_list_heatmap_sc %in% c("Cluster DGE", "Custom DGE")) {
        d$dge_info_heatmap_sc <- "No DE genes."
      } else if(!is.null(input$dge_list_heatmap_sc) && input$dge_list_heatmap_sc == "Manual") {
        d$dge_info_heatmap_sc <- "No genes found in data."
      } 
    } else {
      dims <- dim(heattran1_sc()[[1]])
      geneTxt <- " genes, "
      cellsTxt <- paste0(dims[2], " cells")
      if(dims[1] == 1) geneTxt <- " gene, "
      if(!is.null(input$col_type_sc) && input$col_type_sc && dims[2] > 5000) {
        cellsTxt <- "5000 cells (downsampled)"
      }
      d$dge_info_heatmap_sc <- paste0("Data size: ", dims[1], geneTxt, cellsTxt)
    }
  })
  
  # SC-DGE-GSE - output msgs after running gse
  output$dge_info_heatmap_sc <- renderUI({
    if(!is.null(input$dge_list_heatmap_sc) && input$dge_list_heatmap_sc == "Manual" &&
       (is.null(input$gene_list_heatmap_sc) || input$gene_list_heatmap_sc == "")) return()
    if(is.null(d$dge_info_heatmap_sc) || d$dge_info_heatmap_sc == "") return()
    h4(d$dge_info_heatmap_sc, style = "margin-top: 0px; margin-right: 30px;")
  })
  
  # SC-DGE-GSE - text header for heatmaps
  output$headheat_sc <- renderUI({
    req(heattran2_sc(), !d$newHeat,
        input$submit_list_sc, input$dge_list_heatmap_sc)
    h3("Interactive heatmap (click on cells)")
  })

  # SC-DGE-GSE - select factors for sample labels and legend for heatmap
  output$heatfactorlabel_sc <- renderUI({
    req(d$inD)
    myChoices <- colnames(d$inD@meta.data)[!(colnames(d$inD@meta.data) %in% c("nCount_RNA","nFeature_RNA"))]
    if("seurat_clusters" %in% myChoices) {
      myChoices <- c("seurat_clusters", myChoices[myChoices != "seurat_clusters"])
    }
    mySelected <- myChoices[1]
    if("Cluster" %in% myChoices) mySelected <- "Cluster"
    selectInput(
      inputId = "heatfactorlabel_sc",
      label = "Choose factor(s) for labeling",
      choices = myChoices,
      selected = mySelected,
      multiple = T
    )
  })

  # SC-DGE-GSE - color schemes for heatmaps
  heatCols_sc <-  reactive({
    if(is.null(input$heatColors_sc)) return()
    if (input$heatColors_sc == "Red-white-blue") {
      cols <- colorRampPalette(c('red', 'white', 'blue'))(100)
    } else if (input$heatColors_sc == "Viridis") {
      cols <- viridis(100)
    } else if (input$heatColors_sc == "Green-yellow-red") {
      cols <- rev(brewer.pal(n=11, name = "RdYlGn"))
    } else {
      cols <- rev(brewer.pal(n = 11, name = "RdBu"))
    }
  })
  
  # Returns genes to be used in single-cell GSE
  gsegenes_sc <- reactive({
    req(SubmitData$data_type, input$dge_list_sc, normcounts_sc())
    if(SubmitData$data_type != "Single-cell") return()
    if(input$dge_list_sc == "Manual") {
      if(is.null(input$gene_list_sc)) return()
      genes <- toupper(unlist(strsplit(input$gene_list_sc, split = "\\s+|,\\s?")))
      genes <- genes[genes %in% toupper(rownames(normcounts_sc()))]
      if(length(genes) == 0) return()
      return(genes)
    } else if(input$dge_list_sc == "Cluster DGE"){
      if(is.null(sc_mytab()) || nrow(sc_mytab()) == 0) return()
      dge <- sc_mytab()
      if(nrow(dge) == 0) return()
      if(d$n_genes_sc != "") {
        if(!is.null(input$all_genes_sc) && !(input$all_genes_sc)) {
          if(suppressWarnings(is.na(as.integer(na.omit(d$n_genes_sc))))) return()
          dge <- dge[1:min(as.integer(d$n_genes_sc), nrow(dge)), , drop = F]
        }
      }
      return(rownames(dge)) 
    } else {
      if(is.null(input$gse_customdge_fact) || is.null(input$gse_customdge_contrast)) return()
      if(is.null(d$sc_dge_byfactor[[input$gse_customdge_fact]][[input$gse_customdge_contrast]])) return()
      dge <- d$sc_dge_byfactor[[input$gse_customdge_fact]][[input$gse_customdge_contrast]]
      if(nrow(dge) == 0) return()
      sortCol <- "Log2 fold-change"
      if(!(sortCol %in% colnames(dge))) sortCol <- "Log fold-change"
      dge <- dge[order(dge[[sortCol]], decreasing = T), , drop = F]
      if(is.null(input$sc_gene_list_filtering)) return()
      if(input$sc_gene_list_filtering == 2) {
        dge <- dge[dge[[sortCol]] > 0, , drop = F]
      }
      if(input$sc_gene_list_filtering == 3) {
        dge <- dge[dge[[sortCol]] < 0, , drop = F]
      }
      if(nrow(dge) == 0) return()
      if(d$n_genes_sc != "") {
        if(!is.null(input$all_genes_sc) && !(input$all_genes_sc)) {
          if(suppressWarnings(is.na(as.integer(na.omit(d$n_genes_sc))))) return()
          dge <- dge[1:min(as.integer(d$n_genes_sc), nrow(dge)), , drop = F]
        }
      }
      return(rownames(dge))
    }
  })
  
  # SC-DGE-DGE - data input for dge analyis
  sc_mytab_heatmap <- reactive({
    if(is.null(d$SCV) || is.null(d$res) || is.null(input$sc_contrasts_heatmap) || 
       is.null(d$SCV[[d$res]]) || is.null(input$sc_gene_list_filtering_heatmap)) return()
    if(!(input$sc_contrasts_heatmap %in% names(d$SCV[[d$res]]@DEcombn))) return()
    df <- d$SCV[[d$res]]@DEcombn[names(d$SCV[[d$res]]@DEcombn) %in% input$sc_contrasts_heatmap][[1]]
    if(nrow(df) == 0) return()
    df <- df[order(1/abs(df$logGER), df$pVal, decreasing = F), , drop = F]
    if(input$sc_gene_list_filtering_heatmap == 2) df <- df[df$logGER > 0, , drop = F]
    if(input$sc_gene_list_filtering_heatmap == 3) df <- df[df$logGER < 0, , drop = F]
    if(nrow(df) == 0) return()
    return(df)
  })
  
  # Reactive expression to return raw counts based on current selections
  rawcounts_sc <- reactive({
    req(SubmitData$data_type, d$inD)
    if(SubmitData$data_type != "Single-cell") return()
    GetAssayData(object = d$inD, slot = "counts")
  })
  
  # Reactive expression to return normalized data based on current selections
  normcounts_sc <- reactive({
    req(SubmitData$data_type, d$inD)
    if(SubmitData$data_type != "Single-cell") return()
    GetAssayData(object = d$inD, slot = "data")
  })
  
  # Returns single-cell clusters based on current selections
  clusters_sc <- reactive({
    req(SubmitData$data_type, d$SCV, d$res, input$dge_list_heatmap_sc)
    if(SubmitData$data_type != "Single-cell") return()
    if(!(d$res %in% names(d$SCV))) return()
    Clusters(d$SCV[[d$res]])
  })
  
  # Returns cells to be used in heatmap
  heatmapcells_sc <- reactive({
    req(SubmitData$data_type, d$inD, d$SCV, d$res, input$dge_list_heatmap_sc)
    if(SubmitData$data_type != "Single-cell" || is.null(clusters_sc())) return()
    if(!(d$res %in% names(d$SCV))) return()
    cells <- names(clusters_sc())
    if(!is.null(input$contrast_groups_only_sc) && input$contrast_groups_only_sc &&
       input$dge_list_heatmap_sc %in% c("Cluster DGE", "Custom DGE") && 
       !is.null(input$sc_contrasts_heatmap)) {
      if(input$dge_list_heatmap_sc == "Cluster DGE") {
        if(is.null(input$sc_contrasts_heatmap)) return()
        groups <- strsplit(input$sc_contrasts_heatmap, split = "-")[[1]]
        cells <- cells[as.character(clusters_sc()) %in% groups]
      } else {
        if(is.null(input$heatmap_customdge_contrast)) return()
        groups <- strsplit(input$heatmap_customdge_contrast, split = "-")[[1]]
        if(!("Rest" %in% groups)) {
          cells <- colnames(d$inD)[d$inD@meta.data[[input$heatmap_customdge_fact]] %in% groups]
        }
      }
    }
    return(cells)
  })
  
  # Returns genes to be used in single-cell heatmap
  heatmapgenes_sc <- reactive({
    req(SubmitData$data_type, input$dge_list_heatmap_sc)
    if(SubmitData$data_type != "Single-cell" || is.null(normcounts_sc())) return()
    if(input$dge_list_heatmap_sc == "Manual") {
      if(is.null(input$gene_list_heatmap_sc) || input$gene_list_heatmap_sc == "") return()
      genes <- toupper(unlist(strsplit(input$gene_list_heatmap_sc, split = "\\s+|,\\s?")))
      if(length(genes) == 0) return()
      genes <- genes[genes %in% toupper(rownames(normcounts_sc()))]
      if(length(genes) == 0) return()
      return(genes)
    } else if(input$dge_list_heatmap_sc == "Cluster DGE") {
      if(is.null(sc_mytab_heatmap()) || nrow(sc_mytab_heatmap()) == 0 ) return()
      dge <- sc_mytab_heatmap()
      if(d$n_genes_heatmap_sc != "") {
        if(!(input$all_genes_heatmap_sc)) {
          if(suppressWarnings(is.na(as.integer(na.omit(d$n_genes_heatmap_sc))))) return()
          dge <- dge[1:min(as.integer(d$n_genes_heatmap_sc), nrow(dge)), , drop = F]
        }
      }
      return(rownames(dge))
    } else {
      if(is.null(input$heatmap_customdge_fact) || is.null(input$heatmap_customdge_contrast)) return()
      if(is.null(d$sc_dge_byfactor[[input$heatmap_customdge_fact]][[input$heatmap_customdge_contrast]])) return()
      dge <- d$sc_dge_byfactor[[input$heatmap_customdge_fact]][[input$heatmap_customdge_contrast]]
      if(nrow(dge) == 0) return()
      sortCol <- "Log2 fold-change"
      if(!(sortCol %in% colnames(dge))) sortCol <- "Log fold-change"
      dge <- dge[order(dge[[sortCol]], decreasing = T), , drop = F]
      if(is.null(input$sc_gene_list_filtering_heatmap)) return()
      if(input$sc_gene_list_filtering_heatmap == 2) {
        dge <- dge[dge[[sortCol]] > 0, , drop = F]
      }
      if(input$sc_gene_list_filtering_heatmap == 3) {
        dge <- dge[dge[[sortCol]] < 0, , drop = F]
      }
      if(nrow(dge) == 0) return()
      if(d$n_genes_heatmap_sc != "") {
        if(!is.null(input$all_genes_heatmap_sc) && !(input$all_genes_heatmap_sc)) {
          if(suppressWarnings(is.na(as.integer(na.omit(d$n_genes_heatmap_sc))))) return()
          dge <- dge[1:min(as.integer(d$n_genes_heatmap_sc), nrow(dge)), , drop = F]
        }
      }
      return(rownames(dge))
    }
  })
  
  # SC-DGE-Heatmap - pull dge genes or custom gene list
  heattran1_sc <- reactive({
    req(SubmitData$data_type)
    
    if(SubmitData$data_type != "Single-cell") return()
    if(is.null(heatmapcells_sc()) || length(heatmapcells_sc()) == 0) return()
    if(is.null(heatmapgenes_sc()) || length(heatmapgenes_sc()) == 0) return()
    if(is.null(normcounts_sc()) || nrow(normcounts_sc()) == 0) return()
    
    hmData <- as.matrix(normcounts_sc()[heatmapgenes_sc(), heatmapcells_sc(), drop = F])
    if(nrow(hmData) == 0 || ncol(hmData) == 0) return()
    return(list(hmData))
  })

  # SC-DGE-GSE - store filtered top 50/user selected genes
  # in a reactive exp
  heattran2_sc <- eventReactive(input$submit_list_sc, {
    withProgress(message = "Building heatmap...", value = 0, {
      incProgress(1/4)
      rescaled_mat <- heattran1_sc()[[1]]
      rescaled_mat_unclipped <- rescaled_mat
      keyLabel <- "Normalized expression"
      if(!is.null(input$scale_genes_sc) && input$scale_genes_sc) {
        keyLabel <- "Row z-score"
        rescaled_mat <- t(scale(t(rescaled_mat)))
        rescaled_mat_unclipped <- rescaled_mat
        # Clip values like Seurat
        rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat > 2.5), values = 2.5)
        rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat < -2.5), values = -2.5)
        # re-scale matrix -2 to 2 as per CM's request
        rescaled_mat <- scales::rescale(rescaled_mat, to = c(-2, 2))
        rescaled_mat_unclipped <- scales::rescale(rescaled_mat_unclipped, to = c(-2, 2))
      } else {
        # Clip values like Seurat
        rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat > 6), values = 6)
        rescaled_mat <- replace(x = rescaled_mat, list = (rescaled_mat < -2.5), values = -2.5)
      }
      # replace NAs with 0s
      rescaled_mat[is.na(rescaled_mat)] <- 0
      rescaled_mat_unclipped[is.na(rescaled_mat_unclipped)] <-  0
      # remove rows with all 0s
      rescaled_mat <- rescaled_mat[which(rowSums(rescaled_mat) != 0), , drop = F]
      if(!is.null(input$col_type_sc) && input$col_type_sc && ncol(rescaled_mat) > 5000) {
        set.seed(10)
        rescaled_mat <- rescaled_mat[, sample(x = 1:ncol(rescaled_mat), size = 5000), drop = F]
      }
      rescaled_mat_unclipped <- rescaled_mat_unclipped[rownames(rescaled_mat), colnames(rescaled_mat), drop = F]
      
      cells <- colnames(rescaled_mat)
      if(!is.null(input$contrast_groups_only_sc) && input$contrast_groups_only_sc &&
         input$dge_list_heatmap_sc == "Cluster DGE") {
        sampleLabels <- as.character(Clusters(d$SCV[[d$res]])[cells])
        names(sampleLabels) <- cells
      } else if(!is.null(input$contrast_groups_only_sc) && input$contrast_groups_only_sc &&
                input$dge_list_heatmap_sc == "Custom DGE") {
        sampleLabels <- as.character(d$inD@meta.data[cells, input$heatmap_customdge_fact])
        names(sampleLabels) <- cells
      } else {
        meta <- d$inD@meta.data
        sampleLabelDF <- meta[cells, input$heatfactorlabel_sc, drop = F]
        if(d$resType == "seurat_res" & "Cluster" %in% colnames(sampleLabelDF)) {
          sampleLabelDF$Cluster <- Clusters(d$SCV[[d$res]])[cells]
        }
        sampleLabels <- apply(sampleLabelDF, 1, function(i) paste(i, collapse = "_"))
      }
      sampleLabels <- factor(sampleLabels, levels = mixedsort(unique(sampleLabels), decreasing = T))
      pal <- hue_pal()(length(levels(sampleLabels)))
      names(pal) <- levels(sampleLabels)
      
      legendLabels <- unique(sampleLabels)
      if(length(legendLabels) == 1) {
        legendHeight <- 80
      } else {
        legendHeight <- ifelse(test = (length(legendLabels) == 2), yes = 60, no = 40) * length(legendLabels)
      }
      
      heatmapHeight <- max(round(800 * (nrow(heattran1_sc()[[1]]) / 30)), 267)
      if(nrow(heattran1_sc()[[1]]) > 50) heatmapHeight <- 600
      
      # Keep dendrogram and side colorbar same height when heatmap height changes
      subplotHeights <- c(40) / heatmapHeight
      subplotHeights <- c(subplotHeights, 1 - sum(subplotHeights))
      
      incProgress(1/4)
      
      colDend <- rowDend <- F
      if(input$row_type_sc && nrow(rescaled_mat) >= 2) {
        rowDist <- stats::dist(rescaled_mat_unclipped)
        rowHClust <- hclust(rowDist)
        rowDend <- as.dendrogram(rowHClust)
        if(nrow(rescaled_mat_unclipped) <= 2000) {
          rowDend <- seriate_dendrogram(dend = rowDend, x = rowDist, method = "OLO")
        }
      } 
      
      if(input$col_type_sc && nrow(rescaled_mat) >= 2) {
        # Keep dendrogram and side colorbar same height when heatmap height changes
        subplotHeights <- c(80, 40) / heatmapHeight
        subplotHeights <- c(subplotHeights, 1 - sum(subplotHeights))
        
        colDist <- stats::dist(t(rescaled_mat_unclipped))
        colHClust <- hclust(colDist)
        colDend <- as.dendrogram(colHClust)
        if(ncol(rescaled_mat_unclipped) <= 2000) {
          colDend <- seriate_dendrogram(dend = colDend, x = colDist, method = "OLO")
        }
      } else {
        rescaled_mat <- ReorderColsByFactor(x = rescaled_mat, fact = sampleLabels)
      }
      
      incProgress(1/4)
      
      rescaled_mat_unclipped <- rescaled_mat_unclipped[rownames(rescaled_mat), colnames(rescaled_mat), drop = F]
      
      customHoverMat <- matrix(
        data = rownames(rescaled_mat_unclipped), nrow = nrow(rescaled_mat_unclipped), ncol = ncol(rescaled_mat_unclipped), byrow = F
      )
      
      if(ncol(rescaled_mat_unclipped) <= 100) {
        customHoverMat <- matrix(paste0(
          matrix(paste0("value: ", as.matrix(round(rescaled_mat_unclipped, digits = 5)), "<br>")),
          matrix(paste0("gene: ", rownames(rescaled_mat_unclipped), "<br>"), nrow = nrow(rescaled_mat_unclipped), ncol = ncol(rescaled_mat_unclipped), byrow = F)
        ), nrow = nrow(rescaled_mat_unclipped), ncol = ncol(rescaled_mat_unclipped))
      }
      
      if(nrow(rescaled_mat) > 50) {
        ticks <- c(FALSE, FALSE)
      } else {
        ticks <- c(FALSE, TRUE)
      }
      cols <- heatCols_sc()
      
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
      
      colorbarLabel <- paste(input$heatfactorlabel_sc, collapse = "_")
      if(!is.null(input$contrast_groups_only_sc) && input$contrast_groups_only_sc) {
        colorbarLabel <- "seurat_clusters"
        if(!is.null(d$res) && grepl(pattern = "^Comp:", x = d$res)) {
          colorbarLabel <- "Set"
        }
      }
      
      sampleLabels <- sampleLabels[colnames(rescaled_mat)]
      df <- data.frame("a" = sampleLabels)
      names(df) <- colorbarLabel
      
      # Keep color scale key same height when heatmap height changes
      colorscaleKeyHeight <- 200 / heatmapHeight
      
      # Color scale and legend y position
      legendYPos <- subplotHeights[length(subplotHeights)]
      
      bottomMargin <- max(d$totalFigureHeight_sc - heatmapHeight, 0)
      if(!is.dendrogram(rowDend) || !is.dendrogram(colDend)) bottomMargin <- 0
      
      incProgress(1/4)
    })

    return(list(
      rescaled_mat = rescaled_mat,
      rescaled_mat_unclipped = rescaled_mat_unclipped,
      cols = cols,
      ticks = ticks,
      customHoverMat = customHoverMat,
      df = df,
      pal = pal,
      colorscaleKeyHeight = colorscaleKeyHeight,
      legendYPos = legendYPos,
      keyLabel = keyLabel,
      rowDend = rowDend,
      colDend = colDend,
      heatmapHeight = heatmapHeight,
      subplotHeights = subplotHeights,
      colorbarLabel = colorbarLabel,
      legendHeight = legendHeight,
      label = label,
      bottomMargin = bottomMargin,
      sampleLabels = sampleLabels
    ))
    
  })

  # SC-DGE-GSE - update select gene from heatmap plotly
  observeEvent(input$submit_list_sc, {
    # Reset event_data("plotly_selected")
    runjs("Shiny.setInputValue('plotly_selected-A', null);")
    d$heatmap1Click <- NULL
    d$newHeat <- F
    
    if(is.null(heattran1_sc())) return()
    # Get legend labels and use to determine total heatmap height
    if(!is.null(input$contrast_groups_only_sc) && input$contrast_groups_only_sc &&
       input$dge_list_heatmap_sc == "Cluster DGE") {
      sampleLabels <- as.character(Clusters(d$SCV[[d$res]]))
    } else if(!is.null(input$contrast_groups_only_sc) && input$contrast_groups_only_sc &&
              input$dge_list_heatmap_sc == "Custom DGE") {
      sampleLabels <- d$inD@meta.data[[input$heatmap_customdge_fact]]
    } else {
      meta <- d$inD@meta.data
      sampleLabelDF <- meta[, input$heatfactorlabel_sc, drop = F]
      if(d$resType == "seurat_res" & "Cluster" %in% colnames(sampleLabelDF)) {
        sampleLabelDF$Cluster <- Clusters(d$SCV[[d$res]])
      }
      sampleLabels <- apply(sampleLabelDF, 1, function(i) paste(i, collapse = "_"))
    } 
    legendLabels <- unique(sampleLabels)

    # Legend height
    if(length(legendLabels) == 1) {
      legendHeight <- 80
    } else {
      legendHeight <- ifelse(test = (length(legendLabels) == 2), yes = 60, no = 40) * length(legendLabels)
    }

    heatmapHeight <- max(round(800 * (nrow(heattran1_sc()[[1]]) / 30)), 267)
    if(nrow(heattran1_sc()[[1]]) > 50) heatmapHeight <- 600

    if(input$col_type_sc) {
      if(legendHeight > heatmapHeight) {
        totalFigureHeight <- legendHeight + 120
      } else {
        totalFigureHeight <- heatmapHeight
      }
    } else {
      if(legendHeight > heatmapHeight) {
        totalFigureHeight <- legendHeight + 40
      } else {
        totalFigureHeight <- heatmapHeight
      }
    }
    
    d$totalFigureHeight_sc <- totalFigureHeight
  })

  # SC-DGE-GSE - heatmap of top 50/user selected genes
  output$heatplot1_sc <- renderPlotly({
    
    if(is.null(heattran2_sc())) return()
    
    isolate({
      withProgress(message = "Rendering heatmap...", value = 0, {
        incProgress(1/4)
        p <- heatmaply(
          x = heattran2_sc()$rescaled_mat,
          plot_method = "plotly",
          colors = heattran2_sc()$cols,
          dend_hoverinfo = F,
          custom_hovertext = heattran2_sc()$customHoverMat,
          showticklabels = heattran2_sc()$ticks,
          col_side_colors = heattran2_sc()$df,
          col_side_palette = heattran2_sc()$pal,
          colorbar_len = heattran2_sc()$colorscaleKeyHeight,
          key.title = heattran2_sc()$keyLabel,
          colorbar_xpos = -0.3,
          colorbar_ypos = heattran2_sc()$legendYPos,
          colorbar_yanchor = "top",
          fontsize_row = 14,
          fontsize_col = 14,
          Rowv = heattran2_sc()$rowDend,
          Colv = heattran2_sc()$colDend,
          revC = T,
          width = 1000,
          height = heattran2_sc()$heatmapHeight,
          subplot_heights = heattran2_sc()$subplotHeights
        ) %>%
          colorbar(
            tickfont = list(size = 14, family = "Noto Sans JP"),
            titlefont = list(size = 16, family = "Noto Sans JP"),
            which = 2
          ) %>%
          colorbar(
            title = heattran2_sc()$colorbarLabel,
            len = heattran2_sc()$legendHeight,
            lenmode = "pixels",
            y = heattran2_sc()$legendYPos,
            yanchor = "top",
            tickfont = list(size = 14, family = "Noto Sans JP"),
            titlefont = list(size = 16, family = "Noto Sans JP"),
            which = 1
          ) %>%
          layout(
            font = list(family = "Noto Sans JP", size = 14),
            hoverlabel = heattran2_sc()$label,
            legend = list(font = list(size = 14, family = "Noto Sans JP")),
            margin = list(b = heattran2_sc()$bottomMargin)
          ) 
        incProgress(1/4)
        if(is.dendrogram(heattran2_sc()$colDend) && is.dendrogram(heattran2_sc()$rowDend)) {
          p$x$data[[3]]$hoverinfo <- "none"
        } else if(!is.dendrogram(heattran2_sc()$colDend)) {
          p$x$data[[1]]$text <- "none"
        } else {
          p$x$data[[2]]$text <- "none"
          p$x$layout$yaxis$showticklabels <- F
        }
        incProgress(1/4)
      })
      isolate(return(p))
    })
  })
  
  output$heatplot1_sc_ui <- renderUI({
    plotlyOutput("heatplot1_sc", width = 1000, height = d$totalFigureHeight_sc)
  })

  # SC-DGE-GSE - download button for heatmap (PDF)
  output$dlqcheatplot1pdf_sc <- renderUI({
    req(input$submit_list_sc, !d$newHeat)
    downloadButton("dlqcheatplot1pdfimg_sc", "Download static plot (PDF)")
  })

  # SC-DGE-GSE - download file heatmap (PDF)
  output$dlqcheatplot1pdfimg_sc <- downloadHandler(
    filename =  function() {
      paste("qc-heatmap_sc.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing for download...", value = 0, {
        incProgress(1/2)
        pdf(file, width = 7.5, height = 7.5, onefile = FALSE)
        qcHeatMap(
          heat = heattran2_sc()$rescaled_mat,
          color = heattran2_sc()$cols,
          rows = heattran2_sc()$rowDend,
          col = heattran2_sc()$colDend,
          nam = heattran2_sc()$sampleLabels,
          keyLabel = heattran2_sc()$keyLabel,
          colorLabel = heattran2_sc()$colorbarLabel,
          pal = heattran2_sc()$pal,
          type = "pdf"
        )
        dev.off()
        incProgress(1/2)
      })
    }
  )

  # SC-DGE-GSE - download button for heatmap (PNG)
  output$dlqcheatplot1png_sc <- renderUI({
    req(input$submit_list_sc, !d$newHeat)
    downloadButton("dlqcheatplot1pngimg_sc", "Download static plot (PNG)")
  })

  # SC-DGE-GSE - download file for heatmap (PNG)
  output$dlqcheatplot1pngimg_sc <- downloadHandler(
    filename =  function() {
      paste("qc-heatmap_sc.png")
    },
    content = function(file) {
      withProgress(message = "Preparing for download...", value = 0, {
        incProgress(1/2)
        png(file, width = 900, height = 900)
        qcHeatMap(
          heat = heattran2_sc()$rescaled_mat,
          color = heattran2_sc()$cols,
          rows = heattran2_sc()$rowDend,
          col = heattran2_sc()$colDend,
          nam = heattran2_sc()$sampleLabels,
          keyLabel = heattran2_sc()$keyLabel,
          colorLabel = heattran2_sc()$colorbarLabel,
          pal = heattran2_sc()$pal,
          type = "png"
        )
        dev.off()
        incProgress(1/2)
      })
    }
  )

  # SC-DGE-GSE - choose factor for heatmap
  output$heatfactor_sc <- renderUI({
    req(d$inD)
    if(is.null(input$submit_list_sc) || input$submit_list_sc == 0 ||
       is.null(d$heatmap1Click)) return()
    myChoices <- colnames(d$inD@meta.data)[!(colnames(d$inD@meta.data) %in% c("nCount_RNA", "nFeature_RNA"))]
    if("seurat_clusters" %in% myChoices) {
      myChoices <- c("seurat_clusters", myChoices[myChoices != "seurat_clusters"])
    }
    mySelected <- myChoices[1]
    selectInput(
      inputId = "heatfactor_sc",
      label = "Choose factor",
      choices = myChoices,
      selected = mySelected
    )
  })

  # SC-DGE-GSE - reactive expression to store gene clicked in heatmap
  gene_selection_sc <- eventReactive(input$submit_list_sc, {
    req(input$heatfactor_sc, !d$newHeat, input$submit_list_sc)
    s <- d$heatmap1Click
    return(s)
  })

  # SC-DGE-GSE - boxplot of counts by gene selected from heatmap
  output$heatplot2_sc <- renderPlotly({
    if(is.null(d$resType) || is.null(input$heatfactor_sc) || d$newHeat || 
       is.null(input$submit_list_sc) || input$submit_list_sc == 0 || 
       is.null(d$heatmap1Click) || is.null(input$dge_list_heatmap_sc) ||
       is.null(input$scale_genes_sc)) return()
    if(d$resType == "seurat_res" && is.null(input$overall_res)) return()
    
    s <- d$heatmap1Click
    if(is.null(s)) return()
    
    if(is.dendrogram(heattran2_sc()$rowDend)) {
      gene <- rev(labels(heattran2_sc()$rowDend))[s$y]
    } else {
      gene <- rev(rownames(heattran2_sc()$rescaled_mat))[s$y]
    }

    if (d$resType == "seurat_res") {
      res <- gsub(".*_|\\:.*", "", input$overall_res)
      if(!(res %in% names(seurat_only()))) return(NULL)
      seur <- seurat_only()[grepl(res, names(seurat_only()))]
      seur <- seur[[1]]
    } else {
      seur <- seurat_only()
    }
    rc.data <- as.matrix(GetAssayData(seur, slot = "data"))
    meta <- seur@meta.data
    
    test <- getGenes(
      rc.data = rc.data,
      id = gene,
      coldata = meta,
      type = SubmitData$data_type
    )

    tooltips <- paste0("<b>Sample:</b> ",
                       test$sample,
                       "<br />",
                       "<b>Counts:</b> ",
                       round(test$counts, 3))
    fact <- test[, input$heatfactor_sc]
    
    m <- list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
    )
    
    y <- test[, "counts"]
    x_ord <- factor(fact, levels = mixedsort(unique(fact), decreasing = F))
    pal <- hue_pal()(length(levels(x_ord)))
    names(pal) <- rev(levels(x_ord))
    
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
    
    p <- plot_ly(
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      x = x_ord,
      y = y,
      color = x_ord,
      colors = pal,
      text = tooltips,
      marker = list(size = 7),
      hoverinfo = "text",
      hoverlabel = label
    ) %>%
      layout(
        title = paste(gene, "Counts"),
        titlefont = list(size = 16, family = "Noto Sans JP"),
        xaxis = list(title = paste(input$heatfactor_sc), tickfont = list(size = 14, family = "Noto Sans JP")),
        yaxis = list(title = "Normalized counts", tickfont = list(size = 14, family = "Noto Sans JP"), titlefont = list(size = 16, family = "Noto Sans JP")),
        legend = list(font = list(size = 14, family = "Noto Sans JP")),
        margin = m,
        font = list(family = "Noto Sans JP")
      )
    isolate(p)
  })

  # SC-DGE-GSE - download button for heat counts (PDF)
  output$dlqcheatplot2pdf_sc <- renderUI({
    req(input$submit_list_sc, !d$newHeat, d$heatmap1Click)
    downloadButton("dlqcheatplot2pdfimg_sc", label = "Download static plot (PDF)")
  })

  # SC-DGE-GSE - download file for heat counts (PDF)
  output$dlqcheatplot2pdfimg_sc <- downloadHandler(
    filename =  function() {
      paste("qc-heat-counts_sc.pdf")
    },
    content = function(file) {
      s <- gene_selection_sc()
      
      if(is.dendrogram(heattran2_sc()$rowDend)) {
        gene <- rev(labels(heattran2_sc()$rowDend))[s$y]
      } else {
        gene <- rev(rownames(heattran2_sc()$rescaled_mat))[s$y]
      }
      
      if (d$resType == "seurat_res") {
        res <- gsub(".*_|\\:.*", "", input$overall_res)
        if(!(res %in% names(seurat_only()))) return(NULL)
        seur <- seurat_only()[grepl(res, names(seurat_only()))]
        seur <- seur[[1]]
      } else {
        seur <- seurat_only()
      }
      rc.data <- as.matrix(GetAssayData(seur, slot = "data"))
      meta <- seur@meta.data
      
      test <- getGenes(
        rc.data = rc.data,
        id = gene,
        coldata = meta,
        type = SubmitData$data_type
      )
      
      fact <- test[, input$heatfactor_sc]
      fact <- factor(fact, levels = mixedsort(unique(fact), decreasing = F))
      
      y <- test[, "counts"]

      pdf(file, width = 7, height = 6.5, onefile = FALSE) # open the pdf device
      qcHeatCount(
        data = test,
        fact = fact,
        var = y,
        xaxis = input$heatfactor_sc,
        title = gene,
        type = "sc"
      )
      dev.off()
    }
  )

  # SC-DGE-GSE - download button for heat counts (PNG)
  output$dlqcheatplot2png_sc <- renderUI({
    req(input$submit_list_sc, !d$newHeat,
        d$heatmap1Click)
    downloadButton("dlqcheatplot2pngimg_sc", label = "Download static plot (PNG)")
  })

  # SC-DGE-GSE - download file for heat counts (PNG)
  output$dlqcheatplot2pngimg_sc <- downloadHandler(
    filename =  function() {
      paste("qc-heat-counts_sc.png")
    },
    content = function(file) {
      s <- gene_selection_sc()
      
      if(is.dendrogram(heattran2_sc()$rowDend)) {
        gene <- rev(labels(heattran2_sc()$rowDend))[s$y]
      } else {
        gene <- rev(rownames(heattran2_sc()$rescaled_mat))[s$y]
      }
      
      if (d$resType == "seurat_res") {
        res <- gsub(".*_|\\:.*", "", input$overall_res)
        if(!(res %in% names(seurat_only()))) return(NULL)
        seur <- seurat_only()[grepl(res, names(seurat_only()))]
        seur <- seur[[1]]
      } else {
        seur <- seurat_only()
      }
      rc.data <- as.matrix(GetAssayData(seur, slot = "data"))
      meta <- seur@meta.data
      
      test <- getGenes(
        rc.data = rc.data,
        id = gene,
        coldata = meta,
        type = SubmitData$data_type
      )
      
      fact <- test[, input$heatfactor_sc]
      fact <- factor(fact, levels = mixedsort(unique(fact), decreasing = F))
      
      y <- test[, "counts"]

      png(file, width = 800, height = 750)
      qcHeatCount(
        data = test,
        fact = fact,
        var = y,
        xaxis = input$heatfactor_sc,
        title = gene,
        type = "sc"
      )
      dev.off()
    }
  )

  ######### MANUALLY SELECT CELLS FOR DGE TAB CODE: SCRNA-SEQ#####################

  # SC-DGE-MSC - for manually select cells
  output$manuallyselectcells_header <- renderUI({
    
    h3("Manually select cells for DGE")
  })

  # SC-DGE-MSC - select metadata to overlay onto plot
  output$tsneSelDEcol <- renderUI({
    req(d$inD, d$res)
    tempChoices = colnames(d$inD@meta.data)
    tempChoices = tempChoices[tempChoices != "var"]
    tempSel = tempChoices[1]
    selectInput("tsneSelDEcol","Metadata overlay and filtering",
                choices=tempChoices, selected = tempSel)
  })

  # SC-DGE-MSC - store filtering choices in object d
  observeEvent(input$plusFilt,{
    d$filts <- unique(c(d$filts,input$tsneSelDEcol))
  })

  # SC-DGE-MSC - store unselected filtering choices in object d
  observeEvent(input$minusFilt,{
    d$filts <- d$filts[-which(d$filts == input$tsneSelDEcol)]
  })

  # SC-DGE-MSC - set d$filts to NULL when input$minusFiltAll is
  # selected
  observeEvent(input$minusFiltALL,{ d$filts <- NULL })

  # SC-DGE-MSC - reactive expression for MDpicker cluster selection
  filtValues <- reactive({
    sapply(d$filts,function(MD) {
      if (MD == paste("Clusters:",d$res)) {
        temp_inputSlot <- "MDpicker_clusts"
      } else {
        temp_inputSlot <- paste0("MDpicker_", which(colnames(d$MD) == MD))
      }
      return(input[[temp_inputSlot]])
    },simplify=F)
  })

  # SC-DGE-MSC - reactive exp to store selected cells
  makeMDpicker <- reactive({
    if(is.null(d$filts)) return()
    fluidPage(
      style = "border: 1px solid #eee; padding-left: 10px; padding-bottom: 10px;",
      lapply(d$filts,function(MD) {
        temp_val <- isolate(filtValues()[[MD]])
        if (MD == "") {
        } else if (MD == paste("Clusters:",d$res)) {
          selectInput("MDpicker_clusts",label="Select cluster(s)",
                      choices=levels(Clusters(d$SCV[[d$res]])),
                      multiple=T,selected=temp_val, width = "")
        } else if (is.factor(d$MD[,MD]) | is.character(d$MD[,MD])) {
          selectInput(paste0("MDpicker_",which(colnames(d$MD) == MD)),
                      label=paste0("Select cells by ",MD),
                      choices=levels(as.factor(d$MD[,MD])),
                      multiple=T,selected=temp_val)
        } else {
          if (is.null(temp_val)) { temp_val <- range(as.numeric(d$MD[,MD])) }
          sliderInput(paste0("MDpicker_",which(colnames(d$MD) == MD)),
                      label=paste0("Select cells by ",MD," range"),
                      min=min(as.numeric(d$MD[,MD])),max=max(as.numeric(d$MD[,MD])),value=temp_val)
        }
      })
    )
  })

  # SC-DGE-MSC - reactive exp select cells
  plotlySelected <- reactive({
    plotlySelected <- event_data("plotly_selected")
    if(is.null(plotlySelected)) return()
    plotlySelected <- plotlySelected[!is.na(plotlySelected$key), ]
    return(plotlySelected$key)
  })

  # SC-DGE-MSC - reactive hover select cells
  plotlyHover <- reactive({
    plotlyHover = event_data("plotly_hover")
    if(is.null(plotlyHover)) NULL else plotlyHover$key[1]
  })

  # SC-DGE-MSC - remove selected cells & reset
  output$MDfilts <- renderUI({
    req(d$res, d$inD, d$MD, d$SCV)
    if(is.null(makeMDpicker())) return()
    makeMDpicker()
  })

  # SC-DGE-MSC - action button to remove all filters
  output$MDfiltsRemoveAll <- renderUI({
    req(d$res, d$inD, d$SCV)
    if (length(d$filts) > 0) {
      actionButton("minusFiltALL","Remove all filters",icon("minus"),
                   style="color: #008000; background: #fff")
    }
  })

  # SC-DGE-MSC - action button for selecting set A of cells
  output$addCellsA <- renderUI({
    actionButton("addCellsA", "Set A: Add cells", icon("plus"),
                 style = "color: #fff; border-color: #00bfff; background:#00bfff")
  })

  # SC-DGE-MSC - action button for removing set A cells
  output$removeCellsA <- renderUI({
    actionButton("removeCellsA", "Set A: Remove cells", icon("minus"),
                 style = "color: #fff; border-color: #00bfff")
  })

  # SC-DGE-MSC - action button for clearing set A cells
  output$clearA <- renderUI({
    actionButton("clearA", "Set A: Clear",
                 style = "color: black; background: lightgray; border-color: lightgray")
  })

  # SC-DGE-MSC - reset plotly-select for A cells
  observeEvent(input$clearA, {
    d$a <- NULL
    runjs("Shiny.setInputValue('plotly_selected-A', null);")
  })

  # SC-DGE-MSC - action button for selecting set B of cells
  output$addCellsB <- renderUI({
    actionButton("addCellsB", "Set B: Add cells", icon("plus"),
                 style = "color: #fff; border-color: #00bfff; background:#00bfff")
  })

  # SC-DGE-MSC - action button for removing set B cells
  output$removeCellsB <- renderUI({
    actionButton("removeCellsB", "Set B: Remove cells", icon("minus"),
                 style = "color: #fff; border-color: #00bfff")
  })

  # SC-DGE-MSC - action button for clearing set B cells
  output$clearB <- renderUI({
    actionButton("clearB", "Set B: Clear",
                 style = "color: black; background: lightgray; border-color: lightgray")
  })

  # SC-DGE-MSC - reset plotly-select for B cells
  observeEvent(input$clearB, {
    d$b <- NULL
    runjs("Shiny.setInputValue('plotly_selected-A', null);")
  })

  # SC-DGE-MSC - select dimensionality reduction
  output$SelDE_EmbType <- renderUI({
    req(d$inD)
    temp_embs <- hasEmb(d$inD)
    temp_embs <- temp_embs[
      sapply(temp_embs,function(X) ncol(getEmb(d$inD,X))) >= 2 &
        sapply(temp_embs,function(X) nrow(getEmb(d$inD,X))) == nrow(getMD(d$inD))]
    temp_embs <- toupper(temp_embs)
    temp_embs <- gsub("TSNE", "tSNE", temp_embs)
    selectInput("SelDE_EmbType",label="Embedding",
                choices = temp_embs,
                selected = temp_embs[temp_embs %in% c("tSNE","UMAP")][1])
  })

  # SC-DGE-MSC - select x-axis label
  output$SelDE_EmbDimX <- renderUI({
    req(d$inD, input$SelDE_EmbType)
    selectInput("SelDE_EmbDimX",label="x-axis",
                choices= gsub("_", " ", colnames(getEmb(d$inD,input$SelDE_EmbType))),
                selected= gsub("_", " ", colnames(getEmb(d$inD,input$SelDE_EmbType))[1]))
  })

  # SC-DGE-MSC - select y-axis label
  output$SelDE_EmbDimY <- renderUI({
    req(d$inD, input$SelDE_EmbType)
    selectInput("SelDE_EmbDimY",label="y-axis",
                choices= gsub("_", " ", colnames(getEmb(d$inD,input$SelDE_EmbType))),
                selected= gsub("_", " ", colnames(getEmb(d$inD,input$SelDE_EmbType))[2]))
  })

  # SC-DGE-MSC - warning for manually selecting cells
  output$gse_warning_manuallyselectcells_sc <- renderText(d$gse_warning_msg)

  # SC-DGE-MSC - plot dim reduction based on inputs for manually select cells
  output$tsneSelDE <- renderPlotly({
    req(d$inD, d$SCV, d$res, input$tsneSelDEcol, input$SelDE_EmbType,
        input$SelDE_EmbDimX, input$SelDE_EmbDimY)
    if(!(input$tsneSelDEcol %in% colnames(d$inD@meta.data))) return()
    if(length(d$res) == 0 || !(d$res %in% names(d$SCV))) {
      p <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
        theme(panel.border = element_rect(fill = NA)) +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("Select a clustering resolution in the",
                               "'Cluster/Set comparison of gene statistics' section",
                               "and click 'View' before using this tool.", sep = "\n"))
      yaxis <- list(automargin = TRUE)
      ggplotly(p, width = 750, height = 500) %>%
        layout(
          autosize = F, 
          yaxis = yaxis,
          font = list(family = "noto-sans-jp")
        )
    } else {
      cell <- colnames(d$inD)
      factorTooltipText <- paste0("</br>  <b>", input$tsneSelDEcol, " :</b> ")
      plotDims <- as.integer(c(
        sub("^.*(\\d+)$", replacement = "\\1", x = input$SelDE_EmbDimX),
        sub("^.*(\\d+)$", replacement = "\\1", x = input$SelDE_EmbDimY)
      ))
      xLab <- sub("_", replacement = " ", x = input$SelDE_EmbDimX)
      yLab <- sub("_", replacement = " ", x = input$SelDE_EmbDimY)
      plotType <- "dimplot"
      cols <- NULL
      if(is.numeric(d$inD@meta.data[[input$tsneSelDEcol]]) || input$tsneSelDEcol %in% c("nCount_RNA","nFeature_RNA")) {
        plotType <- "featureplot"
        cols <- rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
      }
      if(input$tsneSelDEcol == "seurat_clusters") {
        if(plotType == "dimplot") {
          p <- DimPlot(
            object = d$inD,
            dims = plotDims,
            reduction = tolower(input$SelDE_EmbType),
            group.by = input$tsneSelDEcol,
            label = F,
            cols = cols
          ) + xlab(xLab) + ylab(yLab)
        } else {
          p <- FeaturePlot(
            object = d$inD,
            features = input$tsneSelDEcol,
            dims = plotDims,
            reduction = tolower(input$SelDE_EmbType),
            label = F,
            cols = cols
          ) + xlab(xLab) + ylab(yLab)
        }
        p$data$cell <- cell
        p <- p + aes(
          key = cell, 
          text = paste(
            factorTooltipText, 
            p$data[[input$tsneSelDEcol]],
            "</br>",
            "<b>", xLab, ": </b>", round(p$data[, 1], 3),
            "</br>",
            "<b>", yLab, ": </b>", round(p$data[, 2], 3),
            "</br>",
            "<b> Cell : </b>", p$data$cell
          )
        )
        yaxis <- list(automargin = TRUE)
        ggplotly(p, width = 750, height = 500, tooltip = "text") %>%
          layout(
            autosize = F, 
            yaxis = yaxis, 
            dragmode = "lasso",
            font = list(family = "Noto Sans JP")
          )
      } else {
        if(plotType == "dimplot") {
          p <- DimPlot(
            object = d$inD,
            dims = plotDims,
            reduction = tolower(input$SelDE_EmbType),
            group.by = input$tsneSelDEcol,
            label = F,
            cols = cols
          ) + xlab(xLab) + ylab(yLab)
        } else {
          p <- FeaturePlot(
            object = d$inD,
            features = input$tsneSelDEcol,
            dims = plotDims,
            reduction = tolower(input$SelDE_EmbType),
            label = F,
            cols = cols
          ) + xlab(xLab) + ylab(yLab)
        }
        p$data$cell <- cell
        p <- p + aes(
          key = cell, 
          text = paste(
            factorTooltipText, p$data[[input$tsneSelDEcol]],
            "</br>",
            "<b>", xLab, ": </b>", round(p$data[, 1], 3),
            "</br>",
            "<b>", yLab, ": </b>", round(p$data[, 2], 3),
            "</br>",
            "<b> Cell : </b>", p$data$cell
          )
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
        
        yaxis <- list(automargin = TRUE)
        ggplotly(p, width = 750, height = 500, tooltip = "text") %>%
          layout(
            autosize = F, 
            yaxis = yaxis, 
            dragmode = "lasso",
            font = list(family = "Noto Sans JP")
          ) %>%
          style(hoverlabel = label)
      }
    }
  })

  # SC-DGE-MSC - reactive exp to store cell selection from filters and/or brush
  currSel <- reactive({
    temp_points <- rownames(brushedPoints(
      as.data.frame(getEmb(d$inD,input$SelDE_EmbType)[,c(gsub(" ", "_", input$SelDE_EmbDimX), gsub(" ", "_", input$SelDE_EmbDimY))]),
      input$tsneSelDEbrush,
      xvar=gsub(" ", "_", input$SelDE_EmbDimX), yvar=gsub(" ", "_", input$SelDE_EmbDimY)
    ))
    temp_picker <- sapply(names(filtValues()),function(X) {
      if (length(filtValues()[[X]]) < 1) {
        rep(T,nrow(d$MD))
      } else {
        if (X == paste("Clusters:",d$res)) {
          Clusters(d$SCV[[d$res]]) %in% filtValues()[[X]]
        } else if (is.factor(d$MD[,X]) | is.character(d$MD[,X])) {
          d$MD[,X] %in% filtValues()[[X]]
        } else {
          d$MD[,X] >= filtValues()[[X]][1] & d$MD[,X] <= filtValues()[[X]][2]
        }
      }
    },simplify=F)
    temp_picker <- as.logical(Reduce("*",temp_picker))
    if (length(temp_points) > 0 & length(temp_picker) > 0) {
      return(rownames(d$MD)[rownames(d$MD) %in% temp_points & temp_picker])
    } else if (length(temp_picker) > 0 & !all(temp_picker)) {
      return(rownames(d$MD)[temp_picker])
    } else if (length(temp_points) > 0) {
      return(temp_points)
    } else { return(character()) }
  })

  # SC-DGE-MSC - text for cells being hovered over in dim plot
  output$cellsHovered <- renderText(
    paste("Hovering over cell(s):", plotlyHover())
  )

  # SC-DGE-MSC - cells selected by lasso
  output$plotlySelected <- renderPrint({
    plotlySelected <- event_data("plotly_selected", source = "A")
    if(is.null(plotlySelected)) "No cells selected" else plotlySelected
  })

  # SC-DGE-MSC - stored data
  output$pData <- renderPrint({
    if(is.null(pData())) "pData is NULL" else pData
  })

  # SC-DGE-MSC - observed cells in pop A
  observeEvent(input$addCellsA,{
    d$a <- unique(c(d$a,currSel()))
    d$a <- unique(c(d$a,plotlySelected()))
  }, ignoreInit = T)

  # SC-DGE-MSC - removed cells in pop A
  observeEvent(input$removeCellsA,{
    d$a <- d$a[!d$a %in% currSel()]
    d$a <- d$a[!d$a %in% plotlySelected()]
  })

  # SC-DGE-MSC - observed cells in pop B
  observeEvent(input$addCellsB,{
    d$b <- unique(c(d$b,currSel()))
    d$b <- unique(c(d$b,plotlySelected()))
  }, ignoreInit = T)

  # SC-DGE-MSC - removed cells in pop B
  observeEvent(input$removeCellsB,{
    d$b <- d$b[!d$b %in% currSel()]
    d$b <- d$b[!d$b %in% plotlySelected()]
  })

  # SC-DGE-MSC - if a and b are null, render no text
  observeEvent(input$go2, {
    d$a <- d$b <- NULL
    output$calcText <- renderText("")
  })

  # SC-DGE-MSC - number of cells in pop A
  output$textSetA <- renderText(paste(length(d$a),"cells in set A."))

  # SC-DGE-MSC - number of cells in pop B
  output$textSetB <- renderText(paste(length(d$b),"cells in set B."))

  # SC-DGE-MSC - overlap message b/t pop A and B
  output$textOverlap <- renderText(
    if (!is.null(d$a) & !is.null(d$b) & length(intersect(d$a,d$b)) > 0) {
      paste(length(intersect(d$a,d$b)),"cells in both sets.",
            "Cells must be assigned to a single set prior to calculation.")
    }
  )

  # SC-DGE-MSC - name the comparison for easier id
  output$DEsetName <- renderUI({
    textInput("DEsetName","Short name for this comparison",
              value = "", placeholder="A-z0-9_ only please")
  })

  # SC-DGE-MSC - action button to run dge
  output$calcDE <- renderUI({
    actionButton("calcDE","Calculate differential gene expression", icon("play"))
  })

  # SC-DGE-MSC - render calculated text value
  output$calcText <- renderText(d$calcTextVal)

  # SC-DGE-MSC - gse for manually selecting cells tab
  sc_manuallyselectcells_func_enrich <- reactive({
    if(is.null(input$comps_manuallyselectcells_sc)) return()
    if(is.null(d$enrichRLib_manuallyselectcells)) {
      d$gse_warning_msg <- "Please choose at least one EnrichR library."
    } else if(length(d$enrichRLib_manuallyselectcells) == 0) {
      d$gse_warning_msg <- "Please choose at least one EnrichR library."
    } else if(d$enrichRLib_manuallyselectcells[1] == "") {
      d$gse_warning_msg <- "Please choose at least one EnrichR library."
    } else if(is.null(input$all_genes_manuallyselectcells_sc)) {
    } else if(!input$all_genes_manuallyselectcells_sc & d$n_genes_manuallyselectcells_sc != "" & suppressWarnings(is.na(as.integer(na.omit(d$n_genes_manuallyselectcells_sc))))) {
      d$gse_warning_msg <- "Please provide an integer in 'Use top n genes'."
    } else if(is.null(d$SCV[[input$comps_manuallyselectcells_sc]]@DEcombn["Set A-Set B"][[1]])) {
      d$gse_warning_msg <- "No DE genes."
    } else {
      d$gse_warning_msg <- ""
      withProgress(message = "Computing functional enrichment gene lists...", value = 0, {
        incProgress(0.3)
        isolate({
          # update input data w/ scClustViz
          dge <- d$SCV[[input$comps_manuallyselectcells_sc]]@DEcombn["Set A-Set B"][[1]]
          dge$id <- rownames(dge)
          if(input$sc_manuallyselectcells_gene_list_filtering == 2) dge <- dge[dge$logGER > 0, ]
          if(input$sc_manuallyselectcells_gene_list_filtering == 3) dge <- dge[dge$logGER < 0, ]
          dge <- dge[order(1/abs(dge$logGER), dge$pVal, decreasing = F), ]
          if(d$n_genes_manuallyselectcells_sc != "") {
            if(!(input$all_genes_manuallyselectcells_sc)) {
              if(suppressWarnings(!is.na(as.integer(na.omit(d$n_genes_manuallyselectcells_sc))))) {
                dge <- dge[1:min(as.integer(d$n_genes_manuallyselectcells_sc), nrow(dge)), ]
              } else {
                return(NULL)
              }
            }
          }

          uni_genes <- rownames(dge)
          d$sc_mytabTemp <- dge
          # load library list
          libs <- d$enrichRLib_manuallyselectcells
          # convert genes to a vector
          genes <- dge %>% pull(id)
          # run enrichR w/selected libraries
          enr <- enrichr(genes, libs)
          enr <- enr[sapply(enr, function(x) nrow(x) > 0)]
          if (length(enr) == 0) {
            return(NULL)
          } else {
            enrCheck <- sapply(enr, function(i) {
              if(ncol(i) == 1) {
                if(colnames(i)[1] == "X.html.") return(FALSE)
              }
              return(TRUE)
            })
            if(all(!enrCheck)) {
              d$gse_warning_msg <- "EnrichR is unavailable. Please try again later."
              return()
            }

            # return significant genes
            enr <- getSigTerms(enr, libs)
            enr <- enr[sapply(enr, function(x) !is.null(x))]
            enr <- rbindlist(enr)
            if (length(enr) == 0) {
              return(NULL)
            } else {
              incProgress(0.6)
              # function for removing redundant go terms
              enr <- enr[order(-score)]
              enr.go <- enr[grepl(pattern = "GO", x = libName) & grepl(pattern = "\\(GO", x = term), , drop = F]
              if (nrow(enr.go) > 0) {
                enr.go[,GOID := tstrsplit(term, "\\(GO")[2]]
                enr.go[,GOID := paste0('GO', GOID)]
                enr.go[,GOID := gsub("\\)", '', GOID)]
                enr.go[,term.short := tstrsplit(term, "\\(")[1]]
                enr.go[,term.short := trimws(term.short, which="right")]
                enr.go <- enr.go %>%
                  mutate(term = paste(term.short, " (", GOID, ")", sep = "")) %>%
                  dplyr::select(-c(term.short, GOID))
                enr.rest <- enr[!grepl('GO', libName),]
                enr <- rbind(enr.rest, enr.go)
              }
              enr <- enr[order(-score)]
              enr <- enr %>%
                mutate(pval = round(pval, digits = 3),
                       adjPval = round(adjPval, digits = 3),
                       Z_score = round(Z_score, digits = 3),
                       score = round(score, digits = 3))
              names(enr) <- c("Library name", "Library rank", "Gene count",
                              "Term", "Overlap", "P-value", "Adjusted p-value",
                              "Old p-value", "Old adjusted p-value", "Z-score",
                              "Score", "Gene list")
              enr <- enr %>%
                arrange(`P-value`)
              enr <- enr[, -c("Old p-value", "Old adjusted p-value"), with = FALSE]
              return(enr)
            }
          }
        })
      })
    }
  })

  # DGE - run dge for manually selected cells
  observeEvent(input$calcDE,{
    runjs("Shiny.setInputValue('plotly_selected-A', null);") # Reset event_data("plotly_selected")
    newRes <- paste0("Comp:",gsub("[^A-Za-z0-9_]","",input$DEsetName))
    if (!is.null(d$a) & !is.null(d$b) & length(intersect(d$a,d$b)) > 0) {
      d$calcTextVal <- "Sets can't overlap (please assign cells to only one set)."
    } else if (any(sapply(list(d$a,d$b),length) < 3)) {
      d$calcTextVal <- "Each set must contain at least 3 cells."
    } else if (nchar(input$DEsetName) < 1) {
      d$calcTextVal <- "Please name this comparison (in text box above)."
    } else if (newRes %in% d$compNames) {
      d$calcTextVal <- "This comparison name has already been used."
    } else if (grepl("^Clusters: Comp:", input$tsneSelDEcol)) {
      d$calcTextVal <- "Please choose a valid resolution from the 'Cluster/Set comparison of gene statistics' tab."
    } else {
      compNames <- c(d$compNames, newRes)
      d$calcTextVal <- ""
      withProgress({
        temp_warn <- options("warn")
        options(warn=-1)
        temp <- rep("Unselected",ncol(getExpr(d$inD,Param(d$SCV[[1]],"assayType")[1])))
        names(temp) <- colnames(getExpr(d$inD,Param(d$SCV[[1]],"assayType")[1]))
        temp[d$a] <- "Set A"
        temp[d$b] <- "Set B"
        d$SCV[[newRes]] <- sCVdata(Clusters=factor(temp,levels=c("Set A","Set B")),
                                   params=Param(d$SCV[[1]]))
        # gene stats per set
        incProgress(amount=1/6, detail="Cluster-wise gene stats")
        ClustGeneStats(d$SCV[[newRes]]) <- fx_calcCGS(nge=getExpr(d$inD,Param(d$SCV[[1]],"assayType")[1]),
                                                      cl=Clusters(d$SCV[[newRes]]),
                                                      exponent=Param(d$SCV[[newRes]], "exponent"),
                                                      pseudocount=Param(d$SCV[[newRes]], "pseudocount"))
        # de per cluster vs all other data
        incProgress(amount=2/6,detail="DE vs tissue logGER calculations")
        deTes <- fx_calcESvsRest(nge=getExpr(d$inD,Param(d$SCV[[1]],"assayType")[1]),
                                 cl=Clusters(d$SCV[[newRes]]),
                                 CGS=ClustGeneStats(d$SCV[[newRes]]),
                                 exponent=Param(d$SCV[[newRes]], "exponent"),
                                 pseudocount=Param(d$SCV[[newRes]], "pseudocount"),
                                 DRthresh=Param(d$SCV[[newRes]], "DRthresh"))
        incProgress(amount = 1/6, detail = "DE vs tissue Wilcoxon rank sum calculations")
        DEvsRest(d$SCV[[newRes]]) <- fx_calcDEvsRest(nge=getExpr(d$inD,Param(d$SCV[[1]],"assayType")[1]),
                                                     cl=Clusters(d$SCV[[newRes]]),
                                                     deTes=deTes)
        # de per cluster vs each other cluster
        incProgress(amount=1/6,detail="Calculating Set A vs Set B")
        deMes <- fx_calcEScombn(cl=Clusters(d$SCV[[newRes]]),
                                CGS=ClustGeneStats(d$SCV[[newRes]]),
                                DRthresh=Param(d$SCV[[newRes]], "DRthresh"))
        DEcombn(d$SCV[[newRes]]) <- fx_calcDEcombn(nge=getExpr(d$inD,Param(d$SCV[[1]],"assayType")[1]),
                                                   cl=Clusters(d$SCV[[newRes]]),
                                                   deMes=deMes)
        incProgress(amount=1/6,detail="Done")
        d$a <- d$b <- NULL
        d$filts <- NULL
        options(warn=temp_warn$warn)
      },message="DE calculations:")
      d$res <- newRes # Automatically update the view to show the calculated results.
    }
  }, ignoreInit = T)

  # SC-DGE-MSC - pull top 50 filtered dge genes or user selected gene list
  heattran1_manuallyselectcells_sc <- reactive({
    req(input$sc_manuallyselectcells_gene_list_filtering, d$SCV)
    if(is.null(input$sc_manuallyselectcells_gene_list_filtering)) return()
    if(d$calcTextVal != "") return()
    if(is.null(input$comps_manuallyselectcells_sc)) return()
    if(is.null(d$SCV)) return()
    if(is.null(d$SCV[[input$comps_manuallyselectcells_sc]])) return()
    if(is.null(d$SCV[[input$comps_manuallyselectcells_sc]]@DEcombn["Set A-Set B"][[1]])) return()
    cells <- names(Clusters(d$SCV[[input$comps_manuallyselectcells_sc]]))[!is.na(Clusters(d$SCV[[input$comps_manuallyselectcells_sc]]))]
    dge <- d$SCV[[input$comps_manuallyselectcells_sc]]@DEcombn["Set A-Set B"][[1]]
    genes <- rownames(dge)
    dge$id <- genes
    dge_filter = input$sc_manuallyselectcells_gene_list_filtering
    logfc_col = "logGER"
    pval_col = "pVal"
    if(dge_filter == 2) dge <- dge[dge[, logfc_col] > 0, ]
    if(dge_filter == 3) dge <- dge[dge[, logfc_col] < 0, ]
    dge <- dge[order(1/abs(dge[, logfc_col]), dge[, pval_col], decreasing = F), ]
    if(d$n_genes_manuallyselectcells_sc != "") {
      if(!(input$all_genes_manuallyselectcells_sc)) {
        if(suppressWarnings(!is.na(as.integer(na.omit(d$n_genes_manuallyselectcells_sc))))) {
          dge <- dge[1:min(as.integer(d$n_genes_manuallyselectcells_sc), nrow(dge)), ]
        } else {
          return(NULL)
        }
      }
    }
    uni_genes <- unique(dge$id)
    if (length(uni_genes) == 0) {
      return(NULL)
    }
    # convert genes to a vector
    genes <- dge %>% pull(id)
    if (length(genes) <= 1) {
      return(NULL)
    }
    dds.counts <- ddsout()[[3]]
    heat.counts <- ddstran()[[1]]
    heat.counts <- assay(heat.counts)
    topID <- order(rowMeans(dds.counts), decreasing = TRUE)
    heat.mat <- heat.counts[topID, ]
    heat.mat <- heat.mat[rownames(heat.mat) %in% genes, cells]
    if (length(genes) < 50 & length(genes) > 1) {
      heat.mat <- heat.mat[1:length(genes), ,drop = FALSE]
    }
    return(list(heat.mat))
  })

  # SC-DGE-MSC - header for downloading results
  output$manuallyselectcells_downloadGSE_header <- renderUI({
    req(input$comps_manuallyselectcells_sc, clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()
    h3("Download gene set enrichment results", style = "margin-top: 0px;")
  })
  
  # SC-DGE-MSC - complete gse on all, up in group 1 or group 2
  output$sc_manuallyselectcells_gene_list_filtering <- renderUI({
    req(input$calcDE, input$comps_manuallyselectcells_sc)
    if(d$calcTextVal != "") return()
    tempChoices = 1:3
    names(tempChoices) = c("All DE genes", "Up in Set A", "Up in Set B")
    awesomeRadio("sc_manuallyselectcells_gene_list_filtering", label = "Gene list filtering",
                 choices = tempChoices,
                 selected = 1, status = "success")
  })

  # SC-DGE-MSC - specify # of top genes to complete gse
  output$top_n_genes_manuallyselectcells_sc <- renderUI({
    req(input$calcDE, input$comps_manuallyselectcells_sc)
    if(d$calcTextVal != "") return()
    textInput("n_genes_manuallyselectcells_sc", label = "Top n genes", value = 100, placeholder = "ALL")
  })

  # SC-DGE-MSC - use all genes passing filter for gse
  output$use_all_genes_manuallyselectcells_sc <- renderUI({
    req(input$calcDE, input$comps_manuallyselectcells_sc)
    if(d$calcTextVal != "") return()
    prettyCheckbox("all_genes_manuallyselectcells_sc", label = "All genes\npassing filters", value = F, status = "default", icon = icon("check"))
  })

  # SC-DGE-MSC - reactive exp for storing data in object d
  observeEvent(input$all_genes_manuallyselectcells_sc, {
    d$fea_all_genes <- input$all_genes_manuallyselectcells_sc
  })

  # SC-DGE-MSC - dge messages on manually selected cells
  output$dge_info_manuallyselectcells_sc <- renderUI({
    req(input$calcDE, input$comps_manuallyselectcells_sc)
    if(is.null(heattran1_manuallyselectcells_sc())) {
      if(d$calcTextVal != "") return()
      d$dge_info_manuallyselectcells_sc <- "No DE genes."
    } else {
      dims = dim(heattran1_manuallyselectcells_sc()[[1]])
      geneTxt <- " genes, "
      if(dims[1] == 1) geneTxt <- " gene, "
      
      d$dge_info_manuallyselectcells_sc <- paste0("Data size: ", dims[1], geneTxt, dims[2], " samples")
    }
    h4(d$dge_info_manuallyselectcells_sc, style = "margin-top: 0px;")
  })

  # SC-DGE-MSC - updated dge messages on manually selected cells
  observeEvent(input$dge_filter_manuallyselectcells_sc, {
    if(is.null(heattran1_manuallyselectcells_sc()[[1]])) {
      d$dge_info_manuallyselectcells_sc <- "No DE genes."
    } else {
      dims = dim(heattran1_manuallyselectcells_sc()[[1]])
      geneTxt <- " genes, "
      if(dims[1] == 1) geneTxt <- " gene, "
      
      d$dge_info_manuallyselectcells_sc <- paste0("Data size: ", dims[1], geneTxt, dims[2], " samples")
    }
  }, ignoreNULL = T)

  # SC-DGE-MSC - warning messages for dge for manually selected cells
  observeEvent(d$n_genes_manuallyselectcells_sc, {
    if(d$n_genes_manuallyselectcells_sc != "" & suppressWarnings(is.na(as.integer(na.omit(d$n_genes_manuallyselectcells_sc))))) {
      if(!is.null(input$all_genes_manuallyselectcells_sc)) {
        if(!input$all_genes_manuallyselectcells_sc) {
          d$gse_warning_msg <- "Please provide an integer in 'Use top n genes'."
        }
      }
    } else {
      d$gse_warning_msg <- ""
    }
    if(is.null(heattran1_manuallyselectcells_sc()[[1]])) {
      d$dge_info_manuallyselectcells_sc <- "No DE genes."
    } else {
      dims = dim(heattran1_manuallyselectcells_sc()[[1]])
      geneTxt <- " genes, "
      if(dims[1] == 1) geneTxt <- " gene, "
      
      d$dge_info_manuallyselectcells_sc <- paste0("Data size: ", dims[1], geneTxt, dims[2], " samples")
    }
  }, ignoreNULL = T)

  # SC-DGE-MSC - warning messages for dge for manually selected cells
  observeEvent(input$all_genes_manuallyselectcells_sc, {
    if(d$n_genes_manuallyselectcells_sc != "" & suppressWarnings(is.na(as.integer(na.omit(d$n_genes_manuallyselectcells_sc))))) {
      if(!input$all_genes_manuallyselectcells_sc) {
        d$gse_warning_msg <- "Please provide an integer in 'Use top n genes'."
      } else {
        d$gse_warning_msg <- ""
      }
    } else {
      d$gse_warning_msg <- ""
    }
    if(is.null(heattran1_manuallyselectcells_sc()[[1]])) {
      d$dge_info_manuallyselectcells_sc <- "No DE genes."
    } else {
      dims = dim(heattran1_manuallyselectcells_sc()[[1]])
      geneTxt <- " genes, "
      if(dims[1] == 1) geneTxt <- " gene, "
      
      d$dge_info_manuallyselectcells_sc <- paste0("Data size: ", dims[1], geneTxt, dims[2], " samples")
    }
    if(input$all_genes_manuallyselectcells_sc) {
      d$gse_warning_msg <- ""
    }
  }, ignoreNULL = T)

  observe({
    if(is.null(input$enrichRLib_manuallyselectcells)) {
      d$enrichRLib_manuallyselectcells <- scan("./data/enrichrlibrary.txt", what = "character", quiet = T)[1:6]
    } else {
      d$enrichRLib_manuallyselectcells <- input$enrichRLib_manuallyselectcells
    }
  })
  
  # SC-DGE-MSC - action button to download gse data
  output$sc_manuallyselectcells_gseActionButton <- renderUI({
    req(input$calcDE, d$enrichRLib_manuallyselectcells, input$comps_manuallyselectcells_sc)
    if(is.null(d$enrichRLib_manuallyselectcells)) return()
    if(length(d$enrichRLib_manuallyselectcells) == 0) return()
    if(d$enrichRLib_manuallyselectcells[1] == "") return()
    if(is.null(heattran1_manuallyselectcells_sc())) return()
    if(d$calcTextVal != "") return()
    actionButton("sc_manuallyselectcells_gseActionButton", label = "Download GSE data")
  })

  # SC-DGE-MSC - warning messages after running gse
  observeEvent(input$sc_manuallyselectcells_gseActionButton, {
    if(!input$all_genes_manuallyselectcells_sc & d$n_genes_manuallyselectcells_sc != "" &
       suppressWarnings(is.na(as.integer(na.omit(d$n_genes_manuallyselectcells_sc))))) {
      d$gse_warning_msg <- "Please provide an integer in 'Use top n genes'."
    } else if(is.null(d$enrichRLib_manuallyselectcells)) {
      d$gse_warning_msg <- "Please provide select at least one EnrichR library."
    } else if (length(d$enrichRLib_manuallyselectcells) == 0 |
               d$enrichRLib_manuallyselectcells[1] == "") {
      d$gse_warning_msg <- "Please provide select at least one EnrichR library."
    } else if(is.null(sc_manuallyselectcells_func_enrich())) {
      if(d$gse_warning_msg == "") {
        d$gse_warning_msg <- "No enrichment."
      }
    } else {
      d$gse_warning_msg <- ""
      shinyjs::runjs("$('#sc_manuallyselectcells_downenrData')[0].click();")
    }
  })

  # SC-DGE-MSC - select enrichr libraries
  output$enrichRLib_manuallyselectcells <- renderUI({
    req(input$calcDE, input$comps_manuallyselectcells_sc)
    if(d$calcTextVal != "") return()
    enrichRLibs <- scan("./data/enrichrlibrary.txt", what = "character", quiet = T)
    selectInput(
      inputId = "enrichRLib_manuallyselectcells", 
      label = "EnrichR libraries",
      choices = enrichRLibs, 
      selected = enrichRLibs[1:6], 
      multiple = T,
      width = "600px"
    )
  })

  # SC-DGE-MSC - download handler for gse data
  output$sc_manuallyselectcells_downenrData <- downloadHandler(
    filename = function() {
      paste0(sub("^Comp:(.+)$", replacement = "\\1", x = input$comps_manuallyselectcells_sc),
             "_ManuallySelected_GSE.csv")
    },
    content = function(file) {
      enrData <- sc_manuallyselectcells_func_enrich() %>%
        mutate_if(is.numeric, round, digits = 4)

      write.csv(
        enrData,
        file,
        row.names = FALSE
      )
    }
  )

  output$manuallyselectcells_downloadGSE_ui <- renderUI({
    # req(input$comps_manuallyselectcells_sc, clustList())
    req(clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()
    fluidPage(
      style = "border: 1px solid #eee; padding: 10px;",
      uiOutput("manuallyselectcells_downloadGSE_header"),
      fluidRow(
        fluidPage(
          div(
            style = "display: inline-block;",
            uiOutput("sc_manuallyselectcells_gene_list_filtering")
          ),
          div(
            style = "display: inline-block; width: 200px; margin-left: 20px;",
            fluidPage(
              uiOutput("top_n_genes_manuallyselectcells_sc"),
              uiOutput("use_all_genes_manuallyselectcells_sc")
            )
          ),
          div(
            style = "display: inline-block;",
            uiOutput("enrichRLib_manuallyselectcells")
          )
        )
      ),
      fluidRow(
        fluidPage(
          div(style = "display: inline-block;", uiOutput("dge_info_manuallyselectcells_sc")),
          div(
            style = "display: inline-block; margin-left: 30px;",
            uiOutput("sc_manuallyselectcells_gseActionButton")
          )
        )
      ),
      div(textOutput("gse_warning_manuallyselectcells_sc"), style = "color: red;"),
      downloadLink("sc_manuallyselectcells_downenrData", label = ""),
      br()
    )
  })
  
  # SC-DGE-MSC - header for downloading results
  output$manuallyselectcells_downloadDGE_header <- renderUI({
    req(input$comps_manuallyselectcells_sc, clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()
    h3("Download differential gene expression results", style = "margin-top: 0px;")
  })

  # SC-DGE-MSC - select previous comparisons from manually selected cells DGE
  output$comps_manuallyselectcells_sc <- renderUI({
    req(clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()
    comps <- sub("^(Comp:.+):.+$", replacement = "\\1", x = comps)
    selectInput("comps_manuallyselectcells_sc", label = "Comparison", choices = comps, selected = comps[1])
  })
  
  observeEvent(input$comps_manuallyselectcells_sc, {
    d$res <- input$comps_manuallyselectcells_sc
  })

  # SC-DGE-MSC - download button for manually select cells
  output$sc_manuallyselectcells_downloaddge <- renderUI({
    # req(input$calcDE)
    req(input$calcDE, input$comps_manuallyselectcells_sc)
    if(d$calcTextVal != "") return()
    downloadButton("sc_manuallyselectcells_downloaddgebutton", "Download DGE data")
  })

  # SC-DGE-MSC - select abs log2fc threshold
  output$manuallyselectcells_lfc <- renderUI({
    req(input$comps_manuallyselectcells_sc, clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()

    numericInput("manuallyselectcells_lfc", label = "Abs. log2 fold-change",
                 value = 0.25, min = 0)
  })

  # SC-DGE-MSC - adj p-val threshold
  output$manuallyselectcells_padj <- renderUI({
    req(input$comps_manuallyselectcells_sc, clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()

    numericInput("manuallyselectcells_padj", label = "Adj. p-value",
                 value = 0.05, min = 0, max = 1)
  })

  # SC-DGE-MSC - download handler dge
  output$sc_manuallyselectcells_downloaddgebutton <- downloadHandler(
    filename = function() {
      paste0(sub("^Comp:(.+)$", replacement = "\\1", x = input$comps_manuallyselectcells_sc),
             "_ManuallySelected_DGE.csv")
    },
    content = function(file) {
      a <- ClustGeneStats(d$SCV[[input$comps_manuallyselectcells_sc]])[["Set A"]][,c("MGE","DR","MDGE")]  %>%
        rownames_to_column(., var = "Gene")
      names(a)[2:4] <- c("MGE_Set A", "DR_Set A", "MDGE_Set A")
      b <- ClustGeneStats(d$SCV[[input$comps_manuallyselectcells_sc]])[["Set B"]][,c("MGE","DR","MDGE")] %>%
        rownames_to_column(., var = "Gene")
      names(b)[2:4] <- c("MGE_Set B", "DR_Set B", "MDGE_Set B")
      temp <- c("Set A-Set B", "Set B-Set A")
      tempName <- temp[temp %in% names(DEcombn(d$SCV[[input$comps_manuallyselectcells_sc]]))]
      out_table <- DEcombn(d$SCV[[input$comps_manuallyselectcells_sc]])[[tempName]] %>%
        rownames_to_column(., var = "Gene")
      # out_table <- out_table[abs(as.numeric(out_table$logGER)) > 0.25 & as.numeric(out_table$pVal) <= 0.001, ]
      out_table <- out_table[abs(as.numeric(out_table$logGER)) > input$manuallyselectcells_lfc &
                               as.numeric(out_table$FDR) <= input$manuallyselectcells_padj, ]
      all <- out_table %>%
        left_join(., a, by = "Gene") %>%
        left_join(., b, by = "Gene")  %>%
        mutate_if(., is.numeric, round, 3)
      names(all) <- gsub("logGER", "logFC", names(all))
      names(all) <- gsub("Wstat", "Wilcoxon statistic", names(all))
      names(all) <- gsub("pVal", "P-value", names(all))
      names(all) <- gsub("dDR", "Difference in detection rate", names(all))
      names(all) <- gsub("RNA\\.", "Avg exp cluster ", names(all))

      nam <- names(all)[grepl("Avg exp cluster ", names(all))]
      nam <- nam[!grepl("Set A", nam)]
      nam <- nam[!grepl("Set B", nam)]
      all <- all %>%
        dplyr::select(-nam)
      for(i in 1:ncol(all)) all[, i] = as.character(all[, i])
      write.csv(
        all,
        file,
        row.names = FALSE,
        quote = T
      )
    }
  )

  output$manuallyselectcells_downloadDGE_ui <- renderUI({
    req(clustList())
    comps <- grep("^Comp:", x = clustList(), value = T)
    if(length(comps) == 0) return()
    
    fluidPage(
      style = "border: 1px solid #eee; padding: 10px;",
      uiOutput("manuallyselectcells_downloadDGE_header"),
      div(
        style = "display: inline-block; width: 200px;", 
        uiOutput("comps_manuallyselectcells_sc")
      ),
      div(
        style = "display: inline-block; width: 200px; margin-left: 10px; vertical-align: top;", 
        uiOutput("manuallyselectcells_lfc")
      ),
      div(
        style = "display: inline-block; width: 200px; margin-left: 10px; vertical-align: top;", 
        uiOutput("manuallyselectcells_padj")
      ),
      div(
        style = "display: inline-block; margin-left: 10px; vertical-align: top; margin-top: 23px;", 
        uiOutput("sc_manuallyselectcells_downloaddge")
      )
    )
  })
  
  # SC-DGE-MSC - action button to add to cell selection
  output$plusFilt <- renderUI({
    req(d$res, d$inD, d$SCV)
    actionButton("plusFilt", "Add", icon("plus"),
                 style = "color: #fff; border-color: #00bfff; background:#00bfff")
  })

  # SC-DGE-MSC - action button to remove cell selection
  output$minusFilt <- renderUI({
    req(d$res, d$inD, d$SCV)
    actionButton("minusFilt", "Remove", icon("minus"),
                 style = "color: #fff; border-color: #00bfff")
  })

  ######### CLUSTERING TAB CODE: SCRNA-SEQ#####################

  # SC-DGE-CLUS - actionbutton for clustering
  output$goclust_sc <- renderUI({
    
    actionButton(
      "goclust_sc",
      "Launch clustering analysis",
      icon = icon("space-shuttle")
    )
  })

  # SC-DGE-CLUS - warning for clustering
  output$wgcna_warning_msg_sc <- renderText(d$wgcna_warning_sc)

  # SC-DGE-CLUS - set newCluster in object d as FALSE
  observeEvent(input$goclust_sc, {
    d$newCluster <- F
  })

  # SC-DGE-CLUS - header for clustering analysis
  output$headclust_sc <- renderUI({
    
    h4("Clustering analysis")
  })

  # SC-DGE-CLUS - set clustering cut-off
  observeEvent(SubmitData$data_type, {
    output$clustvarnumber_sc <- renderUI({
      
      textInput(
        inputId = "clustvarnumber_sc",
        label = div(style = "width: 200px;", "Top variable genes"),
        value = 500
      )
    })
  })

  # SC-DGE-CLUS - choose clustering algorithm
  observeEvent(SubmitData$data_type, {
    output$clustalg_sc <- renderUI({
      
      selectInput(
        inputId = "clustalg_sc",
        label = div(style = "width: 200px;", "Clustering algorithm"),
        choices = c(
          "WGCNA" = "wgcna",
          "K-medoids" = "kmed"
        ),
        selected = "WGCNA"
      )
    })
  })

  # SC-DGE-CLUS - set minimum gene module size
  output$min_module_size_sc <- renderUI({
    req(SubmitData$data_type, input$clustalg_sc)
    if(is.null(SubmitData$data_type)) return()
    if(SubmitData$data_type != "Single-cell") return()
    if(input$clustalg_sc != "wgcna") return()
    div(
      style = "width: 100px; margin-left: 40px;",
      numericInput(
        inputId = "min_module_size_sc",
        label = div(style = "width: 180px;", "Min. module size"),
        value = 30,
        min = 8,
        step = 1
      )
    )
  })

  # SC-DGE-CLUS - update min module size
  observeEvent(input$min_module_size_sc, {
    if(is.null(input$min_module_size_sc)) return()
    if(is.numeric(input$min_module_size_sc)) {
      if(input$min_module_size_sc >= 8) return()
    }
    updateNumericInput(
      session = session,
      inputId = "min_module_size_sc",
      value = 30
    )
  })

  # SC-DGE-CLUS - reactive exp for running clustering
  clustout_sc <- eventReactive(input$goclust_sc, {
    num <- input$clustvarnumber_sc

    if (d$resType == "seurat_res") {
        res <- input$overall_res
        res <- gsub(".*_|\\:.*", "", res)
        if(!(res %in% names(seurat_only()))) return(NULL)
        seur <- seurat_only()[grepl(res, names(seurat_only()))]
        seur <- seur[[1]]
        cts <- as.matrix(GetAssayData(object = seur, slot = "data"))
        tran <- cts
    } else {
      seur <- seurat_only()
      cts <- as.matrix(GetAssayData(object = seur, slot = "data"))
      tran <- cts
    }
    topID <- order(rowVars(cts), decreasing = TRUE)
    cts.var <- tran[topID, ]
    dds_mat <- cts.var[1:num, ]
    if (input$clustalg_sc == "wgcna") {
      withProgress(message = "Running WGCNA...", value = 0, {
        incProgress(1/3)
        enableWGCNAThreads(max(1, parallel::detectCores()-2))
        dds_mat <- as.matrix(dds_mat)
        gene.names <- sort(rownames(dds_mat))
        datExpr <- t(dds_mat)
        # Run this to check if there are gene outliers
        gsg <- goodSamplesGenes(datExpr, verbose = 3)
        gsg$allOK
        # Create an object called "datTraits" that contains your
        # trait data
        datTraits <- seur@meta.data
        datTraits <- datTraits[, grep("nCount_RNA|nFeature_RNA|silWidth|Sample", x = colnames(datTraits), invert = T), drop = F]
        nLevels <- apply(datTraits, 2, function(col) length(unique(col)))
        
        # removed the max bc it makes the subsequent code throw errors
        # if for example we have more than 12 seurat clusters
        datTraits <- datTraits[, nLevels > 1,
                               #& nLevels <= 12, 
                               drop = F]

        # Form a data frame analogous to expression data that will
        # hold the clinical traits.
        # should return TRUE if datasets align correctly, otherwise your
        # names are out of order
        table(rownames(datTraits) == rownames(datExpr))
        # calculates the whole network connectivity
        A <- adjacency(t(datExpr), type="signed")
        k <- as.numeric(apply(A, 2, sum))-1 # standardized connectivity
        Z.k <- scale(k)
        thresholdZ.k <- -2.5 # often -2.5
        outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")
        sampleTree <- flashClust(stats::as.dist(1-A), method = "average")
        # convert traits to color: where red indicates high values
        datTraits$seurat_clusters <- as.numeric(datTraits$seurat_clusters)
        
        if (ncol(datTraits) > 2) {
          datTraits[c(1, 3:ncol(datTraits))] <- lapply(c(1, 3:ncol(datTraits)), function(x) as.character(names(datTraits)[x]))
        } else {
          datTraits$orig.ident <- as.character(as.factor(datTraits$orig.ident))
        }

        # MAY NEED TO UPDATE THIS WITH RAINBOW FOR FACTORS > 12
        traitColors <- datTraits
        for(i in 1:ncol(traitColors)) {
          if(i %% 2 == 0) {
            traitColors[, i] <- labels2colors(traitColors[, i],
                                              colorSeq = brewer.pal(12, "Set3"))
          } else {
            traitColors[, i] <- labels2colors(traitColors[, i],
                                              colorSeq = brewer.pal(12, "Paired"))
          }
        }

        dimnames(traitColors)[[2]] <- paste(names(datTraits))
        datColors <- data.frame(outlier = outlierColor, traitColors)
        incProgress(1/3)
        # TOM analysis - (computationally expensive)
        softPower <- 18
        adjacency <- adjacency(datExpr, power = softPower, type = "signed") #specify network type
        TOM <- TOMsimilarity(adjacency, TOMType = "signed") # specify network type
        dissTOM <- 1 - TOM
        geneTree <- flashClust(stats::as.dist(dissTOM), method="average")
        # sets the minimum number of genes to cluster into a module
        minModuleSize <- input$min_module_size_sc
        dynamicMods <- cutreeDynamic(
          dendro = geneTree,
          distM = dissTOM,
          deepSplit = 2,
          pamRespectsDendro = FALSE,
          minClusterSize = minModuleSize
        )
        if(length(unique(dynamicMods)) == 1) {
          if(minModuleSize == 8) return()
          minModuleSizes <- (minModuleSize-1):8
          for(i in minModuleSizes) {
            dynamicMods <- cutreeDynamic(
              dendro = geneTree,
              distM = dissTOM,
              deepSplit = 2,
              pamRespectsDendro = FALSE,
              minClusterSize = i
            )
            if(length(unique(dynamicMods)) > 1) {
              d$wgcna_minModuleSizeTried_sc <- i
              break
            }
          }
          if(length(unique(dynamicMods)) == 1) {
            d$wgcna_minModuleSizeTried_sc <- 8
            return()
          }
        } else {
          d$wgcna_minModuleSizeTried_sc <- NULL
        }
        dynamicColors <- labels2colors(dynamicMods, colorSeq = brewer.pal(12, "Set3"))
        MEList <- moduleEigengenes(
          datExpr,
          colors = dynamicColors,
          softPower = 18
        )
        MEs <- MEList$eigengenes
        MEDiss <- 1 - cor(MEs)
        METree <- flashClust(stats::as.dist(MEDiss), method = "average")
        # set a threhold for merging modules. In this example we are
        # not merging so MEDissThres=0.0
        MEDissThres <- 0.0
        merge <- mergeCloseModules(
          datExpr,
          dynamicColors,
          cutHeight = MEDissThres,
          verbose = 3
        )
        mergedColors <- merge$colors
        mergedMEs <- merge$newMEs
        # set the diagonal of the dissimilarity to NA
        diag(dissTOM) = NA;
        # export modules to data frame
        module_colors= setdiff(unique(dynamicColors), "grey")
        modlist <- list()
        for (color in module_colors){
          module = gene.names[which(dynamicColors == color)]
          modlist[[color]] <- list(
            gene = module
          )
        }
        moddf <- data.frame(unlist(modlist))
        moddf$module <- gsub("\\..*", "", row.names(moddf))
        colnames(moddf)[1] <- "gene"
        rownames(moddf) <- seq_len(nrow(moddf))
        moddf$gene <- as.character(moddf$gene)
        moddf$module <- as.factor(moddf$module)
        sampleDF <- data.frame(
          sample = sampleTree$labels,
          outlier = datColors$outlier
        )
        disableWGCNAThreads()
        incProgress(1/3)
        return(
          list(
            sampleTree,
            datColors,
            geneTree,
            dynamicColors,
            mergedColors,
            dissTOM,
            moddf,
            sampleDF
          )
        )
      })
    } else if (input$clustalg_sc == "kmed") {
      withProgress(message = "Running K-medoids...", value = 0, {
        incProgress(1/2)
        num <- as.matrix(dds_mat)
        mrwdist <- distNumeric(num, num, method = "mrw")
        # Detect cores of machine
        nclust <- parallel::detectCores() - 1
        message("Using ", nclust, " threads...")
        result <- fastkmed(mrwdist, ncluster = nclust, iterate = 50)
        # a simple and fast k-medoids function for bootstrap evaluation
        parkboot <- function(x, nclust) {
          res <- fastkmed(x, nclust, iterate = 50)
          return(res$cluster)
        }
        fastkmedboot <- clustboot(mrwdist, nclust = nclust, parkboot, nboot = 50)
        # consensus matrix
        wardorder <- function(x, nclust) {
          res <- fastcluster::hclust(x, method = "ward.D2")
          member <- cutree(res, nclust)
          return(member)
        }
        consensusfastkmed <- consensusmatrix(fastkmedboot, nclust = nclust, wardorder)
        # data frame generation
        output <- data.frame(
          gene_id = rownames(num),
          cluster = result$cluster
        )
        rownames(output) <- seq_len(nrow(output))
        output$gene_id <- as.character(output$gene_id)
        output$cluster <- as.factor(output$cluster)
        return(
          list(
            consensusfastkmed,
            output
          )
        )
        incProgress(2/2)
      })
    }
  })

  # SC-DGE-CLUS - output DT for custom annotations for wgcna
  observeEvent(input$goclust_sc, {
    if(input$clustalg_sc != "wgcna") {
      updateTabsetPanel(inputId = "sc_clustering_tabsetPanel", selected = "sc_clustPlotW01")
      d$sc_showClusteringLinks <- T
      return()
    }

    if(is.null(clustout_sc())) {
      d$wgcna_warning_sc <- "No gene modules found after setting 'Min. module size' to 8."
      return()
    }
    if(!is.null(d$wgcna_minModuleSizeTried_sc)) {
      d$wgcna_warning_sc <- paste0("Gene modules found after lowering 'Min. module size' to ",
                                   d$wgcna_minModuleSizeTried_sc, ".")
    } else {
      d$wgcna_warning_sc <- NULL
    }

    datColors <- clustout_sc()[[2]]
    
    if (d$resType == "seurat_res") {
      res <- input$overall_res
      res <- gsub(".*_|\\:.*", "", res)
      if(!(res %in% names(seurat_only()))) return(NULL)
      seur <- seurat_only()[grepl(res, names(seurat_only()))]
      seur <- seur[[1]]
    } else {
      seur <- seurat_only()

    }
    coldata <- seur@meta.data
    rowText <- coldata[clustout_sc()[[1]]$labels,
                       colnames(coldata)[colnames(coldata) %in% colnames(datColors)]]
    headerCallback <- c(
      "function(thead, data, start, end, display){",
      "  $('th', thead).css('border-bottom', '1px solid #ddd');",
      "}"
    )

    for(i in 1:ncol(rowText)) rowText[, i] <- as.character(rowText[, i])
    datColors <- datColors[, colnames(rowText)]
    rownames(datColors) <- rownames(rowText)
    for(colname in colnames(rowText)) {
      local({
        df <- data.frame(Group = unique(rowText[,colname]),
                         Dummy = as.character(letters[1:length(unique(rowText[,colname]))]),
                         Color = datColors[match(unique(rowText[,colname]),
                                                 table = rowText[,colname]),
                                           colname])
        df <- df[order(df$Group), ]
        outputName <- paste0(colname, "_legend_table_sc")
        tableTitle <- colname
        output[[outputName]] <- DT::renderDataTable({
          DT::datatable(
            df[,c("Group","Dummy")],
            rownames = F,
            colnames = rep.int("", times = 2),
            escape = F, selection = "none",
            caption = htmltools::tags$caption(
              style = "text-align: center; color: black; margin-bottom: -1.5em;",
              tableTitle
            ),
            callback = htmlwidgets::JS("$('table.dataTable.no-footer').css('border-bottom', '1.5px solid #ddd');"),
            options = list(lengthChange = F, bFilter = F, bInfo = F, bPaginate = F,
                           ordering = F,
                           headerCallback = htmlwidgets::JS(headerCallback)),
            class = "cell-border") %>%
            formatStyle(columns = c("Group","Dummy"),
                        color = styleEqual(df$Dummy, df$Color),
                        backgroundColor = styleEqual(df$Dummy, df$Color),
                        fontSize = "12px",
                        lineHeight = "50%")
        })
      })
    }
    d$sc_showClusteringLinks <- T
    updateTabsetPanel(inputId = "sc_clustering_tabsetPanel", selected = "sc_clustPlotW02")
  })

  # SC-DGE-CLUS - header for sample dendrogram
  output$headclustplotW01_sc <- renderUI({
    req(input$goclust_sc, clustout_sc(), !d$newCluster)
    isolate({
      if (input$goclust_sc == 0) {
        return()
      } else if (is.null(clustout_sc())) {
        return()
      } else if (input$clustalg_sc == "wgcna") {
        return()
      } else if (input$clustalg_sc == "kmed") {
        if (length(clustout_sc()) == 2) {
          h4(strong("K-Medoids - consensus matrix heatmap"))
        }
      }
    })
  })

  # SC-DGE-CLUS - sample dendrogram plot
  output$clustplotW01_sc <- renderPlot({
    req(input$goclust_sc, clustout_sc())
    if(input$goclust_sc == 0) return()
    if(length(clustout_sc()) != 2) return()

    # K-medoids
    clustheatmap(clustout_sc()[[1]], title = "K-medoids consensus matrix heatmap")
  })
  
  # Use renderUI to avoid blank area for sample dendrogram
  output$clustplotW01_sc_ui <- renderUI({
    req(input$goclust_sc, clustout_sc())
    if(input$goclust_sc == 0) return()
    if(length(clustout_sc()) != 2) return()
    fluidRow(column(12, align = "left", plotOutput("clustplotW01_sc", width = 900, height = 500)))
  })

  # SC-DGE-CLUS - download button for sample dendrogram (PNG)
  output$downloadclustplotW01png_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    isolate({
      if (input$clustalg_sc == "kmed") {
        if (length(clustout_sc()) == 2) {
          downloadButton(
            "downloadclustplotK01pngimg_sc",
            "Download plot (PNG)"
          )
        }
      }
    })
  })

  # SC-DGE-CLUS - download handler for consensus matrix
  output$downloadclustplotK01pngimg_sc <- downloadHandler(
    filename = function() {
      paste("kmed-consensus-matrix-sc.png")
    },
    content = function(file) {
      req(clustout_sc())
      sample <- clustout_sc()[[1]]
      png(file, width = 800, height = 600)
      output <- clustheatmap(sample, "K-medoids consensus matrix heatmap") + 
        theme(text = element_text(family = "noto-sans-jp"))
      print(output)
      dev.off()
    }
  )

  # SC-DGE-CLUS - download button sample dendrogram (PDF)
  output$downloadclustplotW01pdf_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    isolate({
      if (input$clustalg_sc == "kmed") {
        if (length(clustout_sc()) == 2) {
          downloadButton(
            "downloadclustplotK01pdfimg_sc",
            "Download plot (PDF)"
          )
        }
      }
    })
  })

  # SC-DGE-CLUS - download handler for consensus matrix
  output$downloadclustplotK01pdfimg_sc <- downloadHandler(
    filename = function() {
      paste("kmed-consensus-matrix-sc.pdf")
    },
    content = function(file) {
      CairoPDF(file, width = 12, height = 6.5)
      consensusfastkmed <- clustout_sc()[[1]]
      p <- clustheatmap(
        consensusfastkmed,
        "K-medoids consensus matrix heatmap"
      ) + 
      theme(text = element_text(family = "noto-sans-jp"))
      print(p)
      dev.off()
    }
  )

  # SC-DGE-CLUS - header for gene dendrogram
  output$headclustplotW02_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    isolate({
      if (input$clustalg_sc == "wgcna") {
        if (length(clustout_sc()) == 8) {
          h4(strong("WGCNA - gene dendrogram"))
        }
      } else {
        return()
      }
    })
  })

  # SC-DGE-CLUS - gene dendrogram
  output$clustplotW02_sc <- renderPlot({
    validate(need(input$goclust_sc, ""))
    req(input$goclust_sc, clustout_sc(), !d$newCluster)
    isolate({
      withProgress(message = "Creating gene dendrogram...", value = 0, {
        incProgress(1/2)
        geneTree <- clustout_sc()[[3]]
        dynamicColors <- clustout_sc()[[4]]
        mergedColors <- clustout_sc()[[5]]
        datColors <- data.frame("Dynamic tree cut" = dynamicColors, check.names = F)
        incProgress(1/2)
        GGDend(geneTree,
               colorDF = datColors,
               leafLabels = F,
               touchZero = F,
               plotTitle = "Gene dendrogram and module colors")
      })
    })
  })

  # SC-DGE-CLUS - download button gene dendrogram (PNG)
  output$downloadclustplotW02png_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    if(length(clustout_sc()) != 8) return()
    downloadButton("downloadclustplotW02pngimg_sc", label = "Download plot (PNG)")
  })

  # SC-DGE-CLUS - - download file gene dendrogram (PNG)
  output$downloadclustplotW02pngimg_sc <- downloadHandler(
    filename = function() {
      paste("wgcna-gene-dendrogram-sc.png")
    },
    content = function(file) {
      png(file, width = 800, height = 400)
      geneTree <- clustout_sc()[[3]]
      dynamicColors <- clustout_sc()[[4]]
      datColors <- data.frame("Dynamic tree cut" = dynamicColors, check.names = F)
      g <- GGDend(geneTree,
                  colorDF = datColors,
                  leafLabels = F,
                  touchZero = F,
                  plotTitle = "Gene dendrogram and module colors",
                  plotTitleSize = 30)
      print(g)
      dev.off()
      #ggsave(file, plot = g, width = 12, height = 8)
    }
  )

  # SC-DGE-CLUS - download button gene dendrogram (PDF)
  output$downloadclustplotW02pdf_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    if(length(clustout_sc()) != 8) return()
    downloadButton("downloadclustplotW02pdfimg_sc", label = "Download plot (PDF)")
  })

  # SC-DGE-CLUS - download file gene dendrogram (PDF)
  output$downloadclustplotW02pdfimg_sc <- downloadHandler(
    filename = function() {
      paste("wgcna-gene-dendrogram-sc.pdf")
    },
    content = function(file) {
      CairoPDF(file, width = 12, height = 6.5)
      geneTree <- clustout_sc()[[3]]
      dynamicColors <- clustout_sc()[[4]]
      datColors <- data.frame("Dynamic tree cut" = dynamicColors, check.names = F)
      g <- GGDend(geneTree,
                  colorDF = datColors,
                  leafLabels = F,
                  touchZero = F,
                  plotTitle = "Gene dendrogram and module colors",
                  plotTitleSize = 30)
      print(g)
      dev.off()
      #ggsave(file, plot = g, width = 12, height = 8)
    }
  )

  # SC-DGE-CLUS - header for sample dendrogram
  output$headclustplotW03_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    isolate({
      if (input$clustalg_sc == "wgcna") {
        if (length(clustout_sc()) == 8) {
          h4(strong("WGCNA - topological overlap matrix"))
        }
      } else {
        return()
      }
    })
  })

  # SC-DGE-CLUS - TOM plot
  output$clustplotW03_sc <- renderPlot({
    validate(need(input$goclust_sc, ""))
    req(input$goclust_sc, clustout_sc(), !d$newCluster)
    isolate({
      if(input$goclust_sc == 0) return()
      if(input$clustalg_sc != "wgcna") return()
      if(length(clustout_sc()) != 8) return()
      withProgress(message = "Creating TOM plot...", value = 0, {
        incProgress(1/2)
        geneTree <- clustout_sc()[[3]]
        dynamicColors <- clustout_sc()[[4]]
        dissTOM <- clustout_sc()[[6]]
        incProgress(1/2)
        GGTom(geneTree, colors = dynamicColors, distMat = dissTOM, plotTitle = NULL)
      })
    })
  })

  # SC-DGE-CLUS - download button TOM plot (PNG)
  output$downloadclustplotW03png_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    if(length(clustout_sc()) != 8) return()
    downloadButton("downloadclustplotW03pngimg_sc", label = "Download plot (PNG)")
  })

  # SC-DGE-CLUS - download file TOM plot (PNG)
  output$downloadclustplotW03pngimg_sc <- downloadHandler(
    filename = function() {
      paste("wgcna-tom-plot-sc.png")
    },
    content = function(file) {
      png(file, width = 600, height = 400)
      geneTree <- clustout_sc()[[3]]
      dynamicColors <- clustout_sc()[[4]]
      dissTOM <- clustout_sc()[[6]]
      g <- GGTom(geneTree, colors = dynamicColors,
                 distMat = dissTOM)
      print(g)
      dev.off()
    }
  )

  # SC-DGE-CLUS - download button TOM plot (PDF)
  output$downloadclustplotW03pdf_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    if(length(clustout_sc()) != 8) return()
    downloadButton("downloadclustplotW03pdfimg_sc", label = "Download plot (PDF)")
  })

  # SC-DGE-CLUS - download file TOM plot (PDF)
  output$downloadclustplotW03pdfimg_sc <- downloadHandler(
    filename = function() {
      paste("wgcna-tom-plot-sc.pdf")
    },
    content = function(file) {
      CairoPDF(file, width = 8, height = 8)
      geneTree <- clustout_sc()[[3]]
      dynamicColors <- clustout_sc()[[4]]
      dissTOM <- clustout_sc()[[6]]
      g <- GGTom(geneTree, colors = dynamicColors,
                 distMat = dissTOM)
      print(g)
      dev.off()
    }
  )

  # SC-DGE-CLUS - header for download gene modules
  output$headclustmoddown_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    isolate({
      if (input$clustalg_sc == "wgcna") {
        if (length(clustout_sc()) == 8) {
          h4(strong("WGCNA - download modules"))
        }
      } else {
        return()
      }
    })
  })

  # SC-DGE-CLUS - download button for gene modules
  output$downloadclustmod_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    if(length(clustout_sc()) != 8) return()
    downloadButton("downloadclustmod2_sc", label = "Download gene modules (CSV)")
  })

  # SC-DGE-CLUS - download file for gene modules
  output$downloadclustmod2_sc <- downloadHandler(
    filename = function() {
      paste("wgcna-gene-modules-sc.csv")
    },
    content = function(file) {
      moddf <- clustout_sc()[[7]]
      colorPalette <- c("teal","lightyellow","lavender",
                        "salmon","blue","orange","green",
                        "lightpink","gray","violet","mint","yellow")
      names(colorPalette) <- brewer.pal(12, "Set3")
      moddf$module <- colorPalette[moddf$module]
      write.csv(moddf, file, row.names = FALSE)
    }
  )

  # SC-DGE-CLUS - download button for sample module
  output$downloadclustsample_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    isolate({
      if (input$clustalg_sc == "wgcna") {
        if (length(clustout_sc()) == 8) {
          downloadButton(
            "downloadclustsample2_sc",
            "Download sample modules (CSV)"
          )
        }
      }
    })
  })

  # SC-DGE-CLUS - download file for sample module
  output$downloadclustsample2_sc <- downloadHandler(
    filename = function() {
      paste("wgcna-sample-modules-sc.csv")
    },
    content = function(file) {
      sampledf <- clustout_sc()[[8]]
      metaDF <- as.data.frame(colData(ddsout()[[1]]))
      metaDF <- metaDF[sampledf$sample, grep("nCount_RNA|nFeature_RNA|silWidth", x = colnames(metaDF), invert = T)]
      sampledf <- cbind(sampledf, metaDF)
      sampledf <- sampledf[, colnames(sampledf) != "outlier"]
      write.csv(sampledf, file, row.names = FALSE)
    }
  )

  # SC-DGE-CLUS - download button for gene module
  output$downloadclustmodK_sc <- renderUI({
    req(clustout_sc(), !d$newCluster)
    if(length(clustout_sc()) != 2) return()
    downloadButton("downloadclustmodK2_sc", label = "Download clusters (CSV)")
  })

  # SC-DGE-CLUS - download file for gene module
  output$downloadclustmodK2_sc <- downloadHandler(
    filename = function() {
      paste("kmed-gene-clusters-sc.csv")
    },
    content = function(file) {
      clustdf <- clustout_sc()[[2]]
      write.csv(clustdf, file, row.names = FALSE, col.names = TRUE)
    }
  )

  # SC-DGE-CLUS - download handler for gene module
  output$downloadclustmodM2_sc <- downloadHandler(
    filename = function() {
      paste("mcl-gene-clusters-sc.csv")
    },
    content = function(file) {
      clustdf <- clustout_sc()[[2]]
      write.csv(clustdf, file, row.names = FALSE, col.names = TRUE)
    }
  )

  ########### COMBINE CLUSTERS OR CELLS: SCRNA-SEQ ONLY ##############

  # SC-DGE-CCC - header for combining clusters
  output$combineClustersHeader <- renderUI({
    req(d$resType, clustList())
    h4("Combine clusters")
  })

  # SC-DGE-CCC - header for cluster combination
  output$clustercombination_header <- renderUI({
    req(d$SCV)
    h3("Cluster combination")
  })

  # SC-DGE-CCC - header for combining cells
  output$clustercells_header <- renderUI({
    req(d$SCV)
    h3("Combine cells")
  })

  # SC-DGE-CCC - create a drop-down for selecting the
  # clusters you would like to combine
  output$clusters_comb <- renderUI({
    req(d$res)
    if (is.null(seurat_only())) {
    } else if (length(seurat_only()) > 0) {
      if (d$resType == "seurat_res") {
        if (!is.null(input$overall_res)) {
          res <- gsub("^RNA_snn_(.*):.*$", replacement = "\\1", x = input$overall_res)
          if(!(res %in% names(seurat_only()))) return(NULL)
          seur <- seurat_only()[[res]]
          clusters <- mixedsort(unique(seur@meta.data[,"seurat_clusters"]))
          selectizeInput("clusters_comb", "Select clusters to combine",
                         choices = clusters, multiple = T)
        }
      } else {
        meta <- ddsout()[[1]]
        meta <- data.frame(colData(meta))
        meta <- mixedsort(unique(meta$seurat_clusters))
        selectizeInput("clusters_comb", "Select clusters to combine",
                       choices = meta, multiple = T)
      }
    }
  })

  # SC-DGE-CCC - re-create df of number of cells by cluster
  # colored by sample with recombined clusters
  recomb_seur <- reactive({
    if(is.null(d$res) || is.null(d$resType) || is.null(seurat_only()) 
       || length(seurat_only()) == 0) return()
    if(d$resType == "seurat_res" && is.null(input$overall_res)) return()
    
    if(d$resType == "seurat_res") {
      res <- gsub(".*_|\\:.*", "", input$overall_res)
      if(!(res %in% names(seurat_only()))) return(NULL)
      seur <- seurat_only()[grepl(res, names(seurat_only()))]
      seur <- seur[[1]]
      if(!is.null(input$clusters_comb)) {
        clus <- input$clusters_comb
        comb <- ifelse(seur@meta.data$seurat_clusters %in% clus, 100, seur@meta.data$seurat_clusters)
        df <- data.frame("vec" = mixedsort(unique(comb)), "new" = seq(0, length(unique(comb)) - 1, by = 1))
        seur$clus_comb <- df$new[match(comb, df$vec)]
        meta <- seur@meta.data
        cols <- names(meta)[grepl("orig|nCount|nFeature|clus_comb", names(meta))]
        meta <- meta[, names(meta) %in% cols]
        meta$Sample <- rownames(meta)
        grp <- names(meta)[grepl("orig.ident|Sample|clus_comb", names(meta))]
        add <- meta %>%
          dplyr::group_by_at(vars(one_of(grp))) %>%
          dplyr::summarise(`Number of Cells` = n())
        meta <- meta %>%
          full_join(., add)
        return(meta)
      }
    } else {
      seur <- seurat_only()
      if(!is.null(input$clusters_comb)) {
        clus <- input$clusters_comb
        comb <- ifelse(seur@meta.data$seurat_clusters %in% clus, 100, seur@meta.data$seurat_clusters)
        df <- data.frame("vec" = mixedsort(unique(comb)), "new" = seq(0, length(unique(comb)) - 1, by = 1))
        seur$clus_comb <- df$new[match(comb, df$vec)]
        meta <- seur@meta.data
        meta$Sample <- rownames(meta)
        grp <- names(meta)[grepl("orig.ident|Sample|clus_comb", names(meta))]
        add <- meta %>%
          dplyr::group_by_at(vars(one_of(grp))) %>%
          dplyr::summarise(`Number of Cells` = n())
        meta <- meta %>%
          full_join(., add)
        return(meta)
      }
    }
  })

  # SC-DGE-CCC - re-create plot of number of cells by cluster
  # for combined clusters
  output$comb_mdFactor_plotly <- renderPlotly({
    req(d$inD, d$MD, d$res, d$resType)
    if(d$resType == "seurat_res") {
      if(!(d$res %in% colnames(d$MD))) return()
      clusters <- d$MD[, d$res]
    } else {
      if(!("seurat_clusters" %in% colnames(d$MD))) return()
      clusters <- d$MD[, "seurat_clusters"]
    }
    cells <- data.frame(orig_cluster = clusters,
                        clus_comb = clusters,
                        Sample = colnames(d$inD),
                        "Number of Cells" = 1,
                        check.names = F)

    if(!is.null(recomb_seur())) {
      if(ncol(recomb_seur()) > 0) {
        res <- "Number of Cells"
        add <- c("orig.ident", "clus_comb", "Sample")
        cols <- c(add, res)
        cells <- recomb_seur() %>%
          dplyr::select_at(vars(one_of(cols))) %>%
          distinct(.)
        colnames(cells)[colnames(cells) == "orig.ident"] <- "orig_cluster"
        cells$orig_cluster <- clusters
      }
    }
    cells$x <- as.factor(cells$clus_comb)
    cells$y <- cells[, grepl("Number of ", names(cells))]
    cells <- cells %>%
      dplyr::group_by(x) %>%
      dplyr::mutate(count = n()) %>%
      distinct(.)
    len <- length(unique(cells$x))
    gg_color <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cells <- cells %>% arrange(clus_comb)
    cells$clus_comb <- factor(cells$clus_comb, levels = mixedsort(unique(cells$clus_comb)))
    cells$hoverLabels <- paste0("Cluster: ", cells$orig_cluster, ", Cell: ", cells$Sample)

    sum <- cells %>%
      distinct(.) %>%
      arrange(clus_comb) %>%
      group_by(clus_comb) %>%
      summarise(n = n()) %>%
      distinct(.) %>%
      mutate(clus_comb = factor(clus_comb, levels = clus_comb))
    tex <- unlist(sapply(1:nrow(sum), function(y) {
      vec <- c(sum$n[y], rep("", sum$n[y] - 1))
    }))
    cells$tex <- tex

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
    
    plot_ly(cells, textposition = "none") %>%
      add_trace(
        x = ~x,
        y = ~y,
        width = 1,
        type = "bar",
        color = ~as.character(orig_cluster),
        colors = gg_color(len),
        text = ~hoverLabels,
        hovertemplate = "<b>%{text}</b>",
        hoverlabel = label
      ) %>%
      layout(barmode = "stack", legend = list(traceorder = "normal")) %>%
      add_text(x = ~x, y = ~count, text = ~tex, textposition = 'top', showlegend = F,
               textfont=list(size = 14, color="black")) %>%
      layout(yaxis = list(title = "Number of cells", tickfont = list(size = 14), titlefont = list(size = 16)),
             xaxis = list(title = "Cluster", tickfont = list(size = 14), titlefont = list(size = 16)),
             legend = list(font = list(size = 14)),
             font = list(family = "Noto Sans JP")
             )
  })

  # SC-DGE-CCC - create col name for new metadata
  output$username <- renderUI({
    if(is.null(seurat_only()) || is.null(recomb_seur())) return()
    textInput(
      inputId = "username", 
      label = "Username",
      value = ""
    )
  })

  # SC-DGE-CCC - create a note to add to metadata
  output$note <- renderUI({
    if(is.null(seurat_only()) || is.null(recomb_seur())) return()
    textAreaInput(
      inputId = "note", 
      label = "Comment",
      value = "",
      width = "200px"
    )
  })

  # SC-DGE-CCC - specify file output name
  output$table_name <- renderUI({
    req(input$username)
    if(is.null(seurat_only()) || is.null(recomb_seur()) || is.null(d$datasetNames)) return()
    datasets <- paste(d$datasetNames, collapse = "_")
    # need to push info to sc_useradd
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                      dbname = scdb, host = ec_host, port = p)
    dat <- dbReadTable(mydb, "sc_useradd")
    dbDisconnect(mydb)
    idx <- dat %>%
      filter(grepl(paste(datasets, input$username, sep = "-"), table_name))

    if (nrow(idx) == 0) {
      idx <- 1
    } else {
      idx <- idx %>%
        group_by(experiment_id) %>%
        summarise(count = n()) %>%
        pull(count)
      idx <- idx + 1
    }
    datasets <- paste(datasets, sep = "_")
    datasets <- paste(datasets, input$username, idx, sep = "-")
    textInput(
      inputId = "tablename", 
      label = "Updated table name",
      value = datasets 
    )
  })

  # SC-DGE-CCC - create col name for samples used
  samples_used <- reactive({
    req(seurat_only(), recomb_seur())
    if (is.null(seurat_only()) | is.null(recomb_seur())) return()
    if (is.null(input$save_load_selected)) {
      sel_samps <- rownames(SubmitData$selectedData$metadata)
      # specify if the end user chose samples (action button)
    } else if (input$save_load_selected != 0) {
      sel_samps <- rownames(SubmitData$selectedData$metadata)
    } else {
      sel_samps <- rownames(SubmitData$selectedData$metadata)
    }
    return(sel_samps)
  })

  # SC-DGE-CCC - action button for pushing to SQL db
  output$submit_mysql <- renderUI({
    req(input$username, input$tablename)
    if(is.null(recomb_seur())) return()
    actionButton("submit_mysql", label = "Save to database")
  })

  # SC-DGE-CCC - download button for seurat object after combining clusters
  output$download_seurat_comb <- renderUI({
    if(is.null(recomb_seur())) return()
    downloadButton("download_seurat_comb_button", label = "Download Seurat (RDS)")
  })

  # SC-DGE-CCC - download button to download seurat object after combining clusters
  output$download_seurat_comb_button <- downloadHandler(
    filename = function() {
      "SeuratCombined_MergedClusters.rds"
    },
    content=function(file) {
      seuratData <- d$inD
      meta <- recomb_seur()
      rownames(meta) <- meta$Sample
      seuratData@meta.data <- meta[colnames(seuratData), ]
      saveRDS(
        seuratData,
        file = file
      )
    }
  )

  # SC-DGE-CCC - download button for metadata after combining clusters
  output$download_metadata_comb <- renderUI({
    # req(d$inD, recomb_seur(), count_content(), meta_content())
    if(is.null(recomb_seur())) return()
    downloadButton("download_metadata_comb_button", label = "Download metadata (CSV)")
  })

  # SC-DGE-CCC - download file for metadata after combining clusters
  output$download_metadata_comb_button <- downloadHandler(
    filename = function() {
      "Metadata_Combined_MergedClusters.csv"
    },
    content=function(file) {
      meta <- recomb_seur()
      rownames(meta) <- meta$Sample
      write.csv(
        meta,
        file = file,
        quote = F,
        row.names = F
      )
    }
  )

  output$comb_buttons_ui <- renderUI({
    if(is.null(recomb_seur())) return()
    tagList(
      div(style = "display: inline-block;", uiOutput("download_seurat_comb")),
      div(style = "display: inline-block; margin-left: 20px;", uiOutput("download_metadata_comb")),
      div(style = "display: inline-block; margin-left: 20px;", uiOutput("submit_mysql"))
    )
  })
  
  # SC-DGE-CCC - reactive expression to re-create metdata for new column & note
  new_meta <- reactive ({
    dat <- recomb_seur()
    ind <- which(names(dat) %in% "clus_comb")
    names(dat)[5] <- input$column_name
    dat$Note <- input$note
    return(dat)
  })
  
  # SC-DGE-CCC - push updated metadata to sql db
  observeEvent(input$submit_mysql, {
    if(SubmitData$data_type == "Single-cell") {
      if(length(input$clusters_comb) >= 2) {
        withProgress(message = "Pushing recombined cluster data to the MySQL database.", value = 0, {
          incProgress(1/3)

          # need to push info to sc_useradd
          mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sup, password = pwd_sc_sup,
                            dbname = scdb, host = ec_host, port = p)
          dat <- dbReadTable(conn = mydb, name = "sc_useradd")
          if(nrow(dat) == 0) {
            idx <- 1
          } else {
            idx <- dat %>%
              dplyr::select(iterator) %>%
              dplyr::slice(n()) %>%
              mutate(iterator = iterator + 1) %>%
              pull(iterator)
          }
          datasets <- data.frame(
            experiment_id = as.integer(SubmitData$datasets_table$experiment_id[SubmitData$datasets_table$experiment_name %in% d$datasetNames][1])
          )

          df <- dbReadTable(mydb, "sc")
          names(df)[1] <- gsub("^X\\.", "", names(df)[1])
          
          sel_samps <- samples_used()
          sel_samps <- paste(sel_samps, collapse = "; ")

          df$experiment_id <- as.integer(df$experiment_id)
          df <- df %>%
            filter(experiment_id %in% datasets$experiment_id) %>%
            dplyr::select(c(experiment_id, unique_table)) %>%
            mutate(username = input$username,
                   uploaded_date = as.character(gsub(" .*", "", Sys.time())),
                   table_name = input$tablename,
                   samples_used = sel_samps,
                   note = input$note,
                   iterator = idx)
          df <- df %>% dplyr::select(c(iterator, experiment_id,
                                       unique_table, username,
                                       uploaded_date, table_name,
                                       samples_used, note))
          dbWriteTable(conn = mydb, name = "sc_useradd", value = df, append = TRUE, row.names = FALSE)
          # need to push table to scdb
          meta <- recomb_seur() %>%
            dplyr::select(-`Number of Cells`)
          meta$note <- NA
          meta$note[1] <- input$note
          
          # dbWriteTable throws error from Sample and sample in column names because
          # it is case-insensitive, so change sample to sample_DUP. This is corrected
          # when reading in user-created metadata
          if("sample" %in% colnames(meta)) {
            colnames(meta)[colnames(meta) == "sample"] = "sample_DUP"
          }
          
          dbWriteTable(conn = mydb, name = input$tablename, value = meta, overwrite = TRUE, row.names = FALSE)
          dbDisconnect(mydb)
          incProgress(2/3)
        })
      }
    }
  })

  # SC-DGE-CCC - combine cells based on expression of goi
  output$mode = renderUI({
    req(d$inD_orig)
    awesomeRadio("mode", label = "", choices = c("Single gene", "Multiple genes"),
                 selected = "Single gene", inline = T, status = "success")
  })

  # SC-DGE-CCC - choose metric for cell combination of multiple gois
  output$metric <- renderUI({
    req(d$inD_orig, input$gene)
    if(is.null(d$inD_orig) || is.null(input$gene) || length(input$gene) == 1) return()
    div(
      style = "width: 120px; margin-left: 20px;",
      selectInput(
        inputId = "metric", 
        label = "Summary metric",
        choices = c("Mean", "Median", "Sum"), 
        selected = "Mean"
      )
    )
  })
  
  # SC-DGE-CCC - choose genes for cell combination
  output$gene <- renderUI({
    req(d$inD_orig)
    selectizeInput("gene", label = "Select genes", choices = rownames(d$inD_orig),
                   selected = rownames(d$inD_orig)[1], multiple = T,
                   options = list(maxOptions = 1000000))
  })

  output$varRange_ui <- renderUI({
    if(is.null(input$gene)) return()
    div(
      style = "width: 300px; margin-left: 20px;",
      fluidPage(
        # fluidRow(uiOutput("varRangeHeader")),
        fluidRow(
          column(
            style = "width: 140px;",
            width = 5, 
            uiOutput("varRangeText1"),
            uiOutput("inclusive1")
          ),
          column(
            style = "width: 140px; margin-left: 20px;",
            width = 5, 
            offset = 1,
            uiOutput("varRangeText2"),
            uiOutput("inclusive2")
          )
        )
      )
    )
  })
  
  # SC-DGE-CCC - add updated metadata to object d after cell comb
  observeEvent(list(input$gene, input$metric), {
    if(is.null(input$gene)) return()
    
    if(is.null(d$inD_orig)) return()
    if(!(all(input$gene %in% rownames(d$inD_orig)))) return()
    if(length(input$gene) == 1) {
      d$inD_goi = AddMetaData(d$inD_orig, col.name = "var", metadata = FetchData(d$inD_orig, vars = input$gene))
    } else {
      if(is.null(input$metric)) return()
      if(input$metric == "Mean") {
        var = rowMeans(FetchData(d$inD_orig, vars = input$gene))
      } else if(input$metric == "Median") {
        var = apply(FetchData(d$inD_orig, vars = input$gene), 1, median)
      } else if(input$metric == "Sum") {
        var = rowSums(FetchData(d$inD_orig, vars = input$gene))
      }
      d$inD_goi = AddMetaData(d$inD_orig, col.name = "var", metadata = var)
    }
  }, ignoreInit = T)

  # SC-DGE-CCC - replot combined cells based on goi(s)
  output$selectClustersByGOI <- renderPlotly({
    req(d$inD_goi, input$gene)
    input$deleteAll
    if(!("var" %in% colnames(d$inD_goi@meta.data))) return()
    if(length(input$gene) == 1) {
      plotTitle = input$gene
    } else {
      req(input$metric)
      if(input$metric == "Mean") {
        plotTitle = "Mean gene expression"
      } else if(input$metric == "Median") {
        plotTitle = "Median gene expression"
      } else if(input$metric == "Sum") {
        plotTitle = "Sum of gene expression"
      }
    }
    cell = colnames(d$inD_goi)
    p = FeaturePlot(object = d$inD_goi,
                    features = "var",
                    reduction = "umap",
                    label = F) +
      xlab("UMAP 1") + ylab("UMAP 2") +
      labs(title = plotTitle)
    p$data$cell = cell
    
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
    
    p = p + aes(key = cell, text = paste("</br> <b>Cluster:</b> ", p$data$Cluster,
                                         "</br>",
                                         "<b>UMAP 1: </b>", round(p$data[, 1], 3),
                                         "</br>",
                                         "<b>UMAP 2: </b>", round(p$data[, 2], 3),
                                         "</br>",
                                         "<b> Cell: </b>", p$data$cell))
    ggplotly(p, tooltip = "text", width = 540, height = 440) %>%
      layout(
        autosize = F, 
        dragmode = "lasso",
        font = list(family = "Noto Sans JP"),
        margin = list(r = 0)
      ) %>%
      style(hoverlabel = label)
  })

  # SC-DGE-CCC - header for new set name
  output$setNameHeader = renderUI({
    req(d$inD_goi)
    h5("New set name")
  })

  # SC-DGE-CCC - specify name for updated combined cell dataset
  output$setName = renderUI({
    req(d$inD_goi)
    suggestedValue = paste0("set ", length(d$setNames) + 1)
    # textInput("setName", label = NULL, value = suggestedValue)
    textInput("setName", label = "New set name", value = suggestedValue)
  })

  # SC-DGE-CCC - action button for creating a set to add cells to (lassoing)
  output$addSet = renderUI({
    req(d$inD_goi)
    actionButton("addSet", label = "+New set")
  })

  # SC-DGE-CCC - updated info for: action button to add cells to
  # a given set, the total number of cells added to a set
  # and warning msgs
  observeEvent(input$addSet, {
    if(input$setName %in% d$setNames) return()
    d$setNames <- unique(c(d$setNames, input$setName))
    d$addCellsClicks[[input$setName]] <- 0
    setName <- input$setName
    output[[paste0("addTo", setName)]] <- renderUI({
      actionButton(paste0("addTo", setName), label = paste0("Add to ", setName))
    })
    output[[paste0("cellInfo", setName)]] <- renderText({
      cellInfo <- ""
      if(length(d$cellSets[[setName]]) > 0) {
        if(length(d$cellSets[[setName]]) == 1) {
          cellInfo <- paste0(length(d$cellSets[[setName]]), " cell")
        } else {
          cellInfo <- paste0(length(d$cellSets[[setName]]), " cells")
        }
      }
      cellInfo
    })
    output[[paste0("error", setName)]] <- renderText({
      if(is.null(d$addCellsOverlap)) return()
      if(!(setName %in% names(d$addCellsOverlap))) return()
      if(is.null(d$addCellsOverlap[[setName]])) return()
      if(d$addCellsOverlap[[setName]] == 1) {
        paste0(d$addCellsOverlap[[setName]], " cell overlaps other sets.")
      } else {
        paste0(d$addCellsOverlap[[setName]], " cells overlap other sets.")
      }
    })
  })

  # SC-DGE-CCC - ui output for clearing all sets, adding cells
  # to a set and errors for combining cells
  output$addCells <- renderUI({
    req(d$inD_goi, d$setNames)
    c(
      list(
        div(
          style = "display: inline-block; width: 100px;",
          actionButton("clearAll", label = "Clear all sets")
        ),
        div(
          style = "display: inline-block; width: 100px; margin-left: 10px;",
          actionButton("deleteAll", label = "Delete all sets")
        ),
        br(),
        br()
      ),
      lapply(d$setNames, function(i) {
        tagList(
          fluidPage(
            style = "border: 1px solid #eee; padding: 10px;",
            div(
              div(style = "display: inline-block;", uiOutput(paste0("addTo", i))),
              div(
                style = "display: inline-block; margin-left: 5px;", 
                h4(textOutput(paste0("cellInfo", i)))
              )
            ),
            div(textOutput(paste0("error", i)), style = "color: red;")
          ),
          br()
        )
      })
    )
  })

  # SC-DGE-CCC - reactive expression for lasso cell selection
  plotlySelected_GOI <- reactive({
    plotlySelected <- event_data("plotly_selected")
    if(is.null(plotlySelected)) return()
    return(plotlySelected$key[!is.na(plotlySelected$key)])
  })

  # SC-DGE-CCC - cell selection from filters and/or brush
  currSel_GOI <- reactive({
    if(is.null(input$inclusive1) || is.null(input$inclusive2)) return()
    temp_points <- rownames(brushedPoints(
      as.data.frame(getEmb(d$inD_goi, "umap")[,c("UMAP_1", "UMAP_2")]),
      input$selectClustersByGOIbrush,
      xvar = "UMAP_1", yvar = "UMAP_2"
    ))
    if(input$inclusive1 && input$inclusive2) {
      temp_picker = as.numeric(d$inD_goi@meta.data[, "var"]) >= input$varRangeText1 &
        as.numeric(d$inD_goi@meta.data[, "var"]) <= input$varRangeText2
    } else if(input$inclusive1) {
      temp_picker = as.numeric(d$inD_goi@meta.data[, "var"]) >= input$varRangeText1 &
        as.numeric(d$inD_goi@meta.data[, "var"]) < input$varRangeText2
    } else if(input$inclusive2) {
      temp_picker = as.numeric(d$inD_goi@meta.data[, "var"]) > input$varRangeText1 &
        as.numeric(d$inD_goi@meta.data[, "var"]) <= input$varRangeText2
    } else {
      temp_picker = as.numeric(d$inD_goi@meta.data[, "var"]) > input$varRangeText1 &
        as.numeric(d$inD_goi@meta.data[, "var"]) < input$varRangeText2
    }
    if (length(temp_points) > 0 & length(temp_picker) > 0) {
      return(colnames(d$inD_goi)[colnames(d$inD_goi) %in% temp_points & temp_picker])
    } else if (length(temp_picker) > 0 & !all(temp_picker)) {
      return(colnames(d$inD_goi)[temp_picker])
    } else if (length(temp_points) > 0) {
      return(temp_points)
    } else { return(character()) }
  })

  # SC-DGE-CCC - header for selecting expression range
  output$varRangeHeader = renderUI({
    req(d$inD_goi, input$gene, d$inD_goi@meta.data[["var"]])
    h4("Select range", style = "margin-top: 0px;")
  })

  # SC-DGE-CCC - gene expression range min
  output$varRangeText1 = renderUI({
    req(d$inD_goi, input$gene, d$inD_goi@meta.data[["var"]])
    var = as.numeric(d$inD_goi@meta.data[, "var"])
    textInput("varRangeText1", label = "From (min 0)", value = 0)
  })

  # SC-DGE-CCC - gene expression range max
  output$varRangeText2 = renderUI({
    req(d$inD_goi, input$gene, d$inD_goi@meta.data[["var"]])
    var = as.numeric(d$inD_goi@meta.data[, "var"])
    textInput("varRangeText2", label = paste0("To (max ", round(max(var), digits = 2), ")"),
              value = round(max(var), digits = 2))
  })

  # SC-DGE-CCC - include lower limits in selection
  output$inclusive1 = renderUI({
    req(d$inD_goi, input$gene, d$inD_goi@meta.data[["var"]])
    prettyCheckbox("inclusive1", label = "Inclusive",
                  value = T, status = "default", icon = icon("check"))
  })

  # SC-DGE-CCC - include upper limits in selection
  output$inclusive2 = renderUI({
    req(d$inD_goi, input$gene, d$inD_goi@meta.data[["var"]])
    prettyCheckbox("inclusive2", label = "Inclusive",
                  value = T, status = "default", icon = icon("check"))
  })

  # SC-DGE-CCC - updated minimum range
  observeEvent(input$varRangeText1, {
    if(is.null(input$varRangeText1) | is.null(d$inD_goi)) return()
    if(input$varRangeText1 == "") return()
    if(is.null(d$inD_goi@meta.data[["var"]])) return()
    if(as.numeric(input$varRangeText1) < 0 | as.numeric(input$varRangeText1) > max(d$inD_goi@meta.data[, "var"])) {
      updateTextInput(session, "varRangeText1", value = 0)
    }
  }, ignoreNULL = F)

  # SC-DGE-CCC - updated maximum range
  observeEvent(input$varRangeText2, {
    if(is.null(input$varRangeText2) | is.null(d$inD_goi)) return()
    if(input$varRangeText2 == "") return()
    if(is.null(d$inD_goi@meta.data[["var"]])) return()
    if(as.numeric(input$varRangeText2) < 0 | as.numeric(input$varRangeText2) > max(as.numeric(d$inD_goi@meta.data[, "var"]))) {
      updateTextInput(session, "varRangeText2", value = round(max(as.numeric(d$inD_goi@meta.data[, "var"])), digits = 2))
    }
  }, ignoreNULL = F)

  # SC-DGE-CCC - calculates the number of cells being combined
  observe({
    if(length(d$setNames) < 1) return()
    for(setName in d$setNames) {
      button = paste0("addTo", setName)
      if(is.null(input[[button]])) next
      if(input[[button]] == 0) {
        d$addCellsClicks[[setName]] = 0
        next
      }
      if(input[[button]] > d$addCellsClicks[[setName]]) {
        d$addCellsClicks[[setName]] = input[[button]]
        if(length(currSel_GOI()) %in% c(0, ncol(d$inD_goi))) {
          selectedCells = plotlySelected_GOI()
        } else if(length(plotlySelected_GOI()) == 0) {
          selectedCells = currSel_GOI()
        } else {
          selectedCells = plotlySelected_GOI()[plotlySelected_GOI() %in% currSel_GOI()]
        }
        if(length(selectedCells) == 0) return()
        overlap = sum(selectedCells %in% unlist(d$cellSets))
        if(overlap > 0) {
          d$addCellsOverlap[[setName]] = overlap
          return()
        }
        for(i in names(d$addCellsOverlap)) d$addCellsOverlap[[i]] = NULL
        d$cellSets[[setName]] = unique(c(d$cellSets[[setName]], selectedCells))
        d$observeCount = d$observeCount + 1
        return()
      }
    }
  })

  # SC-DGE-CCC - if the cell sets are cleared, set in object d
  observeEvent(input$clearAll, {
    for(setName in d$setNames) {
      d$cellSets[[setName]] = NULL
      d$addCellsOverlap[[setName]] = NULL
    }
  })
  
  observeEvent(input$deleteAll, {
    d$cellSets <- d$addCellsOverlap <- d$addCellsClicks <- list()
    d$setNames <- NULL
    runjs("Shiny.setInputValue('plotly_selected-A', null);")
  })

  # SC-DGE-CCC - create col name for new metadata (cells combined)
  output$username2 <- renderUI({
    if(is.null(d$inD_goi) || is.null(input$gene)) return()
    textInput(
      inputId = "username2", 
      label = "Username",
      value = ""
    )
  })

  # SC-DGE-CCC - create a note to add to metadata (cells combined)
  output$note2 <- renderUI({
    if(is.null(d$inD_goi) || is.null(input$gene)) return()
    textAreaInput(
      inputId = "note2", 
      label = "Comment",
      value = "", 
      width = "200px"
    )
  })

  # SC-DGE-CCC - specify file output name (cells combined)
  output$table_name2 <- renderUI({
    if(is.null(input$username2) || input$username2 == "" || is.null(d$datasetNames) ||
       is.null(d$inD_goi) || is.null(input$gene)) return()
    datasets <- paste(d$datasetNames, collapse = "_")

    # need to push info to sc_useradd
    mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sc, password = pwd_sc,
                      dbname = scdb, host = ec_host, port = p)
    dat <- dbReadTable(mydb, "sc_useradd")
    dbDisconnect(mydb)
    idx <- dat %>%
      filter(grepl(paste(datasets, input$username2, sep = "-"), table_name))

    if (nrow(idx) == 0) {
      idx <- 1
    } else {
      idx <- idx %>%
        group_by(experiment_id) %>%
        summarise(count = n()) %>%
        pull(count)
      idx <- idx + 1
    }

    datasets <- paste(datasets, sep = "_")
    datasets <- paste(datasets, input$username2, idx, sep = "-")

    textInput(
      inputId = "tablename2", 
      label = "Updated table name",
      value = datasets
    )
  })

  # SC-DGE-CCC - reactive expression storing selected data
  # in combined cells
  sc_samples_used <- reactive({
    req(d$inD_goi, input$gene)
    if (is.null(d$inD_goi) || is.null(input$gene)) return()

    if (is.null(input$save_load_selected)) {
      sel_samps <- rownames(d$selectedData$metadata)
      # specify if the end user chose samples (action button)
    } else if (input$save_load_selected != 0) {
      sel_samps <- rownames(d$selectedData$metadata)
    } else {
      sel_samps <- rownames(d$selectedData$metadata)
    }
    return(sel_samps)
  })

  # SC-DGE-CCC - action button to submit combined cells to sql db
  output$submit_mysql2 <- renderUI({
    req(input$username2, input$tablename2)
    if(is.null(d$cellSets) || length(d$cellSets) == 0) return()
    actionButton("submit_mysql2", label = "Save to database")
  })

  # SC-DGE-CCC - push updated metadata for combined cells to sql db
  observeEvent(input$submit_mysql2, {
    if(is.null(input$tablename2) || input$tablename2 == "") return()
    
    selectedCells = unlist(d$cellSets)
    cellSetsVec = rep.int(NA, length(selectedCells))
    names(cellSetsVec) = selectedCells
    for(setName in names(d$cellSets)) {
      if(length(d$cellSets[[setName]]) > 0) cellSetsVec[d$cellSets[[setName]]] = setName
    }

    withProgress(message = "Pushing recombined cluster data to the server.", value = 0, {
      incProgress(1/3)
      # need to push info to sc_useradd
      mydb <- dbConnect(RMariaDB::MariaDB(), user = usr_sup, password = pwd_sc_sup,
                        dbname = scdb, host = ec_host, port = p)
      dat <- dbReadTable(conn = mydb, name = "sc_useradd")
      if(nrow(dat) == 0) {
        idx <- 1
      } else {
        idx <- dat %>%
          dplyr::select(iterator) %>%
          dplyr::slice(n()) %>%
          mutate(iterator = iterator + 1) %>%
          pull(iterator)
      }

      datasets <- data.frame(
        experiment_id = as.integer(SubmitData$datasets_table$experiment_id[SubmitData$datasets_table$experiment_name %in% d$datasetNames][1])
      )
      
      df <- dbReadTable(mydb, "sc")
      names(df)[1] <- gsub("^X\\.", "", names(df)[1])

      sel_samps <- sc_samples_used()
      sel_samps <- paste(sel_samps, collapse = "; ")
      
      df$experiment_id <- as.integer(df$experiment_id)
      df <- df %>%
        filter(experiment_id %in% datasets$experiment_id) %>%
        dplyr::select(c(experiment_id, unique_table)) %>%
        mutate(username = input$username2,
               uploaded_date = as.character(gsub(" .*", "", Sys.time())),
               table_name = input$tablename2,
               samples_used = sel_samps,
               note = input$note2,
               iterator = idx)
      df <- df %>% dplyr::select(c(iterator, experiment_id,
                                   unique_table, username,
                                   uploaded_date, table_name,
                                   samples_used, note))
      dbWriteTable(conn = mydb, name = "sc_useradd", value = df, append = TRUE, row.names = FALSE)
      
      incProgress(2/3)

      seuratData <- d$inD_goi[, selectedCells]
      meta <- seuratData@meta.data
      meta <- meta %>%
        mutate(Set = cellSetsVec,
               note = NA) %>%
        dplyr::select(-var)
      meta$note[1] <- input$note2
      
      meta$Sample <- rownames(meta)
      
      # dbWriteTable throws error from Sample and sample in column names because
      # it is case-insensitive, so change sample to sample_DUP. This is corrected
      # when reading in user-created metadata
      if("sample" %in% colnames(meta)) {
        colnames(meta)[colnames(meta) == "sample"] = "sample_DUP"
      }
      
      dbWriteTable(conn = mydb, name = input$tablename2, value = meta, overwrite = TRUE, row.names = FALSE)
      dbDisconnect(mydb)
    })
  })

  # SC-DGE-CCC - download button for seurat data after combining cells
  output$download_new_dataset_GOI = renderUI({
    if(is.null(d$cellSets) || !(is.list(d$cellSets)) || length(unlist(d$cellSets)) == 0) return()
    downloadButton("download_new_dataset_GOI_button", label = "Download Seurat (RDS)")
  })

  # SC-DGE-CCC - download file for seurat obj after combining cells
  output$download_new_dataset_GOI_button <- downloadHandler(
    filename = function() {
      "SeuratCombined_GOI.rds"
    },
    content=function(file) {
      selectedCells = unlist(d$cellSets)
      cellSetsVec = rep.int(NA, length(selectedCells))
      names(cellSetsVec) = selectedCells
      for(setName in names(d$cellSets)) {
        if(length(d$cellSets[[setName]]) > 0) cellSetsVec[d$cellSets[[setName]]] = setName
      }
      seuratData = d$inD_goi[, selectedCells]
      seuratData = AddMetaData(seuratData, metadata = cellSetsVec, col.name = "Set")
      saveRDS(
        seuratData,
        file = file
      )
    }
  )

  # SC-DGE-CCC - download button for metadata after combining cells
  output$download_new_metadata_GOI = renderUI({
    if(is.null(d$cellSets) || !(is.list(d$cellSets)) || length(unlist(d$cellSets)) == 0) return()
    downloadButton("download_new_metadata_GOI_button", label = "Download metadata (CSV)")
  })

  # SC-DGE-CCC - download file for metadata after combining cells
  output$download_new_metadata_GOI_button <- downloadHandler(
    filename = function() {
      "Metadata_Combined_GOI.csv"
    },
    content=function(file) {
      selectedCells = unlist(d$cellSets)
      cellSetsVec = rep.int(NA, length(selectedCells))
      names(cellSetsVec) = selectedCells
      for(setName in names(d$cellSets)) {
        if(length(d$cellSets[[setName]]) > 0) cellSetsVec[d$cellSets[[setName]]] = setName
      }
      seuratData = d$inD_goi[, selectedCells]
      seuratData = AddMetaData(seuratData, metadata = cellSetsVec, col.name = "Set")
      write.csv(
        data.frame(Cell = rownames(seuratData@meta.data),
                   seuratData@meta.data),
        file = file,
        quote = F,
        row.names = F
      )
    }
  )

  ###################################################################
  ###################################################################
  ### SECTION 05 - MORE
  ###################################################################
  ###################################################################
  
  # MORE - Add session info verbatim
  output$sessinfo <- renderPrint({
    sessionInfo()
  })
  
}