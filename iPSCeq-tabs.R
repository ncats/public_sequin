#---------------------------------------------------------------------
# Title:         NCATS SEQUIN - Tabs
# Author:        Marissa Hirst
# Author2:       Ben Ernest
# Last Modified: 2021-04-01
# --
# Created:       2018-01-26 11:29:39 CDT
#---------------------------------------------------------------------

# # this has to go before load page b/c it 
# is called in load page
# qc: count summary, metadata, qc plots
tab.submit <- tabPanel(
  title = div("Summary", style = "margin-top: -1px;"),
  fluid = TRUE,
  value = "summary",
  sidebarLayout(
    sidebarPanel(
      width = 2,
      uiOutput("showDataSummaryLinks"),
      uiOutput("dataSummaryLinks"),
      uiOutput("showQCLinks"),
      uiOutput("qcLinks")
    ),
    mainPanel(
      width = 10,
      div(
        style = "background-color: #eee;",
        fluidRow( 
          column(width = 3, uiOutput("pre")),
          column(width = 3, uiOutput("post")),
          column(width = 3, uiOutput("samplecount")),
        ),
      ),
      hr(),
      tabsetPanel(
        id = "dataSummaryQCTabset",
        type = "hidden",
        tabPanelBody(
          value = "dataSummary",
          tabsetPanel(
            id = "dataSummary_tabsetPanel",
            type = "hidden",
            tabPanelBody(
              value = "dataSummary_subsetCountData",
              fluidPage(
                div(style = "display: inline-block;", uiOutput("filesummarycts")),
                div(
                  style = "display: inline-block; float: right;",
                  div(
                    style = "display: inline-block;",
                    uiOutput("download_raw_counts")
                  ),
                  div(
                    style = "display: inline-block; margin-left: 10px;",
                    uiOutput("download_normalized_data")
                  )
                ),
                uiOutput("fileoutputcts"),
                br()
              )
            ),
            tabPanelBody(
              value = "dataSummary_metadata",
              fluidPage(
                uiOutput("filesummarycoldata"),
                uiOutput("fileoutputcoldata"),
                br()
              )
            )
          )
        ),
        tabPanelBody(
          value = "qc",
          div(
            style = "display: inline-block; float: right; margin-right: 20px;",
            uiOutput("countsummarylabel")
          ),
          tabsetPanel(
            id = "qc_tabsetPanel",
            type = "hidden",
            tabPanelBody(
              value = "qc_boxAndWhisker",
              fluidPage(
                uiOutput("countbox"),
                uiOutput("boxplot_rotate"),
                withSpinner(
                  plotlyOutput("boxplot"),
                  type = 8, size = 1, color = "black", proxy.height = "300px"
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqcboxplotpdf")
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqcboxplotpng")
                ),
                br(),br()
              )
            ),
            tabPanelBody(
              value = "qc_histogram",
              fluidPage(
                uiOutput("counthist"),
                withSpinner(
                  plotlyOutput("hist"),
                  type = 8, size = 1, color = "black", proxy.height = "300px"
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqchistpdf")
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqchistpng")
                )
              ),
              br()
            ),
            tabPanelBody(
              value = "qc_totalReads",
              fluidPage(
                uiOutput("counttotal"),
                uiOutput("barplot_rotate"),
                withSpinner(
                  plotlyOutput("barplot"),
                  type = 8, size = 1, color = "black", proxy.height = "300px"
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqcbarplotpdf")
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqcbarplotpng")
                )
              ),
              br()
            )
          )
        )
      )
    )
  )
)

# load page: select datasets, choose samples
tab.selectdata <- tabPanel(
  title = "Submit and quality control",
  value = "val1",
  icon = icon("filter"),
  LoadDataUI("load_data")
)
  
# DDA: correlation matrix, distance matrix, dimensional 
# reduction plot
tab.prelim <- tabPanel(
  title = "Data structure ",
  icon = icon("search"),
  fluid = TRUE,
  value = "val2",
  sidebarLayout(
    sidebarPanel(
      width = 2,
      uiOutput("showCorrelationLinks"),
      uiOutput("correlationLinks"),
      uiOutput("showClustering")
    ),
    mainPanel(
      width = 10,
      tabsetPanel(
        id = "dda_tabsetPanel",
        type = "hidden",
        tabPanelBody(
          value = "dda_correlation",
          tabsetPanel(
            id = "dda_correlation_tabsetPanel",
            type = "hidden",
            tabPanelBody(
              value = "dda_correlation_correlation",
              fluidPage(
                ################
                fluidRow(
                  div(style = "display: inline-block;", uiOutput("correlation_options_header")),
                  div(
                    style = "display: inline-block; float: right;", 
                    actionLink("cor_help_button", label = "Help")
                  )
                ),
                fluidRow(
                  column(
                    width = 2,
                    style = "border: 1px solid #eee; padding-left: 10px; width: 180px; height: 180px;",
                    fluidPage(
                      fluidRow(div(style = "width: 100px;", uiOutput("cor_nGenes"))),
                      fluidRow(uiOutput("cor_all_genes")),
                      fluidRow(uiOutput("cor_default_genes"))
                    )
                  ),
                  column(
                    width = 3,
                    style = "border: 1px solid #eee; padding-left: 10px; width: 180px; margin-left: 20px; height: 180px;",
                    fluidPage(
                      fluidRow(
                        div(style = "display: inline-block; width: 100px;", uiOutput("cor_nSamples")),
                        div(style = "display: inline-block; margin-left: 10px;", uiOutput("cor_default_samples"))
                      ),
                      uiOutput("cor_seed")
                    )
                  ),
                  column(
                    width = 4,
                    style = "margin-left: 20px;",
                    fluidRow(div(style = "width: 250px;", uiOutput("plotColors"))),
                    fluidRow(
                      div(style = "color: red; width: 400px;", uiOutput("cor_downsample_msg"))
                    ),
                    br(),
                    fluidRow(
                      uiOutput("build_cor_matrix")
                    )
                  )
                ),
                fluidRow(hr()),
                ################
                uiOutput("headcor"),
                fluidRow(column(12, uiOutput("showLabels"))),
                withSpinner(
                  plotlyOutput("corplot1", width = 1000, height = 600),
                  type = 8, size = 1, color = "black", proxy.height = "300px"
                ),
                br(),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqccorplot1pdf")
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqccorplot1png")
                ),
                br(),
                br(),
                br(),
                uiOutput("corplot2Header"),
                withSpinner(
                  plotlyOutput("corplot2", width = 800, height = 600),
                  type = 8, size = 1, color = "black", proxy.height = "300px"
                ),
                br(),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqccorplot2pdf")
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqccorplot2png")
                ),
                br(), br()
              )
            ),
            tabPanelBody(
              value = "dda_correlation_distance",
              fluidPage(
                fluidRow(
                  div(style = "display: inline-block;", uiOutput("distance_options_header")),
                  div(
                    style = "display: inline-block; float: right;", 
                    actionLink("dist_help_button", label = "Help")
                  )
                ),
                fluidRow(
                  column(
                    width = 2,
                    style = "border: 1px solid #eee; padding-left: 10px; width: 180px; height: 180px;",
                    fluidPage(
                      fluidRow(div(style = "width: 100px;", uiOutput("dist_nGenes"))),
                      fluidRow(uiOutput("dist_all_genes")),
                      fluidRow(uiOutput("dist_default_genes"))
                    )
                  ),
                  column(
                    width = 3,
                    style = "border: 1px solid #eee; padding-left: 10px; width: 180px; margin-left: 20px; height: 180px;",
                    fluidPage(
                      fluidRow(
                        div(style = "display: inline-block; width: 100px;", uiOutput("dist_nSamples")),
                        div(style = "display: inline-block; margin-left: 10px;", uiOutput("dist_default_samples"))
                      ),
                      uiOutput("dist_seed")
                    )
                  ),
                  column(
                    width = 4,
                    style = "margin-left: 20px;",
                    fluidRow(div(style = "width: 250px;", uiOutput("dist_plotColors"))),
                    fluidRow(
                      div(style = "color: red; width: 400px;", uiOutput("dist_downsample_msg"))
                    ),
                    br(),
                    fluidRow(
                      uiOutput("build_dist_matrix")
                    )
                  )
                ),
                fluidRow(hr()),
                uiOutput("headcor2"),
                fluidRow(column(12, uiOutput("dist_labels"))),
                plotlyOutput("corplot3", width = 1000, height = 600),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqcorplot3pdf")
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("dlqcorplot3png")
                ),
                br(), br()
              )
            )
          )
        ),
        tabPanelBody(
          value = "dda_clustering",
          uiOutput("dda_dimred"),
        )
      )
    )
  )
)

# DGE tab
tab.dge <- tabPanel(
  title = "Analysis",
  value = "val3",
  icon = icon("search"),
  fluid = TRUE,
  uiOutput("dge")
)

# Profiler tab
tab.profiler <- tabPanel(
  title = "iPSC profiler",
  value = "val4",
  fluid = TRUE,
  uiOutput("profiler")
)

# Help ----
tab.tutorial <- tabPanel(
  title = "Tutorial",
  fluid = TRUE,
  fluidPage(
    tags$head(tags$style(HTML(
      "
        img{
            max-width: 100%;
        }          
      "
    ))),
    fluidRow(
      column(
        width = 10,
        includeHTML("markdown/ncats-SEQUIN-tutorial.html")
      )
    )
  )
)

# FAQ ----
tab.faq <- tabPanel(
    title = "FAQ",
    fluid = TRUE,
    sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
            includeMarkdown("./markdown/faq.md")
        )
    )
)

# About Us ----
tab.about <- tabPanel(
    title = "About us",
    fluid = TRUE,
    sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
            includeMarkdown("./markdown/about-us.md")
        )
    )
)

# System Info ----
tab.sessinfo <- tabPanel(
    title = "Session info",
    fluid = TRUE,
    sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
            h2("R session info"),
            verbatimTextOutput("sessinfo")
        )
    )
)