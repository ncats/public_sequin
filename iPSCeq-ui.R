#---------------------------------------------------------------------
# Title:         NCATS SEQUIN - User Interface
# Author:        Marissa Hirst
# Author2:       Ben Ernest
# Last Modified: 2021-04-01
# --
# Created:       2018-01-26 11:29:39 CDT
#---------------------------------------------------------------------

iPSCeqUI <- tagList(
  useShinyjs(),
  withMathJax(),
  # .css file to style app
  tags$head(tags$link(
    rel = "stylesheet", type = "text/css", href = "layout.css"
  )),
  tags$head(tags$link(
    rel = "shortcut icon", href = "favicon.ico"
  )),
  use_font(
    id = "lato",
    css_path = "www/css/noto-sans-jp.css"
  ),
  tabsetPanel(
    id = "main_tabsetPanel",
    type = "hidden",
    tabPanelBody(
      value = "splashPage"
    ),
    tabPanelBody(
      value = "mainApp",
      tags$head(tags$style(HTML(
          "
          .short-text { display: none; }
          @media (max-width: 1200px) {
            .short-text { display: inline-block; }
            .full-text { display: none; }
          }
          
          "
      ))),
      
      tagList(
        div(
          style = "background-color: #662E6B; height: 30px;",
          img(
            src = "masthead-hhs-logo.png", 
            height = 18,
            style = "display: inline-block; margin-left: 65px; margin-bottom: 3px;"
          ),
          a(
            div(class = "full-text", "U.S. Department of Health and Human Services"),
            div(class = "short-text", "HHS"),
            href = "https://www.hhs.gov/",
            target = "_blank",
            style = paste0(
              "display: inline-block; color: #fff !important; font-size: 14px; margin-left: 3px; ",
              "margin-top: 5px; font-family: sans-serif;"
            )
          ),
          img(
            src = "masthead-divider.png",
            height = 28,
            style = "display: inline-block; margin-left: 6px;"
          ),
          img(
            src = "masthead-nih-logo.png",
            height = 18,
            style = "display: inline-block; margin-left: 4px; margin-bottom: 3px;"
          ),
          a(
            div(class = "full-text", "National Institutes of Health"),
            div(class = "short-text", "NIH"),
            href = "https://www.nih.gov/",
            target = "_blank",
            style = paste0(
              "display: inline-block; color: #fff !important; font-size: 14px; margin-left: 3px; ",
              "margin-top: 5px; font-family: sans-serif;"
            )
          ),
          img(
            src = "masthead-divider.png",
            height = 28,
            style = "display: inline-block; margin-left: 6px;"
          ),
          img(
            src = "masthead-nih-logo.png",
            height = 18,
            style = "display: inline-block; margin-left: 4px; margin-bottom: 3px;"
          ),
          a(
            div(class = "full-text", "National Center for Advancing Translational Sciences"),
            div(class = "short-text", "NCATS"),
            href = "https://www.ncats.nih.gov/",
            target = "_blank",
            style = paste0(
              "display: inline-block; color: #fff !important; font-size: 14px; margin-left: 3px; ",
              "margin-top: 5px; font-family: sans-serif;"
            )
          ),
          img(
            src = "masthead-divider.png",
            height = 28,
            style = "display: inline-block; margin-left: 6px;"
          )
        ),
        div(
          style = "background-color: #fff; height: 110px;",
          a(
            img(
              src = "ncats-logo.png",
              width = 237,
            ),
            href = "https://www.ncats.nih.gov/",
            target = "_blank",
            style = paste0(
              "display: inline-block; margin-left: 65px; ",
              "margin-top: 29px; padding-right: 15px; border-right: 2px solid #aaa;"
            )
          ),
          div(
            style = paste0(
              "display: inline-block; color: #672e6c; font-size: 40px; ",
              "padding-left: 15px; vertical-align: middle; margin-bottom: 12px; ",
              "font-weight: 350;"
            ),
            "SEQUIN"
          )
        ),
        navbarPage(
          position = "static-top",
          theme = shinytheme("cerulean"),
          windowTitle = "SEQUIN",
          title = uiOutput("placeholderDesc"),
          id = "tab_structure",
          header = uiOutput("topBar"),
          tab.selectdata,
          tab.submit,
          tab.prelim,
          tab.dge,
          tab.profiler,
          navbarMenu(
            title = "More",
            tab.tutorial,
            tab.faq,
            tab.about,
            tab.sessinfo
          )
        )
      )
    )
  )
)
