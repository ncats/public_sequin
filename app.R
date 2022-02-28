#---------------------------------------------------------------------
# Title:         NCATS Complex Shiny Application
# Author:        Brandon Monier
# Author2:       Marissa Hirst
# Author3:       Ben Ernest
# Last Modified: 2020-03-04
# --
# Created:       2018-01-26 11:29:39 CDT
#---------------------------------------------------------------------

# Run as standalone app (and set local data directory) or connect to database
if(is.null(getOption("standalone"))) options(standalone = T)
if(is.null(getOption("localDir"))) options(localDir = "example_data")

# Turn off just-in-time compilation to speed initial load time
compiler::enableJIT(0)

# Load packages ----
source("iPSCeq-pack-load.R")

# Load custom fxns
source("iPSCeq-functions.R")

# source custom modules
source("iPSCeq-modules.R")

# Load tabs used in ui.R
source("iPSCeq-tabs.R")

# Load user interface
source("iPSCeq-ui.R")

# Load server side interface
source("iPSCeq-server.R")

# load credentials required
if(!getOption("standalone")) source("sql.R")

# load embedded fonts
font_add_google(name = "Noto Sans JP", family = "noto-sans-jp")
showtext_auto(enable = TRUE)

# Run shiny ----
shinyApp(iPSCeqUI, iPSCeqServer)
