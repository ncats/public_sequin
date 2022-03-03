# Set local data directory for running as standalone app
options(localDir = "example_data")

################################################################################
# Run as standalone app or connect to existing database
options(standalone = T)

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
