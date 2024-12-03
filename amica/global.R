require(shiny)
require(shinyjs)
require(shinyBS)
require(shinyalert)
require(ggplot2)
require(ggfortify)
require(plotly)
require(heatmaply)
require(DT)
require(tools)
require(DEqMS)
require(gprofiler2)
require(bslib)
require(reshape2)
require(igraph)
require(visNetwork)
require(UpSetR)
require(colourvalues)
require(pheatmap)
require(Rmisc)
require(data.table)
require(dplyr)
require(RColorBrewer)
require(colourpicker)
require(eulerr)
require(cowplot)
require(tidyr)
require(Cairo)
require(rmarkdown)
###
source("R/utils.R")
source("R/ProteomicsData.R")



rawPrefix="^Intensity."
intensityPrefix="LFQ.intensity."
abundancePrefix="iBAQ"
razorUniqueCount="razorUniqueCount"
razorUniqueCountPrefix="razorUniqueCount."
spectraCount="spectraCount"
proteinId="Majority.protein.IDs"
geneName="Gene.names"
contaminantCol="Potential.contaminant"
pvalPrefix="P.Value_"
padjPrefix="adj.P.Val_"
logfcPrefix="logFC_"
avgExprPrefix="AveExpr_"
filterVal="+"

removePlotlyBars <- list(
  'sendDataToCloud',
  'autoScale2d',
  'zoomIn2d',
  'zoomOut2d',
  'toggleSpikelines',
  'hoverClosestCartesian',
  'hoverCompareCartesian'
)

amicaGlobalVars <- new.env()
local <- new.env()

# splash screen
local$defaultSplashScreenEnabled <- FALSE
amicaGlobalVars$splashScreenEnabled <- as.logical(Sys.getenv("AMICA_SPLASHSCREEN_ENABLED", local$defaultSplashScreenEnabled))
local$defaultSplashScreenTitle <- "Welcome to amica"
amicaGlobalVars$splashScreenTitle <- as.character(Sys.getenv("AMICA_SPLASHSCREEN_TITLE", local$defaultSplashScreenTitle))
local$defaultSplashScreenHtmlContentSource <- "mpl/mplSplashScreen.html"
amicaGlobalVars$splashScreenHtmlContentSource <- as.character(Sys.getenv("AMICA_SPLASHSCREEN_CONTENT_SOURCE", local$defaultSplashScreenHtmlContentSource))

# default footer / info
# if using images in the footer and the source is a file, it needs to be added to the www folder under amica
local$defaultOverrideFooter <- FALSE
amicaGlobalVars$overrideFooter <- as.logical(Sys.getenv("AMICA_FOOTER_CUSTOMIZE", local$defaultOverrideFooter))
local$defaultFooterHtmlContentSource <- "default/footer.html"
if (amicaGlobalVars$overrideFooter == TRUE){
  amicaGlobalVars$footerHtmlContentSource <- as.character(Sys.getenv("AMICA_FOOTER_CONTENT_SOURCE", local$defaultFooterHtmlContentSource))
} else {
  amicaGlobalVars$footerHtmlContentSource <-local$defaultFooterHtmlContentSource 
}


# logout user after n seconds/minutes...etc.
# https://stackoverflow.com/questions/33839543/shiny-server-session-time-out-doesnt-work/53207050
local$defaultTimeoutSeconds <- 60*15
local$timeoutSeconds <- strtoi(Sys.getenv("AMICA_SHINY_IDLE_TIMEOUT_SECONDS", local$defaultTimeoutSeconds))

local$defaultAmicaVersion <- "3.0.1"
amicaGlobalVars$amicaVersion <- Sys.getenv("AMICA_VERSION_OVERRIDE", local$defaultAmicaVersion)

local$defaultAmicaSourceHyperlink <- "https://www.github.com/tbaccata/amica"
amicaGlobalVars$amicaSourceHyperlink <- Sys.getenv("AMICA_SOURCE_HYPERLINK", local$defaultAmicaSourceHyperlink)

# will be used to replace every text matching {{placeholder_x}} (incl braces) with the assigned value
amicaGlobalVars$placeholder_replacement_map <- list(
  placeholder_amica_version = amicaGlobalVars$amicaVersion, 
  placeholder_repository_link = amicaGlobalVars$amicaSourceHyperlink
  )


inactivity <- sprintf("function idleTimer() {
var t = setTimeout(logout, %s);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
Shiny.setInputValue('timeOut', '%ss')
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();", local$timeoutSeconds*1000, local$timeoutSeconds, local$timeoutSeconds*1000)



inline = function (x) {
  tags$div(style="display:inline-block;", x)
}


# replace placeholders inside footer if they exist:
local$footerHtmlContent <- amicaUtil$replacePlaceholders(
  paste(readLines(amicaGlobalVars$footerHtmlContentSource)), 
  amicaGlobalVars$placeholder_replacement_map
  )

# used to insert the same, resuable footer to various Shiny panels used throughout amica
# shiny syntax requires it to be provided as a function, not a value, e.g.
# tabPanel(
#  title = 'My Title',
#  value = '',
#  HTML('<p>Text</p>'),
#  footer()
# )
footer = function(x) {
  tags$div(
    style="footer{position: absolute; bottom:5%; left: 33%; padding:5px;}",
    HTML(local$footerHtmlContent))
}

