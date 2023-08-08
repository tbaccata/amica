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

# logout user after n seconds/minutes...etc.
# https://stackoverflow.com/questions/33839543/shiny-server-session-time-out-doesnt-work/53207050
timeoutSeconds <- 60*15

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
idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)



inline = function (x) {
  tags$div(style="display:inline-block;", x)
}

footer = function(x) {
  tags$div(
    style="footer{position: absolute; bottom:5%; left: 33%; padding:5px;}",
    HTML("
    <h5>How to cite us</h5>
         <p>Please cite Didusch, S., Madern, M., Hartl, M., & Baccarini, M. BMC Genomics 23, 817 (2022). 
         amica: an interactive and user-friendly web-platform for the analysis of proteomics data.
         DOI: <a href='https://doi.org/10.1186/s12864-022-09058-7' 
         target='_blank'>https://doi.org/10.1186/s12864-022-09058-7</a>.\n
         </p>
         <p>All code and online documentation can be found on 
         <a href='https://www.github.com/tbaccata/amica' target='_blank'>github</a>.</p>
         <p>amica version 3.0.1
         <img src='maxperutzlabs.jpg' width='100px'>
         </p>
         ")
  )
}
