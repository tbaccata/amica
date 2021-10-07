require(shiny)
library(shinyjs)
require(shinyBS)
require(ggplot2)
require(ggfortify)
require(plotly)
require(heatmaply)
require(DT)
require(tools)
require(DEqMS)
require(gprofiler2)
library(bslib)
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
#require(pryr)
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

