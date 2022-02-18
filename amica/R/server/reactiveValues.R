reacValues <-
  reactiveValues(
    proteinData = NULL,
    expDesign=NULL,
    uploadSuccess=NULL,
    analysisSuccess=NULL,
    inputParameterSummary=NULL,
    filtData = NULL,
    dataLimma = NULL,
    overlapDf = NULL,
    dataLimmaOriginal = NULL,
    dataComp = NULL,
    geneNames = NULL,
    reacConditions = NULL,
    groupFactors = NULL,
    uniqueGroups = NULL,
    selectHeatmapGroups = NULL,
    dataHeatmap = NULL,
    dataDotplot = NULL,
    dotplotGroupsDf = NULL,
    dotplotFactors = NULL,
    GostPlot = NULL,
    clusterRows = NULL,
    clusterCols = NULL,
    scaleHeatmap = NULL,
    combinedData = NULL,
    dataCompAmica = NULL,
    amicaInput = FALSE,
    show_analysis = FALSE,
    networkData = NULL,
    nsubmits = 0,
    show_heatmap = FALSE,
    show_ora = FALSE,
    show_single = FALSE,
    show_multi = FALSE,
    selection = NULL,
    fcCutoff = NULL,
    enrichmentChoice = NULL,
    sigCutoffValue = NULL,
    nwNodes = NULL,
    nwEdges = NULL,
    newMultiNames = NULL,
    compareAmicaSelected = FALSE,
    compareAmicasToggled = FALSE
  )
nclicks <- reactiveVal(0)



### RESET 

observeEvent(input$resetAnalysis,{
  reacValues$proteinData <- NULL
  reacValues$amicaInput <- FALSE
  reacValues$dbTool <- NULL
  reacValues$inputParameterSummary <- NULL
  reacValues$dataHeatmap <- NULL
  reacValues$dotplotGroupsDf <- NULL
  reacValues$dataDotplot <- NULL
  reacValues$dotplotFactors <- NULL
  reacValues$GostPlot <- NULL
  reacValues$expDesign <- NULL
  reacValues$contrastMatrix <- NULL
  reacValues$uploadSuccess <- NULL
  reacValues$analysisSuccess <- NULL
  reacValues$filtData <- NULL
  reacValues$overlapDf <- NULL
  reacValues$dataLimma <- NULL
  reacValues$dataLimmaOriginal <- NULL
  reacValues$dataComp <- NULL
  reacValues$geneNames <- NULL
  reacValues$reacConditions <- NULL
  reacValues$uniqueGroups <- NULL
  reacValues$groupFactors <- NULL
  reacValues$selection <- NULL
  reacValues$newMultiNames <- NULL
  reacValues$dataCompAmica <- NULL
  reacValues$dataGprofiler <- NULL
  reacValues$compareAmicaSelected <- FALSE
  reacValues$compareAmicasToggled <- FALSE
})