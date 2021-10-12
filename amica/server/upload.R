### UPLOAD
observeEvent(input$submitAnalysis, {
  reacValues$inputParameterSummary <- NULL
  ### EXAMPLE
  if (input$source == "example") {
    sourcePath <- "data/PXD0016455/"
    
    tmpData <-
      read.table(
        paste0(sourcePath, "design.txt"),
        header = T,
        stringsAsFactors = F
      )
    reacValues$expDesign <- tmpData
    
    withProgress(message = "Reading in amica file", {
      outData <-
        readInAmicaSumm(paste0(sourcePath, "amica_proteinGroups.tsv"),
                        reacValues$expDesign)
      reacValues$proteinData <- outData$protData
      reacValues$contrastMatrix = outData$contrasts
      reacValues$dataLimma = outData$comparisons
      reacValues$dataLimmaOriginal <- reacValues$dataLimma
      
      ### filtData
      reacValues$filtData <-
        rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData),]
    })
    
    reacValues$amicaInput = TRUE
    reacValues$analysisSuccess <- TRUE
    
    comps <-
      grep(logfcPrefix, colnames(outData$comparisons), value = T)
    reacValues$reacConditions <- gsub(logfcPrefix, "", comps)
  }
  ###
  
  if (input$source != "example" &
      (is.null(input$groupSpecification))) {
    showNotification(paste("Need to upload experimental design."), type = "message")
    return("")
  }
  
  # EXPERIMENTAL DESIGN
  if (input$source != "example" &
      !is.null(input$groupSpecification))  {
    inFile <- input$groupSpecification
    tmpData <- validateFile(inFile, c("groups", "samples"))
    reacValues$expDesign <- tmpData
  }
  
  if (input$source == "amica") {
    if (is.null(input$amicaFile)) {
      showNotification(paste("Need to upload amica output."), type = "error")
      return(NULL)
    }
    
    withProgress(message = "Reading in amica file", {
      outData <-
        readInAmicaSumm(input$amicaFile$datapath, reacValues$expDesign)
      reacValues$proteinData <- outData$protData
      reacValues$contrastMatrix = outData$contrasts
      reacValues$dataLimma = outData$comparisons
      reacValues$dataLimmaOriginal <- reacValues$dataLimma
      
      reacValues$amicaInput = TRUE
      reacValues$analysisSuccess <- TRUE
      
      comps <-
        grep(logfcPrefix, colnames(outData$comparisons), value = T)
      reacValues$reacConditions <- gsub(logfcPrefix, "", comps)
      
      ### filtData
      reacValues$filtData <-
        rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData),]
    })
  }
  
  if (input$source == "custom") {
    if (is.null(input$customFile)) {
      showNotification(paste("Need to upload custom tab-separated file."), type = "error")
      return(NULL)
    }
    
    if (is.null(input$specFile)) {
      showNotification(paste("Need to upload specification file to map relevant columns."),
                       type = "error")
      return(NULL)
    }
    
    specs <- validateFile(input$specFile, c("Variable", "Pattern"))
    reacValues$proteinData <-
      readInCustomSumm(input$customFile$datapath, specs, reacValues$expDesign)
    showNotification(paste("Reading in custom file format ..."), type = "message")
  }
  
  if (input$source == "maxquant") {
    if (is.null(input$maxquantFile)) {
      showNotification(paste("Need to upload MaxQuant or FragPipe output."), type = "error")
      return(NULL)
    }
    
    reacValues$uniqueGroups <- unique(reacValues$expDesign$groups)
    
    header <-
      unlist(strsplit(readLines(input$maxquantFile$datapath, n = 1), "\t"))
    header <- make.names(header)
    mqNames <-
      c("Majority.protein.IDs",
        "Gene.names",
        "Razor...unique.peptides")
    fragNames <-
      c("Indistinguishable.Proteins",
        "Protein.ID",
        "Protein.Probability")
    
    if (all(mqNames %in% header)) {
      reacValues$proteinData <-
        readInMQproteinGroupsSumm(input$maxquantFile$datapath,
                                  reacValues$expDesign)
      reacValues$dbTool <- "maxquant"
      showNotification(paste("Reading in MaxQuant ..."), type = "message")
    } else if (all(fragNames %in% header)) {
      reacValues$proteinData <-
        readInFragPipeProteinGroupsSumm(input$maxquantFile$datapath,
                                        reacValues$expDesign)
      reacValues$dbTool <- "fragpipe"
      showNotification(paste("Reading in FragPipe ..."), type = "message")
    } else {
      showNotification(
        paste0(
          "Unrecognized input format\n",
          "\tallowed input files:\n",
          "\tMaxQuant: proteinGroups.txt\n",
          "\tFragPipe: combined_protein.tsv"
        ),
        type = "error"
      )
      stop("Unrecognized input format")
    }
    
  }
  
  # contrasts
  if (input$source != "example" & input$source != "amica") {
    if (is.null(input$contrastMatrix)) {
      showNotification(paste("Need to upload contrast matrix."), type = "error")
      return(NULL)
    }
    
    inFile <- input$contrastMatrix
    tmpData <- validateFile(inFile, NULL)
    reacValues$contrastMatrix <- tmpData
    
    for (elem in unique(c(tmpData[[1]], tmpData[[2]]))) {
      if (!is.null(reacValues$expDesign) &
          length(grep(elem, reacValues$expDesign)) < 1) {
        shiny:::reactiveStop(showNotification(
          paste(
            "Group",
            elem,
            "in contrasts not in uploaded experimental design"
          ),
          type = "warning",
          duration = 100
        ))
      }
    }
    
  } else {
    reacValues$nsubmits <- reacValues$nsubmits + 1
    if (reacValues$nsubmits < 2) {
      #toggle('ibaq_help')
      toggle(selector = "#navbar li a[data-value=qctab]")
      toggle(selector = "#navbar li a[data-value=quanttab]")
      toggle(selector = "#navbar li a[data-value=comparemicatab]")
      toggle(id = 'hide_before_input', anim = T)
    }
  }
  
  reacValues$uploadSuccess <- TRUE
})