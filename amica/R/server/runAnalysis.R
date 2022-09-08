output$filterValuesInput <- renderUI({
  req(reacValues$uniqueGroups)
  selectizeInput(
    "filterValuesInput",
    "Select groups to be considered for filtering on valid values.
      If none is selected all groups are considered.",
    c("", reacValues$uniqueGroups),
    multiple = T,
    selected = NULL
  )
})

output$intensitySelection <- renderUI({
  req(reacValues$proteinData)
  intensities <- c()
  selected <-  "LFQIntensity"
  
  if (reacValues$dbTool == "maxquant") {
    intensities <- assayNames(reacValues$proteinData)
  } else {
    intensities <- assayNames(reacValues$proteinData)
    
    if ("LFQIntensity" %in% intensities) {
      selected <-  "LFQIntensity"
    } else {
      selected <-
        ifelse("RazorIntensity" %in% intensities,
               "RazorIntensity",
               "Intensity")
    }
  }
  intensities <- intensities[! intensities %in% "ImputedIntensity"]
  
  selectizeInput(
    "quantIntensity",
    "Which intensities should be quantified?",
    intensities,
    multiple = F,
    selected = selected
  )
})

### ANALYSIS
observeEvent(input$runAnalysis, {
  if (is.null(reacValues$proteinData)) {
    showNotification(paste0("No data uploaded."), type = "error")
    return(NULL)
  }
  reacValues$combinedData = NULL
  reacValues$dataCompAmica = NULL
  reacValues$amicaInput = FALSE
  
  if (reacValues$amicaInput == FALSE) {
    ###FILTDATA BEGIN
    reacValues$inputParameterSummary <- NULL
    
    quantIntensity <- "LFQIntensity"
    
    if (!is.null(reacValues$dbTool)) {
      if (is.null(input$quantIntensity) ||
          length(input$quantIntensity) < 2) {
        if (reacValues$dbTool == "fragpipe") {
          quantIntensity <-
            ifelse(
              "RazorIntensity" %in% assayNames(reacValues$proteinData),
              "RazorIntensity",
              "Intensity"
            )
          
          if (!quantIntensity %in% assayNames(reacValues$proteinData)) {
            quantIntensity <- assayNames(reacValues$proteinData)[1]
          }
          
        }
      } else {
        quantIntensity <- input$quantIntensity
      }
      reacValues$inputParameterSummary <- paste0(
        reacValues$inputParameterSummary,
        "DB search tool:\t",
        reacValues$dbTool,
        "\n"
      )
    } else {
      reacValues$inputParameterSummary <- paste0(
        reacValues$inputParameterSummary,
        "DB search tool:\tunknown\n"
      )
    }
    
    reacValues$inputParameterSummary <-
      paste0(
        reacValues$inputParameterSummary,
        "Intensities used for quantification:\t",
        quantIntensity,
        "\n"
      )
    ### filter on values
    impDf <- assay(reacValues$proteinData, quantIntensity)
    reacValues$proteinData <-
      setAssay(x = reacValues$proteinData,
               assay = impDf,
               assayName = "LFQIntensity")
    
    
    rnames <- filterOnMinValuesRnames(
      y = reacValues$proteinData,
      minMSMS = input$minMSMS,
      minRazor = input$minRazor
    )
    
    reacValues$inputParameterSummary <-
      paste0(
        reacValues$inputParameterSummary,
        "Minimum MS/MS counts:\t",
        input$minMSMS,
        "\n",
        "Minimum razor/unique peptides:\t",
        input$minRazor,
        "\n"
      )

    tmp <- tryCatch({
      filterOnValidValues(
        impDf[rnames, ],
        colData(reacValues$proteinData),
        input$filterValuesInput,
        minValue = input$minValidValue,
        method = input$validValuesGroup
      )
    },
    error = function(cond) {
      message(paste(cond))
      shinyalert(
        title = "Filter on valid values failed.",
        text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    },
    warning = function(cond) {
      message(paste(cond))
    }
    ,
    finally = {
      showNotification(paste("Filtering values..."), type = "message")
    })

    reacValues$inputParameterSummary <-
      paste0(
        reacValues$inputParameterSummary,
        "Filter on groups:\t",
        paste(ifelse(is.null(input$filterValuesInput), "all", 
                     paste(input$filterValuesInput, collapse = ',') ),
              collapse = ","),
        "\n",
        "Filter on min. value in group:\t",
        input$minValidValue,
        "\n",
        "Filter on min. value in:\t",
        input$validValuesGroup,
        "\n"
      )
    
    rowsDf <- rowData(reacValues$proteinData)
    rowsDf$quantified <- ""
    rowsDf[tmp, 'quantified'] <- "+"
    
    reacValues$proteinData <-
      setRowData(reacValues$proteinData, rowsDf)
    
    normDf <- impDf[tmp,]
    
    
    # renormalization
    if (input$renormalizationMethod != "None") {
      normDf <- tryCatch({
        renormalizeIntensities(normDf, input$renormalizationMethod)
      },
      error = function(cond) {
        message(paste(cond))
        shinyalert(
          title = "Normalization failed.",
          text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
          size = "s", 
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = FALSE,
          type = "error",
          showConfirmButton = TRUE,
          showCancelButton = FALSE,
          confirmButtonText = "OK",
          confirmButtonCol = "#669966",
          timer = 0,
          imageUrl = "",
          animation = TRUE
        )
      },
      warning = function(cond) {
        message(paste(cond))
      }
      ,
      finally = {
        showNotification(paste("Normalizing intensities..."), type = "message")
      })
    }
    
    reacValues$inputParameterSummary <-
      paste0(
        reacValues$inputParameterSummary,
        "Re-normaliztion method:\t",
        input$renormalizationMethod,
        "\n"
      )
    
    normDf <- tryCatch({
      imputeIntensities(
        normDf,
        method = input$impMethod,
        downshift = input$downshift,
        width = input$width
      )
    },
    error = function(cond) {
      message(paste(cond))
      shinyalert(
        title = "Imputation failed.",
        text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    },
    warning = function(cond) {
      message(paste(cond))
    }
    ,
    finally = {
      showNotification(paste("Imputing intensities..."), type = "message")
      #message("Imputing intensities...", )
    })
    
    reacValues$inputParameterSummary <-
      paste0(
        reacValues$inputParameterSummary,
        "Imputation method:\t",
        input$impMethod,
        "\n"
      )
    
    if (input$impMethod == 'normal' ||
        input$impMethod == 'global') {
      reacValues$inputParameterSummary <-
        paste0(
          reacValues$inputParameterSummary,
          "Downshift:\t",
          input$downshift,
          "\nWidth:\t",
          input$width,
          "\n"
        )
    } else if (input$impMethod == 'min') {
      reacValues$inputParameterSummary <-
        paste0(
          reacValues$inputParameterSummary,
          "Min Intensity:\t",
          min(normDf, na.rm = T),
          "\n"
        )
    }
    
    impDf[row.names(normDf), ] <- normDf
    reacValues$proteinData <-
      setAssay(x = reacValues$proteinData,
               assay = impDf,
               assayName = "ImputedIntensity")
    
    
    # filtered row data as df
    reacValues$filtData <-
      rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData),]
    
    if (is.null(reacValues$filtData)) {
      shinyalert(
        title = "No data left to analyze after set input parameters:",
        text = paste(reacValues$inputParameterSummary, sep = '\n\n',
                     paste0('#of quantified proteins: ', 
                            length(isQuantRnames(reacValues$proteinData)) )),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    }
    
    normDf <- NULL
    tmp <- NULL
    impDf <- NULL
    
    
    pep.count.table <- NULL
    if (input$limmaTrend) {
      countIdx <-
        grep("razorUniqueCount.", colnames(reacValues$filtData))
      pep.count.table = data.frame(
        count = apply(reacValues$filtData[, countIdx], 1, FUN = min),
        row.names = rownames(reacValues$filtData)
      )
      pep.count.table$count = pep.count.table$count + 1
    }
    
    isPilot <-
      ifelse(any(duplicated(colData(
        reacValues$proteinData
      )$groups)), FALSE, TRUE)
    tool <- ifelse(input$limmaTrend, "DEqMS", "limma")
    if (isPilot) tool <- "None"
    reacValues$daTool <- tool
    
    out <- tryCatch({
      tmpOut <- groupComparisons(
        as.matrix(assay(
          reacValues$proteinData, "ImputedIntensity"
        )[isQuantRnames(reacValues$proteinData), ]),
        reacValues$contrastMatrix,
        reacValues$expDesign,
        pep.count.table
      )
    },
    error = function(cond) {
      # showNotification(paste("error in", tool, cond),
      #                  type = "error",
      #                  duration = 100)
      shinyalert(
        title = paste(tool, "failed."),
        text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    }
    #,
    #warning = function(cond) {
      #out <- tmpOut
      #message(paste(cond))
      # showNotification(
      #   paste("warning in", tool, "\nResults may be invalid!\nPlease deselect DEqMS button to continue.\n", cond),
      #   type = "warning",
      #   duration = 100
      # )
      # shiny:::reactiveStop(conditionMessage(cond))
    #}
    )
    
    reacValues$inputParameterSummary <-
      paste0(
        reacValues$inputParameterSummary,
        "Differential abundance statistics:\t",
        tool,
        "\n"
      )
    
    out$Gene.names <- reacValues$filtData[, geneName, drop = T]
    
    tmp <- grep(logfcPrefix, colnames(out), value = T)
    reacValues$reacConditions <- gsub(logfcPrefix, "", tmp)
    reacValues$dataLimma <- out[, c(ncol(out), 1:(ncol(out) - 1))]
    reacValues$dataLimmaOriginal <- reacValues$dataLimma
    
    gc()
  }
  
  reacValues$analysisSuccess <- TRUE
  
  reacValues$nsubmits <- reacValues$nsubmits + 1
  if (reacValues$nsubmits < 2) {
    #toggle('ibaq_help')
    toggle(selector = "#navbar li a[data-value=qctab]")
    toggle(selector = "#navbar li a[data-value=quanttab]")
    toggle(selector = "#navbar li a[data-value=comparemicatab]")
    toggle(id = 'hide_before_input', anim = T)
  }
})