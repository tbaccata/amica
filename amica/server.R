source("global.R")
# library(pryr)

server <- function(input, output, session) {
  # max 75MB file upload restriction
  options(shiny.maxRequestSize = 75 * 1024 ^ 2)
  
  reacValues <-
    reactiveValues(
      proteinData = NULL,
      expDesign=NULL,
      uploadSuccess=NULL,
      analysisSuccess=NULL,
      inputParameterSummary=NULL,
      filtData = NULL,
      dataLimma = NULL,
      dataLimmaOriginal = NULL,
      dataComp = NULL,
      geneNames = NULL,
      reacConditions = NULL,
      uniqueGroups = NULL,
      selectHeatmapGroups = NULL,
      dataHeatmap = NULL,
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
      compareAmicaSelected = FALSE,
      compareAmicasToggled = FALSE
    )
  nclicks <- reactiveVal(0)
  

  observe({
    hide(selector = "#navbar li a[data-value=qctab]")
    hide(selector = "#navbar li a[data-value=quanttab]")
    hide(selector = "#navbar li a[data-value=comparemicatab]")
  })
  
  shinyjs::onclick("toggleAdvancedHeatmap",
                   shinyjs::toggle(id = "advancedHeatmap"))
  
  shinyjs::onclick("heatmapParams",
                   shinyjs::toggle(id = "toggle_heatmap_params"))
  
  shinyjs::onclick("toggleCompareAmicaInput",
                   shinyjs::toggle(id = "compareAmicaInput"))
  
  
  ### PLOT PARAMS ONCLICK
  shinyjs::onclick("pcaParams",
                   shinyjs::toggle(id = "toggle_pca_params"))
  
  shinyjs::onclick("boxplotParams",
                   shinyjs::toggle(id = "toggle_boxplot_params"))
  
  shinyjs::onclick("densityParams",
                   shinyjs::toggle(id = "toggle_density_params"))
  
  shinyjs::onclick("corParams",
                   shinyjs::toggle(id = "toggle_cor_params"))
  
  shinyjs::onclick("cvParams",
                   shinyjs::toggle(id = "toggle_cv_params"))
  
  shinyjs::onclick("contaminantsParams",
                   shinyjs::toggle(id = "toggle_contaminants_params"))
  
  shinyjs::onclick("abundantParams",
                   shinyjs::toggle(id = "toggle_abundant_params"))
  
  shinyjs::onclick("barplotIdParams",
                   shinyjs::toggle(id = "toggle_barplotId_params"))
  
  shinyjs::onclick("barplotMvParams",
                   shinyjs::toggle(id = "toggle_barplotMv_params"))
  
  shinyjs::onclick("scatterParams",
                   shinyjs::toggle(id = "toggle_scatter_params"))
  
  shinyjs::onclick("volcanoParams",
                   shinyjs::toggle(id = "toggle_volcano_params"))
  
  shinyjs::onclick("maParams",
                   shinyjs::toggle(id = "toggle_ma_params"))
  
  shinyjs::onclick("fcParams",
                   shinyjs::toggle(id = "toggle_fc_params"))
  
  shinyjs::onclick("profileParams",
                   shinyjs::toggle(id = "toggle_profile_params"))
  
  
  shinyjs::onclick("oraBarParams",
                   shinyjs::toggle(id = "toggle_oraBar_params"))
  
  shinyjs::onclick("fcamicaParams",
                   shinyjs::toggle(id = "toggle_fcamica_params"))
  
  shinyjs::onclick("scatteramicaParams",
                   shinyjs::toggle(id = "toggle_scatteramica_params"))
  
  shinyjs::onclick("upsetParams",
                   shinyjs::toggle(id = "toggle_upset_params"))
  
  shinyjs::onclick("corAmicasParams",
                   shinyjs::toggle(id = "toggle_corAmicas_params"))

  ### RESET 
  
  observeEvent(input$resetAnalysis,{
    reacValues$proteinData <- NULL
    reacValues$amicaInput <- FALSE
    reacValues$dbTool <- NULL
    reacValues$inputParameterSummary <- NULL
    reacValues$dataHeatmap <- NULL
    reacValues$GostPlot <- NULL
    reacValues$expDesign <- NULL
    reacValues$contrastMatrix <- NULL
    reacValues$uploadSuccess <- NULL
    reacValues$analysisSuccess <- NULL
    reacValues$filtData <- NULL
    reacValues$dataLimma <- NULL
    reacValues$dataLimmaOriginal <- NULL
    reacValues$dataComp <- NULL
    reacValues$geneNames <- NULL
    reacValues$reacConditions <- NULL
    reacValues$uniqueGroups <- NULL
    reacValues$selection <- NULL
    reacValues$dataCompAmica <- NULL
    reacValues$dataGprofiler <- NULL
    reacValues$compareAmicaSelected <- FALSE
    reacValues$compareAmicasToggled <- FALSE
  })
  
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
          rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData), ]
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
          rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData), ]
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
        c(
          "Indistinguishable.Proteins",
          "Protein.ID",
          "Protein.Probability"
        )
      
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
        if (!is.null(reacValues$expDesign) & length(grep(elem, reacValues$expDesign)) < 1) {
          shiny:::reactiveStop(showNotification(
            paste("Group", elem, "in contrasts not in uploaded experimental design"),
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
      selected <-
        ifelse("RazorIntensity" %in% intensities,
               "RazorIntensity",
               "Intensity")
    }
    
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
    
    if (reacValues$amicaInput == FALSE) {
      ###FILTDATA BEGIN
      reacValues$inputParameterSummary <- NULL
      
      quantIntensity <- "LFQIntensity"
      
      if (!is.null(reacValues$dbTool)) {
        if (is.null(input$quantIntensity) | length(input$quantIntensity) < 2) {
          if (reacValues$dbTool == "fragpipe") {
            quantIntensity <- ifelse("RazorIntensity" %in% assayNames(reacValues$proteinData),
                                     "RazorIntensity",
                                     "Intensity")
          }
        } else {
          quantIntensity <- input$quantIntensity
        }
      }

      reacValues$inputParameterSummary <- paste0(reacValues$inputParameterSummary,
                                                 "Intensities used for quantification:\t", quantIntensity, "\n")
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
      
      reacValues$inputParameterSummary <- paste0(reacValues$inputParameterSummary,
                                                 "Minimum MS/MS counts:\t", input$minMSMS, "\n",
                                                 "Minimum razour/unique peptides:\t", input$minRazor, "\n")

      tmp <- tryCatch({
        filterOnValidValues(
          impDf[rnames,],
          colData(reacValues$proteinData),
          input$filterValuesInput,
          minValue = input$minValidValue,
          method = input$validValuesGroup
        )
      },
      error = function(cond) {
        message(cond)
      },
      warning = function(cond) {
        message(cond)
      }
      ,
      finally = {
        showNotification(paste("Filtering values..."), type = "message")
      }
      )
      
      reacValues$inputParameterSummary <- paste0(reacValues$inputParameterSummary,
                                                 "Filter on groups:\t", paste(input$filterValuesInput, collapse = ","), "\n",
                                                 "Filter on min. value in group:\t", input$minValidValue, "\n",
                                                 "Filter on min. value in:\t",input$validValuesGroup, "\n")

      rowsDf <- rowData(reacValues$proteinData)
      rowsDf$quantified <- ""
      rowsDf[tmp, 'quantified'] <- "+" 
      
      reacValues$proteinData <- setRowData(reacValues$proteinData, rowsDf)
      
      normDf <- impDf[tmp, ]
      
      
      
      # renormalization
      if (input$renormalizationMethod != "None") {
        
        normDf <- tryCatch({
          renormalizeIntensities(normDf, input$renormalizationMethod)
        },
        error = function(cond) {
          message(cond)
        },
        warning = function(cond) {
          message(cond)
        }
        ,
        finally = {
          showNotification(paste("Normalizing intensities..."), type = "message")
        }
        )
      }
      
      reacValues$inputParameterSummary <- paste0(reacValues$inputParameterSummary,
                                                 "Re-normaliztion method:\t",
                                                 input$renormalizationMethod, "\n")
      
      normDf <- tryCatch({
        imputeIntensities(
          normDf,
          method = input$impMethod,
          downshift = input$downshift,
          width = input$width
        )
      },
      error = function(cond) {
        message(cond)
      },
      warning = function(cond) {
        message(cond)
      }
      ,
      finally = {
        showNotification(paste("Imputing intensities..."), type = "message")
        #message("Imputing intensities...", )
      }
      )
      
      reacValues$inputParameterSummary <- paste0(reacValues$inputParameterSummary,
                                                 "Imputation method:\t",
                                                 input$impMethod,
                                                 "\nDownshift:\t",input$downshift,
                                                 "\nWidth:\t",input$width,"\n")
      
      impDf[row.names(normDf),] <- normDf
      reacValues$proteinData <- setAssay(x = reacValues$proteinData, assay = impDf, assayName = "ImputedIntensity")


      # filtered row data as df
      reacValues$filtData <-
        rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData), ]
      

      normDf <- NULL
      tmp <- NULL
      impDf <- NULL
      
      
      pep.count.table <- NULL
      if (input$limmaTrend) {
        countIdx <- grep("razorUniqueCount.", colnames(reacValues$filtData ) )
        pep.count.table = data.frame(count = apply(reacValues$filtData[,countIdx],1,FUN=min),
                                     row.names = rownames(reacValues$filtData ) )
        pep.count.table$count = pep.count.table$count+1
      }

      tool <- ifelse(input$limmaTrend, "DEqMS", "limma")
      
      out <- tryCatch({
        tmpOut <- groupComparisons(
              as.matrix(
                assay(reacValues$proteinData, "ImputedIntensity")[isQuantRnames(reacValues$proteinData),]
              ),
              reacValues$contrastMatrix,
              reacValues$expDesign,
              pep.count.table
            )
      },
      error=function(cond) {
        showNotification(paste("error in", tool, cond), type = "error", duration = 100)
      },
      warning=function(cond) {
        showNotification(
          paste("warning in", tool, "\nResults may be invalid!\nPlease deselect DEqMS button to continue.\n", cond),
          type = "warning",
          duration = 100
        )
        shiny:::reactiveStop(conditionMessage(cond))
      },finally = invalidateLater(1)
      )
      
      reacValues$inputParameterSummary <- paste0(reacValues$inputParameterSummary,
                                                 "Differential abundance statistics:\t",
                                                 tool,"\n")
      
      out$Gene.names <- reacValues$filtData[,geneName, drop=T]
      
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
  
  output$brewercols <- renderPlot({
    display.brewer.all()
  })
  
  
  output$brewerOptionsQual <- renderUI({
    req(reacValues$proteinData)
    selectInput(
      "brewerOptionsQual",
      "Choose a color palette for qualitative data:",
      row.names(brewer.pal.info),
      selected = "Set2"
    )
  })
  
  output$brewerOptionsDiv <- renderUI({
    req(reacValues$proteinData)
    selectInput(
      "brewerOptionsDiv",
      "Choose a color palette for diverging data:",
      row.names(brewer.pal.info[brewer.pal.info$category != "qual", ]),
      selected = "RdYlBu"
    )
  })
  
  myColors <- reactive({
    color <-
      ifelse(is.null(input$brewerOptionsQual),
             "Accent",
             input$brewerOptionsQual)
    reverse <- ifelse(input$revQual == "yes", TRUE, FALSE)
    nc <- nrow(colData(reacValues$proteinData))
    
    mycols <- c()
    maxcol <- brewer.pal.info[color, "maxcolors"]
    if (maxcol < nc) {
      mycols <- colorRampPalette(brewer.pal(maxcol, color))(nc)
    } else {
      nc <- ifelse(nc < 3, 3, nc)
      mycols <- brewer.pal(nc, color)
    }
    if (reverse)
      mycols <- rev(mycols)
    mycols
  })
  
  mySelectionColors <- reactive({
    color <- ifelse(is.null(input$brewerOptionsQual), "Set2", input$brewerOptionsQual)
    reverse <- ifelse(input$revQual=="yes", TRUE, FALSE)
    #nc <- length(unique(reacValues$proteinData$groups))
    
    nc <- length(unique(colData(reacValues$proteinData)$groups))
    
    mycols <- c()
    maxcol <- brewer.pal.info[color, "maxcolors"]
    if (maxcol < nc) {
      mycols <- colorRampPalette(brewer.pal(maxcol, color))(nc)
    } else {
      nc <- ifelse(nc<3, 3, nc)
      mycols <- brewer.pal(nc, color)
    }
    
    names(mycols) <- unique(colData(reacValues$proteinData)$groups)
    
    if (reverse) mycols <- rev(mycols)
    mycols
  })
  
  output$myColorPanel <- renderUI({
    req(reacValues$proteinData)
    
    groups <- unique(colData(reacValues$proteinData)$groups)
    
    lapply(seq_along(groups), function(i) {
      colourInput(paste("col", i, sep="_"), paste0(groups[i]), mySelectionColors()[i])
    })
    })
  
  output$volcanoMAColors <- renderUI({
    pal <- brewer.pal(3, "Set2")
    
    lapply(seq_along(1:2), function(i) {
      colourInput(paste("volcanocol", i, sep="_"), paste0("Color ", i, ":"), pal[i])
    })
  })
  
  output$fcPlotColors <- renderUI({
    pal <- myScatterColors()
    
    lapply(seq_along(pal), function(i) {
      colourInput(paste("fcplotcol", i, sep="_"), paste0("Color ", i, ":"), pal[i])
    })
  })
  
  myGroupColors <- reactive({
    groups <- unique(colData(reacValues$proteinData)$groups)
    
    colors <- lapply(seq_along(groups), function(i) {
      input[[paste("col", i, sep="_")]]
    })
    
    if (is.null(input$col_1)) colors <- mySelectionColors()
    
    names(colors) <- groups
    colors
  })

  myScatterColors <- reactive({
    color <- ifelse(is.null(input$brewerOptionsScatter), "Set2", input$brewerOptionsScatter)
    reverse <- ifelse(input$revScatter=="yes", TRUE, FALSE)
    
    pal <- brewer.pal(5, color)
    if (reverse) pal <- rev(pal)
    pal
  })
  
  output$mtcarsBar <- renderPlotly({
    req(reacValues$proteinData)
    tmp <- mtcars
    tmp$car <- row.names(tmp)
    p <-
      ggplot(tmp[1:min(nc <-
                         length(unique(colData(
                           reacValues$proteinData
                         )$groups)), nrow(tmp)), ], aes(x = car, y = mpg, fill = car)) +
      geom_bar(stat = "Identity") + 
      theme_minimal(base_size = 14) +
      scale_fill_manual(values=myColors() ) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 10
      ),legend.title = element_blank()) + xlab("")
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = list(
        'sendDataToCloud',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'toggleSpikelines',
        'hoverClosestCartesian',
        'hoverCompareCartesian'
      ),
      toImageButtonOptions = list(
        format = "svg",
        width = 676,
        height = 676,
        filename = "mtcars_bar_example"
      )
    )
  })
  
  heatColors <- reactive({
    color <- ifelse(is.null(input$brewerOptionsDiv), "RdYlBu", input$brewerOptionsDiv)
    maxcol <- brewer.pal.info[color, "maxcolors"]
    
    reverse <- ifelse(input$revDiv=="yes", TRUE, FALSE)
    
    if (reverse) {
      pal <- colorRampPalette(rev(brewer.pal(n = maxcol, name =
                                               color)))(100)
    } else {
      pal <- colorRampPalette(brewer.pal(n = maxcol, name =
                                           color))(100)
    }
    pal
  })
  
  output$mtcarsHeatmap <- renderPlotly({
    heatmaply(
      t(mtcars),
      scale = "row",
      key.title = 'Z-score',
      row_dend_left = T,
      showticklabels = c(T, T),
      plot_method = "plotly",
      colors = heatColors()
    ) %>%
      config(displaylogo = F,
             modeBarButtonsToRemove = list(
               'sendDataToCloud',
               'autoScale2d',
               'zoomIn2d',
               'zoomOut2d',
               'toggleSpikelines',
               'hoverClosestCartesian',
               'hoverCompareCartesian'
             ),
             toImageButtonOptions = list(format = "svg",
                                         width = 676,
                                         height = 676,
                                         filename = "heatmap_example")
      )
  })
  
  output$irisScatter <- renderPlotly({
    req(reacValues$proteinData)

    p <-
      ggplot(iris, aes(x=Sepal.Length, y=Petal.Length, color=Species)) + geom_point() +
      geom_point() + theme_minimal(base_size = 14) +  scale_color_manual(values=myScatterColors() )
    
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = list(
        'sendDataToCloud',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'toggleSpikelines',
        'hoverClosestCartesian',
        'hoverCompareCartesian'
      ),
      toImageButtonOptions = list(
        format = "svg",
        width = 676,
        height = 676,
        filename = "iris_example"
      )
    )
  })
  
  
  # collapse
  observeEvent(input$p1Button, ({
    updateCollapse(session, "collapseExample", open = "Advanced")
  }))
  
  
  output$expDesignDT <- renderDT({
    req(reacValues$expDesign )
    datatable(
      reacValues$expDesign,
      caption = "Experimental design",
      rownames = F,
      options = list(
        #dom = 't',
        selected = rownames(reacValues$expDesign )[1:nrow(reacValues$expDesign )]
      )
    )
  })
  
  output$amicaInput <- reactive({
    req(reacValues$uploadSuccess)
    reacValues$amicaInput
  })
  outputOptions(output, "amicaInput", suspendWhenHidden = FALSE)
  
  output$uploadSummary <- renderText({
    req(reacValues$uploadSuccess)
    
    if (!is.null(reacValues$analysisSuccess) || reacValues$amicaInput == TRUE ||
        input$source == "example"
        ) return(NULL)
    
    paste0(
      "Successfully uploaded data!\n",
      "Open the 'Advanced' tab in the sidebar to chose parameters."
    )
  })
  
  output$uploadSuccessMsg <- renderUI({
    req(reacValues$uploadSuccess)
    if (!is.null(reacValues$analysisSuccess)) return(NULL)
    
    HTML(
      "<h4>Successfully uploaded data!</h4>
      <p>
      If you have uploaded an amica file format you can inspect your data 
      in the main tab bar. <br>
      Otherwise open the 'Analysis' section on the sidebar and press 'Analyze'.
      </p>"
    )
  })
  
  output$summaryText <- renderText({
    req(reacValues$analysisSuccess)
    
    if (contaminantCol %in% names(rowData(reacValues$proteinData))) {
      nsummary <-
        nrow(rowData(reacValues$proteinData)[rowData(reacValues$proteinData)[[contaminantCol]] != filterVal, ])
      
    } else {
      nsummary <- nrow(rowData(reacValues$proteinData))
    }
    
    
    
    ngroups <- length(unique(reacValues$expDesign$groups))
    ncomps <- nrow(reacValues$contrastMatrix )
    
    nquant <- length(which(rowData(reacValues$proteinData)$quantified=="+"))
    
    paste0(
      "Number of proteins in file (wo. contaminants, reverse proteins, only identified by site):\n",
      nsummary,
      "\nNumber of quantified proteins:\n",
      nquant,
      "\nNumber of conditions: \n",
      ngroups,
      "\nNumber of group comparisons: \n",
      ncomps,
      "\n\n\nYou can now inspect\n",
      "\t1) QC\n",
      "\t2) Quantitative results\n",
      "\t3) Protein-protein interaction networks (only applicable for H.sapiens at the moment)\n",
      "\t4) or upload another amica file to compare experiments!"
    )
  })
  
  
  output$quantSummary <- renderText({
    req(reacValues$proteinData)
    numQuant <- length(isQuantRnames(reacValues$proteinData))
    paste0("Number of quantified proteins in ImputedIntensity: ", numQuant)
  })
  
  output$inputParameterSummary <- renderText({
    req(reacValues$analysisSuccess)
    req(reacValues$inputParameterSummary)
    reacValues$inputParameterSummary
  })
  
  # ------------------------------------------- intensity boxplots
  
  dataAssay <- reactive({
    req(input$assayNames)
    
    object <- 0
    if (input$assayNames=="ImputedIntensity") {
      object <- assay(reacValues$proteinData, "ImputedIntensity")
      
      #object[isQuantRnames(reacValues$proteinData),]
      
      object <- toLongFormat(
        object[isQuantRnames(reacValues$proteinData),],
        reacValues$proteinData,
        addGroup = TRUE,
        addContaminant = FALSE,
        addGeneName = FALSE
      )
    } else {
      object <- toLongFormat(
        assay(reacValues$proteinData, input$assayNames),
        reacValues$proteinData,
        addGroup = TRUE,
        addContaminant = FALSE,
        addGeneName = FALSE
      )
    }
    object
  })
  
  
  # ------------------------------------------- PCA
  
  pcaPlotly <- eventReactive(input$submitPCA,{
    assayNames <- isolate(input$assayNames)
    groupInputs <- isolate(input$pcaSamplesInput)
    
    validate(need(assayNames != "", "Please provide an intensity."))
    validate(need(length(groupInputs) > 1 | is.null(groupInputs), "Please select at least two groups. If none is selected all are considered."))
    
    if (is.null(groupInputs)) groupInputs <- unique(colData(reacValues$proteinData)$groups)
    
    shinyjs::show('hide_before_submit_pca')
    plotData <- 0
    tmpCols <- colData(reacValues$proteinData)
    if (assayNames == "ImputedIntensity") {
      plotData <-
        assay(reacValues$proteinData, "ImputedIntensity")[isQuantRnames(reacValues$proteinData),
                                                          tmpCols$samples[tmpCols$groups %in% groupInputs]]
    } else {
      plotData <- assay(reacValues$proteinData, assayNames)
      
      plotData <- plotData[,tmpCols$samples[tmpCols$groups %in% groupInputs]]
      plotData <- plotData[complete.cases(plotData), ]
    }
    
    withProgress(message = "Plotting PCA ", {
      p <-
        plotPCA(
          plotData,
          reacValues$expDesign,
          myGroupColors(),
          input$pca_base,
          input$pca_legend,
          input$pca_pointsize,
          input$pca_show_label
        )
    })
    
    
    # p
    p
  })
  
  
  output$pcaSamplesInput <- renderUI({
    req(input$assayNames)
    selectizeInput(
      "pcaSamplesInput",
      "Press on 'Plot PCA'.
      Leave the input blank to plot all samples.
      Only select specific groups when you have really have to.",
      unique(reacValues$expDesign$groups),
      multiple = T,
      options = list(minItems = 2)
    )
  })
  
  
  output$pca <- renderPlotly({
    req(input$assayNames)
    req(reacValues$uploadSuccess)
    pcaPlotly()  %>%  config(displaylogo = F,
                             modeBarButtonsToRemove = list(
                               'sendDataToCloud',
                               'autoScale2d',
                               'zoomIn2d',
                               'zoomOut2d',
                               'toggleSpikelines',
                               'hoverClosestCartesian',
                               'hoverCompareCartesian'
                             ),
                             toImageButtonOptions = list(format = input$pca_format,
                                                         width = input$pca_width,
                                                         height = input$pca_height,
                                                         #width = 768,
                                                         #height = 676,
                                                         filename = "pca")
    )
  })
  
  ### ------------------------------------------- BOX PLOT
  
  boxplotPlotly <- eventReactive(input$submitBoxplot,{
    p <- ggplot(dataAssay()[!is.na(dataAssay()$value), ], aes(x = colname, y = value, fill =
                                                                group)) +
      geom_boxplot(outlier.shape = NA,
                   outlier.color = NULL,
                   outlier.fill = NULL) +
      theme_minimal(base_size = input$boxplot_base) +
      scale_fill_manual(values = myGroupColors()) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 10
        ),
        legend.text = element_text(size = input$boxplot_legend),
        legend.title = element_blank()
      ) + ylab("Intensity (log2)") + xlab("")
    
    ggplotly(p)
  })
  
  output$boxPlot <- renderPlotly({
    req(reacValues$uploadSuccess)
    withProgress(message = "Plotting boxplot of intensities", {
      boxplotPlotly()  %>%
        config(
          displaylogo = F,
          modeBarButtonsToRemove = list(
            'sendDataToCloud',
            'autoScale2d',
            'zoomIn2d',
            'zoomOut2d',
            'toggleSpikelines',
            'hoverClosestCartesian',
            'hoverCompareCartesian'
          ),
          toImageButtonOptions = list(
            format = input$boxplot_format,
            width = input$boxplot_width,
            height = input$boxplot_height,
            filename = "boxplot"
          )
        )
    })
  })
  
  ### ------------------------------------------- DENSITY PLOT
  
  densityPlotly <- eventReactive(input$submitDensity,{
    set.seed(123)
    df <- dataAssay()[sample(nrow(dataAssay()), nrow(dataAssay())%/%3  ),]
    
    p <- ggplot(df[!is.na(df$value), ], aes(x = value, color = colname)) +
      # ggplot(dataAssay()[!is.na(dataAssay()$value), ], aes(x = value, color =
      #                                                        colname)) +
      geom_density() +
      scale_color_manual(values = myColors()) +
      theme_minimal(base_size = input$density_base) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 10
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = input$density_legend),
      ) + xlab("intensity (log2)") + ylab("Density")
    
    ggplotly(p) 
  })

  output$densityPlot <- renderPlotly({
    req(reacValues$uploadSuccess)

    withProgress(message = "Plotting density plot", {
      densityPlotly() %>%
        config(
          displaylogo = F,
          modeBarButtonsToRemove = list(
            'sendDataToCloud',
            'autoScale2d',
            'zoomIn2d',
            'zoomOut2d',
            'toggleSpikelines',
            'hoverClosestCartesian',
            'hoverCompareCartesian'
          ),
          toImageButtonOptions = list(
            format = input$density_format,
            width = input$density_width,
            height = input$density_height,
            filename = "density_plot"
          )
        )
    })
  })
  
  # ------------------------------------------- Input
  
  output$samples <- renderUI({
    req(reacValues$proteinData)
    selectInput("samples", "Inspect a sample", c("",colData(reacValues$proteinData)$samples), selected = NULL)
  })
  
  output$assayNames <- renderUI({
    req(reacValues$proteinData)
    selectInput("assayNames",
                "Which intensities do you want to plot?",
                c(assayNames(reacValues$proteinData)),
                selected = "ImputedIntensity")
  })
  
  output$comparisons <- renderUI({
    req(reacValues$dataLimma)
    comps <-
      colnames(reacValues$dataLimma)[grep(logfcPrefix, colnames(reacValues$dataLimma))]
    comps <- gsub(logfcPrefix, "", comps)
    selectInput("comparisons", "Inspect a group comparison", comps)
  })
  
  output$compareComparisons <- renderUI({
    req(reacValues$dataLimma)
    compareComparisons <-
      colnames(reacValues$dataLimma)[grep(logfcPrefix, colnames(reacValues$dataLimma))]
    compareComparisons <- gsub(logfcPrefix, "", compareComparisons)
    
    selectInput(
      "compareComparisons",
      "Compare (multiple) group comparison(s)",
      compareComparisons,
      multiple = T
    )
  })
  
  output$compareScatterPlot1 <- renderUI({
    selectizeInput(
      "selectScatterSample1",
      "Select a sample for the x-axis of the scatter plot",
      c("",colData(reacValues$proteinData)$samples),
      selected = NULL,
      options = list(maxItems = 1)
    )
  })
  
  output$compareScatterPlot2 <- renderUI({
    selectizeInput(
      "selectScatterSample2",
      "Select a sample for the y-axis of the scatter plot",
      c("",colData(reacValues$proteinData)$samples),
      selected = NULL,
      options = list(maxItems = 1)
    )
  })
  
  
  output$assayName1 <- renderUI({
    req(reacValues$proteinData)
    selectInput(
      "assayName1",
      "Which intensities to plot at x-axis?",
      c("", assayNames(reacValues$proteinData)),
      selected = NULL,
      multiple = F
    )
  })
  output$assayName2 <- renderUI({
    req(reacValues$proteinData)
    selectInput(
      "assayName2",
      "Which intensities to plot at y-axis?",
      c("", assayNames(reacValues$proteinData)),
      selected = NULL,
      multiple = F
    )
  })
  
  
  output$scatterPlot <- renderPlotly({
    
    input$submitScatterPlot
    
    req(reacValues$uploadSuccess)
    
    selection1 <- isolate(input$selectScatterSample1)
    selection2 <- isolate(input$selectScatterSample2)
    assay1 <- isolate(input$assayName1)
    assay2 <- isolate(input$assayName2)
    
    plot_width <- isolate(input$scatter_width)
    plot_height <- isolate(input$scatter_height)
    plot_fontsize <- isolate(input$scatter_base)
    plot_legendsize <- isolate(input$scatter_legend)
    scatter_format <- isolate(input$scatter_format)
    showLine <- isolate(input$scatter_showLine)
    
    validate(need(
      !is.null(selection1) & !is.null(selection2),
      "Please enter 2 samples to compare."
    ))
    validate(need(
      !is.null(assay1),
      "Please enter an assay for x-axis."
    ))
    validate(need(
      !is.null(assay2),
      "Please enter an assay for y-axis."
    ))
    
    xLabel <- paste0(assay1, "_", selection1)
    yLabel <- paste0(assay2, "_", selection2)
    plotData <- data.frame()
    
    if (assay1 == assay2) {
      
      validate(need(selection1 != selection2, "Cannot plot the same column on x - and y-axis."))
      
      plotData <- assay(reacValues$proteinData, assay1)[, c(selection1, selection2)]
      names(plotData) <- c("x", "y")
      
    } else {
      plotData <-
        assay(reacValues$proteinData, assay1)[, selection1, drop = F]
      names(plotData) <- c("x")
      plotData[[selection2]] <-
        assay(reacValues$proteinData, assay2)[, selection2, drop = T]
    }
    names(plotData)[2] <- "y"
    
    plotData$Contaminant <- "no"
    
    if (contaminantCol %in% names(rowData(reacValues$proteinData))) {
      plotData$Contaminant <- ifelse(
        rowData(reacValues$proteinData)[[contaminantCol]]=="+",
        "yes",
        "no"
      )
    }
    
    plotData$Gene <- 
      rowData(reacValues$proteinData)[[geneName]]
    
    fit1 <- lm(y~x, data=plotData)
    fit1.intercept <- fit1$coefficients[[1]]
    fit1.slope <- fit1$coefficients[[2]]
    
    title <-
      paste0(signif(fit1$coef[[1]], 5), 
             "x + ", 
             signif(fit1$coef[[2]], 5),
             ", r2 = ",
             signif(summary(fit1)$r.squared, 5))
    
    withProgress(message = "Plotting scatter plot", {
      pu <-
        ggplot(plotData,
               aes(
                 x=x, 
                 y=y,
                 label = Gene,
                 color = Contaminant
               )) + theme_minimal(base_size = plot_fontsize) + labs(x=xLabel, y=yLabel)  + 
        theme(legend.title = element_text(size=plot_legendsize),
              legend.text=element_text(size=plot_legendsize)) +
        geom_point() + scale_color_manual(values=myScatterColors() ) # + ggtitle(title) #scale_color_brewer(palette = "Paired")
      
      intercept <- 0
      slope <- 1
      if (showLine == "linear regression") {
        intercept <- fit1.intercept
        slope <- fit1.slope
        pu <- pu + ggtitle(title)
      } else if (showLine == "straight line") {
        intercept <- 0
        slope <- 1
      }
      
      if (showLine != "none") {
        pu <- pu + geom_abline(
          intercept = intercept,
          slope = slope,
          size = 1,
          alpha = 0.5,
          color = myScatterColors()[3]
        )
      }  
    })
    
    m = list(
      l = 100,
      r = 40,
      b = 100,
      t = 50,
      pad = 0
    )
    
    ggplotly(pu, margin = m) %>%
    #ggplotly(pu, width = plot_width, height = plot_height, margin = m, autosize=T) %>%
      config(
        displaylogo = F,
        modeBarButtonsToRemove = list(
          'sendDataToCloud',
          'autoScale2d',
          'zoomIn2d',
          'zoomOut2d',
          'toggleSpikelines',
          'hoverClosestCartesian',
          'hoverCompareCartesian'
        ),
        toImageButtonOptions = list(
          format = scatter_format,
          width = plot_width,
          height = plot_height,
          filename = "scatter_plot"
        )
      )
  })
  
  output$corrSamplesInput <- renderUI({
    req(input$assayNames)
    selectizeInput(
      "corrSamplesInput",
      "Press on 'Plot Correlation'.
      Leave the input blank to plot all samples.
      Only select specific groups when you have really have to.",
      unique(reacValues$expDesign$groups),
      multiple = T,
      options = list(minItems = 2)
    )
  })
  
  corrBasePlot <- eventReactive(input$submitCor,{
    assayName <- isolate(input$assayNames)
    annotSamples <- isolate(input$cor_annot)
    groupInputs <- isolate(input$corrSamplesInput)
    tmpCols <- colData(reacValues$proteinData)
    
    req(reacValues$uploadSuccess)
    
    validate(need(assayName != "", "Please provide an intensity."))
    
    validate(need(length(groupInputs) > 1 | is.null(groupInputs), "Please select at least two groups. If none is selected all are considered."))
    
    if (is.null(groupInputs)) groupInputs <- unique(colData(reacValues$proteinData)$groups)
    
    df <- assay(reacValues$proteinData, assayName)
    df <- df[,tmpCols$samples[tmpCols$groups %in% groupInputs]]
    
    annot <- colData(reacValues$proteinData)
    row.names(annot) <- annot$samples 
    annot <- annot[names(df),]
    annot$samples <- NULL
    
    corDf <- cor(df, method = "pearson", use = "complete.obs")
    
    if (annotSamples) {
      p <- heatmaply_cor(
        round(corDf, 3),
        xlab = "", 
        ylab = "",
        limits = c(min(corDf)-0.003, 1),
        #row_side_palette = mapping,
        row_side_palette = myGroupColors(),
        row_side_colors = annot,
        plot_method = "plotly",
        key.title = "Pearson Correlation"
      )
    } else {
      p <- heatmaply_cor(
        round(corDf, 3),
        xlab = "", 
        ylab = "",
        limits = c(min(corDf)-0.003, 1),
        plot_method = "plotly",
        key.title = "Pearson Correlation"
      )
    }
    p
  })
  
  output$corrPlotly <- renderPlotly({

    withProgress(message = "Plotting correlation plot ", {
    p <- corrBasePlot()
    })
    p %>%  config(displaylogo = F,
                  modeBarButtonsToRemove = list(
                    'sendDataToCloud',
                    'autoScale2d',
                    'zoomIn2d',
                    'zoomOut2d',
                    'toggleSpikelines',
                    'hoverClosestCartesian',
                    'hoverCompareCartesian'
                  ),
                  toImageButtonOptions = list(format = input$cor_format,
                                              width = input$cor_width,
                                              height = input$cor_height,
                                              filename = "corrplot")
    )
  })
  
  
  ### ------------------------------------------- CONTAMINANTS plots
  
  contaminantsPlotly <- eventReactive(input$submitConts, {
    validate(need(
      abundancePrefix %in% assayNames(reacValues$proteinData),
      "No data to plot (no iBAQ columns available)."
    ))
    
    tmp <- 2 ^ assay(reacValues$proteinData, abundancePrefix)
    allInts <- apply(tmp, 2, sum, na.rm = T)
    contsInts <-
      apply(tmp[rowData(reacValues$proteinData)[[contaminantCol]] == "+", ], 2, sum, na.rm =
              T)
    contsInts <- contsInts / allInts
    
    contsInts <- reshape2::melt(contsInts)
    midx <-
      match(row.names(contsInts), colData(reacValues$proteinData)$samples)
    contsInts$group <- colData(reacValues$proteinData)$groups[midx]
    contsInts$Sample <- row.names(contsInts)
    contsInts$value <- 100 * contsInts$value
    
    p <- ggplot(contsInts, aes(x = Sample, y = value, fill = group)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = input$contaminants_base) +
      scale_fill_manual(values = myGroupColors()) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = input$contaminants_base
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = input$contaminants_legend)
      ) + xlab("") + ylab("% Contaminant")
    
    ggplotly(p)
  })
  
  
  output$contaminants <- renderPlotly({
    req(reacValues$uploadSuccess)
    withProgress(message = "Plotting barplot of contaminants", {
      contaminantsPlotly()  %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = list(
          'sendDataToCloud',
          'autoScale2d',
          'zoomIn2d',
          'zoomOut2d',
          'toggleSpikelines',
          'hoverClosestCartesian',
          'hoverCompareCartesian'
        ),
        toImageButtonOptions = list(
          format = input$contaminants_format,
          width = input$contaminants_width,
          height = input$contaminants_height,
          filename = "contaminants_plot"
        )
      )
    })
  })
  
  
  ### ------------------------------------------- MOST ABUNDANT per sample barplot
  
  output$mostAbundant <- renderPlotly({
    req(input$samples)
    req(reacValues$uploadSuccess)
    validate(need( abundancePrefix %in% assayNames(reacValues$proteinData), "No data to plot (no iBAQ columns available)." ))
    
    tmp <- object <- toLongFormat(
      2^assay(reacValues$proteinData, abundancePrefix),
      reacValues$proteinData,
      addGroup = TRUE,
      addContaminant = TRUE,
      addGeneName = TRUE
    )
    
    tmp <- tmp[tmp$colname==input$samples,]
    tmp <- tmp[order(tmp$value, decreasing = T),]
    
    tmp$value <- 100 * tmp$value/sum(tmp$value,na.rm = T)
    
    withProgress(message = "Plotting most abundant proteins", {
      p <- ggplot(head(tmp,15), aes(x=reorder(Gene.name, -value), y=value)) + 
        geom_bar(stat="identity", fill=myScatterColors()[1]) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 10
        ),
        legend.text = element_text(size = input$abundant_legend),
        legend.title = element_blank() ) + xlab("") + ylab("%Signal of protein")
    })
    
    ggplotly(p) %>% config(displaylogo = F,
                           modeBarButtonsToRemove = list(
                             'sendDataToCloud',
                             'autoScale2d',
                             'zoomIn2d',
                             'zoomOut2d',
                             'toggleSpikelines',
                             'hoverClosestCartesian',
                             'hoverCompareCartesian'
                           ),
                           toImageButtonOptions = list(format = input$abundant_format,
                                                       width = 676,
                                                       height = 676,
                                                       filename = "plot")
    )
  })
  
  ### ------------------------------------------- NUM ID PROTEINS
  
  numIdPlotly <- eventReactive(input$submitNumProts, {
    validate(
      need(
        "Intensity" %in% assayNames(reacValues$proteinData) |
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        "No data to plot (no columns starting with 'LFQIntensity_' or 'Intensity_' available)."
      )
    )
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    
    
    df <-
      as.data.frame(apply(assay(reacValues$proteinData, aname), 2, function(x)
        sum(!is.na(x))))
    names(df) <- "Number"
    df$Sample <- row.names(df)
    
    midx <- match(df$Sample, colData(reacValues$proteinData)$samples)
    df$group <- colData(reacValues$proteinData)$groups[midx]
    
    p <- ggplot(df, aes(x = Sample, y = Number, fill = group)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = input$barplotId_base) +
      scale_fill_manual(values = myGroupColors()) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = input$barplotId_base
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = input$barplotId_legend)
      ) + xlab("") + ylab("#of identified proteins")
    
    ggplotly(p)
  })
  
  output$barplotProteins <- renderPlotly({
    req(reacValues$uploadSuccess)
    withProgress(message = "Plotting Number of inferred protein groups", {
      numIdPlotly()  %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = list(
          'sendDataToCloud',
          'autoScale2d',
          'zoomIn2d',
          'zoomOut2d',
          'toggleSpikelines',
          'hoverClosestCartesian',
          'hoverCompareCartesian'
        ),
        toImageButtonOptions = list(
          format = input$barplotId_format,
          width = input$barplotId_width,
          height = input$barplotId_height,
          filename = "barplot_identified_proteins"
        )
      )
    })
  })
  
  ### ------------------------------------------- % MVs
  
  pctMvsPlotly <- eventReactive(input$submitNumMVs, {
    validate(
      need(
        "Intensity" %in% assayNames(reacValues$proteinData) ||
          "LFQIntensity" %in% assayNames(reacValues$proteinData),
        "No data to plot (no columns starting with 'LFQIntensity_' or 'Intensity_' available)."
      )
    )
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    
    df <-
      as.data.frame(apply(assay(reacValues$proteinData, aname), 2, function(x) {
        100 * sum(is.na(x)) / length(x)
      }))
    names(df) <- "Number"
    df$Sample <- row.names(df)
    
    midx <- match(df$Sample, colData(reacValues$proteinData)$samples)
    df$group <- colData(reacValues$proteinData)$groups[midx]
    
    p <- ggplot(df, aes(x = Sample, y = Number, fill = group)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = myGroupColors()) +
      theme_minimal(base_size = input$barplotMv_base) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = input$barplotMv_base
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = input$barplotMv_legend)
      ) + xlab("") + ylab("missing values (%)")
    
    ggplotly(p) 
  })
  
  output$barplotMissingValues <- renderPlotly({
    req(reacValues$uploadSuccess)
    withProgress(message = "Plotting barplot of missing values ", {
      pctMvsPlotly() %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = list(
          'sendDataToCloud',
          'autoScale2d',
          'zoomIn2d',
          'zoomOut2d',
          'toggleSpikelines',
          'hoverClosestCartesian',
          'hoverCompareCartesian'
        ),
        toImageButtonOptions = list(
          format = input$barplotMv_format,
          width = input$barplotMv_width,
          height = input$barplotMv_height,
          filename = "barplot_missing_values"
        )
      )
    })
  })
  
  ### -------------------------------------------- CVs BOX PLOT
  
  cvsPlotly <- eventReactive(input$submitCVs, {
    calcData <- dataAssay()
    calcData$value <- 2^calcData$value
    cvs <-
      aggregate(value ~ rowname + group, calcData, function(x)
        100 * (sd(x) / mean(x)) )
    
    p <- ggplot(cvs[!is.na(cvs$value),], aes(x=group, y=value, fill=group)) + 
      geom_boxplot(outlier.shape = NA,
                   outlier.color = NULL,
                   outlier.fill = NULL) +
      scale_fill_manual(values = myGroupColors() ) + 
      theme_minimal(base_size = input$cv_base) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = input$cv_base
      ),
      legend.text = element_text(size = input$cv_legend),
      legend.title = element_blank()) + ylab("Coefficient of Variation (%)") + xlab("")
    
    ggplotly(p)
  })
  
  output$boxplotCV <- renderPlotly({
    req(reacValues$uploadSuccess)
    validate(need(any(duplicated(colData(reacValues$proteinData)$groups)), "Cannot output CV plot for a pilot without replicates")) 
    #validate(need(input$assayName != "", "Please provide an intensity."))
    withProgress(message = "Plotting Coefficient of Variations ", {
      cvsPlotly() %>% config(displaylogo = F,
                             modeBarButtonsToRemove = list(
                               'sendDataToCloud',
                               'autoScale2d',
                               'zoomIn2d',
                               'zoomOut2d',
                               'toggleSpikelines',
                               'hoverClosestCartesian',
                               'hoverCompareCartesian'
                             ),
                             toImageButtonOptions = list(format = input$cv_format,
                                                         width = input$cv_width,
                                                         height = input$cv_height,
                                                         filename = "cv_plot")
      )
    })
  })
  
  comparisons <- reactive({input$compareComparisons})
  
  
  observeEvent(input$submitHeatmap,{
    req(reacValues$dataComp)
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    ridx <- rownames(reacValues$dataComp[ridx,])
    
    validate(need(nrow(reacValues$dataComp)>0, "No data to plot"))
    df <- assay(reacValues$proteinData, "ImputedIntensity")
    df <- df[ridx,]
    
    if (!any(duplicated(rowData(reacValues$proteinData)[ridx, geneName]))) {
      rownames(df) <-
        as.character(rowData(reacValues$proteinData)[ridx, geneName] )
    } else {
      rownames(df) <-
        make.names(rowData(reacValues$proteinData)[ridx, geneName] , unique = T)
    }
    
    if (!is.null(input$heatmapSamplesInput) & length(input$heatmapSamplesInput) > 1) {
      idxs <- c()
      for (elem in input$heatmapSamplesInput) {
        relSamples <- reacValues$expDesign$samples[reacValues$expDesign$groups==elem]
        idxs <- c(idxs, grep( paste0(relSamples, collapse = "|"), colnames(df) ) )
      }
      df <- df[, idxs]
    }
    reacValues$clusterRows = input$clusterRows
    reacValues$clusterCols = input$clusterCols
    reacValues$scaleHeatmap = input$scaleHeatmap
    reacValues$dataHeatmap <- df
    
    if (reacValues$show_heatmap == FALSE) {
      reacValues$show_heatmap = TRUE
      toggle(id = 'hide_heatmap_before_submit', anim = T)
    }
  })
  
  
  output$groupComparisonsDT <- renderDT({
    req(reacValues$dataComp)
    matSet <- isolate(enrichedMatrixSet())
    
    validate(need(
      nrow(reacValues$dataComp) >= 1,
      paste0(
        "Need more than ",
        nrow(reacValues$dataComp),
        " proteins selected. Apply less stringent thresholds."
      )
    ))
    
    selection <- gsub(logfcPrefix, "", grep(logfcPrefix, names(reacValues$dataComp), value = T) )
    
    tmp <- matSet[row.names(matSet) %in% 
                                 row.names(reacValues$dataComp),selection, drop = F]
    toAdd <- as.data.frame(ifelse(tmp==1, "yes", "no"))
    names(toAdd) <- paste0("significant_", names(toAdd))
    tmp <- reacValues$dataComp
    
    idx <- ncol(tmp)
    tmp <- cbind(tmp, toAdd)
    tmp <- tmp[, c(1,(idx+1):ncol(tmp), 2:idx)]
    
    hasPadj <- ifelse(length(grep(padjPrefix, names(reacValues$dataComp)))>0, TRUE, FALSE )
      
    dt <- datatable(
      tmp,
      filter = "top",
      rownames = F,
      extensions = c('Buttons'), 
      options = list(
        dom = 'Bfrtip',
        # searching = FALSE,
        pageLength = 10,
        autoWidth = TRUE,
        search = list(regex = TRUE),
        buttons = c('csv')
      )
    ) %>% formatRound(paste0(logfcPrefix,selection), 4) 
    
    if (hasPadj) dt <- dt %>% formatSignif(grep("P.Val", names(tmp)), digits = 4)
    dt
  }, server = F)
  
  
  enrichedMatrixSet <- reactive({
    req(reacValues$dataLimma)
    
    #enrichmentChoice <- isolate(input$enrichmentChoice)
    #sigCutoffValue <- isolate(input$sigCutoffValue)
    #fcCutoff <- isolate(input$fcCutoff)
    
    matrixData <-
      generateEnrichedMatrix(
        reacValues$dataLimma,
        #reacValues$dataLimma[rownames(reacValues$dataComp), ],
        input$enrichmentChoice,
        input$sigCutoffValue,
        input$fcCutoff
      )
    as.data.frame(matrixData)
  })
  
  
  observeEvent(input$submitMultiComp,{
    req(enrichedMatrixSet())
    
    req(input$upset1Sample)

    if (reacValues$show_multi == FALSE) {
      reacValues$show_multi = TRUE
      toggle(id = 'hide_before_multi_submit', anim = T)
    }

    reacValues$dataComp <-
      getComparisonsData2(reacValues$dataLimma,
                          enrichedMatrixSet(),
                          'union',
                          input$upset1Sample)
  })
  
  output$upset1Sample <- renderUI({
    req(reacValues$reacConditions)
    selectInput(
      "upset1Sample",
      "Compare the overlap between comparisons.",
      c("",reacValues$reacConditions),
      selected = NULL,
      multiple = T
    )
  })
  
  
  plotMultiUpset <- function() {
    matrixSet <- isolate(enrichedMatrixSet() )
    samples <- isolate(input$upset1Sample)
    setChoice <- "union"
    fcThresh <- isolate(input$fcCutoff)
    enrichmentChoice <- isolate(input$enrichmentChoice)
    sigCutoffValue <- isolate(input$sigCutoffValue)
    
    reacValues$sigCutoffValue <- sigCutoffValue
    reacValues$selection <- samples
    reacValues$fcCutoff <- fcThresh
    reacValues$enrichmentChoice <- enrichmentChoice
    reacValues$setChoice <- setChoice
    
    upset_ratio <- isolate(input$upset_ratio)
    upset_pointsize <- isolate(input$upset_pointsize)
    upset_sorted <- isolate(input$upset_sorted)
    
    upset_sorted <- ifelse(upset_sorted == 'Degree', 'degree', 'freq')
    
    
    validate(need(!is.null(samples), "Please select at least two comparisons."))
    validate(need(length(colnames(matrixSet) ) > 1, "Need at least two comparisons to render UpSet plot. Only one provided.") )
    if (!is.null(samples)) {
      validate(need(length(samples) > 1, "Please select at least two comparisons.") )
    }
    validate(need(any(matrixSet[,samples]==1),
                  "There are no significant proteins to display."))
    
    mb.ratio <- c(upset_ratio, 1-upset_ratio)
    
    # newNames <- c()
    # for (elem in names(matrixSet)) {
    #   newName <- paste(unlist(strsplit(elem, "__vs__") ), collapse = "/")
    #   newNames <- c(newNames, newName)
    # }
    # names(matrixSet) <- newNames
    
    upset(
      matrixSet,
      sets = samples,
      mb.ratio = mb.ratio,
      #set_size.numbers_size = 8,
      #set_size.show = T,
      order.by = upset_sorted,
      text.scale = c(2, 2,
                     2, 2,
                     2, 3),
      point.size = upset_pointsize
    )
  }
  
  
  output$upsetPlot <- renderPlot({
    input$submitMultiComp
    #input$plotMultiComp
    req(reacValues$dataComp)


    withProgress(message = "Comparing multigroup comparisons", {

      if (reacValues$show_analysis == FALSE) {
        reacValues$show_analysis = TRUE
        toggle(id = 'hide_before_comparisons', anim = T)
      }
      
      print(plotMultiUpset())
    })
  })
  
  output$download_button_upset <- renderUI({
    #req(upsetPlot())
    downloadButton("download_upset", "Download UpSet plot")
  })
  
  output$download_upset <- downloadHandler(
    filename = function() {
      paste("upset_plot.pdf")
    },
    content = function(file) {
      pdf(file, width = input$upset_width, height = input$upset_height, onefile = F)
      print(plotMultiUpset())
      dev.off()
    }
  )
  
  compareHeatmapBase <- eventReactive(input$submitHeatmap,{

    validate(need(reacValues$compareAmicaSelected==FALSE,
                  paste0("Cannot output heatmap for 'multiple_amica_upload'\n
                         Change the data set back to the original data to 
                         plot a heatmap.") ))
    
    
    
    
    
    req(reacValues$dataHeatmap)
    fontsize = isolate(input$heatmap_base)
    plot_width = isolate(input$heatmap_width)
    plot_height = isolate(input$heatmap_height)
    show_labels <- isolate(input$heatmap_labels)
    show_annot <- isolate(input$heatmap_annot)
    format <- isolate(input$heatmap_format)
    
    validate(need(all(names(reacValues$dataHeatmap) %in% colData(reacValues$proteinData)$samples),""))
    
    validate(need(
      nrow(reacValues$dataHeatmap) >= 1,
      paste0(
        "Need more than ",
        nrow(reacValues$dataHeatmap),
        " proteins selected. Apply less stringent thresholds."
      )
    ))
    annot <- colData(reacValues$proteinData)
    row.names(annot) <- annot$samples
    annot <- annot[names(reacValues$dataHeatmap),]
    annot$samples <- NULL
    
    cbarTitle <-
      ifelse(reacValues$scaleHeatmap != "none", "Z-score", "Intensity (log2)")
    
    withProgress(message = "Plotting heatmap ", {
      p <- 0
      
      if (show_annot) {
        p <- heatmaply(
          reacValues$dataHeatmap,
          Rowv = reacValues$clusterRows,
          Colv = reacValues$clusterCols,
          scale = reacValues$scaleHeatmap,
          plot_method = "plotly",
          fontsize_row = fontsize,
          fontsize_col = fontsize,
          showticklabels = c(TRUE, show_labels),
          col_side_palette = myGroupColors(),
          col_side_colors = annot,
          row_dend_left = TRUE,
          column_text_angle = 90,
          key.title = cbarTitle,
          colors = heatColors()
        )
      } else {
        p <- heatmaply(
          reacValues$dataHeatmap,
          Rowv = reacValues$clusterRows,
          Colv = reacValues$clusterCols,
          scale = reacValues$scaleHeatmap,
          plot_method = "plotly",
          fontsize_row = fontsize,
          fontsize_col = fontsize,
          showticklabels = c(TRUE, show_labels),
          row_dend_left = TRUE,
          column_text_angle = 90,
          key.title = cbarTitle,
          colors = heatColors()
        )
      }
      
      calc_height <- min(15 * nrow(reacValues$dataHeatmap), 1200)
      if (calc_height < 600) calc_height <- 600
      
      p %>% layout(height=calc_height)
    })
  })
  
  output$compareHeatmap <- renderPlotly({

    compareHeatmapBase() %>%
      config(displaylogo = F,
             modeBarButtonsToRemove = list(
               'sendDataToCloud',
               'autoScale2d',
               'zoomIn2d',
               'zoomOut2d',
               'toggleSpikelines',
               'hoverClosestCartesian',
               'hoverCompareCartesian'
             ),
             toImageButtonOptions = list(format = input$heatmap_format,
                                         width = input$heatmap_width,
                                         height = input$heatmap_height,
                                         filename = "heatmap")
      )
  })
  
  
  ### FOLD CHANGE PLOT
  output$foldChangeSelection <- renderUI({
    selectizeInput(
      "foldChangeSelection",
      "Select 2 comparisons for a fold change plot",
      c("",reacValues$reacConditions),
      selected = NULL,
      options = list(maxItems = 2)
    )
  })
  
  
  get_fc_data <- reactive({
    event_data("plotly_selected", source = "subset")
  })
  
  
  output$foldChangePlot <- renderPlotly({
    input$sumbitFoldChangePlot
    req(reacValues$dataComp)
    #req(input$foldChangeSelection)
    matrixSet <- isolate(enrichedMatrixSet() )
    
    
    selection <- isolate(input$foldChangeSelection)
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    
    fontsize <- isolate(input$fc_base)
    plot_width <- isolate(input$fc_width)
    plot_height <- isolate(input$fc_height)
    pointsize <- isolate(input$fc_pointsize)
    format <- isolate(input$fc_format)
    showLine <- isolate(input$fc_showLine)
    
    c1 <- isolate(input$fcplotcol_1)
    c2 <- isolate(input$fcplotcol_2)
    c3 <- isolate(input$fcplotcol_3)
    c4 <- isolate(input$fcplotcol_4)
    c5 <- isolate(input$fcplotcol_5)
    
    colors <- c(c1, c2, c3, c4, c5)
    
    if(is.null(c1)) colors <- myScatterColors()
    
    validate(need(length(selection)>1, "Please select two comparisons."))
    validate(need(length(ridx)>0, "No proteins to plot (output table must be empty)."))
    validate(need(nrow(reacValues$dataComp)>0, "No data to plot."))
    rnames <- rownames(reacValues$dataComp[ridx,])
    
    plotData <-
      getFCPlotData(rnames,
                    reacValues$dataLimma,
                    matrixSet,
                    selection)
    
    plotData$show_id <- FALSE
    if (!is.null(get_fc_data() )) {
      plotData[plotData$key %in% get_fc_data()$key, "show_id"] <- TRUE
    }
    
    formula <- as.formula(paste0(paste0(logfcPrefix, selection[2]), ' ~ ',
                                   paste0(logfcPrefix, selection[1])   ) )
    
    fit1 <- lm(formula, data=plotData)
    fit1.intercept <- fit1$coefficients[[1]]
    fit1.slope <- fit1$coefficients[[2]]
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(fit1)[1]), digits = 2),
                          b = format(unname(coef(fit1)[2]), digits = 2),
                          r2 = format(summary(fit1)$r.squared, digits = 3)))
    txt <- as.character(as.expression(eq))
    
    
    title <-
      paste0(signif(fit1$coef[[1]], 5), 
             "x + ", 
             signif(fit1$coef[[2]], 5),
             ", r2 = ",
             signif(summary(fit1)$r.squared, 5))
    
    labelNamesX <- paste(unlist(strsplit(selection[1], "__vs__") ), collapse = "/")
    labelNamesX <- paste0("log2FC(", labelNamesX, ")")
    labelNamesY <- paste(unlist(strsplit(selection[2], "__vs__") ), collapse = "/")
    labelNamesY <- paste0("log2FC(", labelNamesY, ")")
    
    
    plotData$significant <-
      factor(plotData$significant,
             levels = c('both', selection[1], selection[2], 'none'))
    
    withProgress(message = "Plotting Fold Change plot", {
      pu <-
        ggplot(plotData,
               # variable != literal in the R "programming" ""language""
               aes(
                 x = !!sym(paste0(logfcPrefix, selection[1])),
                 y = !!sym(paste0(logfcPrefix, selection[2])),
                 label = Gene.names,
                 key = key,
                 color = significant
               )) + geom_point(size = pointsize) + 
        theme_minimal(base_size = fontsize) + scale_color_manual(values=colors) +
        labs(x = labelNamesX, y = labelNamesY) #, title =  title) 
      
      intercept <- 0
      slope <- 1
      if (showLine == "linear regression") {
        intercept <- fit1.intercept
        slope <- fit1.slope
        pu <- pu + ggtitle(title)
      } else if (showLine == "straight line") {
        intercept <- 0
        slope <- 1
      }
      
      if (showLine != "none") {
        pu <- pu + geom_abline(
          intercept = intercept,
          slope = slope,
          size = 1,
          alpha = 0.5,
          color = colors[5]
        )
      }
      
      pu <- pu + geom_text(
        data = subset(plotData, show_id),
        aes(!!sym(paste0(logfcPrefix, selection[1])),
            !!sym(paste0(logfcPrefix, selection[2])),
            label = Gene.names)
        #hjust=0, vjust=0
        ,position = position_jitter(width=0.25,height=0.25)
      )
    })
    
    p <- ggplotly(pu, source = "subset") %>% layout(dragmode = "select")
    p %>% config(displaylogo = F,
             modeBarButtonsToRemove = list(
               'sendDataToCloud',
               'autoScale2d',
               'zoomIn2d',
               'zoomOut2d',
               'toggleSpikelines',
               'hoverClosestCartesian',
               'hoverCompareCartesian'
             ),
             toImageButtonOptions = list(format = input$fc_format,
                                         width = input$fc_width,
                                         height = input$fc_height,
                                         filename = "foldchange_plot")
    )
  })
  
  output$maVolcanoSample <- renderUI({
    req(reacValues$reacConditions)
    selectInput(
      "maVolcanoSampleSelect",
      "Show Volcano and MA plot for a comparison",
      choices = c("",reacValues$reacConditions),
      multiple = F
    )
  })
  
  
  ## VOLCANO PLOT
  
  volcanoPlotData <- reactive({
    input$maVolcanoSubmit
    #input$maVolcanoPlot
    sample <- isolate(input$maVolcanoSampleSelect)
    fcThresh <- isolate(input$fcCutoff)
    choice <- isolate(input$enrichmentChoice)
    sigCutoffValue <- isolate(input$sigCutoffValue)
    padjY <- isolate(input$volcano_padj_y)
    
    reacValues$selection <- sample
    reacValues$fcCutoff <- fcThresh
    reacValues$enrichmentChoice <- choice
    reacValues$sigCutoffValue <- sigCutoffValue
    
    validate(need(!is.null(sample) & sample!="" & length(sample)>0, "Please select a group comparison."))
    
    padjYBoolean <- ifelse(padjY == "p-values", FALSE, TRUE)
    
    getVolcanoPlotData(
      reacValues$dataLimma,
      sample,
      fcThresh,
      sigCutoffValue,
      choice,
      padjYBoolean
    )
  })
  
  
  get_data <- reactive({
    event_data("plotly_selected", source = "subset")
  })
  
  
  output$volcanoPlot <- renderPlotly({
    input$maVolcanoSubmit
    req(reacValues$dataComp)
    sample <- isolate(input$maVolcanoSampleSelect)
    validate(need(length(grep(sample, names(reacValues$dataComp))) > 0,""))
    
    fontsize <- isolate(input$volcano_base)
    legend_fontsize <- isolate(input$volcano_legend)
    plot_width <- isolate(input$volcano_width)
    plot_height <- isolate(input$volcano_height)
    padjY <- isolate(input$volcano_padj_y)
    pointsize <- isolate(input$volcano_pointsize)
    format <- isolate(input$volcano_format)
    padjYBoolean <- ifelse(padjY == "p-values", FALSE, TRUE)
    
    c1 <- isolate(input$volcanocol_1)
    c2 <- isolate(input$volcanocol_2)
    
    pal <- brewer.pal(3, "Set2")
    if (is.null(c1) | is.null(c2)) {
      c1 <- pal[1]
      c2 <- pal[2]
    }
    pcols <- c(c1,c2,"red")
    setChoice <- "union"#isolate(input$setChoice)
    validate(need(!is.null(sample) & sample!="" & length(sample)>0, "Please select a group comparison."))
    
    pltData <- volcanoPlotData()
    if (!is.null(get_data())) {
      pltData[pltData$key %in% get_data()$key, "show_id"] <- TRUE
    }

    if (reacValues$show_analysis == FALSE) {
      reacValues$show_analysis = TRUE
      toggle(id = 'hide_before_comparisons', anim = T)
    }
    
    validate(need(
      "P.Value" %in% colnames(pltData ),
      paste0(
        "No P-Value available, seems to be a pilot?\nCannot output volcano plot."
      )
    ))
    
    xText <- unlist(strsplit(sample, "__vs__"))
    xText <- paste0("log2FC(", xText[1], "/", xText[2], ")")
    
    withProgress(message = 'Plotting Volcano plot ', {
      p <-
        plotVolcanoPlot(pltData,  pcols, padjYBoolean, pointsize)
      p <-
        p + theme_minimal(base_size = fontsize) + theme(
          legend.text = element_text(size = legend_fontsize),
          legend.title = element_text(size = legend_fontsize)
        )
      p <- p + geom_text(
        data = subset(pltData, show_id),
        aes(logFC, nlog10_pval, label = Gene.names),
        position=position_jitter(width=0.25,height=0.25)
        #,position = position_jitter(seed = 1)#position_nudge(y = -0.1)
      ) + xlab(xText)
      incProgress(amount = 0.75)
      
      p <- ggplotly(p, source = "subset") %>% layout(dragmode = "select")
      incProgress(1/1.01, detail = paste("Finished!"))
      
    })
    print(p %>%
            config(displaylogo = F,
                   modeBarButtonsToRemove = list(
                     'sendDataToCloud',
                     'autoScale2d',
                     'zoomIn2d',
                     'zoomOut2d',
                     'toggleSpikelines',
                     'hoverClosestCartesian',
                     'hoverCompareCartesian'
                   ),
                   toImageButtonOptions = list(format = format,
                                               width = plot_width,
                                               height = plot_height,
                                               filename = "volcano_plot")
            )
    )
  })
  
  
  observeEvent(input$maVolcanoSubmit,{
    if (reacValues$show_single == FALSE) {
      reacValues$show_single = TRUE
      toggle(id = 'hide_before_single_submit', anim = T)
    }
    sigCutoffValue <- isolate(input$sigCutoffValue)
    reacValues$sigCutoffValue <- sigCutoffValue
    sample <- isolate(input$maVolcanoSampleSelect)
    
    reacValues$dataComp <-
      getComparisonsData2(reacValues$dataLimma,
                          enrichedMatrixSet(),
                          'union',
                          input$maVolcanoSampleSelect)
  })
  
  
  ## MA PLOT
  
  output$maplot <- renderPlotly({
    input$maVolcanoSubmit
    req(reacValues$dataComp)
    sample <- isolate(input$maVolcanoSampleSelect)
    fontsize <- isolate(input$volcano_base)
    legend_fontsize <- isolate(input$volcano_legend)
    plot_width <- isolate(input$volcano_width)
    plot_height <- isolate(input$volcano_height)
    pointsize <- isolate(input$volcano_pointsize)
    format <- isolate(input$volcano_format)
    
    pltData <- isolate(volcanoPlotData() )
    c1 <- isolate(input$volcanocol_1)
    c2 <- isolate(input$volcanocol_2)

    pal <- brewer.pal(3, "Set2")
    if (is.null(c1) | is.null(c2)) {
      c1 <- pal[1]
      c2 <- pal[2]
    }
    pcols <- c(c1,c2,"red")
    
    validate(need(!is.null(sample) & sample!="" & length(sample)>0, "Please select a group comparison."))
    
    validate(need(
      "AveExpr" %in% colnames(pltData),
      paste0(
        "No AveExpr_ available.\nCannot output MA plot."
      )
    ))
    
    
    if (!is.null(get_data())) {
      pltData[pltData$key %in% get_data()$key, "show_id"] <- TRUE
    }
    
    withProgress(message = "Plotting MA plot ", {
      p <- plotMAPlot(pltData, pcols, pointsize)
      
      p <-
        p + theme_minimal(base_size = fontsize) + theme(
          legend.text = element_text(size = legend_fontsize),
          legend.title = element_text(size = legend_fontsize)
        )
      
      xText <- unlist(strsplit(sample, "__vs__"))
      xText <- paste0("log2FC(", xText[1], "/", xText[2], ")")
      
      p <- p + geom_text(
        data = subset(pltData, show_id),
        aes(logFC, AveExpr, label = Gene.names),
        position = position_jitter(width = 0.25, height = 0.25)
      ) + xlab(xText)
    })
    
    p <- ggplotly(p, source = "subset") %>% layout(dragmode = "select")
    print(p %>%
            config(displaylogo = F,
                   modeBarButtonsToRemove = list(
                     'sendDataToCloud',
                     'autoScale2d',
                     'zoomIn2d',
                     'zoomOut2d',
                     'toggleSpikelines',
                     'hoverClosestCartesian',
                     'hoverCompareCartesian'
                   ),
                   toImageButtonOptions = list(format = format,
                                               width = plot_width,
                                               height = plot_height,
                                               filename = "ma_plot")
            )
    )
  })
  
  ## PROFILE PLOT
  
  observe({
    updateSelectizeInput(session,
                         server = T,
                         'selectProfilePlotGene',
                         #'selectProfilePlotGene',
                         choices = c("",reacValues$filtData[[geneName]]),
    )
  })
  
  output$profilePlot <- renderPlotly({
    req(reacValues$dataLimma)
    req(input$selectProfilePlotGene)
    
    validate(need(reacValues$compareAmicaSelected==FALSE,
                  paste0("Cannot output profile plot for 'multiple_amica_upload'\n
                         Change the data set back to the original data to 
                         plot a profile plot") ))
    
    rowIdx <- which(rowData(reacValues$proteinData)[[geneName]]==input$selectProfilePlotGene)

    object <- toLongFormat(
      assay(reacValues$proteinData, "ImputedIntensity")[rowIdx, ],
      reacValues$proteinData,
      addGroup = TRUE,
      addContaminant = TRUE,
      addGeneName = TRUE
    )
    
    stats <- Rmisc::summarySE(object, measurevar="value", groupvars=c("rowname","group"))
    stats <- stats[!is.na(stats$value),]
    names(stats)[which(names(stats)=="rowname")] <- "Protein.IDs"
    
    p <- ggplot(stats, aes(x = group, y = value, color = Protein.IDs)) +
      geom_point(size = input$profile_pointsize) + 
      geom_errorbar(aes(ymin = value - se, ymax = value + se), width = .2) +
      xlab("") + ylab("Intensity (log2)") + theme_minimal(base_size = input$profile_base) + ggtitle(input$selectProfilePlotGene) + 
      scale_color_manual(values=myScatterColors() ) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = input$profile_base
      ),legend.title = element_blank())
    
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = list(
        'sendDataToCloud',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'toggleSpikelines',
        'hoverClosestCartesian',
        'hoverCompareCartesian'
      ),
      toImageButtonOptions = list(
        format = input$profile_format,
        width = input$profile_width,
        height = input$profile_height,
        filename = "profile_plot"
      )
    )
  })
  
  
  #### print stat info of protein
  output$geneSummary <- renderDT({
    req(input$selectProfilePlotGene)
    req(reacValues$dataLimma)
    
    idx <- which(reacValues$filtData$Gene.names==input$selectProfilePlotGene)
    hasPadj <- ifelse(length(grep(padjPrefix, names(reacValues$dataLimma)))>0, TRUE, FALSE )
    
    subSelection <-
      reacValues$dataLimma[idx,
                           grep(paste0(c(
                             logfcPrefix,
                             padjPrefix
                           ), collapse = "|"), names(reacValues$dataLimma) )]
    
    subSelection[[proteinId]] <- reacValues$filtData[idx, proteinId]
    
    dt <- NULL
    if (hasPadj) {
      dt <- data.table::melt(setDT(subSelection),
                             measure=patterns(logfcPrefix, padjPrefix),
                             value.name=c("logFC", "adj.P.Val"),
                             variable.factor = T,
                             variable.name="Comparison"
      )
    } else {
      dt <- data.table::melt(setDT(subSelection),
                             measure=patterns(logfcPrefix),
                             value.name=c("logFC"),
                             variable.factor = T,
                             variable.name="Comparison"
      )
    }
    
    x <- as.factor(grep(logfcPrefix, names(subSelection), value = T))
    dt$Comparison <- x[dt$Comparison]
    dt$Comparison <- gsub(logfcPrefix, "",as.character(dt$Comparison))
    
    dt <- datatable(
      dt,
      rownames=F
    ) %>% formatRound("logFC",3) 
    if (hasPadj) dt <- dt %>% formatRound("adj.P.Val",6)
    dt
  })
  
  
  observeEvent(input$submitORA,{
    req(reacValues$dataComp)
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    organism <- isolate(input$species)
    validate(need(nrow(reacValues$dataComp)>0, "No proteins to plot (output table must be empty)."))
    validate(need(length(ridx)>0, "No proteins to plot (output table must be empty)."))
    if ( nrow(reacValues$dataComp) < 2 ) return(NULL)
    
    queryGenes <- gsub(";.*", "", reacValues$dataComp[ridx, geneName])
    #sources <- c("GO:MF","GO:CC","GO:BP","KEGG","REAC","CORUM","WP", "HPA", "TF")
    sources <- c("GO:MF","GO:CC","GO:BP","KEGG","REAC","CORUM","WP")

    withProgress(message = 'Computing ORA ', {
      reacValues$dataGprofiler <- gost(
        query = queryGenes,
        sources = sources,
        evcodes = input$showGenes,
        organism = organism,
        exclude_iea = T,
        significant = input$significantORA
      )
      incProgress(1/1.01, detail = paste("Finished!"))
    })
    
    reacValues$GostPlot <- reacValues$dataGprofiler
    reacValues$dataGprofiler <- reacValues$dataGprofiler$result
    
    cols <- c(
      "source",
      "term_id",
      "term_name",
      "p_value",
      "term_size",
      "intersection_size",
      "intersection"
    )
    if (input$showGenes == F) {
      cols <- cols[1:length(cols) - 1]
    }
    reacValues$dataGprofiler <- reacValues$dataGprofiler[, cols]
    
    if (reacValues$show_ora == FALSE) {
      reacValues$show_ora = TRUE
      toggle(id = 'hide_ora_before_submit', anim = T)
    }
  })
  
  output$gostplot <- renderPlotly({
    req(reacValues$GostPlot )
    validate(
      need(nrow(reacValues$GostPlot$result) >= 1, "No significant over-represented term detected.")
    )
    gostplot(reacValues$GostPlot, capped = F, interactive = T)
    # %>% layout(yaxis = list(title = '-log10(p-adj)'))
  })
  
  oraBarBase <- eventReactive(input$submitORABar,{
    req(reacValues$dataGprofiler)
    orasource <- isolate(input$orasource)
    oracolor <- isolate(input$oraBar_color)
    oraMaxTerm <- isolate(input$oraBar_maxTerms)
    
    oraMaxTerm <- ifelse(oraMaxTerm == 0, 100000, oraMaxTerm)
    
    plotDf <- reacValues$dataGprofiler[reacValues$dataGprofiler$source==orasource,]
    
    validate(need(nrow(plotDf)>0, "No significant over-represented terms detected."))
    
    plotDf$minusLog10p_value <- -log10(plotDf$p_value)
    plotDf <- plotDf[!duplicated(plotDf$term_name),]
    plotDf <- plotDf[order(plotDf$minusLog10p_value),]
    plotDf <- head(plotDf, min(nrow(plotDf), oraMaxTerm))
    plotDf$term_name <- factor(plotDf$term_name, levels = plotDf$term_name)
    
    p <- ggplot(plotDf, aes(x = term_name, y = minusLog10p_value)) +
      geom_bar(stat = "identity", fill=oracolor, width = 0.9)  + xlab("") +
      ylab("-log10(p-value)") +
      theme_minimal( base_size = 14) + coord_flip()
    
    ggplotly(p)
  })
  
  output$oraBarplot <- renderPlotly({
    

    
    oraBarBase() %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = list(
        'sendDataToCloud',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'toggleSpikelines',
        'hoverClosestCartesian',
        'hoverCompareCartesian'
      ),
      toImageButtonOptions = list(
        format = input$oraBar_format,
        width = input$oraBar_width,
        height = input$oraBar_height,
        filename = "ora_barplot"
      )
    )

  })
  
  output$download1 <- downloadHandler(
    filename = function() {
      paste("data_gprofiler2_", Sys.Date(), ".tsv", sep="")
    },
    content = function(file) {
      write.table(reacValues$dataGprofiler, file, row.names = F, sep = "\t", quote = F)
    }
  )
  
  output$download_button_heatmap <- renderUI({
    req(reacValues$dataHeatmap)
    downloadButton("download_heatmap", "Download heatmap")
  })
  
  output$download_button_ora <- renderUI({
    req(reacValues$dataGprofiler)
    downloadButton("download1", "Download results")
  })
  
  
  
  
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("heatmap.pdf")
    },
    content = function(file) {
      
      width <- ifelse(!is.null(input$heatmap_width),
                      input$heatmap_width,
                      7)
      height <- ifelse(!is.null(input$heatmap_height),
                       input$heatmap_height,
                       7)
      
      pdf(file, width = width, height = height, onefile=F)
      
      fontsize = 12
      if (nrow(reacValues$dataHeatmap) > 100) {
        fontsize = 4
      } else if (nrow(reacValues$dataHeatmap) > 74) {
        fontsize = 6
      } else if (nrow(reacValues$dataHeatmap) > 50) {
        fontsize = 8
      } else if (nrow(reacValues$dataHeatmap) > 30) {
        fontsize = 10
      }
      
      clRows <- ifelse(!is.null(reacValues$clusterRows),
                       reacValues$clusterRows,
                       TRUE)
      clCols <- ifelse(!is.null(reacValues$clusterCols),
                       reacValues$clusterCols,
                       TRUE)
      scaleRows <- ifelse(!is.null(reacValues$scaleHeatmap),
                          reacValues$scaleHeatmap,
                          "row")
      print(
        pheatmap(
          reacValues$dataHeatmap,
          cluster_rows = clRows,
          cluster_cols = clCols,
          scale = scaleRows,
          color = heatColors(),
          border_color = NA,
          fontsize_row = fontsize
        )
      )
      dev.off()
    }
  )
  
  output$download_amica <- renderUI({
    req(reacValues$analysisSuccess)
    if (reacValues$amicaInput) return(NULL)
    downloadButton("download2", h4("Download amica output"))
  })
  
  output$download2 <- downloadHandler(
    filename = function() {
      "amica_proteinGroups.tsv"
    },
    content = function(file) {
      outDf <- amicaOutput(reacValues$proteinData, reacValues$dataLimmaOriginal )
      write.table(outDf, file, row.names = F, quote = F, sep = "\t")
    }
  )
  
  output$gprofilerDT <- renderDT({
    req(reacValues$dataGprofiler )
    datatable(
      reacValues$dataGprofiler,
      extensions = 'Buttons',
      filter = "top",
      rownames = F,
      options = list(
        pageLength = 10,
        autoWidth = TRUE
        #dom = 'Bfrtip',
        #buttons = c('csv')
      )
    ) %>% formatSignif(columns = c('p_value'), digits = 3)
  })
  
  # ------------------------------------------- visNetwork
  
  # duplicated code :(
  
  output$networkInput <- renderUI({
    req(reacValues$reacConditions)
    selectInput(
      "networkInput",
      "Choose a network to display.",
      reacValues$reacConditions,
      multiple = F
    )
  })
  
  output$heatmapSamplesInput <- renderUI({
    req(reacValues$dataLimma )
    selectizeInput(
      "heatmapSamplesInput",
      "Select at least 2 groups for the heatmap.",
      reacValues$expDesign$groups,
      multiple = T
    )
  })
  
  values <- reactiveValues(networkData = NULL)
  
  ppi <- reactive({
    simplify(read_graph('data/intact_weighted.edgelist', format="ncol", directed=F))
  })
  
  cellmap <- reactive({
    cellmap <- read.table("data/preys-latest.txt", header = T, sep = "\t", stringsAsFactors = F)
    cellmap <- cellmap[, c("symbol", "MMF.localization", "SAFE.localization")]
    colnames(cellmap) <- c("label", "CellMap MMF localization", "CellMap SAFE localization")
    cellmap
  })
  
  #### PPI NETWORK
  multiNetwork <- function(thresh=0) {
    if (is.na(thresh) || is.null(thresh)) thresh <- 0
    validate(need(nrow(reacValues$dataComp) > 1,"Not enough data to display."))
    
    ridx <- input$groupComparisonsDT_rows_all
    rnames <- rownames(reacValues$dataComp[ridx,])
    
    
    nw.data <-
      getBait2PreyNetwork(reacValues$dataComp[rnames,],
                          enrichedMatrixSet(),
                          reacValues$selection,
                          myScatterColors())
    nodes <- nw.data[[1]]
    edges <- nw.data[[2]]

    idxs <- match(nodes$label, V(ppi() )$name)
    idxs <- idxs[!is.na(idxs)]
    
    if (length(idxs)>1) {
      G <- induced_subgraph(ppi(), idxs)
      G.df <- igraph::as_long_data_frame(G)
      
      names(G.df) <- c("from", "to", "value", "from_label", "to_label")
      
      from_id <- match(G.df$from_label, nodes$label)
      to_id <- match(G.df$to_label, nodes$label)
      
      ppi.edges <-
        data.frame(
          from = from_id,
          to = to_id,
          color = rep(myScatterColors()[2], length(to_id)),
          interaction = rep('PPI', length(to_id)),
          value = G.df$value
        )
      edges <- rbind(edges, ppi.edges[ppi.edges$value >= thresh,])
    }
    
    nodes <-
      merge(nodes, cellmap(), by="label", all.x = T)
    
    reacValues$nwNodes <- nodes
    reacValues$nwEdges <- edges
    
    visNetwork(nodes,
               edges,
               layout = "layout_nicely",
               type = "full") %>%
      visIgraphLayout() %>% visPhysics(stabilization = F) %>%
      visOptions(selectedBy = "CellMap SAFE localization",
                 highlightNearest = TRUE,
                 nodesIdSelection = TRUE) %>%
      visEdges(smooth=F) %>%
      visLegend(
        main = "Legend",
        useGroups = F,
        addNodes = data.frame(
          label = c("Bait", "Prey"),
          shape = c("diamond", "circle"),
          color = c(myScatterColors()[1], myScatterColors()[2])
        ),
        addEdges = data.frame(
          label = c("Bait-Prey", "Protein-Protein"),
          color = c(myScatterColors()[1], myScatterColors()[2])
        ),
        width = 0.3
      ) %>%
      visExport(
        type = "pdf",
        name = "Network",
        float = "left",
        label = "Save network (png)",
        background = "white",
        style = ""
      )
  }

  output$downloadSpecificityNetwork <- downloadHandler(
    filename = function() {
      "amica_specificity_network.gml"
    },
    content = function(file) {
      validate(need(nrow(reacValues$dataComp) > 1,"Not enough data to display."))
      
      outs <- dfs2gml(reacValues$nwNodes, reacValues$nwEdges)
      writeLines(outs, file)
    }
  )


  output$download_button_spec_network <- renderUI({
    req(reacValues$nwEdges)
    downloadButton("downloadSpecificityNetwork", "Download gml")
  })
  # 
  
  singleNetwork <- function(thresh) {
    if (is.na(thresh) || is.null(thresh)) thresh <- 0
    
    req(reacValues$dataComp)
    
    ridx <- input$groupComparisonsDT_rows_all
    rnames <- rownames(reacValues$dataComp[ridx,])
    
    networkData <- toNetworkData(reacValues$dataComp[rnames,], ppi(), cellmap())
    networkData$edges <- networkData$edges[networkData$edges$value >= thresh, ]
    
    msg <- ""
    if (length(rnames) > 1) {
      msg <- "\nCould not map selected proteins to the human PPI network."
    }
    
    validate(need(!is.null(networkData), paste0("#### There are no proteins to display with your selection.", msg)))
    
    colLegend <-
      colour_values(networkData$nodes$log2FC,
                    n_summaries = 6,
                    palette = 'viridis')
    networkData$nodes$color <- colLegend$colours
    
    reacValues$nwNodes <- networkData$nodes
    reacValues$nwEdges <- networkData$edges
    
    lnodes <-
      data.frame(
        id = 1:length(colLegend$summary_values),
        color = colLegend$summary_colours,
        label = colLegend$summary_values,
        shape = "circle"
      )

    visNetwork(networkData$nodes,
               networkData$edges,
               layout = "layout_nicely",
               type = "full") %>%
      visIgraphLayout() %>% visPhysics(stabilization = F) %>%
      visOptions(selectedBy = "CellMap SAFE localization",
                 highlightNearest = TRUE,
                 nodesIdSelection = TRUE) %>%
      visEdges(smooth=F, color=list(color = "grey", highlight = "red")) %>%

      visLegend(
        useGroups = F,
        addNodes = lnodes,
        main = "log2FC",
        width = 0.3,
        ncol = 6,
        stepX = 50
      ) %>%
      visExport(
        type = "png",
        name = "Network",
        float = "left",
        label = "Save network (png)",
        background = "white",
        style = ""
      )
  }
  
  output$network <- renderVisNetwork({
    out <- 0
    if (length(grep(logfcPrefix, names(reacValues$dataComp)) ) > 1) {
      out <- multiNetwork(input$edgeWeightThresh)
    } else {
      out <- singleNetwork(input$edgeWeightThresh)
    }
    out %>% visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}"
    )
  })
  
  #### COMPARE amica EXPERIMENTS
  observeEvent(input$submitAmicaComparisons, {
    
    suffix1 <- isolate(input$suffix1)
    suffix2 <- isolate(input$suffix2)
    delim <- isolate(input$subPattern)
    mergeKey <- isolate(input$mergeKey)
    
    suffix1 <- paste0('_', suffix1)
    suffix2 <- paste0('_', suffix2)
    
    sourcePath <- input$amicaFile2$datapath
    
    amicaData2 <- read.delim(sourcePath, header=T, sep="\t", stringsAsFactors = F)
    if ("Intensity" %in% names(amicaData2)) {
      amicaData2$Intensity <- NULL
    }
    if ("iBAQ" %in% names(amicaData2)) {
      amicaData2$iBAQ <- NULL
    }
    if ("iBAQ.peptides" %in% names(amicaData2)) {
      amicaData2$iBAQ.peptides <- NULL
    }
    
    amicaData2 <- amicaData2[amicaData2$quantified=="+",]
    amicaData2 <- amicaData2[, grep("Majority.protein.IDs|Gene.names|Intensity|iBAQ|__vs__", colnames(amicaData2))]

    names <- assayNames(reacValues$proteinData)
    tmpInts <- assay(reacValues$proteinData, names[1])[isQuantRnames(reacValues$proteinData),]
    names <- names[2:length(names)]
    for (idx in seq_along(names)) {
      name <- names[idx]
      df <- assay(reacValues$proteinData, name)[isQuantRnames(reacValues$proteinData),]
      colnames(df) <- paste0(name, "_", colnames(df))
      tmpInts <- cbind(tmpInts, df)
    }
    
    premerged <-
      cbind(reacValues$dataLimmaOriginal, tmpInts)
    
    reacValues$combinedData <- mergeAmicas(premerged,
                                           amicaData2,
                                           mergeKey,
                                           delim,
                                           suffix1,
                                           suffix2)
    
    if (reacValues$compareAmicasToggled == FALSE) {
      toggle(id = 'compareAmicaInput', anim = T)
      reacValues$compareAmicasToggled <- TRUE
    }
  })
  
  output$multiAmicasInput <- reactive({
    req(reacValues$combinedData)
    TRUE
  })
  outputOptions(output, "multiAmicasInput", suspendWhenHidden = FALSE)

  ####
  
  output$download_merged_amica <- renderUI({
    downloadButton("downloadMerged", h4("Download merged amica files"))
  })
  
  output$downloadMerged <- downloadHandler(
    filename = function() {
      "merged_amica_files.tsv"
    },
    content = function(file) {
      write.table(reacValues$combinedData, file, row.names = F, quote = F, sep = "\t")
    }
  )
  
  ####

  output$summaryMergedAmica <- renderText({
    req(reacValues$combinedData)
    paste0(
      "Successfully uploaded second data set!\n",
      "Number of proteins in merged output:\n",
      nrow(req(reacValues$combinedData)),
      "\nNumber of group comparisons: \n",
      length(grep("logFC", colnames(req(reacValues$combinedData))))
    )
  })
  
  
  observeEvent(input$submitDatasetSelection, {
    if (input$selectedDataset == "original_data") {
      if (reacValues$compareAmicaSelected == TRUE) {
        #reacValues$dataLimmaOriginal <- reacValues$dataLimma
        reacValues$dataLimma <- reacValues$dataLimmaOriginal
        
        comps <-
          grep(logfcPrefix, colnames(reacValues$dataLimma), value = T)
        reacValues$reacConditions <- gsub(logfcPrefix, "", comps)
        
        reacValues$compareAmicaSelected <- FALSE
      }
    } else {
      if (reacValues$compareAmicaSelected == FALSE) {
        reacValues$dataLimma <- reacValues$combinedData

        comps <-
          grep(logfcPrefix, colnames(reacValues$dataLimma), value = T)
        reacValues$reacConditions <- gsub(logfcPrefix, "", comps)
        reacValues$compareAmicaSelected <- TRUE
      }
    }
  })
  
  output$summaryLoadedDataSet <- renderText({
    req(reacValues$combinedData)
    out <- ""
    if (reacValues$compareAmicaSelected) {
      out <- paste0("Current data set: compare multiple amica files")
    } else {
      out <- paste0("Current data set: original upload")
    }
    out
  })

  
  output$assayNamesAmicas <- renderUI({
    req(reacValues$combinedData)
    
    tmp <- grep("Intensity", names(reacValues$combinedData), value = T)
    prefixes <- unique(gsub("Intensity.*", "Intensity", tmp))
    if (length(grep("iBAQ", names(reacValues$combinedData))) > 0 ||
        length(grep("iBAQ", assayNames(reacValues$proteinData))) > 0) {
      prefixes <- c(prefixes, "iBAQ")
    }
    selectInput("assayNamesAmicas",
                "Which intensities do you want to plot?",
                prefixes,
                selected = "ImputedIntensity")
  })
  
  ## AMICA COMPARISON SCATTER PLOTS
  output$compareScatterPlotsAmica <- renderUI({
    req(reacValues$combinedData)
    req(input$assayNamesAmicas)
    assayNameAmicas <- isolate(input$assayNamesAmicas)
    
    samples <- grep(paste0("^",assayNameAmicas),
                    names(reacValues$combinedData),
                    value = T
    )

    samples <- gsub(paste0(assayNameAmicas, "."), "", samples)

    selectizeInput(
      "selectScatterSamplesAmica",
      "Select 2 samples for a scatter plot",
      c("",samples),
      selected = NULL,
      options = list(maxItems = 2)
    )
  })

  
  scatterPlotsAmicaBase <- eventReactive(input$submitScatterAmicas,{
    req(input$selectScatterSamplesAmica)
    req(input$assayNamesAmicas)
    assayNameAmicas <- isolate(input$assayNamesAmicas)
    fontsize <- isolate(input$scatteramica_base)
    
    validate(need(length(input$selectScatterSamplesAmica)==2, ""))
    colors <- isolate(myScatterColors())

    clrs <- colors[1]
    xvar <- input$selectScatterSamplesAmica[1]
    yvar <- input$selectScatterSamplesAmica[2]

    sampleX <- grep(xvar, grep(assayNameAmicas, names(reacValues$combinedData), value = T ), value = T)
    sampleY <- grep(yvar, grep(assayNameAmicas, names(reacValues$combinedData), value = T ), value = T)
    
    formula <- as.formula(paste0(sampleY, ' ~ ', sampleX))
    fit1 <- lm(formula, data=reacValues$combinedData)
    fit1.intercept <- fit1$coefficients[[1]]
    fit1.slope <- fit1$coefficients[[2]]

    title <-
      paste0(signif(fit1$coef[[1]], 5), 
             "x + ", 
             signif(fit1$coef[[2]], 5),
             ", r2 = ",
             signif(summary(fit1)$r.squared, 5)
             )

    pu <-
      ggplot(reacValues$combinedData,
             # variable != literal in the R "programming" ""language""
             aes(
               x = !!sym(paste0(sampleX)),
               y = !!sym(paste0(sampleY)),
               label = Gene.names
             )) + scale_color_brewer(palette="Paired") + theme_minimal(base_size = fontsize)

    pu <- pu + geom_point(color=clrs) + labs(x = paste0(assayNameAmicas, ' ', xvar, ' (log2)'),
                                             y = paste0(assayNameAmicas, ' ', yvar, ' (log2)'))
      
    intercept <- 0
    slope <- 1
    if (input$scatteramica_showLine == "linear regression") {
      intercept <- fit1.intercept
      slope <- fit1.slope
      pu <- pu + ggtitle(title)
    } else if (input$scatteramica_showLine == "straight line") {
      intercept <- 0
      slope <- 1
    }
    
    if (input$scatteramica_showLine != "none") {
      pu <- pu + geom_abline(
        intercept = intercept,
        slope = slope,
        size = 1,
        alpha = 0.5,
        color = colors[5]
      )
    }
    
    p <- ggplotly(pu)
    p
  })
  
  output$scatterPlotsAmica <- renderPlotly({
    
    scatterPlotsAmicaBase() %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = list(
        'sendDataToCloud',
        'autoScale2d',
        'zoomIn2d',
        'zoomOut2d',
        'toggleSpikelines',
        'hoverClosestCartesian',
        'hoverCompareCartesian'
      ),
      toImageButtonOptions = list(
        format = input$scatteramica_format,
        width = input$scatteramica_width,
        height = input$scatteramica_height,
        filename = "scatterplot_experiments"
      )
    )
  })
  
  ### CORR PLOT
  
  output$corrAmicasSamplesInput <- renderUI({
    req(input$assayNamesAmicas)
    
    assayNameAmicas <- isolate(input$assayNamesAmicas)
    samples <- grep(paste0("^",assayNameAmicas),
                    names(reacValues$combinedData),
                    value = T
    )
    samples <- gsub(paste0(assayNameAmicas, "."), "", samples)
    
    selectizeInput(
      "corrAmicasSamplesInput",
      "Press on 'Plot Correlation'.
      Leave the input blank to plot all samples.
      Only select specific groups when you have really have to.",
      samples,
      multiple = T,
      options = list(minItems = 2)
    )
  })
  
  corrBaseAmicasPlot <- eventReactive(input$submitCorAmicas,{
    assayName <- isolate(input$assayNamesAmicas)

    validate(need(assayName != "", "Please provide an intensity."))
    
    groupInputs <- isolate(input$corrAmicasSamplesInput)

    validate(need(length(groupInputs) > 1 | is.null(groupInputs), 
                  "Please select at least two groups. If none is selected all are considered."))
    
    if (is.null(groupInputs)) {
      samples <- grep(paste0("^",assayName),
                      names(reacValues$combinedData),
                      value = T
      )
      groupInputs <- gsub(paste0(assayName, "."), "", samples)
    }

    df <- reacValues$combinedData[, grep(assayName, names(reacValues$combinedData))]
    df <- df[, grep(paste0(groupInputs, collapse = "|"), names(df))]
    names(df) <- gsub(paste0(assayName, "."), "", names(df))
    corDf <- cor(df, method = "pearson", use = "complete.obs")

    heatmaply_cor(
        round(corDf, 3),
        xlab = "", 
        ylab = "",
        limits = c(min(corDf)-0.003, 1),
        plot_method = "plotly",
        key.title = "Pearson Correlation"
      )
  })
  
  output$corrAmicasPlotly <- renderPlotly({
    
    withProgress(message = "Plotting correlation plot ", {
      p <- corrBaseAmicasPlot()
    })
    p %>%  config(displaylogo = F,
                  modeBarButtonsToRemove = list(
                    'sendDataToCloud',
                    'autoScale2d',
                    'zoomIn2d',
                    'zoomOut2d',
                    'toggleSpikelines',
                    'hoverClosestCartesian',
                    'hoverCompareCartesian'
                  ),
                  toImageButtonOptions = list(format = input$corAmicas_format,
                                              width = input$corAmicas_width,
                                              height = input$corAmicas_height,
                                              filename = "corrplot")
    )
  })
  
  ### 

  
  ### headers
  output$compSummary <- renderText({
    req(reacValues$dataComp)
    paste0("There are ", nrow(reacValues$dataComp), " proteins in your selection.")
  })
  
  
  output$parameterSummary <- renderText({
    req(reacValues$dataComp)
    multi = ifelse( length(reacValues$selection)>1, TRUE, FALSE)
    
    paste0(
      "Parameters:\n\tComparison(s): ",
      paste0(reacValues$selection, collapse = ", "),
      ".\n\tFold change threshold: ",
      reacValues$fcCutoff,
      "\n\tP-Value: ", reacValues$sigCutoffValue,
      "\n\tConsidered proteins (fold change): ", reacValues$enrichmentChoice
    )
  })
  
  output$filterDTSummary <- renderText({
    req(reacValues$dataComp)
    req(input$groupComparisonsDT_rows_all)
    
    if (nrow(reacValues$dataComp) == length(input$groupComparisonsDT_rows_all) ) {
      return(NULL)
    } else {
      paste0("There are ", nrow(reacValues$dataComp), " proteins in your selection.\n",
             "After filtering the output table ", length(input$groupComparisonsDT_rows_all),
             " proteins remain for subsequent visualizations.",
             " Remove the filters in the table to visualize all proteins.")
    }
  })
  
  
  output$VolcanoHelpBox <- renderUI({
    if (input$volcanoHelp %% 2){
      helpText("You can highlight proteins by drawing a box or chosing the lasso option.
               This selection stays highlighted when you change to another volcano plot,
               by selecting a different comparison and pressing the 'Submit Analysis' button.")
    } else {
      return()
    }
  })
  
  output$assaysHelpBox <- renderUI({
    if (input$assaysHelp %% 2) {
      
      HTML("
      <p><b>LFQIntensity</b> or <b>Intensity</b> are intensities 
      that still contain missing values and potential contaminants. 
      These intensities are used to calculate the fraction of missing data 
       or the number of identified protein groups in a sample. If the data was quantified 
      with MaxQuant or FragPipe these intensities are already normalized.</p>
      
      <p><b>ImputedIntensity</b> are normalized and imputed intensities 
      that are also used to calculate differential abundance. If the re-normalization 
      option was selected in the input tab the <b>LFQIntensities</b> were normalized 
      after removing potential contaminants, reverse hits, proteins only identified by site 
      and protein groups that had too few valid values per group.</p>
           
           <p><b>iBAQ</b> (intensity-based absolute quantification) values are obtained 
           by dividing protein intensities by the number of theoretically observable tryptic peptides.
            This measure correlates well with protein abundance and is for example used to calculate 
            the percentage of contamination in a sample (if available).
           </p>
           
           <p><b>RawIntensity</b> are non-normalized, summed peptide intensities per protein group.</p>
           ")
    }
  })
  
  output$UpsetHelpBox <- renderUI({
    if (input$upsetHelp %% 2) {
      
      HTML("<p>Set comparison of differentially abundant proteins from selected comparisons 
      under selected thresholds. The dots show which sets are getting compared. 
      A dot not connected to another dot shows the number of proteins specific to that comparisons. 
      The top barplot depicts the number of intersecting proteins, and the barplot on 
      the side shows how many proteins are differentially abundant in the comparison. 
      Change the selected comparisons to your needs.
               <a href='https://jku-vds-lab.at/tools/upset/' target='_blank'>
      <img src='https://jku-vds-lab.at/assets/images/projects/upset/matrix.png' alt='UpSet plot explained'>
      </a>")
    } else {
      return()
    }
  })

  
  output$exampleHelpBox <- renderUI({
    if (!(input$exampleHelp %% 2)){
      HTML(
        
        "
        <p>Teakel, S. L., Ludescher, M., Thejer, B. M., Poschmann, G., Forwood, J. K., Neubauer, H., & Cahill, M. A. 
        (2020). Protein complexes including PGRMC1 and actin-associated proteins are disrupted by AG-205. 
        Biochemical and biophysical research communications, 524(1), 64-69.
        </p>
        <p>PGRMC1 is a protein from the MAPR family with a range of cellular functions. 
                       PGRMC1 has been described to play a role in regulating membrane trafficking and 
                       in cancer progression and response to therapies. To further understand the functions 
                       of PGRMC1 and the mechanism of the small molecule inhibitor of PGRMC1, AG-205, 
                       proteins differentially bound to PGRMC1 were identified following AG-205 
                       treatment of MIA PaCa-2 cells.</p>
                       Data have been re-analyzed from PRIDE identifier PXD0016455.
                       <a href='https://pubmed.ncbi.nlm.nih.gov/31980178'target='_blank'>Further information</a>.</p>
                       "
      )
    }
  })
  
  output$NetworkHelpBox <- renderUI({
    if (input$networkHelp %% 2){
      HTML("<p>PPI (protein-protein interaction) Network from IntAct. All interactions are derived from literature curation or direct user submissions.  
      Edge weights can be further filtered, more information \ can be found <a href='https://www.ebi.ac.uk/intact/' target='_blank'>here</a>.
      When multiple group comparisons are selected two types of edges are created: edges from the group comparison to the proteins 
      and PPI edges from IntAct between the proteins. Networks can be downloaded in GML format enabling the visualization in a network tool.</p>
      <p>Fold changes are color coded onto the nodes, sub-cellular locations can are retrieved from <a href='https://cell-map.org/' target='_blank'>Human Cell Map</a>.</p>
               ")
    } else {
      return()
    }
  })
  
  output$FoldChangePlotHelpBox <- renderUI({
    if (input$FoldChangePlotHelp %% 2){
      helpText("The fold change plot gets plotted for the proteins in your selection.
      Points are colored on their specificity to the comparison. In some cases fold changes 
      seem to be high in both comparisons, but are only significantly differentially abundant in one comparison.
               ")
    }
  })
  
  ### QC HELP
  output$boxplotHelpBox <- renderUI({
    if (input$boxplotHelp %% 2){
      helpText("Box plots show the distribution of selected intensities which
      gives an overview about their similarities.
               ")
    }
  })
  
  output$densityHelpBox <- renderUI({
    if (input$densityHelp %% 2){
      helpText("The density plot shows a smoothed version of a histogram.
      It is especially useful to compare density plots of the intensities before and after imputation.
               ")
    }
  })
  
  output$corHelpBox <- renderUI({
    if (input$corHelp %% 2){
      helpText("
      The Pearson correlation plot visualizes how well samples (e.g replicates)
      correlate.
               ")
    }
  })
  
  output$cvHelpBox <- renderUI({
    if (input$cvHelp %% 2){
      helpText("The Coefficient of Variation (CV) gets calculated by the standard deviation of replicates divided by their mean per protein,
                which gives an estimate on the reproducibility of the experiment.
               ")
    }
  })
  
  output$contaminantHelpBox <- renderUI({
    if (input$contaminantHelp %% 2){
      helpText("The percentage of contaminants is calculated with the iBAQ (Intensity based Absolute Quantification) values and shows
                the percentage of total signal from potential contaminants.
               ")
    }
  })
  
  output$abundantHelpBox <- renderUI({
    if (input$abundantHelp %% 2){
      helpText("The percentage of most abundant proteins from a sample is calculated with the iBAQ (Intensity based Absolute Quantification) values and shows
                the percentage of the 15 most abundant signals from the total signal.
               ")
    }
  })
  
  output$idHelpBox <- renderUI({
    if (input$idHelp %% 2){
      helpText("The number of identified proteins is calculated from the LFQ intensities before imputation.
               ")
    }
  })
  
  output$mvHelpBox <- renderUI({
    if (input$mvHelp %% 2){
      helpText("The percentage of missing values is calculated from the LFQ intensities before imputation.
               ")
    }
  })
  
  output$designTitle <- renderText({
    req(reacValues$proteinData)
    "Experimental Design"
    })
  
  output$analysisSuccessMsg <- renderUI({
    
    req(reacValues$analysisSuccess)
    
    HTML('<h4>Successfully processed data!</h4>\n<p>
         Scroll to the top of the page to visit the analysis tabs.
         </p>')
  })
  
  
  output$exampleFiles <- downloadHandler(
    filename="examples.zip",  # desired file name on client 
    content=function(con) {
      file.copy("data/examples.zip", con)
    }
  )
  
  
  # observeEvent(input$showMemory,{
  #   output$printMemory <- renderText({
  #     paste0(
  #       #"  output: ", capture.output(print(object_size(reacValues))), 
  #       ", mem used: ", capture.output(print(mem_used())), 
  #       "\n")
  #   })
  # })
  
  
  output$isPilot <- reactive({
    req(reacValues$proteinData)
    any(duplicated(colData(reacValues$proteinData)$groups))
  })
  outputOptions(output, "isPilot", suspendWhenHidden = FALSE)
  
  # 
  # conditionalPanel(
  #   condition = 'output.isPilot',
  #   verbatimTextOutput("color_pr")
  #   )
  
  output$exampleAmicaFile <- renderDT({
    varName <-
      c(
        "Protein ID",
        "Gene name",
        "LFQ intensity prefix",
        "Imputed Intensity prefix",
        "razor unique count",
        "razor unique prefix",
        "p-value Prefix",
        "adj. p-value prefix",
        "fold change prefix",
        "Avg. expression prefix",
        "comparison infix",
        "quantified column",
        "Potential contaminant column"
      )
    colName <- c("Majority.protein.IDs",
                 "Gene.names",
                 "LFQIntensity_",
                 "ImputedIntensity_",
                 "razorUniqueCount",
                 "razorUniqueCount",
                 "P.Value_",
                 "adj.P.Val_",
                 "logFC_",
                 "AveExpr_",
                 "__vs__",
                 "quantified",
                 "Potential.contaminant"
                 )
    comment <- c("",
                 "",
                 "MaxQuants (MQs) 'LFQ intensity' columns",
                 "Imputed and re-normalized intensities",
                 "MQs 'razor+unique count' column",
                 "MQs 'razor+unique count' column per sample",
                 "e.g P.Value_group1__vs__group2",
                 "e.g adj.P.Val_group1__vs__group2",
                 "e.g logFC_group1__vs__group2",
                 "e.g AveExp_group1__vs__group2",
                 "see below",
                 "see below",
                 "MQs Potential.contaminants column"
                 )
    df <- data.frame(varName, colName, comment)
    names(df) <- c("Variable name",
                   "Column name or prefix",
                   "Comment")
    
    datatable(
      df,
      filter = "none",
      rownames = F,
      extensions = c('Buttons'),
      options = list(
        searching = FALSE,
        dom = 'Bfrtip',
        pageLength = 15,
        autoWidth = TRUE,
        buttons = list(
          # 'csv',
          list(
            extend = 'csv',
            fieldSeparator = '\t',
            fieldBoundary = ''
          )
        )
      )
    )
  }, server = F)
  
  output$exampleDesign <- renderDT({
    tmpData <-
      read.table(
        "data/PXD0016455/design.txt",
        header = T,
        stringsAsFactors = F
      )
    

    datatable(
      tmpData,
      filter = "none",
      rownames = F,
      extensions = c('Buttons'),
      options = list(
        searching = FALSE,
        dom = 'Bfrtip',
        pageLength = 10,
        autoWidth = TRUE,
        buttons = list(
          # 'csv',
          list(
            extend = 'csv',
            fieldSeparator = '\t',
            fieldBoundary = ''
          )
        )
      )
    )
  }, server = F)
  
  output$exampleContrasts <- renderDT({
    tmpData <-
      read.table(
        "data/PXD0016455/contrasts.txt",
        header = T,
        stringsAsFactors = F
      )
    
    datatable(
      tmpData,
      filter = "none",
      rownames = F,
      extensions = c('Buttons'),
      options = list(
        searching = FALSE,
        dom = 'Bfrtip',
        pageLength = 10,
        autoWidth = TRUE,
        buttons = list(
          #'csv',
          list(
            extend = 'csv',
            fieldSeparator = '\t',
            fieldBoundary = ''
          )
        )
      )
    )
  }, server = F)
  
  
  output$specificationsExplanation <- renderDT({
    Variable <- c("proteinId",
                  "geneName",
                  "intensityPrefix",
                  "abundancePrefix",
                  "razorUniqueCount",
                  "razorUniquePrefix",
                  "spectraCount",
                  "contaminantCol")
    Pattern <- rep("...", length(Variable))
    Mandatory <- c("yes",
                   "yes",
                   "yes",
                   "no",
                   "no",
                   "no",
                   "no",
                   "no")
    df <- data.frame(Variable, Pattern, Mandatory)
    datatable(
      df,
      filter = "none",
      rownames = F,
      extensions = c('Buttons'),
      options = list(
        searching = FALSE,
        dom = 'Bfrtip',
        pageLength = 10,
        autoWidth = TRUE,
        buttons = list(
          #'csv',
          list(
            extend = 'csv',
            fieldSeparator = '\t',
            fieldBoundary = ''
          )
        )
      )
    )
  }, server = F)
  
  output$exampleSpecifications <- renderDT({
    tmpData <-
      read.table(
        "data/PXD0016455/specification.txt",
        header = T,
        stringsAsFactors = F
      )
    
    datatable(
      tmpData,
      filter = "none",
      rownames = F,
      extensions = c('Buttons'),
      options = list(
        searching = FALSE,
        dom = 'Bfrtip',
        pageLength = 10,
        autoWidth = TRUE,
        buttons = list(
          #'csv',
          list(
            extend = 'csv',
            fieldSeparator = '\t',
            fieldBoundary = ''
          )
        )
      )
    )
  }, server = F)
  
  
  observeEvent(input$showFileInput, {
    showModal(modalDialog(
      title = "Accepted File input",
      
      HTML("
      <h4>Inspecting previously analyzed output:</h4>
      <p>
      If you just want to try out the software select the 'Load in example' option 
      and press 'Upload'. After the upload the QC -, Diff. abundance - and Compare 
      amica data sets - tabs open in the main tab bar.
      </p>
      <p>
      If you want to (re-)inspect previously analyzed amica output you can just upload 
      the amica_protein_groups.tsv file together with the experimental design and
      after pressing the 'Upload' button you can 
      </p>
       <p>
           One property that can be changed in the 'Input' tab are the colors for 
           all upcoming plots. Press on 'Choose colors' to toggle the color selection 
           that is separated into different categories, depending on which plots 
           are outputted. The selected colors propagate to all subsequent visualizations 
           (e.g groups from the exp. design will always have the same color in a plot legend.)
           <center><img src='input_tutorial/colors.png' width='100%'></center>
           </p>
      <p>
      <h4>Analzing a data set</h4>
      Automatically recognized input files are MaxQuant's 'proteinGroups.txt' file and 
FragPipe's 'combined_proteins.txt' file, however you can also upload a generic 
tab-separated format which can easily be processed by the addition of a file 
that maps relevant search engine-specific columns to a standard format.
           </p>
           
           <p>
           After uploading all required files and pressing the 'Upload' button 
           a new field becomes accessible which allows for changing analysis parameters:
           <center><img src='input_tutorial/advanced.png' width='100%'></center>
           </p>
          
           "),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$showAmicaInput, {
    showModal(modalDialog(
      title = "amica file input",
      
      HTML("<p>
           amica’s tab-separated protein groups file has following columns:
           </p><br>"), 
      DTOutput("exampleAmicaFile"),
      
      HTML('
      <p>
        <ul style="list-style-type:square">
        <li>IntensityPrefix, ImputedIntensityPrefix and abundancePrefix columns 
        are log2 transformed, all 0s need to be converted to NANs. 
        No INF values allowed. amica searches for all Intensity prefixes in 
        the column names, if you want to provide more than the dafalt intensities.
        However, all intensity prefixes must have the same number of samples in 
        order to get processed.</li>
        <li>ImputedIntensityPrefix should only contain filtered, 
        imputed and normalized values.</li>
        <li>quantCol: All proteins passing spectraCount and 
        razorUniqueCount thresholds that have been quantified are set to "+" in 
        this column. Otherwise no value ("") is written in the
column</li>
        <li>
        comparisonInfix: The infix is important to retrieve the group ids 
        from a group comparison (e.g for downstream visualizations like heatmaps). 
        The groups before and after the "__vs__" infix should match with groups 
        defined in the uploaded experimental design.
        </li>
        <li>
        razorUniqueCount is a column, razorUniquePrefix is the prefix to the 
        count per sample, but they may very well have the same value 
        (just like in MaxQuant’s proteinGroups.txt)
        </li>
        </ul>
      </p><br>
           <p>
           Proteins inferred from reverse hits and peptides ”only identified by 
           site modifications” are not to be written into amica’s output. 
           Additional columns can be added in the future but are at the
moment not considered when uploaded.
           </p>'),
      size = "l",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$showDesign, {
    showModal(modalDialog(
      title = "Example experimental design",
      
      HTML("<p>
           The design file has two columns: <b>samples</b> and <b>groups</b>. The sample names in the samples column need to match 
           the column names of the input file in the order of the input file.</p>"),
      
      DTOutput("exampleDesign"),

      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$showContrasts, {
    showModal(modalDialog(
      title = "Example contrast matrix",
      
      HTML("<p>
           The contrast matrix tells amica which group comparisons to perform. The column names of
this file can be freely chosen, but column names must be provided. For each row in this file the
comparison group1-group2 is performed. If one wants to change the sign of the fold changes the
position of the groups needs to be switched in the file (e.g group2-group1 instead of group1-group2</p>"),
      
      
      DTOutput("exampleContrasts"),

      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$showSpecifications, {
    showModal(modalDialog(
      title = "Example contrast matrix",
      
      HTML("<p>
           Following variables can be parsed:
           </p>"),
      
      DTOutput("specificationsExplanation"),
      HTML("<br>The proteinId column must only contain unique entries.
      If razorUnique count is missing some functionality will be lost (DEqMS).<br>
      It is important that the provided intensities are not log2-transformed.
           <br>An example format is provided in the examples.zip file<br>"),
      DTOutput("exampleSpecifications"),
      
      HTML("<p>
           The specification file needs to be uploaded if a custom tab-delimited file is analyzed.
           The file has two columns, Variable and Pattern, these are used to change the prefixes (or post-
fixes) to identify the relevant columns in your data.</p>"),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$showTutorial, {
    showModal(modalDialog(
      title = "Example use cases",

      HTML('
      <h3>Use case 1: Single group comparison</h3>
      <p>
      These small examples have been produced with the provided example data set. 
      As the example data contains AP-MS data we set the selection parameter to 
      "enriched" to retrieve differentially abundant proteins against the control.
      <center><img src="da_tutorial/global_param.png" width="100%"></center>
      Then we select a single group comparison (PRGMC1 bait vs MIA PaCa-2 cell 
      background) and press "Submit Analysis".
      <center><img src="da_tutorial/single_comp.png" width="100%"></center>
      
      When you hover over the volcano - or MA - plot you can see features to 
      manipulate the plot. When we utilize the select box or lasso tool we can 
      annotate the highest enriched proteins, as seen in the figure above.
      </p>
      
      <p>
      Further information on most plots can be acquired when you press the "info"
      icon. Plot parameters can be changed when pressing the "wrench" icon and 
      can be saved upon hovering over the plot and clicking on the "camera" icon.
      </p>
      
      <h3>Use case 2: Multiple group comparisons</h3>
      <p>
      
      When we select the "Analyze multiple comparisons" tab pill we can compare 
      how the prey proteins change with and without the small molecule inhibitor 
      AG-205, just select the two comparisons of the bait versus the controls: 
      
      <center><img src="da_tutorial/multi_comp.png" width="100%"></center>
      </p>
      
      <p>
      Scrolling down we can figure out the quantitive changes of prey proteins 
      with and without AG-205 treatment in a fold change plot:
      <center><img src="da_tutorial/fcplot.png" width="100%"></center>
      </p>
           
           '),
      size = "l",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  
  
  observeEvent(input$showQueryTutorial, {
    showModal(modalDialog(
      title = "Advanced queries",
      
      HTML('
      <h3>Use case: Visualize proteins from over-represented functional term</h3>
      <p>
      This small examples have been produced with the provided example data set 
      (using the same global paramters as in the previous tutorial). 
      An over-representation analysis was conducted utilizing the 46 enriched 
      proteins from the comparison PGRMC1__vs__MIAPACA. The "Show genes in 
      functional enrichment?" button was selected.

      <center><img src="query_tutorial/ora.png" width="100%"></center>
      </p>
      
      <p>
      Sorting the output table from most significant p-value to least significant 
      we find the term "actin binding" on top of the list. 15 proteins from the 
      enriched proteins are annotated with this term.
      <center><img src="query_tutorial/oratable.png" width="100%"></center>
      </p>
      
      <p>
      
      All visualizations (heatmap, fold change plot and PPI network) work 
      only on the proteins selected in the above output table we can filter that 
      table to only show proteins annotated with our term of interest. The output 
      table can parse "regular expressions", so all we need to do is to copy paste 
      the comma-delimited gene names into a text editor (or text processing tool 
      like MS Word) and replace all commas with viertical line symbols ("|" which 
      is the logical "or" operator) with the "Find and Replace" tool:
      
      <center><img src="query_tutorial/edit.png" width="100%"></center>
      </p>
      
      <p>
      We can now paste the the vertical line delimited gene names into the 
      "Gene.names" search bar in the output table and we could successfully subset 
      the data table to only show proteins of our interest:
      <center><img src="query_tutorial/filtereddt.png" width="100%"></center>
      
      Below the table there is now a text message telling us that the original 
      table has been filtered and that only the remaining proteins are used in 
      subsequent visualizations. As an example you can now observe how the selected 
      proteins compare across different group comparisons in a heatmap or fold 
      change plot.
      </p>
           '),
      size = "l",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$timeOut, { 
    print(paste0("Session (", session$token, ") timed out at: ", Sys.time()))
    showModal(modalDialog(
      title = "Timeout",
      paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
      footer = NULL
    ))
    session$close()
  })
  
}