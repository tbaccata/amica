source("global.R")

server <- function(input, output, session) {
  # max 75MB file upload restriction
  options(shiny.maxRequestSize = 75 * 1024 ^ 2)
  
  source('R/server/reactiveValues.R', local = TRUE)
  source('R/server/shinyjs.R', local = TRUE)
  source('R/server/colorsFactors.R', local = TRUE)
  source('R/server/upload.R', local = TRUE)
  source('R/server/runAnalysis.R', local = TRUE)
  source('R/server/exampleData.R', local = TRUE)
  source('R/server/modals.R', local = TRUE)
  source('R/server/qcPlots.R', local = TRUE)
  source('R/server/diffAbundancePlots.R', local = TRUE)
  
  ### INPUT LOGIC AND PLOTS
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
      theme_cowplot(font_size = 14) + 
      background_grid() +
      scale_fill_manual(values=myColors() ) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 10
      ),legend.title = element_blank()) + xlab("")
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = "svg",
        width = 676,
        height = 676,
        filename = "mtcars_bar_example"
      )
    )
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
             modeBarButtonsToRemove = removePlotlyBars,
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
      geom_point() + theme_cowplot(font_size = 14) +  
      background_grid() + 
      scale_color_manual(values=myScatterColors() )
    
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
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
    req(reacValues$uploadSuccess)
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
  
  output$uploadSuccess <- reactive({
    reacValues$uploadSuccess
  })
  outputOptions(output, "uploadSuccess", suspendWhenHidden = FALSE)
  
  output$analysisSuccess <- reactive({
    reacValues$analysisSuccess
  })
  outputOptions(output, "analysisSuccess", suspendWhenHidden = FALSE)
  
  output$uploadSummary <- renderText({
    req(reacValues$uploadSuccess)
    
    if (!is.null(reacValues$analysisSuccess) || !reacValues$analysisSuccess ||
        reacValues$amicaInput == TRUE ||
        input$source == "example"
        ) return(NULL)
    
    paste0(
      "Successfully uploaded data!\n",
      "Open the 'Advanced' tab in the sidebar to choose parameters."
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
      "\t1) QC (generate QC report)\n",
      "\t2) Quantitative results (generate Diff. Abundance report)\n",
      "\t3) Protein-protein interaction networks (only applicable for H.sapiens at the moment)\n",
      "\t4) or upload another amica file to compare experiments!"
    )
  })
  
  # output$inputParameterSummary <- renderText({
  #   req(reacValues$analysisSuccess)
  #   req(reacValues$inputParameterSummary)
  #   reacValues$inputParameterSummary
  # })
  output$analysisSuccessMsg <- renderUI({
    req(reacValues$analysisSuccess)
    if (reacValues$amicaInput) return(NULL)
    
    HTML(
      "<h4>Successfully analyzed data!</h4>
      <p>
      Scroll to the top of the page to visit the analysis tabs.
      </p>
      <p>
      The table below describes parameters used in the analysis. <br>
      Please download this table in order to reproduce the analysis.
      </p>"
    )
  })
  
  analysisParams <- reactive({
    req(reacValues$inputParameterSummary)
    dat <- read.table(
      text = reacValues$inputParameterSummary,
      sep = '\t',
      header = F
    )
    names(dat) <- c("Analysis option", "Parameter")
    dat$`Analysis option` <- gsub(":", "", dat$`Analysis option`)
    dat
  })
  
  output$inputParamDT <- renderDT({
    req(reacValues$uploadSuccess)
    req(reacValues$inputParameterSummary)
    
    datatable(
      analysisParams(),
      caption = "Summary of analysis parameters used",
      rownames = F,
      extensions = 'Buttons',
      
      options = list(
        pageLength = 12,
        autoWidth = TRUE,
        dom = 'Bfrtip',
        buttons = list(list(
          extend = 'csv',
          filename = paste0('parameter_summary_', Sys.Date())
        ))
      )
    )
  })
  
  
  ### QC LOGIC AND PLOTS 
  
  output$quantSummary <- renderText({
    req(reacValues$proteinData)
    numQuant <- length(isQuantRnames(reacValues$proteinData))
    paste0("Number of quantified proteins in ImputedIntensity: ", numQuant)
  })
  
  # ------------------------------------------- intensity boxplots
  
  dataAssay <- reactive({
    req(input$assayNames)
    req(reacValues$proteinData)
    
    object <- 0
    if (input$assayNames == "ImputedIntensity") {
      object <- assay(reacValues$proteinData, "ImputedIntensity")
      
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
    
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(object$group))) {
      object <- object[object$group %in% reacValues$groupFactors,]
      object$group <-
        factor(object$group, levels = reacValues$groupFactors)
    }
    
    object
  })
  
  
  # ------------------------------------------- PCA
  
  get_pca_data <- reactive({
    event_data("plotly_selected", source = "subset")
  })
  

    output$pca <- renderPlotly({
    validate(need(
      input$assayNames != "",
      "Please press submit or provide an intensity prefix."
    ))
    validate(
      need(
        length(input$pcaSamplesInput) > 1 |
          is.null(input$pcaSamplesInput),
        "Please select at least two groups. If none is selected all are considered."
      )
    )
    
    plotData <- assay(reacValues$proteinData, input$assayNames)
    
    withProgress(message = "Plotting PCA ", {
      olist <-
        computePCA(
          plotData,
          colData(reacValues$proteinData),
          reacValues$groupFactors,
          input$pcaSamplesInput
        )
      
      df_out <- olist[[1]]
      percentage <- olist[[2]]
      
      df_out$show_id <- FALSE
      if (!is.null(get_pca_data())) {
        df_out[df_out$key %in% get_pca_data()$key, "show_id"] <- TRUE
      }
      
      p <- plotPCA(
        df_out,
        percentage,
        myGroupColors(),
        input$assayNames,
        input$pca_base,
        input$pca_legend,
        input$pca_pointsize
      )
      
    })
    ggplotly(p, source = "subset") %>%
      layout(dragmode = "select") %>%
      config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
        toImageButtonOptions = list(
          format = input$pca_format,
          width = input$pca_width,
          height = input$pca_height,
          filename = paste0("pca_", input$assayNames)
        )
      )
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
  
  
  ### ------------------------------------------- INTENSITY BOX PLOT
  
  output$boxPlot <- renderPlotly({

    withProgress(message = "Plotting boxplot of intensities", {
      ggplotly(plotBoxPlot(
        dataAssay(),
        input$assayNames,
        myGroupColors(),
        input$boxplot_base,
        input$boxplot_legend
      )) %>% config(
          displaylogo = F,
          modeBarButtonsToRemove = removePlotlyBars,
          toImageButtonOptions = list(
            format = input$boxplot_format,
            width = input$boxplot_width,
            height = input$boxplot_height,
            filename = paste0("boxplot_", input$assayNames)
          )
        )
    })
  })
  
  ### ------------------------------------------- DENSITY PLOT
  
  output$densityPlot <- renderPlotly({
    
    withProgress(message = "Plotting density plot", {
      p <- plotDensityPlot(
        dataAssay(),
        input$assayNames,
        myColors(),
        base_size = input$density_base,
        legend_size = input$density_legend
      )
      ggplotly(p) %>%
        config(
          displaylogo = F,
          modeBarButtonsToRemove = removePlotlyBars,
          toImageButtonOptions = list(
            format = input$density_format,
            width = input$density_width,
            height = input$density_height,
            filename = paste0("density_plot_", input$assayNames)
          )
        )
    })
  })
  
  # ------------------------------------------- MOST ABUNDANT PROTEINS IN SAMPLE Input
  
  output$samples <- renderUI({
    req(reacValues$proteinData)
    selectInput("samples",
                "Inspect a sample",
                c("", colData(reacValues$proteinData)$samples),
                selected = NULL)
  })
  
  output$assayNames <- renderUI({
    req(reacValues$proteinData)
    selectInput("assayNames",
                "Which intensities do you want to plot?",
                c(assayNames(reacValues$proteinData)),
                selected = "ImputedIntensity")
  })
  
  
  output$compareScatterPlot1 <- renderUI({
    selectizeInput(
      "selectScatterSample1",
      "Select a sample for the x-axis of the scatter plot",
      c("", colData(reacValues$proteinData)$samples),
      selected = NULL,
      options = list(maxItems = 1)
    )
  })
  
  output$compareScatterPlot2 <- renderUI({
    selectizeInput(
      "selectScatterSample2",
      "Select a sample for the y-axis of the scatter plot",
      c("", colData(reacValues$proteinData)$samples),
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
  
  get_scatter_data <- reactive({
    event_data("plotly_selected", source = "subset")
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
    validate(need(!is.null(assay1),
                  "Please enter an assay for x-axis."))
    validate(need(!is.null(assay2),
                  "Please enter an assay for y-axis."))
    
    if (assay1 == assay2) {
      validate(need(
        selection1 != selection2,
        "Cannot plot the same column on x - and y-axis."
      ))
    }
    
    pu <- plotScatterPlot(
      reacValues$proteinData,
      assay1,
      assay2,
      selection1,
      selection2,
      myScatterColors(),
      plot_fontsize,
      plot_legendsize,
      showLine,
      get_scatter_data()
    )
    
    xLabel <- paste0(assay1, "_", selection1)
    yLabel <- paste0(assay2, "_", selection2)

    m = list(
      l = 100,
      r = 40,
      b = 100,
      t = 50,
      pad = 0
    )
    
    ggplotly(pu, margin = m, source = "subset") %>% 
      layout(dragmode = "select") %>%
      config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
        toImageButtonOptions = list(
          format = scatter_format,
          width = plot_width,
          height = plot_height,
          filename = paste0("scatter_plot_", xLabel, "_", yLabel)
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
  
  output$corrPlotly <- renderPlotly({
    withProgress(message = "Plotting correlation plot ", {
      
      p <- plotCorrPlotly(
        assay(reacValues$proteinData, input$assayNames),
        colData(reacValues$proteinData),
        myGroupColors(),
        heatColors(),
        input$assayNames,
        input$cor_annot,
        reacValues$groupFactors,
        input$corrSamplesInput
      )
    })
    p %>%  config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = input$cor_format,
        width = input$cor_width,
        height = input$cor_height,
        filename = paste0("corrplot_", input$assayNames)
      )
    ) %>% layout(
      font = list(
        family = "",
        size = 16,
        color = "black"
      ))
  })
  
  output$overlapHeatmapPlotly <- renderPlotly({
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    
    lfqs <- assay(reacValues$proteinData, aname)
    
    outp <- plotOverlaply(
      lfqs,
      colData(reacValues$proteinData),
      myGroupColors(),
      heatColors(),
      annotSamples = input$overlapHeatmap_annot,
      metric = input$overlapHeatmap_metric,
      groupFactors = reacValues$groupFactors
    )
    
    reacValues$overlapDf <- outp$df
    outp$plot
    
    
    withProgress(message = "Plotting overlap plot ", {
      p <- outp$plot
    })
    p %>%  config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = input$overlapHeatmap_format,
        width = input$overlapHeatmap_width,
        height = input$overlapHeatmap_height,
        filename = paste0("overlap_heatmap_", input$overlapHeatmap_metric)
      )
    )
  })
  
  output$overlapSummaryDT <- renderDT({
    req(reacValues$overlapDf )
    df <- reacValues$overlapDf
    names(df) <- c(
      "Sample 1",
      "Sample 2",
      "# Sample 1",
      "# Sample 2",
      "# shared",
      "Overlap Coeff.",
      "Jaccard index"
    )
    datatable(
      df,
      extensions = 'Buttons',
      filter = "top",
      rownames = F,
      options = list(
        pageLength = 10,
        autoWidth = TRUE
      )
    )
  })
  
  
  ### ------------------------------------------- CONTAMINANTS plots
  
  contaminantsPlotly <- reactive({
    validate(need(
      abundancePrefix %in% assayNames(reacValues$proteinData),
      "No data to plot (no iBAQ columns available)."
    ))

    plotContaminants(
      reacValues$proteinData,
      colData(reacValues$proteinData),
      myGroupColors(),
      groupFactors = reacValues$groupFactors,
      contaminants_base = input$contaminants_base,
      contaminants_legend = input$contaminants_legend
    )
  })
  
  
  output$contaminants <- renderPlotly({
    validate(need(
      abundancePrefix %in% assayNames(reacValues$proteinData),
      "No data to plot (no iBAQ columns available)."
    ))
    withProgress(message = "Plotting barplot of contaminants", {
      ggplotly(contaminantsPlotly())  %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
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
    validate(need(
      abundancePrefix %in% assayNames(reacValues$proteinData),
      "No data to plot (no iBAQ columns available)."
    ))
    
    withProgress(message = "Plotting most abundant proteins", {
      
      tmp <- toLongFormat(
        2 ^ assay(reacValues$proteinData, abundancePrefix),
        reacValues$proteinData,
        addGroup = TRUE,
        addContaminant = TRUE,
        addGeneName = TRUE
      )
      
      ggplotly(plotMostAbundantProteinsInSample(
        dataAssay=tmp,
        sample=input$samples,
        color = myScatterColors()[1],
        abundant_base = input$abundant_base
      ))  %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
        toImageButtonOptions = list(
          format = input$abundant_format,
          width = 676,
          height = 676,
          filename = paste0("most_abundant_proteins_", input$samples)
        )
      )
    })
  })
  
  ### ------------------------------------------- NUM ID PROTEINS
  
  numIdPlotly <- reactive({
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
    
    plotNumberIdentifiedProteins(
      assay(reacValues$proteinData, aname),
      colData(reacValues$proteinData),
      myGroupColors(),
      groupFactors = reacValues$groupFactors,
      barplotId_base = input$barplotId_base,
      barplotId_legend = input$barplotId_legend
    )
  })
  
  output$barplotProteins <- renderPlotly({
    validate(
      need(
        "Intensity" %in% assayNames(reacValues$proteinData) |
          "LFQIntensity" %in% assayNames(reacValues$proteinData),
        "No data to plot (no columns starting with 'LFQIntensity_' or 'Intensity_' available)."
      )
    )
    withProgress(message = "Plotting Number of inferred protein groups", {
      ggplotly(numIdPlotly())  %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
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
  
  pctMvsPlotly <- reactive({
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    
    plotMissingValues(
      assay(reacValues$proteinData, aname),
      colData(reacValues$proteinData),
      myGroupColors(),
      groupFactors = reacValues$groupFactors,
      barplotMv_base = input$barplotMv_base,
      barplotMv_legend = input$barplotMv_legend
    )
  })
  
  output$barplotMissingValues <- renderPlotly({
    validate(
      need(
        "Intensity" %in% assayNames(reacValues$proteinData) ||
          "LFQIntensity" %in% assayNames(reacValues$proteinData),
        "No data to plot (no columns starting with 'LFQIntensity_' or 'Intensity_' available)."
      )
    )
    
    withProgress(message = "Plotting barplot of missing values ", {
      ggplotly(pctMvsPlotly()) %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
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
  
  output$boxplotCV <- renderPlotly({
    req(reacValues$uploadSuccess)
    validate(need(
      any(duplicated(
        colData(reacValues$proteinData)$groups
      )),
      "Cannot output CV plot for a pilot without replicates"
    ))

    withProgress(message = "Plotting Coefficient of Variations ", {
      ggplotly(plotCoeffVarPlot(dataAssay(),
                       myGroupColors(),
                       input$assayNames,
                       input$cv_base,
                       input$cv_legend)) %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
        toImageButtonOptions = list(
          format = input$cv_format,
          width = input$cv_width,
          height = input$cv_height,
          filename = paste0("cv_plot_", input$assayNames)
        )
      )
    })
  })
  
  
  ### DIFFERENTIAL ABUNDANCE LOGIC AND PLOTS 

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
  
  output$dotplotGroupComparisons <- renderUI({
    req(reacValues$dataLimma)
    
    isNotPilot <-
      ifelse(length(grep(
        paste0(padjPrefix, collapse = '|'),
        names(reacValues$dataLimma),
        value = T
      )) > 0,
      TRUE,
      FALSE)
    
    validate(need(
      isNotPilot==TRUE,
      "Error! Cannot output Dotplot for pilots, Dotplot needs p-values."
    ))
    
    validate(need(
      reacValues$compareAmicaSelected == FALSE,
      paste0(
        "Cannot output Dotplot for 'multiple_amica_upload'\n
                         Change the data set back to the original data to
                         plot a Dotplot"
      )
    ))
    
    comps <-
      gsub(logfcPrefix, "", grep(logfcPrefix, names(reacValues$dataLimma), value = T))
    
    selectInput(
      "dotplotGroupComparisons",
      "Choose log2FC and padj statistics from group comparisons:",
      comps,
      multiple = T
    )
  })
  
  output$dotplotGroups <- renderUI({
    input$dotplotSelection
    req(input$dotplotGroupComparisons)
    selected <- isolate(input$dotplotGroupComparisons)
    numVars <- length(selected)
    
    validate(need(
      numVars > 1,
      "Need at least two groups for the Dotplot."
    ))

    C = sapply(1:numVars, function(i){paste0("dotplot_comp_",i)})
    L = sapply(1:numVars, function(i){paste0("dotplot_group_",i)})

    output = tagList()

    for(i in seq_along(1:numVars)){
      output[[i]] = tagList()
      output[[i]][[1]] = selectInput(C[i], "Group comparison:", selected[i])
      groups <- unlist(strsplit(selected[i], '__vs__'))
      output[[i]][[2]] = selectInput(L[i], label = "Which group to show?", groups, multiple = F)
    } ## for loop
    output
  })
  
  output$boolUniqueGroups <- renderUI({
    req(reacValues$dotplotGroupsDf)
    group <- reacValues$dotplotGroupsDf$group
    
    if (any(duplicated(group))) {
      HTML('<h4 style="color:red;">Error!</h4><p style="color:red;">
         Cannot plot Dotplot with duplicated groups.<br>
         Please select unique groups.
         </p>')
    } else {
      return(NULL)
    }
  })

  observeEvent(input$dotplotSelection, {
    req(input$dotplotGroupComparisons)
    selected <- isolate(input$dotplotGroupComparisons)
    comparison <- sapply(1:length(selected), function(i) {input[[paste0("dotplot_comp_",i)]]})
    group <- sapply(1:length(selected), function(i) {input[[paste0("dotplot_group_",i)]]})
    reacValues$dotplotGroupsDf <- data.frame(comparison, group)
  })
  
  # output$submitDotplot <- renderUI({
  #   req(reacValues$dotplotGroupsDf)
  #   actionButton("submitDotplot", label = "Plot Dotplot")
  
  observeEvent(
    c(
      reacValues$dataComp,
      reacValues$dotplotGroupsDf,
      input$groupComparisonsDT_rows_all
    ),
    {
      
    req(reacValues$dotplotGroupsDf)
    #req(reacValues$dataComp)
    
    validate(need(!any(duplicated(reacValues$dotplotGroupsDf$group)), "Error! Please provide 
                  unique groups"  )  )
    
    #group2comps <- reacValues$dotplotGroupsDf
    ridx <- input$groupComparisonsDT_rows_all
    ridx <- rownames(reacValues$dataComp[ridx,])
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    validate(need(aname %in% assayNames(reacValues$proteinData), 
                  "No LFQIntensity or Intensity columns found." ) )
    
    pattern <- padjPrefix
    if (reacValues$sigCutoffValue == "p-value") pattern <- pvalPrefix
    
    selection <- grep(pattern, names(reacValues$dataLimma), value = T)
    widestats <- reacValues$dataLimma[ridx, selection]
    
    selection <- grep(logfcPrefix, names(reacValues$dataLimma), value = T)
    widefcs <- reacValues$dataLimma[ridx, selection]
    
    df <- assay(reacValues$proteinData, aname)
    df <- df[ridx,]
    
    if (reacValues$show_dotplot == FALSE) {
      reacValues$show_dotplot = TRUE
      toggle(id = 'hide_dotplot_before_submit', anim = T)
    }
    
    reacValues$dataDotplot <- getDotplotData(widestats,
               widefcs,
               df,
               reacValues$filtData,
               reacValues$dotplotGroupsDf,
               reacValues$expDesign,
               reacValues$sigCutoffValue)
  })
  
  output$dotplot_color_gradient <- renderUI({
    req(reacValues$dataDotplot)
    minVal <- round(min(reacValues$dataDotplot$log2FC, na.rm = T), 3)
    maxVal <- round(max(reacValues$dataDotplot$log2FC, na.rm = T), 3)
    sliderInput(
      'dotplot_color_gradient',
      'Define Dotplot color gradient',
      min = minVal,
      max = maxVal,
      value = c(
        ifelse(
          is.null(input$dotplot_color_gradient[1]),
          minVal,
          input$dotplot_color_gradient[1]
        ),
        ifelse(
          is.null(input$dotplot_color_gradient[2]),
          maxVal,
          input$dotplot_color_gradient[2]
        )
      )
    )
  })
  
  output$dotplot_size_gradient <- renderUI({
    req(reacValues$dataDotplot)
    sliderInput(
      'dotplot_size_gradient',
      'Define Dotplot size gradient',
      1,
      8,
      value = c(
        ifelse(
          is.null(input$dotplot_size_gradient[1]),
          1,
          input$dotplot_size_gradient[1]
        ),
        ifelse(
          is.null(input$dotplot_size_gradient[2]),
          6,
          input$dotplot_size_gradient[2]
        )
      )
    )
  })
  
  output$dotplot_clustering_option <- renderUI({
    req(reacValues$dataDotplot)
    
    options <- c("log2FC", "AvgIntensity")

    selectInput("dotplot_clustering_option", 
                "How to cluster data in Dotplot?",
                choices = options,
                selected = input$dotplot_clustering_option,
                multiple = F)
  })
  
  output$dotplot_ctrl_substraction <- renderUI({
    req(input$dotplot_clustering_option)
    
    if (input$dotplot_clustering_option != "log2FC") return(NULL)
    
    checkboxInput("dotplot_ctrl_substraction",
                  "Control substraction: Only show proteins with log2FC > 0?",
                  value = T)
  })
    
  
  numberOfDotplotPoints <- reactive({
    req(reacValues$dataDotplot)
    max(450, 20 * length(unique(reacValues$dataDotplot$Gene)))
  } )
  
  numberOfDotplotBaits <- reactive({
    req(reacValues$dataDotplot)
    numCols <- length(unique(reacValues$dataDotplot$Group))
    
    value <- 300 + numCols * 20

    value
  })
  
  
  output$dotplot_cluster_columns <- renderUI({
    req(reacValues$dataDotplot)
    
    checkboxInput("dotplot_cluster_columns",
                  "Cluster columns?",
                  value = ifelse(!is.null(input$dotplot_cluster_columns),
                                 input$dotplot_cluster_columns,
                                 T))
  })
  
  dotplot <- reactive({
    validate(need(
      reacValues$compareAmicaSelected == FALSE,
      paste0(
        "Cannot output Dotplot for 'multiple_amica_upload'\n
                         Change the data set back to the original data to
                         plot a Dotplot"
      )
    ))
    req(reacValues$dataDotplot)
    req(reacValues$dotplotGroupsDf)

    validate(need(length(input$groupComparisonsDT_rows_all) > 1,
                  "Cannot output Dotplot with less than two proteins." ) )
    
    minColorGradient <- input$dotplot_color_gradient[1]
    maxColorGradient <- input$dotplot_color_gradient[2]
    
    minSizeGradient <- input$dotplot_size_gradient[1]
    maxSizeGradient <- input$dotplot_size_gradient[2]
    
    dotplotColors <- dotplotColorPalette(input$dotplot_palette,
                                         input$dotplot_rev_colors)
    
    
    clusteringMetric <- "AvgIntensity"
    if (!is.null(input$dotplot_clustering_option)) {
      clusteringMetric <- input$dotplot_clustering_option
    }
    
    plotDotplot(reacValues$dataDotplot,
                reacValues$dotplotGroupsDf,
                dotplotColors,
                minColorGradient,
                maxColorGradient,
                minSizeGradient,
                maxSizeGradient,
                clusteringMetric,
                input$dotplot_distance_metric,
                input$dotplot_clustering_method,
                input$dotplot_cluster_columns,
                reacValues$sigCutoffValue,
                input$dotplot_ctrl_substraction)
  })

  output$dotplot <- renderPlot({
    dotplot()
  })
  
  output$plot.ui <- renderUI({
    plotOutput('dotplot', height = numberOfDotplotPoints(),
               width = numberOfDotplotBaits(),
               hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")
               )
  })
  
  output$downloadDotPlot <- downloadHandler(
    filename = function(){paste("dotplot",'.pdf',sep='')},
    content = function(file){
      cowplot::ggsave2(file,plot=dotplot(),
                       width = numberOfDotplotBaits() * 1/72, 
                       height = numberOfDotplotPoints() * 1/72,
                       limitsize = F,
                       device = cairo_pdf)
    })
  
  output$downloadDotPlotData <- downloadHandler(
    filename = function(){paste("dotplot_data",'.tsv',sep='')},
    content = function(file){
      write.table(reacValues$dataDotplot[-which(names(reacValues$dataDotplot) == "significant")], 
                  file, row.names = F, sep = '\t', quote = F)
    })

  output$hover_info <- renderUI({
    req(reacValues$dataDotplot)
    hover <- input$plot_hover
    point <- nearPoints(reacValues$dataDotplot, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    left_px <- hover$coords_css$x
    top_px <- hover$coords_css$y
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px;")
    
    text <- paste0("<b> Gene: </b>", point$Gene, "<br/>",
                   "<b> logFC: </b>", round(point$log2FC,2), "<br/>",
                   "<b>", reacValues$sigCutoffValue , ": </b>", scales::scientific(point$padj, digits = 3), "<br/>",
                   "<b> Comparison: </b>", point$Comparison, "<br/>",
                   "<b> Group: </b>", point$Group, "<br/>",
                   "<b> AvgIntensity: </b>", round(point$AvgIntensity, 2), "<br/>"
                   )
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(text))
    )
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
        buttons = list(list(
          extend = 'csv',
          filename = paste0(
            reacValues$enrichmentChoice,
            '_proteins_log2fcThresh_',
            reacValues$fcCutoff, 
            '_signcutoff_',
            reacValues$sigCutoffValue
          )
        ))
      )
    ) %>% formatRound(paste0(logfcPrefix,selection), 4) 
    
    if (hasPadj) dt <- dt %>% formatSignif(grep("P.Val", names(tmp)), digits = 4)
    dt
  }, server = F)
  
  
  enrichedMatrixSet <- reactive({
    req(reacValues$dataLimma)
    
    matrixData <-
      generateEnrichedMatrix(
        reacValues$dataLimma,
        input$enrichmentChoice,
        input$sigCutoffValue,
        input$fcCutoff,
        input$pvalCutoff
      )
    as.data.frame(matrixData)
  })
  
  observeEvent(c(input$submitMultiComp, reacValues$nsubmits),{
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
  
  numOfComps <- reactive({
    input$upset1Sample
  })
  
  output$multiCompLabelsInput <- renderUI({
    req(input$upset1Sample)
    numVars <- length(input$upset1Sample)
    
    C = sapply(1:numVars, function(i){paste0("multi_comp_cols_",i)})
    L = sapply(1:numVars, function(i){paste0("multi_label_",i)})
    
    output = tagList()
    
    for(i in seq_along(1:numVars)){
      output[[i]] = tagList()
      output[[i]][[1]] = selectInput(C[i], "Variable to summarize:", numOfComps()[i])
      output[[i]][[2]] = textInput(L[i], label = "Label for variable:", 
                                   value = "")
    } ## for loop
    
    output
  })
  
  observeEvent(input$changeMultiCompNames, {
    oldNames <- sapply(1:length(numOfComps()), function(i) {input[[paste0("multi_comp_cols_",i)]]})
    newNames <- sapply(1:length(numOfComps()), function(i) {input[[paste0("multi_label_",i)]]})
    reacValues$newMultiNames <- list(old=oldNames, new=newNames)
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
  
  eulerData <- function() {
    comparisons <- isolate(input$upset1Sample)
    binMat <- isolate(enrichedMatrixSet())
    showQuant <- isolate(input$euler_quant)
    bool <- isolate(input$euler_line)
    showLegend <- isolate(input$euler_legend)
    
    validate(need(!is.null(comparisons), "Please select at least two comparisons."))
    validate(need(length(comparisons) > 1, 
                  "Need at least two comparisons to render UpSet plot."))
    validate(need(length(comparisons) < 6, "Cannot output Euler plot for more than 5 sets."))

    cols <- myScatterColors()[1:length(comparisons)]
    if (!is.null(input$eulercol_1)) {
      tmp <- sapply(1:length(comparisons), function(i) {input[[paste0("eulercol_",i)]]})
      if (all(!is.null(tmp))) cols <- tmp
    }
    
    plotEulerDiagram(comparisons,
               binMat,
               showQuant,
               bool,
               showLegend,
               reacValues$newMultiNames,
               cols)
  }
  
  output$eulerrPlot <- renderPlot({
    input$submitMultiComp
    reacValues$nsubmits
    print(eulerData())
  })
  
  output$download_button_eulerr <- renderUI({
    downloadButton("download_eulerr", "Download Euler plot", icon = icon("camera"))
  })
  
  output$download_eulerr <- downloadHandler(
    filename = function() {
      paste("eulerr_plot.pdf")
    },
    content = function(file) {
      pdf(file, width = input$euler_width, 
          height = input$euler_width, onefile = F)
      print(eulerData())
      dev.off()
    }
  )

  plotMultiUpset <- function() {
    matrixSet <- isolate(enrichedMatrixSet() )
    samples <- isolate(input$upset1Sample)
    scale <- isolate(input$upset_scale)
    fcThresh <- isolate(input$fcCutoff)
    enrichmentChoice <- isolate(input$enrichmentChoice)
    sigCutoffValue <- isolate(input$sigCutoffValue)
    
    reacValues$sigCutoffValue <- sigCutoffValue
    reacValues$selection <- samples
    reacValues$fcCutoff <- fcThresh
    reacValues$enrichmentChoice <- enrichmentChoice

    upset_ratio <- isolate(input$upset_ratio)
    upset_pointsize <- isolate(input$upset_pointsize)
    upset_sorted <- isolate(input$upset_sorted)

    validate(need(!is.null(samples), "Please select at least two comparisons."))
    validate(need(length(colnames(matrixSet) ) > 1, "Need at least two comparisons to render UpSet plot. Only one provided.") )
    if (!is.null(samples)) {
      validate(need(length(samples) > 1, "Please select at least two comparisons.") )
    }
    validate(need(any(matrixSet[,samples]==1),
                  "There are no significant proteins to display."))
    
    plotUpsetPlot(
      matrixSet = matrixSet,
      samples = samples,
      scale=scale,
      upset_ratio = upset_ratio,
      upset_pointsize = upset_pointsize,
      upset_sorted = upset_sorted,
      newMultiNames = reacValues$newMultiNames
    )
  }
  
  output$upsetPlot <- renderPlot({
    input$submitMultiComp
    reacValues$nsubmits
    req(reacValues$dataComp)

    withProgress(message = "Comparing multigroup comparisons", {
      print(plotMultiUpset())
      if (reacValues$show_analysis == FALSE) {
        reacValues$show_analysis = TRUE
        toggle(id = 'hide_before_comparisons', anim = T)
      }
    })
  })
  
  output$download_button_upset <- renderUI({
    #req(upsetPlot())
    downloadButton("download_upset", "Download UpSet plot", icon = icon("camera"))
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
  
  observeEvent(input$submitHeatmap, {
    if (reacValues$show_heatmap == FALSE) {
      reacValues$show_heatmap = TRUE
      toggle(id = 'hide_heatmap_before_submit', anim = T)
    }
  })
  
  compareHeatmapBase <- reactive({
    input$submitHeatmap
    validate(need(reacValues$compareAmicaSelected==FALSE,
                  paste0("Cannot output heatmap for 'multiple_amica_upload'\n
                         Change the data set back to the original data to 
                         plot a heatmap.") ))
    
    req(reacValues$dataComp)
    #ridx <- isolate(input$groupComparisonsDT_rows_all)
    ridx <- input$groupComparisonsDT_rows_all
    ridx <- rownames(reacValues$dataComp[ridx,])
    
    validate(need(nrow(reacValues$dataComp)>1, "No data to plot"))
    validate(need(length(ridx)>1, "No data to plot. Please remove filters in datatable."))
    
    heatmapSamplesInput <- isolate(input$heatmapSamplesInput)
    
    df <- assay(reacValues$proteinData, "ImputedIntensity")
    df <- df[ridx,]
    
    reacValues$dataHeatmap <- getHeatmaplyData(df,
                                               rowData(reacValues$proteinData)[ridx, geneName],
                                               reacValues$expDesign,
                                               heatmapSamplesInput
    )

    #req(reacValues$dataHeatmap)
    fontsize = isolate(input$heatmap_base)
    plot_width = isolate(input$heatmap_width)
    plot_height = isolate(input$heatmap_height)
    show_row_labels <- isolate(input$heatmap_row_labels)
    show_col_labels <- isolate(input$heatmap_col_labels)
    show_annot <- isolate(input$heatmap_annot)
    clusterRows = isolate(input$clusterRows)
    clusterCols = isolate(input$clusterCols)
    scaleHeatmap = isolate(input$scaleHeatmap)
    
    validate(need(all(names(reacValues$dataHeatmap) %in% colData(reacValues$proteinData)$samples),""))
    
    validate(need(
      nrow(reacValues$dataHeatmap) >= 1,
      paste0(
        "Need more than ",
        nrow(reacValues$dataHeatmap),
        " proteins selected. Apply less stringent thresholds."
      )
    ))
    
    withProgress(message = "Plotting heatmap ", {
      annot <- colData(reacValues$proteinData)
      p <- plotHeatmaply(dataHeatmap=reacValues$dataHeatmap,
                 annot=annot,
                 myGroupColors=myGroupColors(),
                 heatColors=heatColors(),
                 fontsize=fontsize,
                 scaleHeatmap = scaleHeatmap,
                 show_annot=show_annot,
                 show_col_labels = show_col_labels,
                 show_row_labels = show_row_labels,
                 clusterRows = clusterRows,
                 clusterCols = clusterCols,
                 plot_method = "plotly"
                 )
      calc_height <- min(15 * nrow(reacValues$dataHeatmap), 1200)
      if (calc_height < 600)
        calc_height <- 800
      
      p %>% layout(height = calc_height)
    })
  })
  
  output$compareHeatmap <- renderPlotly({
    input$submitHeatmap
    validate(need(length(input$groupComparisonsDT_rows_all)>1, 
                  "No data to plot. Please remove filters in datatable."))
    
    if (reacValues$show_heatmap == FALSE) {
      reacValues$show_heatmap = TRUE
      toggle(id = 'hide_heatmap_before_submit', anim = T)
    }
    
    compareHeatmapBase() %>%
      config(displaylogo = F,
             modeBarButtonsToRemove = removePlotlyBars,
             toImageButtonOptions = list(format = input$heatmap_format,
                                         width = input$heatmap_width,
                                         height = input$heatmap_height,
                                         filename = "heatmap")
      )
  })
  
  output$foldChangeSelection <- renderUI({
    selectizeInput(
      "foldChangeSelection",
      "Select 2 comparisons for a fold change plot",
      c("",reacValues$reacConditions),
      selected = NULL,
      options = list(maxItems = 2)
    )
  })
  
  output$fcPlotLabelsInput <- renderUI({
    req(input$foldChangeSelection)
    numVars <- length(input$foldChangeSelection)
    
    C = sapply(1:numVars, function(i){paste0("fcplot_cols_",i)})
    L = sapply(1:numVars, function(i){paste0("fcplot_label_",i)})
    
    output = tagList()
    
    for(i in seq_along(1:numVars)){
      output[[i]] = tagList()
      output[[i]][[1]] = selectInput(C[i], "Variable to summarize:", input$foldChangeSelection[i])
      output[[i]][[2]] = textInput(L[i], label = "Label for variable:", 
                                   value = "")
    }
    output
  })
  
  fcPlotLabels <- reactive({
    unlist(sapply(1:2, function(i) {input[[paste0("fcplot_label_",i)]]}))
  })
  
  get_fc_data <- reactive({
    event_data("plotly_selected", source = "subset")
  })
  
  output$foldChangePlot <- renderPlotly({
    input$sumbitFoldChangePlot
    req(reacValues$dataComp)
    #req(input$foldChangeSelection)
    matrixSet <- isolate(enrichedMatrixSet() )
    
    labels <- isolate(fcPlotLabels())
    selection <- isolate(input$foldChangeSelection)
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    
    fontsize <- isolate(input$fc_base)
    legendFontSize <- isolate(input$fc_legend)
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
                    selection,
                    labels)
    
    plotData$show_id <- FALSE
    if (!is.null(get_fc_data() )) {
      plotData[plotData$key %in% get_fc_data()$key, "show_id"] <- TRUE
    }
    
    withProgress(message = "Plotting Fold Change plot", {
    pu <- plotFoldChangePlot(plotData,
                       selection,
                       colors,
                       labels,
                       showLine,
                       fontsize,
                       legendFontSize,
                       pointsize)
    })

    p <- ggplotly(pu, source = "subset") %>% layout(dragmode = "select")
    p %>% config(displaylogo = F,
             modeBarButtonsToRemove = removePlotlyBars,
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
  
  volcanoPlotData <- reactive({
    input$maVolcanoSubmit
    reacValues$nsubmits
    sample <- isolate(input$maVolcanoSampleSelect)
    fcThresh <- isolate(input$fcCutoff)
    choice <- isolate(input$enrichmentChoice)
    sigCutoffValue <- isolate(input$sigCutoffValue)
    #padjY <- isolate(input$volcano_padj_y)
    pvalThresh <- isolate(input$pvalCutoff)
    
    reacValues$selection <- sample
    reacValues$fcCutoff <- fcThresh
    reacValues$enrichmentChoice <- choice
    reacValues$sigCutoffValue <- sigCutoffValue
    
    validate(need(!is.null(sample) & sample!="" & length(sample)>0, "Please select a group comparison."))
    padjYBoolean <- ifelse(sigCutoffValue == "p-value", FALSE, TRUE)
    #padjYBoolean <- ifelse(padjY == "p-values", FALSE, TRUE)
    
    getVolcanoPlotData(
      reacValues$dataLimma,
      sample,
      fcThresh,
      sigCutoffValue,
      choice,
      padjYBoolean,
      pvalThresh
    )
  })
  
  get_data <- reactive({
    event_data("plotly_selected", source = "subset")
  })
  
  observeEvent(input$maVolcanoSubmit, {
    if (reacValues$show_single == FALSE) {
      reacValues$show_single = TRUE
      toggle(id = 'hide_before_single_submit', anim = T)
    }
  })
  
  output$volcanoPlot <- renderPlotly({
    input$maVolcanoSubmit
    req(reacValues$dataComp)
    sample <- isolate(input$maVolcanoSampleSelect)
    validate(need(!is.null(sample) & sample!="" & length(sample)>0,
                  "Please select a group comparison."))
    validate(need(length(grep(sample, names(reacValues$dataComp))) > 0,
                  "Please select a group comparison."))
    
    fontsize <- isolate(input$volcano_base)
    legend_fontsize <- isolate(input$volcano_legend)
    plot_width <- isolate(input$volcano_width)
    plot_height <- isolate(input$volcano_height)
    padjY <- isolate(input$sigCutoffValue) #padjY <- isolate(input$volcano_padj_y)
    pointsize <- isolate(input$volcano_pointsize)
    format <- isolate(input$volcano_format)
    padjYBoolean <- ifelse(padjY == "p-value", FALSE, TRUE)
    
    c1 <- isolate(input$volcanocol_1)
    c2 <- isolate(input$volcanocol_2)
    
    pal <- myScatterColors()[1:3]
    if (is.null(c1) | is.null(c2)) {
      c1 <- pal[1]
      c2 <- pal[2]
    }
    pcols <- c(c1,c2,"red")
    setChoice <- "union"
    
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
      p <- plotVolcanoPlot(pltData,
                           xText,
                           padjYBoolean,
                           pcols,
                           fontsize,
                           legend_fontsize,
                           pointsize)

      incProgress(amount = 0.75)
      
      p <- ggplotly(p, source = "subset") %>% layout(dragmode = "select")
      incProgress(1/1.01, detail = paste("Finished!"))
      
    })
    print(p %>%
            config(displaylogo = F,
                   modeBarButtonsToRemove = removePlotlyBars,
                   toImageButtonOptions = list(format = format,
                                               width = plot_width,
                                               height = plot_height,
                                               filename = paste0("volcano_plot_", sample)
            )
    ))
  })
  
  observeEvent(c(input$maVolcanoSubmit, reacValues$nsubmits),{
    # if (reacValues$show_single == FALSE) {
    #   reacValues$show_single = TRUE
    #   toggle(id = 'hide_before_single_submit', anim = T)
    # }
    req(input$maVolcanoSampleSelect)
    sigCutoffValue <- isolate(input$sigCutoffValue)
    reacValues$sigCutoffValue <- sigCutoffValue
    sample <- isolate(input$maVolcanoSampleSelect)
    reacValues$nsubmits
    reacValues$dataComp <-
      getComparisonsData2(reacValues$dataLimma,
                          enrichedMatrixSet(),
                          'union',
                          input$maVolcanoSampleSelect)
  })
  

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
    
    validate(need(!is.null(sample) & sample!="",
                  "Please select a group comparison."))
    validate(need(length(grep(sample, names(reacValues$dataComp))) > 0,
                  "Please select a group comparison."))
    
    pltData <- volcanoPlotData()
    c1 <- isolate(input$volcanocol_1)
    c2 <- isolate(input$volcanocol_2)

    pal <- myScatterColors()[1:3]
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
    
    xText <- unlist(strsplit(sample, "__vs__"))
    xText <- paste0("log2FC(", xText[1], "/", xText[2], ")")
    
    withProgress(message = "Plotting MA plot ", {
      p <- plotMAPlot(pltData,
                      xText,
                      pcols,
                      fontsize,
                      legend_fontsize,
                      pointsize)
    })
    
    p <- ggplotly(p, source = "subset") %>% layout(dragmode = "select")
    print(p %>%
            config(displaylogo = F,
                   modeBarButtonsToRemove = removePlotlyBars,
                   toImageButtonOptions = list(format = format,
                                               width = plot_width,
                                               height = plot_height,
                                               filename = paste0("ma_plot_", sample))
            )
    )
  })
  

  observe({
    updateSelectizeInput(session,
                         server = T,
                         'selectProfilePlotGene',
                         choices = c("",reacValues$filtData[[geneName]]),
    )
  })
  
  output$profilePlot <- renderPlotly({
    req(reacValues$dataLimma)
    req(input$selectProfilePlotGene)
    plot_type <- input$profile_plot_type
    
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
    
    if (!is.null(reacValues$groupFactors) && 
        all(reacValues$groupFactors %in% unique(object$group))) {
      object <- object[object$group %in% reacValues$groupFactors,]
      object$group <- factor(object$group, levels = reacValues$groupFactors)
    }
    names(object)[which(names(object)=="rowname")] <- "Protein.IDs"
    
    isPilot <- ifelse(any(duplicated(colData(reacValues$proteinData)$groups )), FALSE, TRUE)
    if(isPilot) plot_type <- "data_points"
    
    p <- 0
    if (plot_type == "error_bars") {
      stats <- Rmisc::summarySE(object, measurevar="value", groupvars=c("Protein.IDs","group"), na.rm = T)
      stats <- stats[!is.na(stats$value),]
      
      p <- ggplot(stats, aes(x = group, y = value, color = Protein.IDs)) +
        geom_point(size = input$profile_pointsize) + 
        geom_errorbar(aes(ymin = value - se, ymax = value + se), width = .2)
    } else if (plot_type == "violin") {
      
      p <- ggplot(object, aes(x = group, y = value, fill = Protein.IDs)) + 
        geom_violin(width = 0.5, fill = 'white', trim = T) + 
        stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                     geom = "crossbar", 
                     width = 0.25
                     )
    } else if (plot_type == "data_points") {
      p <- ggplot(object, aes(x = group, y = value, color = Protein.IDs)) + 
        geom_jitter(
        data = object,
        mapping = aes(x = group, y = value, color = Protein.IDs),
        width = 0.1,
        height = 0,
        size = input$profile_jittersize
      )
    }
  
    p <- p +
      xlab("") + ylab("Imputed Intensities (log2)") + 
      theme_cowplot(font_size = input$profile_base) + 
      background_grid() + 
      ggtitle(input$selectProfilePlotGene) + 
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = input$profile_base
      ),legend.title = element_blank())
    
    if (input$profile_add_points && plot_type != "data_points") {
      if (input$profile_plot_type == " violin") {
        p <- p + geom_jitter(
          width = 0.1,
          height = 0,
          shape=23,
          size = input$profile_jittersize
        )
      } else {
        p <- p + geom_jitter(
          data = object,
          mapping = aes(x = group, y = value, color = Protein.IDs),
          width = 0.1,
          height = 0,
          size = input$profile_jittersize
        )
      }
      
    }
    p <- p + scale_color_manual(values=myScatterColors()) +
      scale_fill_manual(values=myScatterColors())
    
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = input$profile_format,
        width = input$profile_width,
        height = input$profile_height,
        filename = paste0("profile_plot_", input$selectProfilePlotGene)
      )
    )
  })
  
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
  
  gprofilerOrganisms <- reactive({
    df <- read.delim('data/gprofiler_sources_2021_10.txt', header=T, sep='\t')
    df
  })
  
  observe({
    updateSelectizeInput(session,
                         server = T,
                         'gprofilerOrganism',
                         choices = gprofilerOrganisms()$id,
                         selected = "hsapiens"
    )
  })
  
  observeEvent(c(input$maVolcanoSubmit, input$submitMultiComp), {
    if (reacValues$show_ora == TRUE) {
      reacValues$show_ora = FALSE
      toggle(id = 'hide_ora_before_submit', anim = T)
    }
    reacValues$dataGprofiler <- NULL
    reacValues$GostPlot <- NULL
    
  })
  
  output$organismSources <- renderText({
    req(input$gprofilerOrganism)
    sources <- gprofilerOrganisms()[gprofilerOrganisms()$id == input$gprofilerOrganism, "sources"]
    paste0("Available sources: ", sources)
  })
  
  observeEvent(input$submitORA, {

    req(reacValues$dataComp)
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    organism <- isolate(input$gprofilerOrganism)
    sources <- isolate(input$oraSources)
    customBg <- isolate(input$oraCustom)
    excludeIea <- isolate(input$oraExcludeIea)

    validate(need(!is.null(organism), "Please select an organism."))
    validate(need(!is.null(sources), "Please select at least one source."))
    validate(need(nrow(reacValues$dataComp)>0, "No proteins to plot (output table must be empty)."))
    validate(need(length(ridx)>0, "No proteins to plot (output table must be empty)."))
    if ( nrow(reacValues$dataComp) < 2 ) return(NULL)

    queryGenes <- gsub(";.*", "", reacValues$dataComp[ridx, geneName])

    withProgress(message = 'Computing ORA ', {

      customBgList <- NULL
      if (customBg) {
        customBgList <- gsub(";.*", "", reacValues$filtData[[geneName]])
      }

      reacValues$dataGprofiler <-
        gost(
          query = queryGenes,
          sources = sources,
          evcodes = input$showGenes,
          organism = organism,
          custom_bg = customBgList,
          exclude_iea = excludeIea,
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
    req(reacValues$GostPlot)
    validate(
      need(nrow(reacValues$GostPlot$result) >= 1, 
           paste("No results to show\n",
           "Please make sure that the organism",
           "is correct or deselect Only show significant terms"))
    )
    # req(reacValues$GostPlot )
    gostplot(reacValues$GostPlot, capped = F, interactive = T)
    # %>% layout(yaxis = list(title = '-log10(p-adj)'))
  })

  oraBarBase <-
    eventReactive(c(input$submitORABar, reacValues$dataGprofiler), {
      req(reacValues$dataGprofiler)
      orasource <- isolate(input$orasource)
      plotDf <-
        reacValues$dataGprofiler[reacValues$dataGprofiler$source == orasource, ]
      print(plotDf)
      validate(need(
        !is.null(reacValues$dataGprofiler) && nrow(plotDf) > 0,
        paste(
          "No results to show for this source\n",
          "Please make sure that the organism",
          "is correct or deselect Only show significant terms"
        )
      ))
      
      oracolor <- isolate(input$oraBar_color)
      oraMaxTerm <- isolate(input$oraBar_maxTerms)
      baseSize <- isolate(input$oraBar_base)
      
      p <- plotGprofilerBarplot(
        plotDf = plotDf,
        orasource = orasource,
        oracolor = oracolor,
        oraMaxTerm = oraMaxTerm,
        baseSize = baseSize
      )
      
      ggplotly(p)
    })
  
  output$oraBarplot <- renderPlotly({

    oraBarBase() %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
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
      paste(
        "data_gprofiler2_",
        reacValues$enrichmentChoice,
        '_proteins_log2fcThresh_',
        reacValues$fcCutoff,
        '_signcutoff_',
        reacValues$sigCutoffValue,
        ".tsv",
        sep = ""
      )
    },
    content = function(file) {
      write.table(reacValues$dataGprofiler, file, row.names = F, sep = "\t", quote = F)
    }
  )
  
  output$download_button_heatmap <- renderUI({
    req(reacValues$dataHeatmap)
    downloadButton("download_heatmap_data", "Download heatmap data")
  })
  
  output$download_heatmap_data <- downloadHandler(
    filename = function() {
      "heatmap_data.tsv"
    },
    content = function(file) {
      df <-  reacValues$dataHeatmap
      df[[geneName]] <- row.names(df)
      df <- df[,c(ncol(df),1:ncol(df)-1)]

      write.table(df, file, row.names = F, quote = F, sep = "\t")
    }
  )
  
  output$download_button_ora <- renderUI({
    req(reacValues$dataGprofiler)
    downloadButton("download1", "Download results")
  })

  output$download_amica <- renderUI({
    req(reacValues$analysisSuccess)
    if (reacValues$amicaInput) return(NULL)
    downloadButton("download2", h4("Download amica output"))
  })
  
  output$download2 <- downloadHandler(
    filename = function() {
      "amica_protein_groups.tsv"
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

  ppi <- reactive({
    simplify(read_graph('data/intact_weighted.edgelist', format="ncol", directed=F))
  })
  
  cellmap <- reactive({
    cellmap <-
      read.table(
        "data/preys-latest.txt",
        header = T,
        sep = "\t",
        stringsAsFactors = F
      )
    cellmap <-
      cellmap[, c("symbol", "MMF.localization", "SAFE.localization")]
    
    cellmap$localization <- cellmap$MMF.localization

    mito <-
      c(
        "mitochondrial matrix",
        "mitochondrial inner membrane, mitochondrial intermembrane space",
        "mitochondrial outer membrane, peroxisome"
      )
    
    er <-
      c("ER lumen",
        "ER membrane",
        "nuclear outer membrane-ER membrane network")
    
    endoLyso <- c("early endosome, recycling endosome",
                  "endosome, lysosome")
    
    nucleus <-
      c(
        "nuclear body",
        "nucleolus",
        "nucleoplasm"
      )
    
    cytoskeleton <-
      c("actin cytoskeleton, cytosol",
        "microtubule cytoskeleton")
    
    plasma <- 
      c("cell junction",
        "plasma membrane")
    
    cellmap[cellmap$localization %in% mito, 'localization'] <- "mitochondrion"
    cellmap[cellmap$localization %in% nucleus, 'localization'] <- "nucleus"
    cellmap[cellmap$localization %in% endoLyso, 'localization'] <-
      "endosome, lysosome"
    cellmap[cellmap$localization %in% er, 'localization'] <- "ER"
    cellmap[cellmap$localization %in% cytoskeleton, 'localization'] <-
      "cytoskeleton"
    cellmap[cellmap$localization %in% plasma, 'localization'] <-
      "plasma membrane, cell junction"
    colnames(cellmap) <-
      c("label",
        "CellMap NMF localization",
        "CellMap SAFE localization",
        "Subcell_localization")
    cellmap
  })

  multiNetwork <- function(thresh = 0) {
    if (is.na(thresh) || is.null(thresh))
      thresh <- 0
    validate(need(nrow(reacValues$dataComp) > 1, "Not enough data to display."))
    
    ridx <- input$groupComparisonsDT_rows_all
    rnames <- rownames(reacValues$dataComp[ridx, ])
    validate(need(length(rnames) > 1, "Not enough data to display."))
    
    nw.data <-
      getBait2PreyNetwork(
        reacValues$dataComp[rnames, ],
        enrichedMatrixSet(),
        reacValues$selection,
        myScatterColors()
      )
    
    plotMultiNetwork(nw.data,
                     ppi(),
                     cellmap(),
                     myScatterColors(),
                     thresh,
                     input$enableNetworkZoom)
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
  
  singleNetwork <- function(thresh) {
    if (is.na(thresh) || is.null(thresh))
      thresh <- 0
    
    req(reacValues$dataComp)
    validate(need(nrow(reacValues$dataComp) > 1, "Not enough data to display."))
    
    ridx <- input$groupComparisonsDT_rows_all
    rnames <- rownames(reacValues$dataComp[ridx, ])
    validate(need(length(rnames) > 1, "Not enough data to display."))
    
    networkData <-
      toNetworkData(reacValues$dataComp[rnames, ], ppi(), cellmap())
    
    msg <- ""
    if (length(rnames) > 1) {
      msg <- "\nCould not map selected proteins to the human PPI network."
    }
    
    validate(need(
      !is.null(networkData),
      paste0("#### There are no proteins to display with your selection.", msg)
    ))
    
    plotSingleNetwork(
      networkData = networkData,
      thresh = thresh,
      enableNetworkZoom = input$enableNetworkZoom
    )
  }
  
  
  ppiNet <- function() {
  # output$network <- renderVisNetwork({
    out <- 0
    if (length(grep(logfcPrefix, names(reacValues$dataComp)) ) > 1) {
      out <- multiNetwork(input$edgeWeightThresh)
    } else {
      out <- singleNetwork(input$edgeWeightThresh)
    }
    out
  }
  
  output$network <- renderVisNetwork({
    ppiNet()
  })
  
  output$networkDT <- renderDT({
    req(reacValues$nwNodes)
    req(ppiNet())
    
    datatable(
      reacValues$nwNodes[, grep("id|color|shape|key", names(reacValues$nwNodes), invert = T)],
      extensions = 'Buttons',
      filter = "top",
      rownames = F,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        dom = 'Bfrtip',
        buttons = c('csv')
      )
    )
  }, server = F)
  
  output$download_network_button <- renderUI({
    downloadButton("downloadNetworkHtml", "Download html")
  })
  
  output$downloadNetworkHtml <- downloadHandler(
    filename = function() {
      paste('network-', Sys.Date(), '.html', sep='')
    },
    content = function(con) {
      ppiNet() %>% visSave(con)
    }
  )
  
  ### COMPARE amica EXPERIMENTS LOGIC AND PLOTS
  observeEvent(input$submitAmicaComparisons, {
    
    suffix1 <- isolate(input$suffix1)
    suffix2 <- isolate(input$suffix2)
    
    tmpSuffixes <- make.names(c(suffix1, suffix2), unique = T)
    suffix1 <- tmpSuffixes[1]
    suffix2 <- tmpSuffixes[2]
    
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
    
    if ("quantified" %in% names(amicaData2)) {
      invisible()
    } else {
      rows <- which(complete.cases( amicaData2[, grep("__vs__|^ImputedIntensity",
                                                    colnames(amicaData2))  ] ))
      amicaData2$quantified <- ""
      amicaData2[rows, 'quantified'] <- "+"
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

  output$download_merged_amica <- renderUI({
    req(reacValues$combinedData)
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
  
  get_scatter_amicas_data <- reactive({
    event_data("plotly_selected", source = "subset")
  })

  output$scatterPlotsAmica <- renderPlotly({
    input$submitScatterAmicas
    req(reacValues$combinedData)
    assayNameAmicas <- isolate(input$assayNamesAmicas)
    fontsize <- isolate(input$scatteramica_base)
    pformat <- isolate(input$scatteramica_format)
    pwidth <- isolate(input$scatteramica_width)
    pheight <- isolate(input$scatteramica_height)
    showLine <- isolate(input$scatteramica_showLine)
    
    sampleSelection <- isolate(input$selectScatterSamplesAmica)
    validate(need(length(sampleSelection)==2, ""))
    colors <- isolate(myScatterColors())
    
    clrs <- colors[1]
    xvar <- isolate(sampleSelection[1])
    yvar <- isolate(sampleSelection[2])
    
    sampleX <- grep(xvar, grep(paste0("^", assayNameAmicas), names(reacValues$combinedData), value = T ), value = T)
    sampleY <- grep(yvar, grep(paste0("^", assayNameAmicas), names(reacValues$combinedData), value = T ), value = T)
    
    plotData <- reacValues$combinedData[, c("Gene.names", sampleX, sampleY)]
    plotData$key <- plotData$Gene.names
    
    plotData$show_id <- FALSE
    if (!is.null(get_scatter_amicas_data() )) {
      plotData[plotData$key %in% get_scatter_amicas_data()$key, "show_id"] <- TRUE
    }
    
    formula <- as.formula(paste0(sampleY, ' ~ ', sampleX))
    fit1 <- lm(formula, data=plotData)
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
      ggplot(plotData,
             aes(
               x = !!sym(paste0(sampleX)),
               y = !!sym(paste0(sampleY)),
               label = Gene.names,
               key = key
             )) + geom_point(color = clrs) +
      labs(
        x = paste0(assayNameAmicas, ' ', xvar, ' (log2)'),
        y = paste0(assayNameAmicas, ' ', yvar, ' (log2)')
      ) +
      scale_color_manual(values = colors) + theme_cowplot(font_size = fontsize) +
      background_grid()
    
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
      aes(!!sym(paste0(sampleX)),
          !!sym(paste0(sampleY)),
          label = Gene.names)
      ,position = position_jitter(width=0.25,height=0.25)
    )
    ggplotly(pu, source = "subset") %>% 
      layout(dragmode = "select") %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = pformat,
        width = pwidth,
        height = pheight,
        filename = "scatterplot_experiments"
      )
    )
  })
  
  output$corrAmicasSamplesInput <- renderUI({
    req(input$assayNamesAmicas)
    req(reacValues$combinedData)
    
    assayNameAmicas <- isolate(input$assayNamesAmicas)
    samples <- grep(paste0("^",assayNameAmicas),
                    names(reacValues$combinedData),
                    value = T
    )
    samples <- gsub(paste0("^", assayNameAmicas, "."), "", samples)
    
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
    req(reacValues$combinedData)
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
      groupInputs <- gsub(paste0("^", assayName, "."), "", samples)
    }

    df <- reacValues$combinedData[, grep(paste0("^", assayName), names(reacValues$combinedData))]
    df <- df[, grep(paste0(groupInputs, collapse = "|"), names(df))]
    names(df) <- gsub(paste0(assayName, "."), "", names(df))
    corDf <- cor(df, method = "pearson", use = "complete.obs")
    limits <- c(min(corDf) - 0.003, 1)
    diag(corDf) <- NA

    heatmaply_cor(
        round(corDf, 3),
        xlab = "", 
        ylab = "",
        limits = limits,
        colors = heatColors(),
        plot_method = "plotly",
        key.title = "Pearson Correlation"
      )
  })
  
  output$corrAmicasPlotly <- renderPlotly({
    req(reacValues$combinedData)
    withProgress(message = "Plotting correlation plot ", {
      p <- corrBaseAmicasPlot()
    })
    p %>%  config(displaylogo = F,
                  modeBarButtonsToRemove = removePlotlyBars,
                  toImageButtonOptions = list(format = input$corAmicas_format,
                                              width = input$corAmicas_width,
                                              height = input$corAmicas_height,
                                              filename = "corrplot")
    )
  })
  
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
  
  output$designTitle <- renderText({
    req(reacValues$uploadSuccess)
    req(reacValues$proteinData)
    "Experimental Design"
    })
  
  output$isPilot <- reactive({
    req(reacValues$proteinData)
    any(duplicated(colData(reacValues$proteinData)$groups))
  })
  outputOptions(output, "isPilot", suspendWhenHidden = FALSE)
  
  observeEvent(input$timeOut, { 
    print(paste0("Session (", session$token, ") timed out at: ", Sys.time()))
    showModal(modalDialog(
      title = "Timeout",
      paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
      footer = NULL
    ))
    session$close()
  })
  
  output$qcReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "qcReport.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "qcReport.Rmd")
      file.copy("reports/qcReport.Rmd", tempReport, overwrite = TRUE)
      
      withProgress(message = "Generating QC report", value = 0, {
      # Set up parameters to pass to Rmd document

      nsummary <- ifelse(
        contaminantCol %in% names(rowData(reacValues$proteinData)),
        nrow(rowData(reacValues$proteinData)[rowData(reacValues$proteinData)[[contaminantCol]] != filterVal, ]),
        nrow(rowData(reacValues$proteinData))
      )
      ngroups <- length(unique(reacValues$expDesign$groups))
      ncomps <- nrow(reacValues$contrastMatrix )
      nquant <- length(which(rowData(reacValues$proteinData)$quantified=="+"))
      
      # prepare output plots
      
      lfqPCA <- NA
      lfqDensityPlot <- NA
      lfqCvBox <- NA
      lfqCorrPlot <- NA
      
      isPilot <-
        ifelse(any(duplicated(colData(
          reacValues$proteinData
        )$groups)), FALSE, TRUE)
      
      lfqAvail <- TRUE
      if (
        "Intensity" %in% assayNames(reacValues$proteinData) ||
        "LFQIntensity" %in% assayNames(reacValues$proteinData)
      ) {
        aname <-
          ifelse(
            "LFQIntensity" %in% assayNames(reacValues$proteinData),
            'LFQIntensity',
            'Intensity'
          )
        
        lfqData <- getAssayData(reacValues$proteinData,
                                assayNames = aname
        )
        
        # PCA
        lfqOlist <-
          computePCA(
            assay(reacValues$proteinData, aname),
            colData(reacValues$proteinData),
            reacValues$groupFactors
          )
        
        lfqPCADf <- lfqOlist[[1]]
        lfqPct <- lfqOlist[[2]]
        
        lfqPCA <- plotPCA(
          lfqPCADf,
          lfqPct,
          myGroupColors(),
          aname
        )
        
        # box
        lfqBoxplot <- plotBoxPlot(lfqData,
                                  aname,
                                  myGroupColors()
        )
        
        # density
        lfqDensityPlot <- plotDensityPlot(
          lfqData,
          aname,
          myColors()
        )
        
        if (!isPilot) {
          lfqCvBox <- plotCoeffVarPlot(lfqData,
                                       myGroupColors(),
                                       aname)
        }
        
        # corr plot
        lfqCorrPlot <- plotCorrPlotly(
          assay(reacValues$proteinData, aname),
          colData(reacValues$proteinData),
          myGroupColors(),
          heatColors(),
          aname,
          annotSamples=T,
          reacValues$groupFactors
        )
        
      } else {
        lfqAvail <- FALSE
        }

      incProgress(3/10)
      impData <- getAssayData(reacValues$proteinData,
                              assayNames = "ImputedIntensity"
      )

      # PCA
      impOlist <-
        computePCA(
          assay(reacValues$proteinData, "ImputedIntensity"),
          colData(reacValues$proteinData),
          reacValues$groupFactors
        )
      
      impPCADf <- impOlist[[1]]
      impPct <- impOlist[[2]]
      
      impPCA <- plotPCA(
        impPCADf,
        impPct,
        myGroupColors(),
        "ImputedIntensity"
      )
      
      # intensity boxplot
      incProgress(4/10)
      impBoxplot <- plotBoxPlot(impData,
                                "ImputedIntensity",
                                myGroupColors()
      )
      incProgress(5/10)
      # density plots
      incProgress(6/10)
      impDensityPlot <- plotDensityPlot(
        impData,
        "ImputedIntensity",
        myColors()
      )
      incProgress(7/10)
      
      # CV-boxplots
      
      impCvBox <- NA
      if (!isPilot) {
        impCvBox <- plotCoeffVarPlot(impData,
                                     myGroupColors(),
                                     "ImputedIntensity")
      }
      
      
      incProgress(8 / 10)
      
      # corrplot
      impCorrPlot <- plotCorrPlotly(
        assay(reacValues$proteinData, "ImputedIntensity"),
        colData(reacValues$proteinData),
        myGroupColors(),
        heatColors(),
        "ImputedIntensity",
        annotSamples=T,
        reacValues$groupFactors
      )
      incProgress(9/10)
      
      ibaqAvail <- ifelse(
        abundancePrefix %in% assayNames(reacValues$proteinData),
        TRUE,
        FALSE
      )
      
      params <- list(isPilot=isPilot,
                     numberOfProteins=nsummary,
                     quantifiedProteins=nquant,
                     numberOfComparisons=ncomps,
                     numberOfGroups=ngroups,
                     lfqAvail=lfqAvail,
                     ibaqAvail=ibaqAvail,
                     lfqPCA=lfqPCA,
                     impPCA=impPCA,
                     lfqBoxplot=lfqBoxplot,
                     impBoxplot=impBoxplot,
                     lfqDensityPlot=lfqDensityPlot,
                     impDensityPlot=impDensityPlot,
                     lfqCvBox=lfqCvBox,
                     impCvBox=impCvBox,
                     lfqCorrPlot=lfqCorrPlot,
                     impCorrPlot=impCorrPlot
                     )
      
      amicaInput <- TRUE
      if (!is.null(reacValues$inputParameterSummary)) {
        params[['analysisParams']] <- analysisParams()
        amicaInput <- FALSE
        params[['daTool']] <- reacValues$daTool
      }
      params[['amicaInput']] <- amicaInput
      
      if (lfqAvail) {
        params[['numid']] <- numIdPlotly()
        params[['missingvals']] <- pctMvsPlotly()
      }
      if (ibaqAvail) {
        params[['contaminants']] <- contaminantsPlotly()
      }

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
      })
    }
  )
  
  output$diffAbundanceReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "diffAbudanceReport.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tmpDir <- tempdir()
      tempReport <- file.path(tmpDir, "diffAbundanceReport.Rmd")
      tempHelpers <- file.path(tmpDir, "diffAbundancePlots.R")
      file.copy("reports/diffAbundanceReport.Rmd", tempReport, overwrite = TRUE)
      file.copy("R/server/diffAbundancePlots.R", tempHelpers, overwrite = TRUE)

      withProgress(message = "Generating Diff. abundance report", value = 0, {
        # Set up parameters to pass to Rmd document
        
        isPilot <-
          ifelse(any(duplicated(colData(
            reacValues$proteinData
          )$groups)), FALSE, TRUE)
        
        enrichedUpSet <- NA
        nMultiEnriched <- 0
        nMultiReduced <- 0
        enrichedUpSet <- NA
        reducedUpSet <- NA
        
        if (length(reacValues$reacConditions) > 1) {
          if (input$enrichmentChoice == "enriched" ||
              input$enrichmentChoice == "absolute") {
            enrMatrix <- generateEnrichedMatrix(
              reacValues$dataLimma,
              "enriched",
              input$sigCutoffValue,
              input$fcCutoff,
              input$pvalCutoff
            )
            enrMatrix <- as.data.frame(enrMatrix)
            if (nrow(enrMatrix>0))
            enrichedUpSet <- upset(enrMatrix, order.by = "freq")
            nMultiEnriched <- nrow(enrMatrix)
          }
          
          if (input$enrichmentChoice == "reduced" ||
              input$enrichmentChoice == "absolute") {
            redMatrix <- generateEnrichedMatrix(
              reacValues$dataLimma,
              "reduced",
              input$sigCutoffValue,
              input$fcCutoff,
              input$pvalCutoff
            )
            redMatrix <- as.data.frame(redMatrix)
            if (nrow(redMatrix>0))
            reducedUpSet <- upset(redMatrix, order.by = "freq")
            nMultiReduced <- nrow(redMatrix)
          }
        }
        
        incProgress(5/10)
        params <- list(isPilot=isPilot,
                       fcCutoff=input$fcCutoff,
                       sigCutoffValue=input$sigCutoffValue,
                       pvalCutoff=input$pvalCutoff,
                       enrichmentChoice=input$enrichmentChoice,
                       comparisons=reacValues$reacConditions,
                       dataLimma=reacValues$dataLimma,
                       myScatterColors=myScatterColors(),
                       enrichedUpSet=enrichedUpSet,
                       reducedUpSet=reducedUpSet,
                       nMultiEnriched=nMultiEnriched,
                       nMultiReduced=nMultiReduced
                       )
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      })
    }
  )
  
}