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
      geom_point() + theme_minimal(base_size = 14) +  scale_color_manual(values=myScatterColors() )
    
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
  
  output$inputParameterSummary <- renderText({
    req(reacValues$analysisSuccess)
    req(reacValues$inputParameterSummary)
    reacValues$inputParameterSummary
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
    
    object <- 0
    if (input$assayNames == "ImputedIntensity") {
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
  
  #pcaPlotly <- eventReactive(input$submitPCA, {
  output$pca <- renderPlotly({
    input$submitPCA
    assayNames <- isolate(input$assayNames)
    groupInputs <- isolate(input$pcaSamplesInput)
    pca_legend <- isolate(input$pca_legend)
    pca_base <- isolate(input$pca_base)
    pca_pointsize <- isolate(input$pca_pointsize)
    pca_show_label <- isolate(input$pca_show_label)
    
    validate(need(assayNames != "", "Please press submit or provide an intensity prefix."))
    validate(
      need(
        length(groupInputs) > 1 |
          is.null(groupInputs),
        "Please select at least two groups. If none is selected all are considered."
      )
    )
    
    if (is.null(groupInputs))
      groupInputs <- unique(colData(reacValues$proteinData)$groups)

    tmpCols <- colData(reacValues$proteinData)
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(tmpCols$groups))) {
      tmpCols$groups <-
        factor(tmpCols$groups, levels = c(reacValues$groupFactors, 
                                          setdiff(groupInputs, reacValues$groupFactors)))
    }
    
    plotData <- 0
    
    if (assayNames == "ImputedIntensity") {
      plotData <-
        assay(reacValues$proteinData, "ImputedIntensity")[isQuantRnames(reacValues$proteinData),
                                                          tmpCols$samples[tmpCols$groups %in% groupInputs]]
    } else {
      plotData <- assay(reacValues$proteinData, assayNames)
      
      plotData <-
        plotData[, tmpCols$samples[tmpCols$groups %in% groupInputs]]
      plotData <- plotData[complete.cases(plotData), ]
    }
    
    withProgress(message = "Plotting PCA ", {

      # idxs <-
      #   grep(paste0(colnames(plotData), collapse = "|"), tmpCols$samples)
      # groups <- tmpCols$groups[idxs]
      
      pca <- prcomp(as.data.frame(t(plotData)))
      df_out <- as.data.frame(pca$x)
      
      df_out$group <- tmpCols$groups[tmpCols$samples %in% row.names(df_out)]
      df_out$sample <- rownames(df_out)
      df_out$key <- df_out$sample
      
      df_out$show_id <- FALSE
      if (!is.null(get_pca_data() )) {
        df_out[df_out$key %in% get_pca_data()$key, "show_id"] <- TRUE
      }

      tmp <- summary(pca)
      percentage <- round(100 * tmp$importance['Proportion of Variance', 1:2], 2)
      percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
      
      p <- ggplot(df_out,aes(x=PC1,y=PC2,color=group, label=sample,  key=key ))
      p <- p + geom_point(size=pca_pointsize) + xlab(percentage[1]) + ylab(percentage[2])
      p <- p + scale_color_manual(values=myGroupColors()) + theme_minimal(base_size = pca_base) +
          theme(legend.text = element_text(size = pca_legend))
      
      # if (pca_show_label) {
      #   p <- p + geom_text(data=df_out,
      #                      aes(x=PC1,y=PC2, label=sample, key=key),
      #                      position=position_jitter(width=0.25,height=0.25)
      #                      )
      # }
      # 
      p <- p + geom_text(
        data = subset(df_out, show_id),
        aes(PC1, PC2, label = sample),
        position=position_jitter(width=0.25,height=0.5)
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
              filename = "pca"
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
  
  
  # output$pca <- renderPlotly({
  #   req(input$assayNames)
  #   req(reacValues$uploadSuccess)
  #   pcaPlotly()  %>%  config(
  #     displaylogo = F,
  #     modeBarButtonsToRemove = removePlotlyBars,
  #     toImageButtonOptions = list(
  #       format = input$pca_format,
  #       width = input$pca_width,
  #       height = input$pca_height,
  #       #width = 768,
  #       #height = 676,
  #       filename = "pca"
  #     )
  #   )
  # })
  
  ### ------------------------------------------- BOX PLOT
  
  boxplotPlotly <- eventReactive(input$submitBoxplot, {
    p <-
      ggplot(dataAssay()[!is.na(dataAssay()$value), ], aes(x = colname, y = value, fill =
                                                             group)) +
      geom_boxplot(
        outlier.shape = NA,
        outlier.color = NULL,
        outlier.fill = NULL
      ) +
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
          modeBarButtonsToRemove = removePlotlyBars,
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
  
  densityPlotly <- eventReactive(input$submitDensity, {
    p <-
      ggplot(dataAssay()[!is.na(dataAssay()$value), ], aes(x = value, color = colname)) +
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
      ) + xlab("Intensity (log2)") + ylab("Density")
    ggplotly(p)
  })
  
  output$densityPlot <- renderPlotly({
    req(reacValues$uploadSuccess)
    
    withProgress(message = "Plotting density plot", {
      densityPlotly() %>%
        config(
          displaylogo = F,
          modeBarButtonsToRemove = removePlotlyBars,
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
    
    xLabel <- paste0(assay1, "_", selection1)
    yLabel <- paste0(assay2, "_", selection2)
    plotData <- data.frame()
    
    if (assay1 == assay2) {
      validate(need(
        selection1 != selection2,
        "Cannot plot the same column on x - and y-axis."
      ))
      
      plotData <-
        assay(reacValues$proteinData, assay1)[, c(selection1, selection2)]
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
      plotData$Contaminant <- ifelse(rowData(reacValues$proteinData)[[contaminantCol]] ==
                                       "+",
                                     "yes",
                                     "no")
    }
    
    plotData$Gene <-
      rowData(reacValues$proteinData)[[geneName]]
    plotData$key <- plotData$Gene
    
    plotData$show_id <- FALSE
    if (!is.null(get_scatter_data() )) {
      plotData[plotData$key %in% get_scatter_data()$key, "show_id"] <- TRUE
    }
    
    fit1 <- lm(y ~ x, data = plotData)
    fit1.intercept <- fit1$coefficients[[1]]
    fit1.slope <- fit1$coefficients[[2]]
    
    title <-
      paste0(
        signif(fit1$coef[[1]], 5),
        "x + ",
        signif(fit1$coef[[2]], 5),
        ", r2 = ",
        signif(summary(fit1)$r.squared, 5)
      )
    
    withProgress(message = "Plotting scatter plot", {
      pu <-
        ggplot(plotData,
               aes(
                 x = x,
                 y = y,
                 label = Gene,
                 key = key,
                 color = Contaminant
               )) + theme_minimal(base_size = plot_fontsize) + labs(x = xLabel, y =
                                                                      yLabel)  +
        theme(
          legend.title = element_text(size = plot_legendsize),
          legend.text = element_text(size = plot_legendsize)
        ) +
        geom_point() + scale_color_manual(values = myScatterColors()) # + ggtitle(title) #scale_color_brewer(palette = "Paired")
      
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
    
    pu <- pu + geom_text(
      data = subset(plotData, show_id),
      aes(x,
          y,
          label = Gene)
      #hjust=0, vjust=0
      ,position = position_jitter(width=0.25,height=0.25)
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
  
  corrBasePlot <- eventReactive(input$submitCor, {
    assayName <- isolate(input$assayNames)
    annotSamples <- isolate(input$cor_annot)
    groupInputs <- isolate(input$corrSamplesInput)
    tmpCols <- colData(reacValues$proteinData)
    
    req(reacValues$uploadSuccess)
    
    validate(need(assayName != "", "Please provide an intensity."))
    
    validate(
      need(
        length(groupInputs) > 1 |
          is.null(groupInputs),
        "Please select at least two groups. If none is selected all are considered."
      )
    )
    
    if (is.null(groupInputs))
      groupInputs <- unique(colData(reacValues$proteinData)$groups)
    
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(tmpCols$groups))) {
      tmpCols$groups <-
        factor(tmpCols$groups, levels = c(reacValues$groupFactors, 
                                          setdiff(groupInputs, reacValues$groupFactors)))
    }
    
    df <- assay(reacValues$proteinData, assayName)
    df <- df[, tmpCols$samples[tmpCols$groups %in% groupInputs]]
    
    annot <- colData(reacValues$proteinData)
    row.names(annot) <- annot$samples
    annot <- annot[names(df), ]
    annot$samples <- NULL
    
    corDf <- cor(df, method = "pearson", use = "complete.obs")
    limits <- c(min(corDf) - 0.003, 1)
    diag(corDf) <- NA
    
    if (annotSamples) {
      p <- heatmaply_cor(
        round(corDf, 3),
        xlab = "",
        ylab = "",
        limits = limits,
        colors = heatColors(),
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
        colors = heatColors(),
        limits = limits,
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
    p %>%  config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = input$cor_format,
        width = input$cor_width,
        height = input$cor_height,
        filename = "corrplot"
      )
    )
  })
  
  
  ### OVERLAP HEATMAP
  overlapBasePlot <- eventReactive(input$submitOverlapHeatmap, {
    req(reacValues$uploadSuccess)
    tmpCols <- colData(reacValues$proteinData)
    
    annotSamples <- isolate(input$overlapHeatmap_annot)
    metric <- isolate(input$overlapHeatmap_metric)
    
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(tmpCols$groups))) {
      tmpCols <- tmpCols[tmpCols$groups %in% reacValues$groupFactors, ]
      tmpCols$groups <-
        factor(tmpCols$groups, levels = c(reacValues$groupFactors))
    }
    groupInputs <- unique(tmpCols$groups)
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    
    lfqs <- assay(reacValues$proteinData, aname)
    lfqs <- lfqs[, tmpCols$samples[tmpCols$groups %in% groupInputs]]
    
    annot <- tmpCols
    row.names(annot) <- annot$samples
    annot$samples <- NULL
    
    idProts <- list()
    for (col in names(lfqs)) {
      tmp <- row.names(lfqs)[!is.na(lfqs[, col])]
      idProts[[col]] <- tmp
    }
    
    reacValues$overlapDf <- calc_pairwise_overlaps(idProts)
    idProts <- NULL
    
    validate(need(ncol(lfqs)>2, "Need more than 2 samples."))
    
    limits <- c(0, 1)
    title <- "Overlap Coefficient"
    df <- head(reacValues$overlapDf)
    df1 <- head(reacValues$overlapDf)
    
    if (metric == 'num_shared') {
      df <- reacValues$overlapDf[, c("sample1", "sample2", "num_shared")]
      df1 <- data.frame(sample1 = df$sample2, sample2 = df$sample1, overlap = df$num_shared)
      limits <- c(min(df[,3], na.rm = T) - 0.05, max(df[,3], na.rm = T))
      title <- "#shared proteins"
    } else if (metric == 'overlap_coefficient') {
      df <- reacValues$overlapDf[, c("sample1", "sample2", "overlap")]
      df1 <- data.frame(sample1 = df$sample2, sample2 = df$sample1, overlap = df$overlap)
      limits <- c(min(df[,3], na.rm = T) - 0.05, 1)
    } else if (metric == 'jaccard_index') {
      df <- reacValues$overlapDf[, c("sample1", "sample2", "jaccard")]
      df1 <- data.frame(sample1 = df$sample2, sample2 = df$sample1, overlap = df$jaccard)
      limits <- c(min(df[,3], na.rm = T) - 0.05, 1)
      title <- "Jaccard Coefficient"
    }
    names(df) <- c("sample1", "sample2", "overlap")
    df <- rbind(df, df1)
    overlap_matrix <- acast(df, sample1 ~ sample2, value.var = "overlap")

    if (annotSamples) {
      p <- heatmaply_cor(
        overlap_matrix,
        xlab = "",
        ylab = "",
        limits = limits,
        row_side_palette = myGroupColors(),
        row_side_colors = annot,
        plot_method = "plotly",
        key.title = title,
        colors = heatColors()
      )
    } else {
      p <- heatmaply_cor(
        overlap_matrix,
        xlab = "",
        ylab = "",
        limits = limits,
        plot_method = "plotly",
        key.title = title,
        colors = heatColors()
      )
    }
    p
  })
  
  output$overlapHeatmapPlotly <- renderPlotly({
    withProgress(message = "Plotting overlap plot ", {
      p <- overlapBasePlot()
    })
    p %>%  config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = input$overlapHeatmap_format,
        width = input$overlapHeatmap_width,
        height = input$overlapHeatmap_height,
        filename = "overlap_heatmap"
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
        #dom = 'Bfrtip',
        #buttons = c('csv')
      )
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
      apply(tmp[rowData(reacValues$proteinData)[[contaminantCol]] == "+",], 2, sum, na.rm =
              T)
    contsInts <- contsInts / allInts
    
    contsInts <- reshape2::melt(contsInts)
    midx <-
      match(row.names(contsInts),
            colData(reacValues$proteinData)$samples)
    contsInts$group <- colData(reacValues$proteinData)$groups[midx]
    contsInts$Sample <- row.names(contsInts)
    contsInts$value <- 100 * contsInts$value
    
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(contsInts$group))) {
      contsInts <-
        contsInts[contsInts$group %in% reacValues$groupFactors, ]
      contsInts$group <-
        factor(contsInts$group, levels = reacValues$groupFactors)
    }
    
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
    req(reacValues$uploadSuccess)
    validate(need(
      abundancePrefix %in% assayNames(reacValues$proteinData),
      "No data to plot (no iBAQ columns available)."
    ))
    
    tmp <- object <- toLongFormat(
      2 ^ assay(reacValues$proteinData, abundancePrefix),
      reacValues$proteinData,
      addGroup = TRUE,
      addContaminant = TRUE,
      addGeneName = TRUE
    )
    
    tmp <- tmp[tmp$colname == input$samples, ]
    tmp <- tmp[order(tmp$value, decreasing = T), ]
    
    tmp$value <- 100 * tmp$value / sum(tmp$value, na.rm = T)
    
    withProgress(message = "Plotting most abundant proteins", {
      p <-
        ggplot(head(tmp, 15), aes(x = reorder(Gene.name,-value), y = value)) +
        geom_bar(stat = "identity", fill = myScatterColors()[1]) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 10
          ),
          legend.text = element_text(size = input$abundant_legend),
          legend.title = element_blank()
        ) + xlab("") + ylab("%Signal of protein")
    })
    
    ggplotly(p) %>% config(
      displaylogo = F,
      modeBarButtonsToRemove = removePlotlyBars,
      toImageButtonOptions = list(
        format = input$abundant_format,
        width = 676,
        height = 676,
        filename = "plot"
      )
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
    
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(df$group))) {
      df <- df[df$group %in% reacValues$groupFactors, ]
      df$group <- factor(df$group, levels = reacValues$groupFactors)
    }
    
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
    
    if (!is.null(reacValues$groupFactors) &&
        all(reacValues$groupFactors %in% unique(df$group))) {
      df <- df[df$group %in% reacValues$groupFactors, ]
      df$group <- factor(df$group, levels = reacValues$groupFactors)
    }
    
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
  
  cvsPlotly <- eventReactive(input$submitCVs, {
    calcData <- dataAssay()
    calcData$value <- 2 ^ calcData$value
    cvs <-
      aggregate(value ~ rowname + group, calcData, function(x)
        100 * (sd(x) / mean(x)))
    
    p <-
      ggplot(cvs[!is.na(cvs$value), ], aes(x = group, y = value, fill = group)) +
      geom_boxplot(
        outlier.shape = NA,
        outlier.color = NULL,
        outlier.fill = NULL
      ) +
      scale_fill_manual(values = myGroupColors()) +
      theme_minimal(base_size = input$cv_base) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = input$cv_base
        ),
        legend.text = element_text(size = input$cv_legend),
        legend.title = element_blank()
      ) + ylab("Coefficient of Variation (%)") + xlab("")
    
    ggplotly(p)
  })
  
  output$boxplotCV <- renderPlotly({
    req(reacValues$uploadSuccess)
    validate(need(
      any(duplicated(
        colData(reacValues$proteinData)$groups
      )),
      "Cannot output CV plot for a pilot without replicates"
    ))
    #validate(need(input$assayName != "", "Please provide an intensity."))
    withProgress(message = "Plotting Coefficient of Variations ", {
      cvsPlotly() %>% config(
        displaylogo = F,
        modeBarButtonsToRemove = removePlotlyBars,
        toImageButtonOptions = list(
          format = input$cv_format,
          width = input$cv_width,
          height = input$cv_height,
          filename = "cv_plot"
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
  
  output$submitDotplot <- renderUI({
    req(reacValues$dotplotGroupsDf)
    actionButton("submitDotplot", label = "Plot Dotplot")
  })
  
  observeEvent(input$submitDotplot,{
    
    req(reacValues$dotplotGroupsDf)
    req(reacValues$dataComp)
    
    validate(need(!any(duplicated(reacValues$dotplotGroupsDf$group)), "Error! Please provide 
                  unique groups"  )  )

    group2comps <- reacValues$dotplotGroupsDf
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    ridx <- rownames(reacValues$dataComp[ridx,])

    pattern <- padjPrefix
    if (reacValues$sigCutoffValue == "p-value") pattern <- pvalPrefix
    
    selection <- grep(pattern, names(reacValues$dataLimma), value = T)
    
    widestats <- reacValues$dataLimma[ridx, selection]

    longStats <- stats::reshape(widestats, idvar = "rowname",
                                ids = rownames(widestats), times = names(widestats),
                                timevar = "colname", varying = list(names(widestats)),
                                direction = "long", v.names = "value")

    longStats$colname <- gsub(pattern, "", longStats$colname)
    names(longStats) <- c('Comparison', 'padj', 'ProteinID')

    longStats <- longStats[longStats$Comparison %in% group2comps$comparison, ]
    midx <- match(longStats$Comparison, group2comps$comparison)
    longStats$Group <- group2comps$group[midx]

    longStats$Comparison <- NULL

    selection <- grep(logfcPrefix, names(reacValues$dataLimma), value = T)
    widefcs <- reacValues$dataLimma[ridx, selection]

    longFCs <- stats::reshape(widefcs, idvar = "rowname",
                              ids = rownames(widefcs), times = names(widefcs),
                              timevar = "colname", varying = list(names(widefcs)),
                              direction = "long", v.names = "value")

    longFCs$colname <- gsub("logFC_", "", longFCs$colname)
    names(longFCs) <- c('Comparison', 'log2FC', 'ProteinID')

    longFCs <- longFCs[longFCs$Comparison %in% group2comps$comparison, ]
    midx <- match(longFCs$Comparison, group2comps$comparison)
    longFCs$Group <- group2comps$group[midx]

    longData <- merge(longFCs, longStats, by = c('ProteinID', 'Group'))
    
    aname <-
      ifelse(
        "LFQIntensity" %in% assayNames(reacValues$proteinData),
        'LFQIntensity',
        'Intensity'
      )
    df <- assay(reacValues$proteinData, aname)
    df <- df[ridx,]
    
    intData <- stats::reshape(df, idvar = "rowname",
                               ids = rownames(df), times = names(df),
                               timevar = "colname", varying = list(names(df)),
                               direction = "long", v.names = "value")
    midx <- match(intData$colname, reacValues$expDesign$samples)
    intData$group <- reacValues$expDesign$groups[midx]
    
    intData <- intData[intData$group %in% group2comps$group,]
    intData <-
      aggregate(value ~ rowname + group, intData, mean)
    names(intData) <- c('ProteinID', 'Group', 'AvgIntensity')
    
    longData <- merge(longData, intData, by = c('ProteinID', 'Group'))
    
    longData$Gene <- reacValues$filtData[longData$ProteinID, "Gene.names"]
    longData$Gene <- gsub("^([^;]*;[^;]*);.*", "\\1", longData$Gene)

    if (reacValues$show_dotplot == FALSE) {
      reacValues$show_dotplot = TRUE
      toggle(id = 'hide_dotplot_before_submit', anim = T)
    }
    longData$significant <- factor(ifelse(longData$padj <= 0.05, 1.5, 0))
    reacValues$dataDotplot <- longData
  })
  
  output$dotplot_color_gradient <- renderUI({
    req(reacValues$dataDotplot)
    minVal <- round(min(reacValues$dataDotplot$log2FC), 3)
    maxVal <- round(max(reacValues$dataDotplot$log2FC), 3)
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
          2,
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
    
    options <- c("AvgIntensity", "log2FC")

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
    max(600, 22 * length(unique(reacValues$dataDotplot$Gene)))
  } )
  
  numberOfDotplotBaits <- reactive({
    req(reacValues$dataDotplot)
    maxNamelength <- max(nchar(as.character(reacValues$dataDotplot$Gene)))
    numCols <- length(unique(reacValues$dataDotplot$Group))
    value <- 60 * numCols + 60
    if (numCols < 3) {
      if (maxNamelength >= 60) {
        value <- value + 500
      } else if (maxNamelength >= 20) {
        value <- value + max(300, 6*maxNamelength)
      }
    }
    max(400, value)
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
    req(reacValues$dataDotplot)
    
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
    
    mat <- reacValues$dataDotplot %>%
      select(Gene, Group, all_of(clusteringMetric)) %>%  # drop unused columns to faciliate widening
      pivot_wider(names_from = Group, values_from = all_of(clusteringMetric)) %>%
      data.frame() # make df as tibbles -> matrix annoying
    
    row.names(mat) <- mat$Gene  # put gene in `row`
    mat$Gene <- NULL #drop gene column as now in rows
    
    tmat <- mat %>% as.matrix()
    tmat[is.na(tmat)] <- 0
    
    dist_meth <- ifelse(!is.null(input$dotplot_distance_metric),
                        input$dotplot_distance_metric, 
                        "canberra")
    clst_meth <- ifelse(!is.null(input$dotplot_clustering_method),
                        input$dotplot_clustering_method,
                        "complete")
    
    #rddr <- reorder(as.dendrogram(rclust), rowMeans(mat))
    #cddr <- reorder(as.dendrogram(cclust), colMeans(mat))
    
    if (!is.null(input$dotplot_cluster_columns) && input$dotplot_cluster_columns) { # order by input
      cclust <- hclust(dist(t(tmat), method = dist_meth),
                       method = clst_meth) # hclust with distance matrix
      reacValues$dataDotplot$Group <-
        factor(reacValues$dataDotplot$Group, levels = names(mat)[as.hclust(cclust)$order])
    } else { # order by clustering
      reacValues$dataDotplot$Group <-
        factor(reacValues$dataDotplot$Group, levels = reacValues$dotplotGroupsDf$group)
    }

    rclust <- hclust(dist(tmat, 
                          method = dist_meth), 
                     method = clst_meth) # hclust with distance matrix
    reacValues$dataDotplot$Gene <-
      factor(reacValues$dataDotplot$Gene, levels = rownames(mat)[as.hclust(rclust)$order])

    
    signficantTitle <- "adj.p-value"
    if (reacValues$sigCutoffValue == "p-value")
      signficantTitle <- "p-value"
    
    filterValues <- TRUE
    if (clusteringMetric == "log2FC" &&
        !is.null(input$dotplot_ctrl_substraction)) {
      filterValues <- input$dotplot_ctrl_substraction
    }

    preFiltered <- reacValues$dataDotplot
    if (filterValues) {
      preFiltered <- reacValues$dataDotplot %>% filter(!!rlang::sym(clusteringMetric) > 0)
    }
    
    x <- preFiltered[[clusteringMetric]]
    
    if (clusteringMetric == "log2FC") {
      preFiltered$relativeAbundance <- (x-min(x,na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
    } else {
      preFiltered$relativeAbundance <- x/max(x, na.rm = T)
    }
    
    p <- preFiltered %>%
      ggplot(aes(
        x = Group,
        y = Gene,
        fill = log2FC,
        color = significant,
        size = relativeAbundance,
        stroke = 1
      )) +
      geom_point(shape = 21) +
      cowplot::theme_cowplot() +
      theme(axis.line  = element_blank()) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      ylab('') + xlab('') +
      theme(axis.ticks = element_blank()) +
      scale_size_continuous(
        range = c(minSizeGradient, maxSizeGradient),
        name = 'Relative Abundance',
        breaks = c(
          min(preFiltered$relativeAbundance ) + 0.01,
          max(preFiltered$relativeAbundance )
        ),
        labels = c('Less', 'More')
        ) +
      scale_fill_gradientn(
        # colours = heatColors(),
        colours = dotplotColors,
        limits = c(
          ifelse(
            minColorGradient < min(preFiltered$log2FC),
            min(preFiltered$log2FC),
            minColorGradient
          ),
          maxColorGradient
        ),
        oob = scales::squish,
        name = 'log2FC',
        guide = guide_colorbar(order = 1)
      ) +
      scale_color_manual(
        signficantTitle,
        values = c('skyblue', 'black'),
        limits = c('0', '1.5'),
        labels = c('> 0.05', '<= 0.05')
        )
    p
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
                       width = numberOfDotplotBaits() * 1/80, 
                       height = numberOfDotplotPoints() * 1/80,
                       device = cairo_pdf)
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
    
    if ("AvgSpec" %in% names(reacValues$dataDotplot))
      text <- paste0(text, 
                     "<b> AvgSpec: </b>", round(point$AvgSpec, 3), "<br/>")
    
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
        buttons = c('csv')
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
    lty <- ifelse(bool, 1, 0)
    
    validate(need(!is.null(comparisons), "Please select at least two comparisons."))
    validate(need(length(comparisons) > 1, 
                  "Need at least two comparisons to render UpSet plot."))
    validate(need(length(comparisons) < 6, "Cannot output Euler plot for more than 5 sets."))
    
    if (!is.null(reacValues$newMultiNames) && all(reacValues$newMultiNames$new != "") &&
        length(reacValues$newMultiNames$new) == length(comparisons)) {
      for (idx in seq_along(reacValues$newMultiNames$old)) {
        names(binMat)[which(names(binMat) == reacValues$newMultiNames$old[idx])] <- reacValues$newMultiNames$new[idx]
        comparisons[idx] <- reacValues$newMultiNames$new[idx]
      }
    }

    cols <- myScatterColors()[1:length(comparisons)]
    if (!is.null(input$eulercol_1)) {
      tmp <- sapply(1:length(comparisons), function(i) {input[[paste0("eulercol_",i)]]})
      if (all(!is.null(tmp))) cols <- tmp
    }

    fit <- euler(binMat[, comparisons])
    comps <- grep("&", names(fit$original), invert = T, value = T)
    numComps <- length(comps)
    plot(fit, quantities=showQuant, 
         fills = list(fill = cols),
         legend = ifelse(showLegend, list(labels = comps), F), 
         lty = lty
    )
  }
  
  output$eulerrPlot <- renderPlot({
    input$submitMultiComp
    print(eulerData())
  })
  
  output$download_button_eulerr <- renderUI({
    downloadButton("download_eulerr", "Download Euler plot")
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

    if (!is.null(reacValues$newMultiNames) && all(reacValues$newMultiNames$new != "") && 
        length(reacValues$newMultiNames$new) == length(samples)) {
      for (idx in seq_along(reacValues$newMultiNames$old)) {
        names(matrixSet)[which(names(matrixSet) == reacValues$newMultiNames$old[idx])] <- reacValues$newMultiNames$new[idx]
        samples[idx] <- reacValues$newMultiNames$new[idx]
      }
    }
    
    upset(
      matrixSet,
      sets = samples,
      mb.ratio = mb.ratio,
      #set_size.numbers_size = 8,
      #set_size.show = T,
      order.by = upset_sorted,
      #cutoff = cutoff,
      text.scale = scale,
      # text.scale = c(2, 2,
      #                2, 2,
      #                2, 3),
      point.size = upset_pointsize
    )
  }
  
  output$upsetPlot <- renderPlot({
    input$submitMultiComp
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
    show_row_labels <- isolate(input$heatmap_row_labels)
    show_col_labels <- isolate(input$heatmap_col_labels)
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
          showticklabels = c(show_col_labels, show_row_labels),
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
          showticklabels = c(show_col_labels, show_row_labels),
          row_dend_left = TRUE,
          column_text_angle = 90,
          key.title = cbarTitle,
          colors = heatColors()
        )
      }
      
      calc_height <- min(15 * nrow(reacValues$dataHeatmap), 1200)
      if (calc_height < 600) calc_height <- 800
      
      p %>% layout(height=calc_height)
    })
  })
  
  output$compareHeatmap <- renderPlotly({

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
    
    levels <- c('both', selection[1], selection[2], 'none')
    if(!is.null(labels) & all(labels != "")) {
      levels[2] <- labels[1]
      levels[3] <- labels[2]
    }
    
    plotData$significant <-
      factor(plotData$significant,
             levels = levels)
    
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
        theme_minimal(base_size = fontsize) + theme(
          legend.text = element_text(size = legendFontSize),
          legend.title = element_text(size = legendFontSize)
        ) + scale_color_manual(values=colors) +
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
    sample <- isolate(input$maVolcanoSampleSelect)
    fcThresh <- isolate(input$fcCutoff)
    choice <- isolate(input$enrichmentChoice)
    sigCutoffValue <- isolate(input$sigCutoffValue)
    padjY <- isolate(input$volcano_padj_y)
    pvalThresh <- isolate(input$pvalCutoff)
    
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
      padjYBoolean,
      pvalThresh
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
    
    pal <- myScatterColors()[1:3]
    if (is.null(c1) | is.null(c2)) {
      c1 <- pal[1]
      c2 <- pal[2]
    }
    pcols <- c(c1,c2,"red")
    setChoice <- "union"
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
                   modeBarButtonsToRemove = removePlotlyBars,
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
                   modeBarButtonsToRemove = removePlotlyBars,
                   toImageButtonOptions = list(format = format,
                                               width = plot_width,
                                               height = plot_height,
                                               filename = "ma_plot")
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
      stats <- Rmisc::summarySE(object, measurevar="value", groupvars=c("Protein.IDs","group"))
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
      xlab("") + ylab("Intensity (log2)") + 
      theme_minimal(base_size = input$profile_base) + 
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
        filename = "profile_plot"
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
  
  output$organismSources <- renderText({
    req(input$gprofilerOrganism)
    sources <- gprofilerOrganisms()[gprofilerOrganisms()$id == input$gprofilerOrganism, "sources"]
    paste0("Available sources: ", sources)
  })
  
  observeEvent(input$submitORA,{
    req(reacValues$dataComp)
    ridx <- isolate(input$groupComparisonsDT_rows_all)
    organism <- isolate(input$gprofilerOrganism)
    sources <- isolate(input$oraSources)
    
    validate(need(!is.null(organism), "Please select an organism."))
    validate(need(!is.null(sources), "Please select at least one source."))
    validate(need(nrow(reacValues$dataComp)>0, "No proteins to plot (output table must be empty)."))
    validate(need(length(ridx)>0, "No proteins to plot (output table must be empty)."))
    if ( nrow(reacValues$dataComp) < 2 ) return(NULL)
    
    queryGenes <- gsub(";.*", "", reacValues$dataComp[ridx, geneName])

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
    baseSize <- isolate(input$oraBar_base)
    
    oraMaxTerm <- ifelse(oraMaxTerm == 0, 100000, oraMaxTerm)
    
    plotDf <- reacValues$dataGprofiler[reacValues$dataGprofiler$source==orasource,]
    
    validate(need(nrow(plotDf)>0, "No significant over-represented terms detected."))
    
    plotDf$minusLog10p_value <- -log10(plotDf$p_value)
    plotDf <- plotDf[!duplicated(plotDf$term_name),]
    plotDf <- plotDf[order(plotDf$minusLog10p_value),]
    plotDf <- tail(plotDf, min(nrow(plotDf), oraMaxTerm))
    plotDf$term_name <- factor(plotDf$term_name, levels = plotDf$term_name)
    
    p <- ggplot(plotDf, aes(x = term_name, y = minusLog10p_value)) +
      geom_bar(stat = "identity", fill=oracolor, width = 0.9)  + xlab("") +
      ylab("-log10(p-value)") +
      theme_minimal( base_size = baseSize) + coord_flip()
    
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
    
    cellmap[cellmap$localization %in% mito, 'localization'] <- "mitchondrion"
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
        "Subcell. localization")
    cellmap
  })

  multiNetwork <- function(thresh=0) {
    if (is.na(thresh) || is.null(thresh)) thresh <- 0
    validate(need(nrow(reacValues$dataComp) > 1,"Not enough data to display."))
    
    ridx <- input$groupComparisonsDT_rows_all
    rnames <- rownames(reacValues$dataComp[ridx,])
    validate(need(length(rnames) > 1,"Not enough data to display."))

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
               type = "full",
               height = "1200px",
               width = "2000px") %>%
      visIgraphLayout() %>% visPhysics(stabilization = F) %>%
      visOptions(selectedBy = "Subcell. localization",
                 highlightNearest = TRUE,
                 nodesIdSelection = TRUE) %>%
      visNodes(font = list(size = 20, strokeWidth = 2)) %>%
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
    validate(need(nrow(reacValues$dataComp) > 1,"Not enough data to display."))
    
    ridx <- input$groupComparisonsDT_rows_all
    rnames <- rownames(reacValues$dataComp[ridx,])
    validate(need(length(rnames) > 1,"Not enough data to display."))
    
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
        shape = "circle",
        font = list(size = 16, color = "white")
      )

    visNetwork(networkData$nodes,
               networkData$edges,
               layout = "layout_nicely",
               type = "full",
               height = "1200px",
               width = "2000px") %>%
      visIgraphLayout() %>% visPhysics(stabilization = F) %>%
      visOptions(selectedBy = "Subcell. localization",
                 highlightNearest = TRUE,
                 nodesIdSelection = TRUE) %>%
      visNodes(font = list(size = 20, strokeWidth = 2)) %>%
      visEdges(smooth=F, color=list(color = "grey", highlight = "red")) %>%

      visLegend(
        useGroups = F,
        addNodes = lnodes,
        main = "log2FC",
        width = 0.3,
        ncol = 6,
        stepX = 50
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
      scale_color_manual(values = colors) + theme_minimal(base_size = fontsize)
    
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
    req(reacValues$proteinData)
    "Experimental Design"
    })
  
  output$analysisSuccessMsg <- renderUI({
    req(reacValues$analysisSuccess)
    HTML('<h4>Successfully processed data!</h4>\n<p>
         Scroll to the top of the page to visit the analysis tabs.
         </p>')
  })
  
  
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