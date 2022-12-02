#' Plot density plot
#' 
#' @param proteinData ProteomicsData object.
#' @param assayNames char 
#' @param groupFactors vector of groups to be ordered.
#' @return df in long format.
#' @examples
getAssayData <-
  function(proteinData,
           assayNames = "ImputedIntensity",
           groupFactors = NULL) {
    object <- 0
    if (assayNames == "ImputedIntensity") {
      object <- assay(proteinData, "ImputedIntensity")
      
      object <- toLongFormat(
        object[isQuantRnames(proteinData), ],
        proteinData,
        addGroup = TRUE,
        addContaminant = FALSE,
        addGeneName = FALSE
      )
    } else {
      object <- toLongFormat(
        assay(proteinData, assayNames),
        proteinData,
        addGroup = TRUE,
        addContaminant = FALSE,
        addGeneName = FALSE
      )
    }
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(object$group))) {
      object <- object[object$group %in% groupFactors, ]
      object$group <-
        factor(object$group, levels = groupFactors)
    }
    object
  }

#' Plot density plot
#' 
#' @param dataAssay data.frame in long format.
#' @param myColors vector of colors.
#' @param base_size int describing the base plot font size.
#' @param legend_size int describing the legend size.
#' @param assayNames char 
#' @return ggplot object.
#' @examples
plotDensityPlot <-
  function(dataAssay,
           assayNames,
           myColors,
           base_size = 14,
           legend_size = 10) {
    p <-
      ggplot(dataAssay[!is.na(dataAssay$value),], aes(x = value, color = colname)) +
      # ggplot(dataAssay[!is.na(dataAssay$value), ], aes(x = value, color =
      #                                                        colname)) +
      geom_density() +
      scale_color_manual(values = myColors) +
      theme_cowplot(font_size = base_size) +
      background_grid() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 10
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_size),
      ) + xlab(paste0(gsub("Intensity", " intensities", assayNames), " (log2)")) + ylab("Density") +
      ggtitle(paste0(
        'Density plot (',
        gsub("Intensity", " intensities", assayNames),
        ')'
      ))
  }



#' Plot box plot
#' 
#' @param dataAssay data.frame in long format.
#' @param myGroupColors vector of colors.
#' @param boxplot_base int describing the base plot font size.
#' @param boxplot_legend int describing the legend size.
#' @param assayNames char 
#' @return ggplot object.
#' @examples
plotBoxPlot <-
  function(dataAssay,
          assayNames,
          myGroupColors,
          boxplot_base = 14,
          boxplot_legend = 12) {
    p <-
      ggplot(dataAssay[!is.na(dataAssay$value),], aes(x = colname, y = value, fill =
                                                        group)) +
      geom_boxplot(
        outlier.shape = NA,
        outlier.color = NULL,
        outlier.fill = NULL
      ) +
      theme_cowplot(font_size = boxplot_base) +
      background_grid() +
      scale_fill_manual(values = myGroupColors) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 10
        ),
        legend.text = element_text(size = boxplot_legend),
        legend.title = element_blank()
      ) + ylab(paste0(gsub("Intensity", " intensities", assayNames),
                      " (log2)")) + xlab("") +
      ggtitle(paste0(
        'Box plot (',
        gsub("Intensity", " intensities", assayNames),
        ')'
      ))
    
    p
  }

#' Compute pca plot
#' 
#' @param plotData data.frame of intensities in wide format.
#' @param myGroupColors vector of colors.
#' @param boxplot_base int describing the base plot font size.
#' @param boxplot_legend int describing the legend size.
#' @param assayNames char 
#' @param groupFactors vector of groups to be ordered.
#' @param groupInputs vector of groups to be plotted.
#' @return list of df_out containing coordinates and percentages char vector of % variance.
#' @examples
computePCA <-
  function(plotData,
           expDesign,
           groupFactors = NULL,
           groupInputs = NULL) {
    # filter only by selected groups
    if (is.null(groupInputs))
      groupInputs <- unique(expDesign$groups)
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(expDesign$groups))) {
      expDesign$groups <-
        factor(expDesign$groups, levels = c(groupFactors,
                                            setdiff(groupInputs, groupFactors)))
    }
    
    plotData <-
      plotData[, expDesign$samples[expDesign$groups %in% groupInputs]]
    plotData <- plotData[complete.cases(plotData),]
    
    # compute PCA
    pca <- prcomp(as.data.frame(t(plotData)))
    df_out <- as.data.frame(pca$x)
    
    df_out$group <-
      expDesign$groups[expDesign$samples %in% row.names(df_out)]
    df_out$sample <- rownames(df_out)
    df_out$key <- df_out$sample
    df_out$show_id <- FALSE
    
    tmp <- summary(pca)
    percentage <-
      round(100 * tmp$importance['Proportion of Variance', 1:2], 2)
    percentage <-
      paste(colnames(df_out), "(", paste(as.character(percentage), "%", ")", sep =
                                           ""))
    return(list(df_out = df_out, percentage = percentage))
  }

#' Compute pca plot
#' 
#' @param df_out df, output from computePCA.
#' @param percentage char, output from computePCA; vector of %var of 1st and 2nd PC.
#' @param myGroupColors vector of colors.
#' @param pca_base int describing the base plot font size.
#' @param pca_legend int describing the legend size.
#' @param pca_pointsize int describing ggplot point size.
#' @param assayNames char 
#' @param groupFactors vector of groups to be ordered.
#' @param groupInputs vector of groups to be plotted.
#' @return list of df_out containing coordinates and percentages char vector of % variance.
#' @examples
plotPCA <-
  function(df_out,
           percentage,
           myGroupColors,
           assayNames = 'LFQIntensity',
           pca_base = 14,
           pca_legend = 12,
           pca_pointsize = 2) {
    p <-
      ggplot(df_out, aes(
        x = PC1,
        y = PC2,
        color = group,
        label = sample,
        key = key
      ))
    p <-
      p + geom_point(size = pca_pointsize) + xlab(percentage[1]) + ylab(percentage[2])
    p <-
      p + scale_color_manual(values = myGroupColors) + theme_cowplot(font_size = pca_base) +
      background_grid() +
      theme(legend.text = element_text(size = pca_legend))
    p <-
      p + ggtitle(paste0('PCA (', gsub("Intensity", " intensities", assayNames), ')'))
    
    p <- p + geom_text(
      data = subset(df_out, show_id),
      aes(PC1, PC2, label = sample),
      position = position_jitter(width = 0.25, height = 0.5)
    )
    
    p
  }

#' Plot correlation plot (plotly)
#' 
#' @param df data.frame with intensity values.
#' @param myGroupColors vector of colors.
#' @param heatColors vector of colors.
#' @param annotSamples boolean to annotate samples in heatmap.
#' @param assayName name of assay to be plotted.
#' @param groupFactors vector of groups to be ordered.
#' @param groupInputs vector of groups to be plotted.
#' @param assayNames char 
#' @return plotly object.
#' @examples
plotCorrPlotly <-
  function(df,
           expDesign,
           myGroupColors,
           heatColors,
           assayName,
           annotSamples = T,
           groupFactors = NULL,
           groupInputs = NULL) {
    if (is.null(groupInputs))
      groupInputs <- unique(expDesign$groups)
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(expDesign$groups))) {
      expDesign$groups <-
        factor(expDesign$groups, levels = c(groupFactors,
                                            setdiff(groupInputs, groupFactors)))
    }
    
    #df <- assay(proteinData, assayName)
    df <- df[, expDesign$samples[expDesign$groups %in% groupInputs]]
    
    annot <- expDesign
    row.names(annot) <- annot$samples
    annot <- annot[names(df),]
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
        main = paste0(
          "<b>Correlation plot (",
          gsub("Intensity", " intensities", assayName),
          ")</b>"
        ),
        colors = heatColors,
        row_side_palette = myGroupColors,
        row_side_colors = annot,
        plot_method = "plotly",
        key.title = "Pearson Correlation"
      )
    } else {
      p <- heatmaply_cor(
        round(corDf, 3),
        xlab = "",
        ylab = "",
        colors = heatColors,
        limits = limits,
        main = paste0(
          "<b>Correlation plot (",
          gsub("Intensity", " intensities", assayName),
          ")</b>"
        ),
        plot_method = "plotly",
        key.title = "Pearson Correlation"
      )
    }
    p
  }

#' Plot CV boxplot plot (plotly)
#' 
#' @param dataAssay df in long format with intensity values.
#' @param myGroupColors vector of colors.
#' @param cv_base int describing the base plot font size.
#' @param cv_legend int describing the legend size.
#' @param assayNames char 
#' @return plotly object.
#' @examples
plotCoeffVarPlot <-
  function(dataAssay,
           myGroupColors,
           assayNames = 'LFQIntensity',
           cv_base = 14,
           cv_legend = 12) {
    dataAssay$value <- 2 ^ dataAssay$value
    cvs <-
      aggregate(value ~ rowname + group, dataAssay, function(x)
        100 * (sd(x) / mean(x)))
    
    p <-
      ggplot(cvs[!is.na(cvs$value),], aes(x = group, y = value, fill = group)) +
      geom_boxplot(
        outlier.shape = NA,
        outlier.color = NULL,
        outlier.fill = NULL
      ) +
      scale_fill_manual(values = myGroupColors) +
      theme_cowplot(font_size = cv_base) +
      background_grid() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = cv_base
        ),
        legend.text = element_text(size = cv_legend),
        legend.title = element_blank()
      ) + ylab("Coefficient of Variation (%)") + xlab("") +
      ggtitle(paste0('CV plot (', gsub("Intensity", " intensities", assayNames), ')'))
    
    p
  }


#' Plot most abundant proteins per sample barplot
#' 
#' @param dataAssay df in long format with intensity values.
#' @param sample choosen sample.
#' @param color int describing the base plot font size.
#' @param n int describing top n most abundant proteins.
#' @param abundant_base int describing the base font size.
#' @return ggplot object.
#' @examples
plotMostAbundantProteinsInSample <-
  function(dataAssay,
           samples,
           color = "skyblue",
           n = 15,
           abundant_base = 14) {
    dataAssay <- dataAssay[dataAssay$colname == samples,]
    dataAssay <- dataAssay[order(dataAssay$value, decreasing = T),]
    dataAssay$value <- 100 * dataAssay$value / sum(dataAssay$value, na.rm = T)
    
    plotdf <- head(dataAssay, n)
    plotdf$Gene.name <- gsub(";.*", "", plotdf$Gene.name)
    
    p <-
      ggplot(plotdf, aes(x = reorder(Gene.name, value), y = value)) +
      geom_bar(stat = "identity", fill = color) +
      coord_flip() +
      theme_cowplot(font_size = abundant_base) +
      background_grid() +
      theme(legend.title = element_blank()) + xlab("") + ylab(paste0("iBAQ intensities (%)")) +
      ggtitle(label = paste0(n, " most abundant proteins in ", samples))
    p
  }


#' Plot number of identified proteins barplot
#' 
#' @param dataAssay data.frame, in long format with intensity values.
#' @param expDesign data.frame, experimetnal design mapping samples to groups.
#' @param groupFactors vector of groups to be ordered.
#' @param barplotId_base int describing the base font size.
#' @param barplotId_legend int describing the legend font size.
#' @return ggplot object.
#' @examples
plotNumberIdentifiedProteins <-
  function(dataAssay,
           expDesign,
           myGroupColors,
           groupFactors = NULL,
           barplotId_base = 14,
           barplotId_legend = 12) {
    df <-
      as.data.frame(apply(dataAssay, 2, function(x)
        sum(!is.na(x))))
    names(df) <- "Number"
    df$Sample <- row.names(df)
    
    midx <- match(df$Sample, expDesign$samples)
    df$group <- expDesign$groups[midx]
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(df$group))) {
      df <- df[df$group %in% groupFactors,]
      df$group <- factor(df$group, levels = groupFactors)
    }
    
    p <- ggplot(df, aes(x = Sample, y = Number, fill = group)) +
      geom_bar(stat = "identity") +
      theme_cowplot(font_size = barplotId_base) +
      background_grid() +
      scale_fill_manual(values = myGroupColors) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = barplotId_base
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = barplotId_legend)
      ) + xlab("") + ylab("#of identified proteins")
    
    p
  }


#' Plot number of missing values per sample barplot
#' 
#' @param dataAssay data.frame, in long format with intensity values.
#' @param expDesign data.frame, experimetnal design mapping samples to groups.
#' @param groupFactors vector of groups to be ordered.
#' @param barplotMv_base int describing the base font size.
#' @param barplotMv_legend int describing the legend font size.
#' @return ggplot object.
#' @examples
plotMissingValues <-
  function(dataAssay,
           expDesign,
           myGroupColors,
           groupFactors = NULL,
           barplotMv_base = 14,
           barplotMv_legend = 12) {
    df <-
      as.data.frame(apply(dataAssay, 2, function(x) {
        100 * sum(is.na(x)) / length(x)
      }))
    names(df) <- "Number"
    df$Sample <- row.names(df)
    
    midx <- match(df$Sample, expDesign$samples)
    df$group <- expDesign$groups[midx]
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(df$group))) {
      df <- df[df$group %in% groupFactors, ]
      df$group <- factor(df$group, levels = groupFactors)
    }
    
    p <- ggplot(df, aes(x = Sample, y = Number, fill = group)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = myGroupColors) +
      theme_cowplot(font_size = barplotMv_base) +
      background_grid() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = barplotMv_base
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = barplotMv_legend)
      ) + xlab("") + ylab("missing values (%)")
    
    p
  }


#' Plot percentage of contaminants per sample 
#' 
#' @param proteinData ProteomicsData object.
#' @param expDesign data.frame, experimetnal design mapping samples to groups.
#' @param myGroupColors vector of colors.
#' @param groupFactors vector of groups to be ordered.
#' @param contaminants_base int describing the base font size.
#' @param contaminants_legend int describing the legend font size.
#' @return ggplot object.
#' @examples
plotContaminants <-
  function(proteinData,
           expDesign,
           myGroupColors,
           groupFactors = NULL,
           contaminants_base = 14,
           contaminants_legend = 12) {
    tmp <- 2 ^ assay(proteinData, "iBAQ")
    
    allInts <-
      apply(2 ^ assay(proteinData, "iBAQ"), 2, sum, na.rm = T)
    contsInts <-
      apply(tmp[rowData(proteinData)[[contaminantCol]] == "+", ], 2, sum, na.rm =
              T)
    contsInts <- contsInts / allInts
    
    contsInts <- reshape2::melt(contsInts)
    
    midx <-
      match(row.names(contsInts),
            expDesign$samples)
    contsInts$group <- expDesign$groups[midx]
    contsInts$Sample <- row.names(contsInts)
    contsInts$value <- 100 * contsInts$value
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(contsInts$group))) {
      contsInts <-
        contsInts[contsInts$group %in% groupFactors, ]
      contsInts$group <-
        factor(contsInts$group, levels = groupFactors)
    }
    
    p <-
      ggplot(contsInts, aes(x = Sample, y = value, fill = group)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = myGroupColors) +
      theme_cowplot(font_size = contaminants_base) +
      background_grid() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = contaminants_base
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = contaminants_legend)
      ) + xlab("") + ylab("iBAQ intensities (%)") +
      ggtitle(label = "Relative Amount of Contaminants")
    
    p
  }


#' Plot id protein overlap plot (plotly)
#' 
#' @param dataAssay data.frame with intensity values.
#' @param myGroupColors vector of colors.
#' @param heatColors vector of colors.
#' @param annotSamples boolean to annotate samples in heatmap.
#' @param assayName name of assay to be plotted.
#' @param metric char, one of "jaccard_index", "overlap_coeff" or "num_shared"
#' @param groupFactors vector of groups to be ordered.
#' @return plotly object.
#' @examples
plotOverlaply <-
  function(dataAssay,
           expDesign,
           myGroupColors,
           heatColors,
           annotSamples = TRUE,
           metric = "jaccard_index",
           groupFactors = NULL) {

    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(expDesign$groups))) {
      expDesign <- expDesign[expDesign$groups %in% groupFactors,]
      expDesign$groups <-
        factor(expDesign$groups, levels = c(groupFactors))
    }
    groupInputs <- unique(expDesign$groups)
    dataAssay <-
      dataAssay[, expDesign$samples[expDesign$groups %in% groupInputs]]
    
    annot <- expDesign
    row.names(annot) <- annot$samples
    annot$samples <- NULL
    
    idProts <- list()
    for (col in names(dataAssay)) {
      tmp <- row.names(dataAssay)[!is.na(dataAssay[, col])]
      idProts[[col]] <- tmp
    }
    
    overlapDf <- calc_pairwise_overlaps(idProts)
    idProts <- NULL
    
    limits <- c(0, 1)
    title <- "Overlap Coefficient"
    df <- head(overlapDf)
    df1 <- head(overlapDf)
    
    if (metric == 'num_shared') {
      df <- overlapDf[, c("sample1", "sample2", "num_shared")]
      df1 <-
        data.frame(
          sample1 = df$sample2,
          sample2 = df$sample1,
          overlap = df$num_shared
        )
      limits <-
        c(min(df[, 3], na.rm = T) - 0.05, max(df[, 3], na.rm = T))
      title <- "#shared proteins"
    } else if (metric == 'overlap_coefficient') {
      df <- overlapDf[, c("sample1", "sample2", "overlap")]
      df1 <-
        data.frame(
          sample1 = df$sample2,
          sample2 = df$sample1,
          overlap = df$overlap
        )
      limits <- c(min(df[, 3], na.rm = T) - 0.05, 1)
    } else if (metric == 'jaccard_index') {
      df <- overlapDf[, c("sample1", "sample2", "jaccard")]
      df1 <-
        data.frame(
          sample1 = df$sample2,
          sample2 = df$sample1,
          overlap = df$jaccard
        )
      limits <- c(min(df[, 3], na.rm = T) - 0.05, 1)
      title <- "Jaccard Coefficient"
    }
    names(df) <- c("sample1", "sample2", "overlap")
    df <- rbind(df, df1)
    overlap_matrix <-
      acast(df, sample1 ~ sample2, value.var = "overlap")
    
    overlap_matrix <- overlap_matrix[row.names(annot),]

    if (annotSamples) {
      p <- heatmaply_cor(
        overlap_matrix,
        xlab = "",
        ylab = "",
        limits = limits,
        row_side_palette = myGroupColors,
        row_side_colors = annot,
        fontsize_row=12,
        fontsize_col=12,
        column_text_angle = -270,
        plot_method = "plotly",
        key.title = title,
        colors = heatColors
      )
    } else {
      p <- heatmaply_cor(
        overlap_matrix,
        xlab = "",
        ylab = "",
        limits = limits,
        fontsize_row=12,
        fontsize_col=12,
        column_text_angle = -270,
        plot_method = "plotly",
        key.title = title,
        colors = heatColors
      )
    }
    return(list(df=overlapDf, plot=p))
  }


#' Scatter plot
#' 
#' @param proteinData ProteomicsData object.
#' @param assay1 char, valid assay in @param proteinData.
#' @param assay2 char, valid assay in @param proteinData.
#' @param selection1 char, valid sample in @param proteinData.
#' @param selection2 char, valid sample in @param proteinData.
#' @param myScatterColors vector of colors.
#' @param plot_fontsize int describing the base font size.
#' @param plot_legendsize int describing the legend font size.
#' @param showLine one of "straight line", "linear regression" or "none".
#' @param get_scatter_data plotly event data.
#' @return ggplot object.
#' @examples
plotScatterPlot <-
  function(proteinData,
           assay1,
           assay2,
           selection1,
           selection2,
           myScatterColors,
           plot_fontsize = 14,
           plot_legendsize = 10,
           showLine = "straight line",
           get_scatter_data = NULL) {
    xLabel <- paste0(assay1, "_", selection1)
    yLabel <- paste0(assay2, "_", selection2)
    plotData <- data.frame()
    
    if (assay1 == assay2) {
      plotData <-
        assay(proteinData, assay1)[, c(selection1, selection2)]
      names(plotData) <- c("x", "y")
      
    } else {
      plotData <-
        assay(proteinData, assay1)[, selection1, drop = F]
      names(plotData) <- c("x")
      plotData[[selection2]] <-
        assay(proteinData, assay2)[, selection2, drop = T]
    }
    names(plotData)[2] <- "y"
    
    plotData$Contaminant <- "no"
    
    if (contaminantCol %in% names(rowData(proteinData))) {
      plotData$Contaminant <-
        ifelse(rowData(proteinData)[[contaminantCol]] ==
                 "+",
               "yes",
               "no")
    }
    
    plotData$Gene <-
      rowData(proteinData)[[geneName]]
    plotData$key <- plotData$Gene
    
    plotData$show_id <- FALSE
    if (!is.null(get_scatter_data)) {
      plotData[plotData$key %in% get_scatter_data$key, "show_id"] <- TRUE
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
    
    pu <-
      ggplot(plotData,
             aes(
               x = x,
               y = y,
               label = Gene,
               key = key,
               color = Contaminant
             )) + theme_cowplot(font_size = plot_fontsize) +
      background_grid() +
      labs(x = xLabel, y = yLabel)  +
      theme(
        legend.title = element_text(size = plot_legendsize),
        legend.text = element_text(size = plot_legendsize)
      ) +
      geom_point() + scale_color_manual(values = myScatterColors) # + ggtitle(title) #scale_color_brewer(palette = "Paired")
    
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
        color = myScatterColors[3]
      )
    }
    
    pu <- pu + geom_text(
      data = subset(plotData, show_id),
      aes(x,
          y,
          label = Gene),
      position = position_jitter(width = 0.25, height = 0.25)
    )
    pu
  }
