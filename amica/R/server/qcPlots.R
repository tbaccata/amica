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
    
    #df <- assay(reacValues$proteinData, assayName)
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
