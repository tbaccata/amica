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

