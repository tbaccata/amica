#' prepare MA plot data
#' 
#' @param data data.frame, processed output from limma.
#' @param comparison char describing the group comparison (column pattern in @param data).
#' @param fcCutoff log2FC threshold.
#' @return data.frame containing relevant MA plot columns.
#' @examples
getMAPlotData <- function(data, comparison, fcCutoff) {
  check <- FALSE
  comp <- data[, c(geneName, paste0(logfcPrefix, comparison), paste0(avgExprPrefix, comparison ) )]
  
  if(length(grep(padjPrefix, colnames(data) )) > 0 ) {
    comp$pval <- data[[paste0(padjPrefix, comparison)]]
    check <- TRUE
  }
  colnames(comp) <-
    gsub(paste0("_", comparison) ,
         "",
         colnames(comp))
  
  comp$significant <- "no"
  if (check) {
    comp$significant[abs(comp$logFC) >= fcCutoff & comp$pval <= 0.05] <-
      "yes"
  } else {
    comp$significant[abs(comp$logFC) >= fcCutoff] <-
      "yes"
  }
  
  return(comp)
}

#' MA plot
#' 
#' @param pltData data.frame containing relevant columns.
#' @param colPalette vector of colors.
#' @param fontsize int describing the base plot font size.
#' @param pca_legend int describing the legend size.
#' @param legend_fontsize int describing ggplot point size.
#' @param xText charof xLabel
#' @param pointsize int describing pointsize in scatterplot.
#' @return ggplot object.
#' @examples
plotMAPlot <- function(pltData,
                       xText,
                       colPalette = NA,
                       fontsize = 14,
                       legend_fontsize = 10,
                       pointsize = 2) {
  
  colPalette <- ifelse(is.na(colPalette),
                       brewer.pal(5, "Set2"),
                       colPalette)
  print(colPalette)
  p <-
    ggplot(pltData,
           aes(
             x = logFC,
             y = AveExpr,
             color = significant,
             label = Gene.names,
             key = key
           )) + geom_point(size = pointsize) + coord_cartesian()  + scale_color_manual(values =
                                                                                         colPalette)
  
  p <-
    p + xlab("logFC") + ylab("Avg. Intensity")
  
  p <-
    p + theme_cowplot(font_size = fontsize) + background_grid() +
    theme(
      legend.text = element_text(size = legend_fontsize),
      legend.title = element_text(size = legend_fontsize)
    )
  
  p <- p + geom_text(
    data = subset(pltData, show_id),
    aes(logFC, AveExpr, label = Gene.names),
    position = position_jitter(width = 0.25, height = 0.25)
  ) + xlab(xText)
  
  return(p)
}

#' prepare volcano plot data
#' 
#' @param data data.frame, processed output from limma.
#' @param comparison char describing the group comparison (column pattern in @param data).
#' @param fcCutoff numeric, log2FC threshold.
#' @param sigCutoffValue char, choice between "adj.P.Val" and "p-value".
#' @param selectionChoice char, choice between "absolute", "enriched" and "reduced".
#' @param padjY bool, whether to plot -log10 padj on y-axis instead of -log10 p-value.
#' @param pvalCutoff numeric, filter on (adj.) p-value.
#' @return data.frame containing relevant volcano plot columns.
#' @examples
getVolcanoPlotData <-
  function(data,
           comparison,
           fcCutoff,
           sigCutoffValue,
           selectionChoice = "absolute",
           padjY = FALSE,
           pvalCutoff = NULL) {
    geneIdx <- which(colnames(data)==geneName)
    
    relIdx <-
      grep(
        paste0(
          logfcPrefix,
          comparison,
          "|",
          padjPrefix,
          comparison,
          "|",
          pvalPrefix,
          comparison,
          "|",
          avgExprPrefix,
          comparison
        ),
        colnames(data)
      )
    
    if (length(relIdx) < 1) return(NULL)
    
    comp <- data[, c(geneIdx, relIdx)]
    colnames(comp) <-
      gsub(paste0("_", comparison) ,
           "",
           colnames(comp))
    comp$key <- row.names(data)
    comp$show_id <- FALSE
    
    pvalColumn <- "adj.P.Val"
    if (sigCutoffValue!="none") {
      pvalColumn <- ifelse(sigCutoffValue == "adj.p-value", pvalColumn, "P.Value") 
    }
    
    if ("P.Value" %in% colnames(comp)) {
      
      if (padjY) {
        comp$nlog10_pval <- -log10(comp$adj.P.Val)
      } else {
        comp$nlog10_pval <- -log10(comp$P.Value)
      }
      
      pvalThresh <- ifelse(!is.null(pvalCutoff), pvalCutoff, 0.05)
      
      comp$significant <- "no"
      
      if (selectionChoice=="absolute") {
        comp$significant[comp[[pvalColumn]] <= pvalThresh & abs(comp$logFC) >= fcCutoff] <-
          "yes"
      } else if (selectionChoice=="enriched") {
        comp$significant[comp[[pvalColumn]] <= pvalThresh & comp$logFC >= fcCutoff] <-
          "yes"
      } else if (selectionChoice=="reduced") {
        comp$significant[comp[[pvalColumn]] <= pvalThresh & comp$logFC <= -fcCutoff] <-
          "yes"
      }
      
      
    } else {
      comp$significant <- "no"
      
      if (selectionChoice=="absolute") {
        comp$significant[abs(comp$logFC) >= fcCutoff] <-
          "yes"
      } else if (selectionChoice=="enriched") {
        comp$significant[comp$logFC >= fcCutoff] <-
          "yes"
      } else if (selectionChoice=="reduced") {
        comp$significant[comp$logFC <= -fcCutoff] <-
          "yes"
      }
      
    }
    
    if (sigCutoffValue == "none") {
      comp$significant <- "no"
      
      if (selectionChoice=="absolute") {
        comp$significant[abs(comp$logFC) >= fcCutoff] <-
          "yes"
      } else if (selectionChoice=="enriched") {
        comp$significant[comp$logFC >= fcCutoff] <-
          "yes"
      } else if (selectionChoice=="reduced") {
        comp$significant[comp$logFC <= -fcCutoff] <-
          "yes"
      }
    }
    
    return(comp)
  }

#' volcano plot
#' 
#' @param pltData data.frame containing relevant columns.
#' @param colPalette vector of colors.
#' @param padjYBoolean boolean, whether to plot -log10 padj on y-axis.
#' @param xText charof xLabel
#' @param fontsize int describing the base plot font size.
#' @param pca_legend int describing the legend size.
#' @param legend_fontsize int describing ggplot point size.
#' @param pointsize int describing pointsize in scatterplot.
#' @return ggplot object.
#' @examples
plotVolcanoPlot <-
  function(pltData,
           xText,
           padjYBoolean,
           colPalette = NA,
           fontsize = 14,
           legend_fontsize = 10,
           pointsize = 2) {
    colPalette <- ifelse(is.na(colPalette),
                         brewer.pal(5, "Set2"),
                         colPalette)
    p <-
      ggplot(pltData,
             aes(
               x = logFC,
               y = nlog10_pval,
               color = significant,
               label = Gene.names,
               key = key
             )) + geom_point(size = pointsize) + scale_color_manual(values =
                                                                      colPalette) +
      coord_cartesian()
    
    ylabel <-
      ifelse(padjYBoolean, "-log10(adj. p-value)", "-log10(p-value)")
    
    p <-
      p + xlab("logFC") + ylab(ylabel)
    
    p <-
      p + theme_cowplot(font_size = fontsize) + background_grid() +
      theme(
        legend.text = element_text(size = legend_fontsize),
        legend.title = element_text(size = legend_fontsize)
      )
    p <- p + geom_text(
      data = subset(pltData, show_id),
      aes(logFC, nlog10_pval, label = Gene.names),
      position = position_jitter(width = 0.25, height = 0.25)
      #,position = position_jitter(seed = 1)#position_nudge(y = -0.1)
    ) + xlab(xText)
    
    return(p)
  }


#' prepare fc-plot data
#' 
#' @param rnames vector, row.names of @param dataLimma.
#' @param dataLimma data.frame, processed output from limma.
#' @param enrichedMatrixSet data.frame of 1s and 0s denoting significance in group comparison, output from generateEnrichedMatrix.
#' @param foldChangeSelection char vector of length 2, comaprison pattern column.names of @param dataLimma.
#' @param lables char vector of length 2, rename x - and y-label.
#' @return data.frame containing relevant fc-plot columns.
#' @examples
getFCPlotData <-
  function(rnames,
           dataLimma,
           enrichedMatrixSet,
           foldChangeSelection,
           labels = NULL) {
    keepCols <-
      c(
        geneName,
        paste0(logfcPrefix, foldChangeSelection),
        paste0(padjPrefix, foldChangeSelection),
        paste0(pvalPrefix, foldChangeSelection)
      )
    
    plotData <- dataLimma[rnames, grep(paste0(keepCols, collapse = "|"),
                                       colnames(dataLimma))]
    
    newNames <- c(foldChangeSelection[1], foldChangeSelection[2])
    if (!is.null(labels) & all(labels != "")) {
      newNames[1] <- labels[1]
      newNames[2] <- labels[2]
    }
    
    plotData$significant <- ifelse(
      enrichedMatrixSet[rnames, foldChangeSelection[1]]==1 &
        enrichedMatrixSet[rnames, foldChangeSelection[2]]==1,
      "both",
      ifelse(
        enrichedMatrixSet[rnames, foldChangeSelection[1]]==1,
        newNames[1],
        ifelse(
          enrichedMatrixSet[rnames, foldChangeSelection[2]]==1,
          newNames[2],
          "none"
        )
      )
    )
    plotData$key <- plotData[[geneName]]
    return(plotData)
  }
