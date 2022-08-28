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
#' sourcePath <- "data/PXD0016455/"
#' 
#' expDesign <-
#'   read.table(
#'     paste0(sourcePath, "design.txt"),
#'     header = T,
#'     stringsAsFactors = F
#'   )
#' 
#'   outData <-
#'     readInAmicaSumm(paste0(sourcePath, "amica_proteinGroups.tsv"),
#'                     expDesign)
#'   proteinData <- outData$protData
#'   contrastMatrix = outData$contrasts
#'   dataLimma = outData$comparisons
#'
#' comparison <- "PGRMC1__vs__MIAPACA" 
#' pltData <- getVolcanoPlotData(
#' dataLimma,
#' sample = comparison,
#' fcThresh = 1.5,
#' sigCutoffValue = "adj.p-value"
#' )
#' 
#' xText <- unlist(strsplit(comparison, "__vs__"))
#' xText <- paste0("log2FC(", xText[1], "/", xText[2], ")")
#' 
#' p <- plotMAPlot(pltData, xText)
#' 
#' 
plotMAPlot <- function(pltData,
                       xText,
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
#' 
#' sourcePath <- "data/PXD0016455/"
#' 
#' expDesign <-
#'   read.table(
#'     paste0(sourcePath, "design.txt"),
#'     header = T,
#'     stringsAsFactors = F
#'   )
#' 
#'   outData <-
#'     readInAmicaSumm(paste0(sourcePath, "amica_proteinGroups.tsv"),
#'                     expDesign)
#'   proteinData <- outData$protData
#'   contrastMatrix = outData$contrasts
#'   dataLimma = outData$comparisons
#'
#' comparison <- "PGRMC1__vs__MIAPACA" 
#' pltData <- getVolcanoPlotData(
#' dataLimma,
#' sample = comparison,
#' fcThresh = 1.5,
#' sigCutoffValue = "adj.p-value"
#' )
#' 
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
#' sourcePath <- "data/PXD0016455/"
#' 
#' expDesign <-
#'   read.table(
#'     paste0(sourcePath, "design.txt"),
#'     header = T,
#'     stringsAsFactors = F
#'   )
#' 
#'   outData <-
#'     readInAmicaSumm(paste0(sourcePath, "amica_proteinGroups.tsv"),
#'                     expDesign)
#'   proteinData <- outData$protData
#'   contrastMatrix = outData$contrasts
#'   dataLimma = outData$comparisons
#'
#' comparison <- "PGRMC1__vs__MIAPACA" 
#' pltData <- getVolcanoPlotData(
#' dataLimma,
#' sample = comparison,
#' fcThresh = 1.5,
#' sigCutoffValue = "adj.p-value"
#' )
#' 
#' xText <- unlist(strsplit(comparison, "__vs__"))
#' xText <- paste0("log2FC(", xText[1], "/", xText[2], ")")
#' 
#' p <- plotVolcanoPlot(pltData, xText)
#' 
#' 
plotVolcanoPlot <-
  function(pltData,
           xText,
           padjYBoolean = TRUE,
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
    plotData$show_id <- FALSE
    return(plotData)
  }

#' fold-change plot
#' 
#' @param plotData data.frame containing relevant columns.
#' @param selection char vector of length 2, containing group comparison names.
#' @param colors vector of colors.
#' @param labels char vector of length 2, containing new labels corresponding to @param selection.
#' @param showLine one of "straight line", "linear regression" or "none".
#' @param fontsize int describing the base plot font size.
#' @param pca_legend int describing the legend size.
#' @param legendFontSize int describing ggplot point size.
#' @param pointsize int describing pointsize in scatterplot.
#' @return ggplot object.
#' @examples
plotFoldChangePlot <- function(plotData, 
                               selection,
                               colors = NA,
                               labels = NULL,
                               showLine = "straight line",
                               fontsize = 14,
                               legendFontSize = 12,
                               pointsize = 2
                               ) {
  
  colors <- ifelse(is.na(colors),
                       brewer.pal(5, "Set2"),
                       colors)
  
  formula <- as.formula(paste0(paste0("logFC_", selection[2]), ' ~ ',
                               paste0("logFC_", selection[1])   ) )
  
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
      theme_cowplot(font_size = fontsize) + background_grid() +
      theme(
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
}

#' Euler diagram
#' 
#' @param comparisons vector of group comparisons.
#' @param binMat binary data.frame of significant proteins (rows) in comparisons (cols).
#' @param showQuant boolean, displays number of elements in sets.
#' @param bool boolean, whether to plot lines around Euler ellipses..
#' @param showLegend boolean, whether to displahy legend.
#' @param newMultiNames named list, with old and new names of comparisons
#' @param colPalette vector, colors.
#' @return Euler plot.
#' @examples
plotEulerDiagram <-
  function(comparisons,
           binMat,
           showQuant = TRUE,
           bool = TRUE,
           showLegend = TRUE,
           newMultiNames = NULL,
           colPalette = NA,
           minSizeGradient=1, 
           maxSizeGradient=7) {
    if (length(comparisons) < 2) {
      stop("Need at least two comparisons to render Euler plot.")
    } else if (length(comparisons) > 5) {
      stop("Cannot output Euler plot for more than 5 sets.")
    }
    
    colPalette <- ifelse(is.na(colPalette),
                         brewer.pal(length(comparisons), "Set2"),
                         colPalette)
    
    if (!is.null(newMultiNames) && all(newMultiNames$new != "") &&
        length(newMultiNames$new) == length(comparisons)) {
      for (idx in seq_along(newMultiNames$old)) {
        names(binMat)[which(names(binMat) == newMultiNames$old[idx])] <-
          newMultiNames$new[idx]
        comparisons[idx] <- newMultiNames$new[idx]
      }
    }
    
    fit <- euler(binMat[, comparisons])
    comps <- grep("&",
                  names(fit$original),
                  invert = T,
                  value = T)
    numComps <- length(comps)
    lty <- ifelse(bool, 1, 0)
    plot(
      fit,
      quantities = showQuant,
      fills = list(fill = colPalette),
      legend = ifelse(showLegend, list(labels = comps), F),
      lty = lty
    )
  }

#' dotplot data
#' 
#' @param widestats data.frame, containing (adj) p-values (rows) in group comparisons (columns).
#' @param widefcs data.frame, containing log2FCs (rows) in group comparisons (columns).
#' @param df data.frame, containing intensities (rows) in samples (columns).
#' @param filtData data.frame, quantified rowData from proteinData.
#' @param group2comps data.frame, containing group comparisons and corresponding group for dotplot.
#' @param expDesign named list, with old and new names of comparisons
#' @param sigCutoffValue char, choice between "adj.P.Val" and "p-value".
#' @return dotplot plot in long format.
#' @examples
getDotplotData <-
  function(widestats,
           widefcs,
           df,
           filtData,
           group2comps,
           expDesign,
           sigCutoffValue = "p-value") {
    
  pattern <- "adj.P.Val_"
  if (sigCutoffValue == "p-value") pattern <- "P.Value_"

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
  
  intData <- stats::reshape(df, idvar = "rowname",
                            ids = rownames(df), times = names(df),
                            timevar = "colname", varying = list(names(df)),
                            direction = "long", v.names = "value")
  midx <- match(intData$colname, expDesign$samples)
  intData$group <- expDesign$groups[midx]
  
  intData <- intData[intData$group %in% group2comps$group,]
  intData <-
    aggregate(value ~ rowname + group, intData, mean, na.rm=TRUE, na.action=NULL)
  names(intData) <- c('ProteinID', 'Group', 'AvgIntensity')
  
  longData <- merge(longData, intData, by = c('ProteinID', 'Group'))
  
  longData$Gene <- filtData[longData$ProteinID, "Gene.names"]
  longData$Gene <- gsub(";.*", "", longData$Gene)

  longData$significant <- factor(ifelse(longData$padj <= 0.05, 1.5, 0))
  longData
}

#' dotplot
#' 
#' @param dataDotplot data.frame, output from getDotplotData.
#' @param dotplotGroupsDf data.frame, containing group comparisons and corresponding group for dotplot.
#' @param dotplotColors vector of colors.
#' @param minColorGradient numeric.
#' @param maxColorGradient numeric.
#' @param minSizeGradient int.
#' @param maxSizeGradient int.
#' @param clusteringMetric char, one of "AvgIntensity" and "log2FC".
#' @param dotplot_distance_metric char, valid option from hclust.
#' @param dotplot_clustering_method char, valid option from hclust.
#' @param dotplot_cluster_columns boolean, performs clustering on columns.
#' @param sigCutoffValue char, choice between "adj.P.Val" and "p-value".
#' @param dotplot_ctrl_substraction boolean, ignores proteins with negative log2FCs.
#' @return ggplot2 object.
#' @examples
plotDotplot <-
  function(dataDotplot,
           dotplotGroupsDf,
           dotplotColors=NULL,
           minColorGradient=NULL,
           maxColorGradient=NULL,
           minSizeGradient=1,
           maxSizeGradient=7,
           clusteringMetric = "AvgIntensity",
           dotplot_distance_metric = "canberra",
           dotplot_clustering_method = "complete",
           dotplot_cluster_columns=TRUE,
           sigCutoffValue="p-value",
           dotplot_ctrl_substraction=TRUE) {
    
    if (is.null(dotplotColors))
      dotplotColors <- viridis(20, option = "viridis")
    
  mat <- dataDotplot %>%
    dplyr::select(Gene, Group, all_of(clusteringMetric)) %>%  # drop unused columns to faciliate widening
    pivot_wider(names_from = Group, values_from = all_of(clusteringMetric)) %>%
    data.frame() # make df as tibbles -> matrix annoying
  
  row.names(mat) <- mat$Gene  # put gene in `row`
  mat$Gene <- NULL #drop gene column as now in rows
  
  tmat <- mat %>% as.matrix()
  tmat[is.na(tmat)] <- 0
  
  dist_meth <- ifelse(!is.null(dotplot_distance_metric),
                      dotplot_distance_metric, 
                      "canberra")
  clst_meth <- ifelse(!is.null(dotplot_clustering_method),
                      dotplot_clustering_method,
                      "complete")
  
  if (!is.null(dotplot_cluster_columns) && dotplot_cluster_columns) { # order by clustering
    cclust <- hclust(dist(t(tmat), method = dist_meth),
                     method = clst_meth) # hclust with distance matrix
    cclust <- reorder(as.dendrogram(cclust), colMeans(mat))
    dataDotplot$Group <-
      factor(dataDotplot$Group, levels = names(mat)[as.hclust(cclust)$order])
  } else { # order by input
    dataDotplot$Group <-
      factor(dataDotplot$Group, levels = dotplotGroupsDf$group)
  }
  
  rclust <- hclust(dist(tmat, 
                        method = dist_meth), 
                   method = clst_meth) # hclust with distance matrix
  rclust <- reorder(as.dendrogram(rclust), rowMeans(mat))
  dataDotplot$Gene <-
    factor(dataDotplot$Gene, levels = rownames(mat)[as.hclust(rclust)$order])
  
  
  signficantTitle <- "adj.p-value"
  if (sigCutoffValue == "p-value")
    signficantTitle <- "p-value"
  
  filterValues <- TRUE
  if (clusteringMetric == "log2FC" &&
      !is.null(dotplot_ctrl_substraction)) {
    filterValues <- dotplot_ctrl_substraction
  }
  
  preFiltered <- dataDotplot
  if (filterValues) {
    preFiltered <- dataDotplot %>% filter(!!rlang::sym(clusteringMetric) > 0)
  }
  
  x <- preFiltered[[clusteringMetric]]
  
  legName <- "Relative AvgIntensity"
  if (clusteringMetric == "log2FC") {
    preFiltered$relativeAbundance <- (x-min(x,na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
    legName <- "Relative Fold Change"
  } else {
    preFiltered$relativeAbundance <- x/max(x, na.rm = T)
  }
  
  if (is.null(minColorGradient)){
    minColorGradient < min(preFiltered$log2FC, na.rm = T)
  }
  if (is.null(maxColorGradient)){
    maxColorGradient < max(preFiltered$log2FC, na.rm = T)
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
    )) +  theme(legend.justification = "top") +
    ylab('') + xlab('') +
    scale_x_discrete(position = "top") +
    scale_size_continuous(
      range = c(minSizeGradient, maxSizeGradient),
      name = legName, #paste("Relative", clusteringMetric),
      guide = guide_legend(order=2, override.aes = list(shape = 19)),
      breaks = c(
        min(preFiltered$relativeAbundance ) + 0.01,
        max(preFiltered$relativeAbundance )
      ),
      labels = c('', '')
    ) +
    scale_fill_gradientn(
      colours = dotplotColors,
      limits = c(
        ifelse(
          minColorGradient < min(preFiltered$log2FC, na.rm = T),
          min(preFiltered$log2FC, na.rm = T),
          minColorGradient
        ),
        maxColorGradient
      ),
      oob = scales::squish,
      name = 'Log2FC',
      guide = guide_colorbar(order = 1)
    ) +
    scale_color_manual(
      signficantTitle,
      values = c('skyblue', 'black'),
      limits = c('0', '1.5'),
      labels = c('> 0.05', '\u2264 0.05')
    )
  p
}