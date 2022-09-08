filterOnValidValues <-
  function(data,
           mappings,
           groupsToConsider = NULL,
           minValue = 2,
           method = "in_one_group") {
    if (!any(duplicated(as.character(mappings$groups))))
      return(row.names(data))
    i <- 1
    
    data <- data[apply(data, 1, function(x) !all(is.na(x)) ),]
    
    
    relevantGroups <- unique(mappings$groups)
    if (!is.null(groupsToConsider)) {
      relevantGroups <- groupsToConsider
    }
    
    boolMat <- matrix(rep(TRUE,times=nrow(data)*length(relevantGroups)), ncol=length(relevantGroups))
    row.names(boolMat) <- row.names(data)
    
    dfBool <- !is.na(data)
    for (i in seq_along(relevantGroups)) {
      group <- relevantGroups[i]
      
      tmpSamples <-
        mappings$samples[mappings$groups == group]
      minValue <- min(length(tmpSamples), minValue)
      
      idxs <-
        grep(paste0(tmpSamples, collapse = "|"), colnames(dfBool))
      
      sumValidValues <- apply(dfBool[,idxs],1, sum)
      validVals <- sumValidValues >= minValue
      boolMat[,i] <- validVals
    }
    
    outBool <- c()
    if (method == "in_one_group") {
      outBool <- apply(boolMat, 1, any)
    } else if (method == "in_each_group") {
      outBool <- apply(boolMat, 1, all)
    }

    return(row.names(data[outBool,])) 
  }


imputeIntensities <- function(df_int,
                              method = "normal",
                              downshift = 1.8,
                              width = 0.3, seed = 12345) {
  set.seed(seed)

  if (method == "normal") {
    df_int[] <-
      lapply(df_int, function(x)
        replace(x, is.na(x), rnorm(
          mean = median(x, na.rm = T) - downshift * sd(x, na.rm = T),
          sd = sd(x, na.rm = T) * width,
          n = sum(is.na(x)) )
          )
        )
  }
  if (method == "min") {
    const <- min(df_int, na.rm = T)
    df_int[is.na(df_int)] <- const
  }
  if (method == "global") {
    medians <- apply(df_int, MARGIN = 2, FUN=median, na.rm=TRUE)
    median_global <- median(medians)
    sd_global <- sd(as.numeric(as.matrix(df_int)), na.rm=TRUE)

    mu_imputed <- median_global - downshift*sd_global
    sd_imputed <- width*sd_global

    df_int[] <-
      lapply(df_int, function(x)
        replace(x, is.na(x), rnorm(
          mean = mu_imputed,
          sd = sd_imputed,
          n = sum(is.na(x)) )
        )
      )
  }
  return(df_int)
}


renormalizeIntensities <- function(data, method = "None") {
  if (method == "Quantile") {
    data <- normalizeBetweenArrays(data, method = "quantile")
  } else if (method == "VSN") {
    data <- normalizeVSN(2 ^ data)
  } else if (method == "Median.centering") {
    data[] <- lapply(data, function(x) {
      gMedian = median(x, na.rm = TRUE)
      x - gMedian
    })
  }
  return(as.data.frame(data))
}


groupComparisons <-
  function(imp_df_int, comparisons, expDesign, pep.count.table = NULL) {
    df <- data.frame()
    idx <- 1

    
    for (i in 1:nrow(comparisons)) {
      group1 <- comparisons[i, 1]
      group2 <- comparisons[i, 2]
      
      group1names <- expDesign[expDesign$groups==group1, "samples"]
      group2names <- expDesign[expDesign$groups==group2, "samples"]
      
      group1_idx <- grep(paste0("^",group1names,"$", collapse = "|"), colnames(imp_df_int))
      group2_idx <- grep(paste0("^",group2names,"$", collapse = "|"), colnames(imp_df_int))
      
      rel_idxs <- c(group2_idx, group1_idx)
      relevant_group_names <- c(rep(group2, length(group2_idx)), rep(group1, length(group1_idx)) )
      comparison <- subset(imp_df_int, select = rel_idxs)
      
      limmaResults <- NULL #data.frame(row.names = rownames(imp_df_int) )

      if (length(group1_idx) < 2 | length(group2_idx) < 2) {
        logFC <- subset(imp_df_int, select = group1_idx) - subset(imp_df_int, select = group2_idx)
        
        AveExpr <-
          (subset(imp_df_int, select = group1_idx) + 
             subset(imp_df_int, select = group2_idx)) / 2
        
        tmp <- data.frame(logFC, AveExpr)

        colnames(tmp) <- c("logFC", "AveExpr")
        limmaResults <- rbind(limmaResults, tmp)
      } else {
        
        class = factor(relevant_group_names, levels=c(group2, group1))
        
        
        design = model.matrix(~0+class) # fitting without intercept
        colnames(design) <- c(group2, group1)

        fit1 = lmFit(comparison, design = design)
        contrastNames <- c(paste0(group1,"-",group2))
        cont <- makeContrasts(contrasts=contrastNames, levels = colnames(design))
        #cont <- makeContrasts(contrasts=eval(paste0(group1,"-",group2)), levels = colnames(design))
        fit2 = contrasts.fit(fit1,contrasts = cont)
        fit3 <- eBayes(fit2)
        
        limmaResults <- NULL
        
        if (!is.null(pep.count.table) ){
          fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
          fit4 = spectraCounteBayes(fit3)
          limmaResults = outputResult(fit4,coef_col = 1)
          #limmaResults <- limmaResults[order(as.numeric( rownames(limmaResults) )), ]
          limmaResults <- limmaResults[rownames(imp_df_int),]
          limmaResults$P.Value <- limmaResults$sca.P.Value
          limmaResults$adj.P.Val <- limmaResults$sca.adj.pval
          limmaResults$sca.P.Value <- NULL
          limmaResults$sca.adj.pval <- NULL

        } else {
          limmaResults <- topTable(
            fit3,
            coef = 1,
            number = Inf,
            adjust = "BH",
            sort.by = "none"
          )
        }
      }
      
      colnames(limmaResults) <-
        paste0(colnames(limmaResults), "_", group1, "__vs__", group2)
      
      if (idx < 2) {
        df <- limmaResults
      } else {
        df <- cbind(df, limmaResults)
      }
      idx <- idx + 1
    }
    
    rownames(df) <- rownames(imp_df_int)
    return(df)
  }


getPepCountTable <- function(filtData) {
  countIdx <- grep(razorUniqueCount, colnames(filtData ) )
  pep.count.table = data.frame(count = apply(filtData[,countIdx],1,FUN=min),
                               row.names = rownames(filtData) )
  pep.count.table$count = pep.count.table$count+1
  return(pep.count.table)
}

generateEnrichedMatrix <-
  function(data,
           enrichmentChoice,
           sigCutoffValue = padjPrefix,
           fcTreshold = 1,
           pvalThresh = 0.05) {
    if (sigCutoffValue == "p-value") {
    #pvalVar <- "P.Value"
    pvalVar <- pvalPrefix
  } else if (sigCutoffValue == "adj.p-value") {
    #pvalVar <- "adj.P.Val"
    pvalVar <- padjPrefix
  } else {
    pvalVar <- 'none'
  }
  
    
  pvalThresh <- ifelse(is.null(pvalThresh), 0.05, pvalThresh)
  fcIdx <- grep(logfcPrefix, colnames(data) )
  pvalIdx <- grep(pvalVar, colnames(data) )
  
  pilot <- ifelse(length(pvalIdx) < 1, TRUE, FALSE)
  
  if (pvalVar == 'none') pilot <- TRUE
  
  comps <- colnames(data)[fcIdx]
  # comps <- gsub("logFC_","", comps)
  comps <- gsub(logfcPrefix,"", comps)
  
  out <- data.frame(dummy = row.names(data))
  
  for (i in seq(1, length(fcIdx) ) ) {
    fc <- fcIdx[i]
    pval <- pvalIdx[i]
    comp <- as.character(comps[i])
    data$tmp <- 0

    try(if (pilot) {
      if (enrichmentChoice == "enriched") {
        data[data[fc] >= fcTreshold & !is.na(data[fc]) ,]$tmp <- 1
      } else if (enrichmentChoice == "absolute") {
        data[abs(data[fc]) >= fcTreshold & !is.na(data[fc]),]$tmp <- 1
      } else if (enrichmentChoice == "reduced") {
        data[data[fc] <= -fcTreshold & !is.na(data[fc]),]$tmp <- 1
      }
    } else {
      if (enrichmentChoice == "enriched") {
        data[data[fc] >= fcTreshold &
               !is.na(data[fc]) &
               data[pval] <= pvalThresh & !is.na(data[pval]), ]$tmp <- 1
      } else if (enrichmentChoice == "absolute") {
        data[abs(data[fc]) >= fcTreshold &
               !is.na(data[fc]) &
               data[pval] <= pvalThresh & !is.na(data[pval]), ]$tmp <- 1
      } else if (enrichmentChoice == "reduced") {
        data[data[fc] <= -fcTreshold &
               !is.na(data[fc]) &
               data[pval] <= pvalThresh & !is.na(data[pval]), ]$tmp <- 1
      }
    }
    , silent = T
    
    )

    out[[comp]] <- data$tmp
  }
  data$tmp <- NULL
  out$dummy <- NULL
  row.names(out) <- row.names(data)
  out <- out[apply(out, 1, function(x) sum(x, na.rm = T)>0),, drop=F]
  return(out)
}


getComparisonsData2 <-
  function(dataLimma,
          enrichmentMatrix,
           setChoice,
           comparisons) {

    if (length(comparisons) < 2) {
       df <-
         enrichmentMatrix[enrichmentMatrix[[comparisons]] == 1,, drop=FALSE]

    } else {
      
      idx <- grep(paste0(comparisons, collapse = "|"),
                  colnames(enrichmentMatrix))
      tmpMat <-
        enrichmentMatrix[, idx] #, drop=F]
      if (setChoice == "union") {
        df <-
          tmpMat[apply(tmpMat, 1, function(x)
            any(x == 1)),]
      } else if (setChoice == "intersection") {
        df <-
          tmpMat[apply(tmpMat, 1, function(x)
            all(x == 1)),]
      } else if (setChoice == "specific" & !is.null(comparisons)) {
        specificRows <- rownames(tmpMat[apply(tmpMat, 1, function(x)
          all(x == 1)),])
        
        tmpMat <- enrichmentMatrix[specificRows, -idx]
        tmpMat <- tmpMat[apply(tmpMat, 1, function(x)
          all(x == 0)),]
        df <- enrichmentMatrix[rownames(tmpMat), idx]
      }
    }
    
    keepCols <-
      c(
        paste0('^', geneName, '$'),
        paste0(logfcPrefix, comparisons),
        paste0(padjPrefix, comparisons),
        paste0(pvalPrefix, comparisons)
      )
    
    
    return(dataLimma[rownames(df),
                     grep(paste0(keepCols, collapse = "|"),
                          colnames(dataLimma))])
  }

validateFile <- function(inFile, columns) {
  if (is.null(inFile)) {
    stop(paste0("Uploaded file ", inFile$name, " is empty."))
  }
    
  isTabDelim <- "\t" %in% strsplit(readLines(inFile$datapath, n = 1), split='')[[1]]
  
  if (!isTabDelim) {
    stop(paste0("Uploaded file ", inFile$name, "\n",
                "is not tab-delimeted.\n",
                "Please provide a tab-delimited file."))
  }
  
  colNames <- strsplit(readLines(inFile$datapath, n = 1), split='\t')[[1]]
  
  wrongColNames <- ""
  for (column in columns) {
    if (length(grep(column, colNames)) < 1) {
      wrongColNames <- paste0(
        wrongColNames,
        paste0(
          "No column named \"",
          column,
          "\" found in ",
          inFile$name,
          "\nPlease provide a column named ",
          "\"",
          column,
          "\"",
          " in ",
          inFile$name,
          ".\n"
        )
      )
    }
  }
  if (wrongColNames != "") {
    stop(wrongColNames)
  }
  
  data <-
    read.delim(
      inFile$datapath,
      stringsAsFactors = F,
      header = TRUE,
      sep = "\t"
    )
  return(data)
}


amicaOutput <- function(protData, comparisons) {
    df <- as.data.frame(rowData(protData))
    snames <- grep("Gene.names", names(comparisons), invert = T, value = T)
    
    df <- merge(df, comparisons[,snames], by=0, all.x=TRUE)
    row.names(df) <- df$Row.names
    df <- subset(df,select=-c(Row.names))
    
    tmp <- data.frame(row.names = row.names(rowData(protData)))
    for (elem in assayNames(protData)) {
      vals <- assay(protData, elem)
      names(vals) <- paste0(elem,"_",names(vals))
      tmp <- cbind(tmp, vals)
    }

    df <- merge(df, tmp, by=0)
    row.names(df) <- df$Row.names
    df <- subset(df,select=-c(Row.names))
    return(df)
  }


mergeAmicas <-
  function(exp1,
           exp2,
           mergeKey,
           splitDelim = NULL,
           exp1.prefix = "_exp1",
           exp2.prefix = "_exp2") {
    key1 <- 0
    key2 <- 0
    if (mergeKey == 'protein') {
      key1 <- row.names(exp1)
      key2 <- exp2$Majority.protein.IDs
    } else {
      key1 <- toupper(exp1$Gene.names)
      key2 <- toupper(exp2$Gene.names)
    }
    
    if (!is.null(splitDelim) && splitDelim != "") {
      key1 <- gsub(paste0(splitDelim,".*"),"",key1)
      key2 <- gsub(paste0(splitDelim,".*"),"",key2)
    }
    names(exp1) <- paste0(names(exp1), exp1.prefix)
    names(exp2) <- paste0(names(exp2), exp2.prefix)
    
    exp1$merge.key <- key1
    exp2$merge.key <- key2
    
    merged.exp <- merge(exp1, exp2, by='merge.key', all=T)
    
    genes.exp1 <- paste0("Gene.names", exp1.prefix)
    genes.exp2 <- paste0("Gene.names", exp2.prefix)
    
    merged.exp$Gene.names <- ifelse(!is.na(merged.exp[[genes.exp1]]),
                                    merged.exp[[genes.exp1]],
                                    merged.exp[[genes.exp2]])
    
    merged.exp <- merged.exp[,c(ncol(merged.exp),1:(ncol(merged.exp)-1))]
    
    return(merged.exp)
  }


#### visNetwork input
toNetworkData <- function(df.genes, ppi, cellmap) {
  
  if(is.null(df.genes) | nrow(df.genes) < 1) {
    return(NULL)
  }
  
  df.genes$Gene.names <- toupper(df.genes$Gene.names)
  df.genes$key <- df.genes$Gene.names
  df.genes$Gene.names <- gsub(";.*", "", df.genes$Gene.names)
  df.genes <- df.genes[, grep("Gene|logFC|key", colnames(df.genes), value = T)]
  df.genes <- df.genes[order(df.genes[[geneName]], decreasing=TRUE),]
  df.genes <- df.genes[!duplicated(df.genes[[geneName]]),]
  
  #idxs <- match(df.genes[[geneName]], V(ppi)$label)
  idxs <- match(df.genes[[geneName]], V(ppi)$name)
  
  idxs <- idxs[!is.na(idxs)]
  
  if(length(idxs) < 1) {
    return(NULL)
  }
  
  g <- induced_subgraph(ppi, idxs)
  df <- toVisNetworkData(g, idToLabel = T)
  df$nodes <-
    merge(df$nodes, df.genes,
          by.x = "label", by.y = geneName)
  colnames(df$nodes) <- c("label", "id", "log2FC", "key")
  colnames(df$edges) <- c("from", "to", "value")
  
  df$nodes <-
    merge(df$nodes, cellmap, by="label", all.x = T)
  
  # convert id to integer for cytoscape  
  df$nodes$id <- 1:nrow(df$nodes)
  tmpFrom <- c()
  tmpTo <- c()
  for (idx in 1:nrow(df$edges)) {
    from <- df$edges[idx, 'from']
    to <- df$edges[idx, 'to']
    fromId <- df$nodes[df$nodes$label == from, 'id']
    toId <- df$nodes[df$nodes$label == to, 'id']
    
    tmpFrom[idx] <- fromId
    tmpTo[idx] <- toId
  }
  
  if (!is.na(tmpFrom)) {
    df$edges$from <- tmpFrom
    df$edges$to <- tmpTo
  }
  
  return(df)
}

readInAmicaSumm <- function(inFilePath, design) {
  protData <- read.delim(inFilePath, header=T, sep="\t", stringsAsFactors = F)
  
  if ("Intensity" %in% names(protData)) {
    protData$Intensity <- NULL
  } 
  if ("iBAQ" %in% names(protData)) {
    protData$iBAQ <- NULL
  } 
  if ("iBAQ.peptides" %in% names(protData)) {
    protData$iBAQ.peptides <- NULL
  }
  if (!(proteinId %in% names(protData))) {
    stop(paste0("Error: ", 
                proteinId, " is not a column name of the uploaded input data."))
  }
  if (!(geneName %in% names(protData))) {
    stop(paste0("Error: ", 
                geneName, " is not a column name of the uploaded input data."))
  }
  if (length(grep("ImputedIntensity", names(protData)) ) < 1 ) {
    stop(paste0("Error: There are no ImputedIntensity columns ", 
                " in the uploaded input data."))
  }
  if (length(grep("__vs__", names(protData)) ) < 1 ) {
    stop(paste0("Error: There are no '__vs__' - infix columns ", 
                " in the uploaded input data."))
  }
  if (length(grep("logFC", names(protData)) ) < 1 ) {
    stop(paste0("Error: There are no 'logFC' - columns ", 
                " in the uploaded input data."))
  }
  
  rownames(protData) <- protData[[proteinId]]
  
  comps <- grep("logFC", colnames(protData), value=T)
  comps <- gsub("logFC_", "", comps)
  
  group1 <- c()
  group2 <- c()
  
  for (elem in strsplit(comps, "__vs__")) {
    group1 <- c(group1,elem[1])
    group2 <- c(group2,elem[2]) 
  }
  contrasts <- data.frame(group1, group2)
  
  if ("Intensity" %in% names(protData)) {
    idx <- which(names(protData)=="Intensity")
    protData[idx] <- NULL
  }
  
  tmp <- grep("Intensity", names(protData), value = T)
  prefixes <- unique(gsub("Intensity.*", "Intensity", tmp))
  
  assayList <- list()
  intIdxs <- c()
  for (i in seq_along(prefixes)) {
    gpref <- paste0("^",prefixes[i],"_" )
    filtIdxs <- grep(gpref, names(protData))
    intIdxs <- c(intIdxs, filtIdxs)
    tmpDf <- protData[,filtIdxs]
    names(tmpDf) <- gsub(gpref, "", names(tmpDf) )
    assayList[[i]] <- tmpDf
  }
  names(assayList) <- prefixes
  
  if (length(grep("iBAQ", names(protData)) ) > 0 ) {
    if ("iBAQ.peptides" %in% names(protData)) protData$iBAQ.peptides <- NULL
    if ("iBAQ" %in% names(protData)) protData$iBAQ <- NULL
    idxs <- grep("iBAQ", names(protData))
    intIdxs <- c(intIdxs, idxs)
    ibaqs <- protData[, idxs]
    names(ibaqs) <- gsub("iBAQ.", "", names(ibaqs))
    assayList[[length(prefixes)+1]] <- ibaqs
    names(assayList)[length(prefixes)+1] <- "iBAQ"
  }
  
  dropIdxs <- c(intIdxs, grep("__vs__", colnames(protData)))
  
  if ("quantified" %in% names(protData)) {
    invisible()
  } else {
    rows <- which(complete.cases( protData[, grep("__vs__|^ImputedIntensity",
                                                  colnames(protData))  ] ))
    protData$quantified <- ""
    protData[rows, 'quantified'] <- "+"
  }

  se <- ProteomicsData(assays = assayList,
                       rowData = protData[,-dropIdxs],
                       colData = design)
  
  comparisons <-
    protData[protData$quantified == "+",
             c(geneName, grep("__vs__", colnames(protData), value = T))]
  
  return(
    list(
      protData = se,
      comparisons = comparisons,
      contrasts = contrasts
    )
  )
}


readInMQproteinGroupsSumm <- function(mqFile, design) {
  protData <- read.delim(mqFile, sep = "\t", header = T, stringsAsFactors = F)
  
  if ("Intensity" %in% names(protData)) {
    protData$Intensity <- NULL
  } 
  if ("iBAQ" %in% names(protData)) {
    protData$iBAQ <- NULL
  } 
  if ("iBAQ.peptides" %in% names(protData)) {
    protData$iBAQ.peptides <- NULL
  }
  if (!("Majority.protein.IDs" %in% names(protData))) {
    stop(paste0("Error: ", 
                "Majority.protein.IDs", " is not a column name of the uploaded input data."))
  }
  if (!(geneName %in% names(protData))) {
    protData[[geneName]] <- ""
  }

  row.names(protData) <- protData$Majority.protein.IDs
  
  if ("iBAQ.peptides" %in% names(protData)) protData$iBAQ.peptides <- NULL
  if ("iBAQ" %in% names(protData)) protData$iBAQ <- NULL
  if ("Intensity" %in% names(protData)) protData$Intensity <- NULL
  
  
  ##
  names(protData)[which(names(protData)=="MS.MS.count")] <- "spectraCount"
  names(protData)[which(names(protData)=="Razor...unique.peptides")] <- "razorUniqueCount"
  idx <-  grep("Razor...unique.peptides.", names(protData))
  names(protData)[idx] <- gsub("Razor...unique.peptides.", "", names(protData)[idx])
  names(protData)[idx] <- paste0("razorUniqueCount.", names(protData)[idx])
  ##
  
  
  rBool <- FALSE
  sBool <- FALSE
  if ("Reverse" %in% names(protData)) rBool <- TRUE
  if ("Only.identified.by.site" %in% names(protData)) sBool <- TRUE
  
  if(rBool & sBool) {
    protData <-
      protData[protData$Reverse != "+" &
                 protData$Only.identified.by.site != "+", ]
  } else if(rBool) {
    protData <-
      protData[protData$Reverse != "+", ]
  } else if (sBool) {
    protData <-
      protData[protData$Only.identified.by.site != "+", ]
  }
  protData$Reverse <- NULL
  protData$Only.identified.by.site <- NULL
  
  
  
  protData[[geneName]] <- ifelse(
    is.na(protData[[geneName]]) |protData[[geneName]] == "",
    as.character(protData[[proteinId]]),
    as.character(protData[[geneName]])
  )
  
  # int_idx <- grep(intensityPrefix, colnames(protData))
  # ibaq_idx <- grep(abundancePrefix, colnames(protData))
  int_idx <- grep("^LFQ.intensity", colnames(protData))
  ibaq_idx <- grep("^iBAQ.", colnames(protData))
  raw_idx <- grep("^Intensity", names(protData))
  
  assayList <- list()
  assayIdx <- 1
  
  raws <- protData[,raw_idx]
  names(raws) <- gsub("Intensity.", "", names(raws))
  raws[raws==0] <- NA
  raws <- log2(raws)
  
  assayList[[assayIdx]] <- raws
  names(assayList) <- c("RawIntensity")
  assayIdx <- assayIdx + 1
  
  if (length(ibaq_idx) > 0) {
    ibaqs <- protData[,ibaq_idx]
    names(ibaqs) <- gsub("iBAQ.", "", names(ibaqs))
    ibaqs[ibaqs==0] <- NA
    ibaqs <- log2(ibaqs)
    
    assayList[[assayIdx]] <- ibaqs
    names(assayList) <- c("RawIntensity", "iBAQ")
    assayIdx <- assayIdx + 1
  }
  
  if (length(int_idx) > 0) {
    lfqs <- protData[, int_idx]
    names(lfqs) <- gsub("LFQ.intensity.", "", names(lfqs))
    lfqs[lfqs == 0] <- NA
    lfqs <- log2(lfqs)
    
    assayList[[assayIdx]] <- lfqs
    names(assayList)[assayIdx] <- "LFQIntensity"
    assayIdx <- assayIdx + 1
    #assayList[[assayIdx]] <- lfqs
    #names(assayList)[assayIdx] <- "ImputedIntensity"
  }
  
  dropIdx <-
    grep(
      'Majority|Gene.names|Potential|spectraCount|razorUniqueCount',
      colnames(protData),
      invert = T
    )
  protData <- protData[,-dropIdx]


  se <- ProteomicsData(assays = assayList,
                             rowData = protData,
                             colData = design)
  return(se)
}




readInFragPipeProteinGroupsSumm <- function(mqFile, design) {
  protData <- read.delim(mqFile, sep = "\t", header = T, stringsAsFactors = F)
  row.names(protData) <- protData$Protein

  spCount <- "Summarized.Razor.Spectral.Count"
  pepCount <- "Unique.Stripped.Peptides"
  sampleSpCountPostfix <- ".Razor.Spectral.Count"
  gName <- ifelse("Gene.Names" %in% names(protData), "Gene.Names", "Gene")
  
  if (length(grep("Summarized", names(protData))) < 1) {
    spCount <- "Combined.Spectral.Count"
    pepCount <- "Combined.Total.Peptides"
    sampleSpCountPostfix <- ".Spectral.Count"
  }

  protData[[contaminantCol]] <- ""
  
  if (!("Protein" %in% names(protData))) {
    stop(paste0("Error: ", 
                "Protein", " is not a column name of the uploaded input data."))
  }
  
  names(protData)[which(names(protData)=="Protein")] <- "Majority.protein.IDs"
  names(protData)[which(names(protData)==gName)] <- "Gene.names"
  names(protData)[which(names(protData)==spCount)] <- "spectraCount"
  names(protData)[which(names(protData)==pepCount)] <- "razorUniqueCount"

  protData$Gene.names <- ifelse(
    protData$Gene.names == "",
    protData$Majority.protein.IDs,
    protData$Gene.names
  )
  
  query <- paste0(design$samples, sampleSpCountPostfix, collapse = "|")
  idx <-  grep(query, names(protData))

  names(protData)[idx] <- gsub(sampleSpCountPostfix, "", names(protData)[idx])
  names(protData)[idx] <- paste0("razorUniqueCount.", names(protData)[idx])
  
  samples <- make.names(design$samples)
  
  
  unique_idx <- grep(paste0(samples, ".Unique.Intensity", collapse = "|"), colnames(protData))
  tot_idx <- grep(paste0(samples, ".Total.Intensity", collapse = "|"), colnames(protData))
  razor_idx <- grep(paste0(samples, ".Razor.Intensity", collapse = "|"), names(protData))
  ints_idx <- grep(paste0(samples, ".Intensity", collapse = "|"), names(protData))
  lfq_idx <- grep(paste0(samples, ".MaxLFQ.Intensity", collapse = "|"), names(protData))
  lfq_unique_idx <- grep(paste0(samples, ".MaxLFQ.Unique.Intensity", collapse = "|"), names(protData))
  lfq_unique_tot_idx <- grep(paste0(samples, ".MaxLFQ.Total.Intensity", collapse = "|"), names(protData))
  
  uniqs <- NULL
  tots <- NULL
  razor <- NULL
  ints <- NULL
  lfqInts <- NULL
  lfqUniqueInts <- NULL
  lfqTotalInts <- NULL
  
  assayList <- list()
  
  if (length(tot_idx) > 0) {
    tots <- protData[,tot_idx]
    names(tots) <- gsub(".Total.Intensity", "", names(tots))
    tots[tots==0] <- NA
    tots <- log2(tots)
    assayList[['TotalIntensity']] <- tots
  }
  
  if (length(unique_idx) > 0) {
    uniqs <- protData[,unique_idx]
    names(uniqs) <- gsub(".Unique.Intensity", "", names(uniqs))
    uniqs[uniqs==0] <- NA
    uniqs <- log2(uniqs)
    assayList[['UniqueIntensity']] <- uniqs
  }
  
  if (length(razor_idx) > 0) {
    razor <- protData[,razor_idx]
    names(razor) <- gsub(".Razor.Intensity", "", names(razor))
    razor[razor==0] <- NA
    razor <- log2(razor)
  }
  
  if (length(ints_idx) > 0) {
    ints <- protData[,ints_idx]
    names(ints) <- gsub(".Intensity", "", names(ints))
    ints[ints==0] <- NA
    ints <- log2(ints)
  }
  
  # MaxLFQ intensities
  if (length(lfq_idx) > 0) {
    lfqInts <- protData[,lfq_idx]
    names(lfqInts) <- gsub(".MaxLFQ.Intensity", "", names(lfqInts))
    lfqInts[lfqInts==0] <- NA
    lfqInts <- log2(lfqInts)
  }
  
  if (length(lfq_unique_idx) > 0) {
    lfqUniqueInts <- protData[,lfq_unique_idx]
    names(lfqUniqueInts) <- gsub(".MaxLFQ.Unique.Intensity", "", names(lfqUniqueInts))
    lfqUniqueInts[lfqUniqueInts==0] <- NA
    lfqUniqueInts <- log2(lfqUniqueInts)
  }
  
  if (length(lfq_unique_tot_idx) > 0) {
    lfqTotalInts <- protData[,lfq_unique_tot_idx]
    names(lfqTotalInts) <- gsub(".MaxLFQ.Total.Intensity", "", names(lfqTotalInts))
    lfqTotalInts[lfqTotalInts==0] <- NA
    lfqTotalInts <- log2(lfqTotalInts)
  }
  
  dropIdx <-
    c(
      unique_idx,
      razor_idx,
      tot_idx,
      ints_idx,
      lfq_idx,
      lfq_unique_idx,
      lfq_unique_tot_idx
    )
  dropIdx <- c(dropIdx, grep("Spectral.Count", names(protData)) )

  if (!is.null(razor)) {
    assayList[["RazorIntensity"]] <- razor
  }
  if (!is.null(ints)) {
    assayList[["Intensity"]] <- ints
  }
  if (!is.null(lfqInts)) {
    assayList[["LFQIntensity"]] <- lfqInts
  }
  if (!is.null(lfqUniqueInts)) {
    assayList[["LFQUniqueIntensity"]] <- lfqUniqueInts
  }
  if (!is.null(lfqTotalInts)) {
    assayList[["LFQTotalIntensity"]] <- lfqTotalInts
  }
  
  if (length(assayList) < 1) {
    stop(paste0("Error: ", 
                "There are intensities in the uploaded FragPipe file."))
  }
  
  se <- ProteomicsData(
    assays = assayList,
    rowData = protData[,-dropIdx],
    colData = design
  )
  return(se)
}


#### specific input

readInCustomSumm <- function(mqFile, specs, design, logtransform=TRUE) {
  protData <- read.delim(mqFile, sep = "\t", header = T, stringsAsFactors = F)
  
  specs$Pattern <- make.names(specs$Pattern, unique = FALSE)
  
  proteinId <- specs$Pattern[specs$Variable=="proteinId"]
  geneName <- specs$Pattern[specs$Variable=="geneName"]
  intensityPrefix <- specs$Pattern[specs$Variable=="intensityPrefix"]
  abundancePrefix <- specs$Pattern[specs$Variable=="abundancePrefix"]
  razorUniqueCount <- specs$Pattern[specs$Variable=="razorUniqueCount"]
  razorUniqueCountPrefix <- specs$Pattern[specs$Variable=="razorUniqueCountPrefix"]
  spectraCount <- specs$Pattern[specs$Variable=="spectraCount"]
  pot.contaminant <- specs$Pattern[specs$Variable=="contaminantCol"]
  
  if (!isTruthy(proteinId) || !(proteinId %in% names(protData))) {
    stop(paste0("Error: Provided proteinId in specification file = ", 
                proteinId, " is not a column name of the uploaded input data."))
  }
  
  if (!isTruthy(geneName) || !(geneName %in% names(protData))) {
    stop(paste0("Error: Provided geneName in specification file = ", 
                geneName, " is not a column name of the uploaded input data."))
  }
  
  if (!isTruthy(intensityPrefix) || length(grep(intensityPrefix, names(protData))) < 1 ) {
    stop(paste0("Error: Provided intensityPrefix in specification file = 
                ", intensityPrefix, " is not a column name of the uploaded input data."))
  }
  
  names(protData)[which(names(protData)==proteinId)] <- "Majority.protein.IDs"
  names(protData)[which(names(protData)==geneName)] <- "Gene.names"
  
  if (isTruthy(spectraCount) && spectraCount %in% names(protData)) {
    names(protData)[which(names(protData)==spectraCount)] <- "spectraCount"
  }
  if (isTruthy(razorUniqueCount) && razorUniqueCount %in% names(protData)) {
    names(protData)[which(names(protData)==razorUniqueCount)] <- "razorUniqueCount"
  }
  
  if (isTruthy(razorUniqueCountPrefix) && length(grep(razorUniqueCountPrefix, names(protData))) > 1) {
    idx <-  grep(razorUniqueCountPrefix, names(protData))
    
    names(protData)[idx] <- gsub(razorUniqueCountPrefix, "", names(protData)[idx])
    names(protData)[idx] <- paste0("razorUniqueCount.", names(protData)[idx])
  }
  
  if (isTruthy(pot.contaminant) && pot.contaminant  %in% names(protData)) {
    names(protData)[which(names(protData)==pot.contaminant)] <- "Potential.contaminant"
  } else {
      protData$Potential.contaminant <- NULL
  }

  row.names(protData) <- protData$Majority.protein.IDs
  
  protData$Gene.names <- ifelse(
    protData$Gene.names == "",
    protData$Majority.protein.IDs,
    protData$Gene.names
  )
  
  int_idx <- grep(intensityPrefix, names(protData))
  razor <- protData[,int_idx]
  names(razor) <- gsub(intensityPrefix, "", names(razor))
  razor[razor==0] <- NA
  if (logtransform) {
    razor <- log2(razor)
  }

  intsList <- list(LFQIntensity=razor, ImputedIntensity=razor)
  dropIdx <- int_idx
  
  if (length(abundancePrefix) > 0) {
    tot_idx <- grep(abundancePrefix, colnames(protData))
    
    tots <- protData[,tot_idx]
    names(tots) <- gsub(abundancePrefix, "", names(tots))
    tots[tots==0] <- NA
    if (logtransform) {
      tots <- log2(tots)
    }
    intsList[[3]] <- tots
    names(intsList)[3] <- "iBAQ"
    dropIdx <- c(dropIdx, tot_idx)
  }
  
  data <- protData[, -dropIdx]
  names(data) <- gsub('Intensity', 'intensity', names(data)) 
  se <- ProteomicsData(assays = intsList,
                       rowData = data,
                       #rowData = protData[, -dropIdx],
                       colData = design)
  
  return(se)
}


toLongFormat <- function(object,
                         se,
                         addGroup = TRUE,
                         addContaminant = FALSE,
                         addGeneName = FALSE
                         ) {
  object <- stats::reshape(object, idvar = "rowname",
                           ids = rownames(object), times = names(object),
                           timevar = "colname", varying = list(names(object)),
                           direction = "long", v.names = "value")
  if (addGroup) {
    midx <- match(object$colname, colData(se)$samples)
    object$group <- colData(se)$groups[midx]
  }
  if (addContaminant) {
    rnames <- which(rowData(se)$Potential.contaminant == "+")
    #which(rowData(se)[rowData(se)$Potential.contaminant == "+", ])
    
    object$Contaminant <- "no"
    object[object$rowname %in% rnames, "Contaminant"] <- "yes"
  }
  if (addGeneName) {
    object$Gene.name <- rowData(se)[object$rowname, geneName]
    # object$Gene.name <- rowData(se[object$rowname,])[[geneName]]
  }
  
  return(object)
}


#### VIS NETWORK BAIT TO PREY

getBait2PreyNetwork <- function(dataComp, matrixSet, samples, palette) {
  tmp <- matrixSet[row.names(matrixSet) %in% 
                     row.names(dataComp),samples, drop = F]
  
  tmp$Gene.names <- dataComp[row.names(tmp), "Gene.names"]
  key <- tmp$Gene.names
  tmp$Gene.names <- toupper(tmp$Gene.names)
  tmp$Gene.names <- gsub(";.*", "", tmp$Gene.names)
  
  tmp <- tmp[, c(ncol(tmp), 1:(ncol(tmp)-1))]
  
  nodes <-
    data.frame(
      id = 1:(ncol(tmp) - 1),
      label = names(tmp)[2:ncol(tmp)],
      key = rep('', ncol(tmp) - 1),
      group = rep('Comparison', ncol(tmp) - 1),
      shape = rep('diamond', ncol(tmp) - 1),
      color = rep(palette[1], ncol(tmp) - 1)
    )
  df <-
    data.frame(
      id = 1:nrow(tmp),
      label = tmp$Gene.names,
      key = key,
      group = rep('Protein', nrow(tmp)),
      #shape = rep('circle', nrow(tmp)),
      shape = rep(NA, nrow(tmp)),
      color = rep(palette[2], nrow(tmp))
    )
  df$id <- df$id + max(nodes$id)
  nodes.df <- rbind(nodes,  df)
  
  edges.df <- data.frame()
  for(bait in names(tmp)[2:ncol(tmp)]) {
    node.id <- nodes$id[nodes$label==bait]
    vec <- ifelse(tmp[[bait]]==1, nodes.df$id[nodes.df$label %in% tmp$Gene.names], NA)
    vec <- vec[!is.na(vec)]
    bait.df <-
      data.frame(
        from = rep(node.id, length(vec)),
        to = vec,
        color = rep(palette[1], length(vec)),
        interaction = rep('experimental', length(vec))
      )
    edges.df <- rbind(edges.df, bait.df)
  }
  edges.df$value <- 0.33
  return(list(nodes.df, edges.df))
}


dfs2gml <- function(nodes, edges) {
  outstring <- 'graph [\n  directed 0\n'
  for (idx in 1:nrow(nodes)) {
    t <- nodes[idx,]
    outstring <- paste0(outstring, '  node [\n')
    
    for (elem in names(nodes)) {
      attribute <- elem
      if (length(grep("NMF localization", elem))>0) {
        attribute <- "CellMap_NMF"
      } else if (length(grep("SAFE localization", elem))>0) {
        attribute <- "CellMap_SAFE"
      } else if (length(grep("Subcell. localization", elem))>0) {
        attribute <- "CellMap_NMF_SUMMARIZED"
      }
      
      if (is(t[[elem]] ,'character') | is.na(t[[elem]]) ) {
        outstring <- paste0(outstring, '\t', attribute, ' "', t[[elem]], '"\n')
      } else {
        outstring <- paste0(outstring, '\t', attribute, ' ', t[[elem]], '\n')
      }
    }
    outstring <- paste0(outstring, '\n  ]\n')

  }
  for (idx in 1:nrow(edges)) {
    t <- edges[idx,]
    outstring <- paste0(outstring, '  edge [\n')
    
    for (elem in names(edges)) {
      
      attribute <- elem
      if (elem == "from") {
        attribute <- "source"
      } else if (elem == "to") {
        attribute <- "target"
      }
      
      if (is(t[[elem]] ,'character') | is.na(t[[elem]]) ) {
        outstring <- paste0(outstring, '\t', attribute, ' "', t[[elem]], '"\n')
      } else {
        outstring <- paste0(outstring, '\t', attribute, ' ', t[[elem]], '\n')
      }
    }
    outstring <- paste0(outstring, '\n  ]\n')

  }
  outstring <- paste0(outstring, ']')
  return(outstring)
}


### ID PROTEIN OVERLAP
# https://blog.jdblischak.com/posts/pairwise-overlaps/
calc_pairwise_overlaps <- function(sets) {
  # Ensure that all sets are unique character vectors
  sets_are_vectors <- vapply(sets, is.vector, logical(1))
  if (any(!sets_are_vectors)) {
    stop("Sets must be vectors")
  }
  sets_are_atomic <- vapply(sets, is.atomic, logical(1))
  if (any(!sets_are_atomic)) {
    stop("Sets must be atomic vectors, i.e. not lists")
  }
  sets <- lapply(sets, as.character)
  is_unique <- function(x) length(unique(x)) == length(x)
  sets_are_unique <- vapply(sets, is_unique, logical(1))
  if (any(!sets_are_unique)) {
    stop("Sets must be unique, i.e. no duplicated elements")
  }
  
  n_sets <- length(sets)
  set_names <- names(sets)
  n_overlaps <- choose(n = n_sets, k = 2)
  
  vec_name1 <- character(length = n_overlaps)
  vec_name2 <- character(length = n_overlaps)
  vec_num_sample1 <- integer(length = n_overlaps)
  vec_num_sample2 <- integer(length = n_overlaps)
  vec_num_shared <- integer(length = n_overlaps)
  vec_overlap <- numeric(length = n_overlaps)
  vec_jaccard <- numeric(length = n_overlaps)
  overlaps_index <- 1
  
  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      name2 <- set_names[j]
      set2 <- sets[[j]]
      
      set_intersect <- set1[match(set2, set1, 0L)]
      set_union <- .Internal(unique(c(set1, set2), incomparables = FALSE,
                                    fromLast = FALSE, nmax = NA))
      num_shared <- length(set_intersect)
      overlap <- num_shared / min(length(set1), length(set2))
      jaccard <- num_shared / length(set_union)
      
      vec_name1[overlaps_index] <- name1
      vec_name2[overlaps_index] <- name2
      vec_num_sample1[overlaps_index] <- length(set1)
      vec_num_sample2[overlaps_index] <- length(set2)
      vec_num_shared[overlaps_index] <- num_shared
      vec_overlap[overlaps_index] <- overlap
      vec_jaccard[overlaps_index] <- jaccard
      
      overlaps_index <- overlaps_index + 1
    }
  }
  
  result <- data.frame(sample1 = vec_name1,
                       sample2 = vec_name2,
                       num_sample1 = vec_num_sample1,
                       num_sample2 = vec_num_sample2,
                       num_shared = vec_num_shared,
                       overlap = round(vec_overlap, 3),
                       jaccard = round(vec_jaccard, 3)
                       )
  return(result)
}


