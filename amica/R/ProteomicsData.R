### This class emulates SummarizedExperiment (SE) because
### SE takes too much RAM (>350MB)
### and takes a longer time to load

setClass(
  "ProteomicsData",
  representation(
    assays = "list",
    rowData = "data.frame",
    colData = "data.frame"
    )
  )

.valid.assays_nrow <- function(x, y) {
  if (length(x) < 2 )
    return(NULL)
  
  lens <- unlist(lapply(x, nrow))
  variance <- var(lens)
  
  if (variance != 0) {
    txt <- sprintf(
      "\n  nb of rows in 'assay' must be equal %s",
      paste(lens, collapse = ", ") )
    return(txt)
  }
  
  assays_nrow <- lens[1]
  
  rowData_nrow <- nrow(y)
  if (assays_nrow != rowData_nrow) {
    txt <- sprintf(
      "\n  nb of rows in 'assay' (%d) must equal nb of rows in 'rowData' (%d)",
      assays_nrow, rowData_nrow)
    return(txt)
  }
  NULL
}

.valid.assays_ncol <- function(x, y) {
  if (length(x) < 2)
    return(NULL)
  
  lens <- unlist(lapply(x, ncol))
  variance <- var(lens)
  
  if (variance != 0) {
    txt <- sprintf(
      "\n  nb of cols in 'assay' must be equal %s",
      paste(lens, collapse = ", ") )
    return(txt)
  }
  
  assays_ncol <- lens[1]
  ###
  
  colData_ncol <- nrow(y)
  
  if (assays_ncol != colData_ncol) {
    txt <- sprintf(
      "\n  nb of cols in 'assay' (%d) must equal nb of rows in 'colData' (%d)",
      assays_ncol, colData_ncol)
    return(txt)
  }
  NULL
}

.valid.assays_dim <- function(x,y,z) {
  c(.valid.assays_nrow(x,y),
    .valid.assays_ncol(x,z))
}

# .valid.ProteomicsData <- function(x) {
#   .valid.ProteomicsData.assays_dim(x)
# }
# 
# setValidity("ProteomicsData", .valid.ProteomicsData)



ProteomicsData <- function(assays, rowData, colData) {
  
  out <- .valid.assays_dim(assays, rowData, colData)
  
  if (!is.null(out)) {
    stop(out)
  }

  rnames <- row.names(rowData)
  cnames <- colData$samples
  
  assays <- lapply(assays, function(x){ row.names(x)<-rnames;colnames(x) <- cnames; x})
  
  new("ProteomicsData",
      rowData=rowData,
      colData=colData,
      assays=assays)
  }


#######


setMethod("show",
          "ProteomicsData",
          function(object) {
            cat("class:", class(object), "\n")
            cat("dim:", dim(object@assays[[1]]), "\n")
            cat("assays:", names(object@assays), "\n")
            cat("rowData:", row.names(object@rowData)[1:min(4, nrow(object@rowData))], "...\n")
            cat("colData samples:",object@colData$samples[1:min(4, nrow(object@colData))], "...\n")
            cat("colData groups:",object@colData$groups[1:min(4, nrow(object@colData))], "...\n")
          }
)

rowData = function (y) {
  if (is(y, 'ProteomicsData')) y@rowData
  else y
}    

colData = function (y) {
  if (is(y, 'ProteomicsData')) y@colData
  else y
}  

assays = function (y) {
  if (is(y, 'ProteomicsData')) y@assays
  else y
}

assayNames = function (y) {
  if (is(y, 'ProteomicsData')) names(y@assays)
  else y
}

assay = function (y, x) {
  if (is(y, 'ProteomicsData')) {
    if (is(x, 'numeric')) {
      y@assays[[x]]
    } else if (is(x, 'character')) {
      idx <- which(assayNames(y)==x)
      y@assays[[idx]]
    }
  }
}

setGeneric("setAssay", function(x, assay, assayName){standardGeneric("setAssay")})
setMethod("setAssay", 
          c(x = "ProteomicsData", assay = "data.frame", assayName = "character"), 
          function(x, assay, assayName){
            x@assays[[assayName]] <- assay
            return(x)})

setGeneric("setRowData", function(x, rowData){standardGeneric("setRowData")})
setMethod("setRowData", 
          c(x = "ProteomicsData", rowData = "data.frame"), 
          function(x, rowData){
            x@rowData <- rowData
            return(x)})


isQuantRnames = function(y) {
  if (is(y, 'ProteomicsData')) {
    if ('quantified' %in% names(y@rowData)) {
      return(row.names(y@rowData)[which(y@rowData$quantified == '+')])
    } else {
      return(NULL)
    }
  }
}

filterOnMinValuesRnames = function(y, minMSMS, minRazor) {
  rnames <- row.names(y@rowData)

  if (contaminantCol %in% names(y@rowData) &&
       spectraCount %in% names(y@rowData) &&
      razorUniqueCount %in% names(y@rowData)) {
    
    rnames <-
      row.names(y@rowData[
        (is.na(y@rowData[[contaminantCol]]) |
        y@rowData[[contaminantCol]] !=
                             "+" ) &
                            y@rowData[[razorUniqueCount]] >=
                             minRazor &
                            y@rowData[[spectraCount]] >=
                             minMSMS,])
  } else if ( spectraCount %in% names(y@rowData) &
             contaminantCol %in% names(y@rowData)) {
    
    rnames <-
      row.names(y@rowData[y@rowData[[contaminantCol]] !=
                             "+" &
                            y@rowData[[spectraCount]] >=
                             minMSMS,])
  } else if (razorUniqueCount %in% names(y@rowData) &&
             contaminantCol %in% names(y@rowData)) {
    
    rnames <-
      row.names(y@rowData[y@rowData[[contaminantCol]] !=
                             "+" &
                            y@rowData[[razorUniqueCount]] >=
                             minRazor,])
  } else if (contaminantCol %in% names(y@rowData) ) {

    rnames <-
      row.names(y@rowData[y@rowData[[contaminantCol]] !=
                            "+",])
  } else if (razorUniqueCount %in% names(y@rowData) ) {
    
    rnames <-
      row.names(y@rowData[y@rowData[[razorUniqueCount]] >=
                  minRazor,])
  } else if (spectraCount %in% names(y@rowData)) {
    
    rnames <-
      row.names(y@rowData[y@rowData[[spectraCount]] >=
                            minMSMS,])
  }
  return(rnames)
}
