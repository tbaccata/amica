### UPLOAD
observeEvent(input$submitAnalysis, {
  reacValues$inputParameterSummary <- NULL
  allFilesUploaded <- TRUE
  
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
        rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData),]
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
    allFilesUploaded <- FALSE
    showNotification(paste("Need to upload experimental design."), type = "message")
    return("")
  }
  
  # EXPERIMENTAL DESIGN
  if (input$source != "example" &
      !is.null(input$groupSpecification))  {
    ###
    tryCatch({
      tmpData <- validateFile(input$groupSpecification, c("groups", "samples"))
    },
    error = function(cond) {
      #message('Upload of experimental design file failed. ', cond)
      allFilesUploaded <- FALSE
      
      shinyalert(
        title = "Upload of experimental design file failed.",
        text = gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
      
      # showNotification(
      #   paste(cond),
      #   duration = 100,
      #   closeButton = T,
      #   type = "error"
      # )
      #return("")
    },
    warning = function(cond) {
      showNotification(paste(cond))
    }
    #, finally = {
    #  showNotification(paste("Succesfully uploaded experimental design."), type = "message")
    #}
    )
    ###
    
    if (exists("tmpData")) {
      if (class(tmpData$samples) == 'integer' ||
          class(tmpData$samples) == 'numeric') {
        tmpData$samples <- as.character(tmpData$samples)
      } else {
        tmpData$samples <- make.names(tmpData$samples)
      }
      
      for (group in unique(tmpData$groups)) {
        if (make.names(group) != group) {
          out <- paste0(
            "The group ",
            group,
            "is not a valid group name.\n",
            "A syntactically valid name does not contain whitespace, consists of letters, ",
            "numbers, dots or underline characters and starts with a letter.\n\n",
            "Please enter only valid group names in your experimental design."
          )
          
          shinyalert(
            title = "Please enter only valid group names in your experimental design.",
            text = paste(out),
            size = "s", 
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = FALSE,
            type = "error",
            showConfirmButton = TRUE,
            showCancelButton = FALSE,
            confirmButtonText = "OK",
            confirmButtonCol = "#669966",
            timer = 0,
            imageUrl = "",
            animation = TRUE
          )
          
          #showNotification(out, duration = NULL, type = "error")
          allFilesUploaded <- FALSE
          return("")
        }
      }
      reacValues$expDesign <- tmpData
    } else {
      return("")
    }
  }
  
  if (input$source == "amica") {
    if (is.null(input$amicaFile)) {
      showNotification(paste("Need to upload amica output."), type = "message")
      allFilesUploaded <- FALSE
      return("")
    }
    
    withProgress(message = "Reading in amica file", {
      #
      tryCatch({
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
          rowData(reacValues$proteinData)[isQuantRnames(reacValues$proteinData),]
        
      },
      error = function(cond) {
        #message('amica upload failed. ', cond)

        shinyalert(
          title = "amica upload failed.",
          text = gsub("Error in .*):|.*Error: ", "Error: ",  cond),
          size = "s", 
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = FALSE,
          type = "error",
          showConfirmButton = TRUE,
          showCancelButton = FALSE,
          confirmButtonText = "OK",
          confirmButtonCol = "#669966",
          timer = 0,
          imageUrl = "",
          animation = TRUE
        )
      },
      warning = function(cond) {
        message(paste(cond))
      }
      # , finally = {
      #   showNotification(paste("Successfully processed amica file"), type = "message")
      # }
      )
      #
    })
  }
  
  if (input$source == "custom") {
    if (is.null(input$customFile)) {
      allFilesUploaded <- FALSE
      showNotification(paste("Need to upload custom tab-separated file."), type = "message")
      return("")
    }
    
    if (is.null(input$specFile)) {
      allFilesUploaded <- FALSE
      showNotification(paste("Need to upload specification file to map relevant columns."),
                       type = "error")
      return("")
    }
    
    ###
    tryCatch({
      specs <- validateFile(input$specFile, c("Variable", "Pattern"))
    },
    error = function(cond) {
      message('Upload of experimental design file failed. ', cond)
      allFilesUploaded <- FALSE
      # showNotification(
      #   paste(cond),
      #   duration = 100,
      #   closeButton = T,
      #   type = "error"
      # )
      shinyalert(
        title = "Upload of specification file failed.",
        text = gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
      #return("")
    },
    warning = function(cond) {
      showNotification(paste(cond))
    }
    #, finally = {
    #  showNotification(paste("Succesfully uploaded experimental design."), type = "message")
    #}
    )
    ###
    
    tryCatch({
      reacValues$proteinData <-
        readInCustomSumm(input$customFile$datapath, specs, reacValues$expDesign, input$customDataLogTransform)
    },
    error = function(cond) {
      message('Custom upload failed. ', cond)
      allFilesUploaded <- FALSE
      
      shinyalert(
        title = "Upload of custom file failed.",
        text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
      
      },
    warning = function(cond) {
      showNotification(paste(cond))
    }, finally = {
      showNotification(paste("Successfully processed custom file"), type = "message")
    }
    )
  }
  
  if (input$source == "maxquant") {
    if (is.null(input$maxquantFile)) {
      allFilesUploaded <- FALSE
      showNotification(paste("Need to upload MaxQuant or FragPipe output."), type = "message")
      return("")
    }
    
    reacValues$uniqueGroups <- unique(reacValues$expDesign$groups)
    
    header <-
      unlist(strsplit(readLines(input$maxquantFile$datapath, n = 1), "\t"))
    header <- make.names(header)
    mqNames <-
      c("Majority.protein.IDs",
        "Razor...unique.peptides")
    fragNames <-
      c("Indistinguishable.Proteins",
        "Protein.ID",
        "Protein.Probability")
    
    if (all(mqNames %in% header &&
            length(grep("^Intensity|^iBAQ.|^LFQ.intensity", header)) > 0
            )) {
      tryCatch({
        reacValues$proteinData <-
          readInMQproteinGroupsSumm(input$maxquantFile$datapath,
                                    reacValues$expDesign)
      },
      error = function(cond) {
        allFilesUploaded <- FALSE
        message('MaxQuant upload failed. ', cond)
        
        shinyalert(
          title = "Upload of MaxQuant file failed.",
          text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
          size = "s", 
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = FALSE,
          type = "error",
          showConfirmButton = TRUE,
          showCancelButton = FALSE,
          confirmButtonText = "OK",
          confirmButtonCol = "#669966",
          timer = 0,
          imageUrl = "",
          animation = TRUE
        )
        return("")
      },
      warning = function(cond) {
        message(paste(cond))
      }, finally = {
        showNotification(paste("Reading in MaxQuant ..."), type = "message")
      }
      )
      reacValues$dbTool <- "maxquant"
      
    } else if (all(fragNames %in% header)) {
      #
      tryCatch({
        reacValues$proteinData <-
          readInFragPipeProteinGroupsSumm(input$maxquantFile$datapath,
                                          reacValues$expDesign)
      },
      error = function(cond) {
        allFilesUploaded <- FALSE
        message('FragPipe upload failed. ', cond)
        
        shinyalert(
          title = "Upload of FragPipe file failed.",
          text =  gsub("Error in .*):|.*Error: ", "Error: ",  cond),
          size = "s", 
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = FALSE,
          type = "error",
          showConfirmButton = TRUE,
          showCancelButton = FALSE,
          confirmButtonText = "OK",
          confirmButtonCol = "#669966",
          timer = 0,
          imageUrl = "",
          animation = TRUE
        )
        
        # showNotification(
        #   paste(cond),
        #   duration = 100,
        #   closeButton = T,
        #   type = "error"
        # )
      },
      warning = function(cond) {
        message(paste(cond))
      }, finally = {
        showNotification(paste("Reading in FragPipe ..."), type = "message")
      }
      )
      #
      reacValues$dbTool <- "fragpipe"
    } else {
      allFilesUploaded <- FALSE
      
      shinyalert(
        title = "Unrecognized input format.",
        text = paste0(
          "Unrecognized input format\n",
          "\tallowed input files:\n",
          "\tMaxQuant: proteinGroups.txt\n",
          "\tFragPipe: combined_protein.tsv"
        ),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
      
      # showNotification(
      #   paste0(
      #     "Unrecognized input format\n",
      #     "\tallowed input files:\n",
      #     "\tMaxQuant: proteinGroups.txt\n",
      #     "\tFragPipe: combined_protein.tsv"
      #   ),
      #   type = "error"
      # )
    }
    
  }
  
  # contrasts
  if (input$source != "example" && input$source != "amica") {
    if (is.null(input$contrastMatrix)) {
      showNotification(paste("Need to upload contrast matrix."), type = "message")
      return("")
    }

    ###
    tryCatch({
      contrastData <- validateFile(input$contrastMatrix, NULL)
    },
    error = function(cond) {
      message('Upload of contrast file failed. ', cond)
      allFilesUploaded <- FALSE

      shinyalert(
        title = "Unrecognized input format.",
        text = gsub("Error in .*):|.*Error: ", "Error: ",  cond),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "error",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#669966",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    },
    warning = function(cond) {
      showNotification(paste(cond))
    }
    #, finally = {
    #  showNotification(paste("Succesfully uploaded experimental design."), type = "message")
    #}
    )
    ###
    if (exists("contrastData")) {
      reacValues$contrastMatrix <- contrastData
      
      for (elem in unique(c(contrastData[[1]], contrastData[[2]]))) {
        if (!is.null(reacValues$expDesign) &&
            length(grep(elem, reacValues$expDesign$groups)) < 1) {
          shiny:::reactiveStop(showNotification(
            paste(
              "Group",
              elem,
              "in contrasts not in uploaded experimental design"
            ),
            type = "warning",
            duration = 100
          ))
        }
      }
    }
  } else {
    reacValues$nsubmits <- reacValues$nsubmits + 1
    if (reacValues$nsubmits < 2 && 
        !is.null(reacValues$proteinData) &&
        !is.null(reacValues$expDesign) &&
        !is.null(reacValues$contrastMatrix)) {
      #toggle('ibaq_help')
      toggle(selector = "#navbar li a[data-value=qctab]")
      toggle(selector = "#navbar li a[data-value=quanttab]")
      toggle(selector = "#navbar li a[data-value=comparemicatab]")
      toggle(id = 'hide_before_input', anim = T)
    }
  }
  if (!is.null(reacValues$proteinData) &&
      !is.null(reacValues$expDesign) &&
      !is.null(reacValues$contrastMatrix)) {
    reacValues$uploadSuccess <- TRUE
  }
})