output$exampleFiles <- downloadHandler(
  filename = "examples.zip",
  # desired file name on client
  content = function(con) {
    file.copy("data/examples.zip", con)
  }
)

output$exampleAmicaFile <- renderDT({
  varName <-
    c(
      "Protein ID",
      "Gene name",
      "LFQ intensity prefix",
      "Imputed Intensity prefix",
      "razor unique count",
      "razor unique prefix",
      "p-value Prefix",
      "adj. p-value prefix",
      "fold change prefix",
      "Avg. expression prefix",
      "comparison infix",
      "quantified column",
      "Potential contaminant column"
    )
  colName <- c(
    "Majority.protein.IDs",
    "Gene.names",
    "LFQIntensity_",
    "ImputedIntensity_",
    "razorUniqueCount",
    "razorUniqueCount",
    "P.Value_",
    "adj.P.Val_",
    "logFC_",
    "AveExpr_",
    "__vs__",
    "quantified",
    "Potential.contaminant"
  )
  comment <- c(
    "",
    "",
    "MaxQuants (MQs) 'LFQ intensity' columns",
    "Imputed and re-normalized intensities",
    "MQs 'razor+unique count' column",
    "MQs 'razor+unique count' column per sample",
    "e.g P.Value_group1__vs__group2",
    "e.g adj.P.Val_group1__vs__group2",
    "e.g logFC_group1__vs__group2",
    "e.g AveExp_group1__vs__group2",
    "see below",
    "see below",
    "MQs Potential.contaminants column"
  )
  df <- data.frame(varName, colName, comment)
  names(df) <- c("Variable name",
                 "Column name or prefix",
                 "Comment")
  
  datatable(
    df,
    filter = "none",
    rownames = F,
    extensions = c('Buttons'),
    options = list(
      searching = FALSE,
      dom = 'Bfrtip',
      pageLength = 15,
      autoWidth = TRUE,
      buttons = list(# 'csv',
        list(
          extend = 'csv',
          fieldSeparator = '\t',
          fieldBoundary = ''
        ))
    )
  )
}, server = F)

output$exampleDesign <- renderDT({
  tmpData <-
    read.table("data/PXD0016455/design.txt",
               header = T,
               stringsAsFactors = F)
  
  
  datatable(
    tmpData,
    filter = "none",
    rownames = F,
    extensions = c('Buttons'),
    options = list(
      searching = FALSE,
      dom = 'Bfrtip',
      pageLength = 10,
      autoWidth = TRUE,
      buttons = list(# 'csv',
        list(
          extend = 'csv',
          fieldSeparator = '\t',
          fieldBoundary = ''
        ))
    )
  )
}, server = F)

output$exampleContrasts <- renderDT({
  tmpData <-
    read.table(
      "data/PXD0016455/contrasts.txt",
      header = T,
      stringsAsFactors = F
    )
  
  datatable(
    tmpData,
    filter = "none",
    rownames = F,
    extensions = c('Buttons'),
    options = list(
      searching = FALSE,
      dom = 'Bfrtip',
      pageLength = 10,
      autoWidth = TRUE,
      buttons = list(#'csv',
        list(
          extend = 'csv',
          fieldSeparator = '\t',
          fieldBoundary = ''
        ))
    )
  )
}, server = F)


output$specificationsExplanation <- renderDT({
  Variable <- c(
    "proteinId",
    "geneName",
    "intensityPrefix",
    "abundancePrefix",
    "razorUniqueCount",
    "razorUniquePrefix",
    "spectraCount",
    "contaminantCol"
  )
  Pattern <- rep("...", length(Variable))
  Mandatory <- c("yes",
                 "yes",
                 "yes",
                 "no",
                 "no",
                 "no",
                 "no",
                 "no")
  df <- data.frame(Variable, Pattern, Mandatory)
  datatable(
    df,
    filter = "none",
    rownames = F,
    extensions = c('Buttons'),
    options = list(
      searching = FALSE,
      dom = 'Bfrtip',
      pageLength = 10,
      autoWidth = TRUE,
      buttons = list(#'csv',
        list(
          extend = 'csv',
          fieldSeparator = '\t',
          fieldBoundary = ''
        ))
    )
  )
}, server = F)

output$exampleSpecifications <- renderDT({
  tmpData <-
    read.table(
      "data/PXD0016455/specification.txt",
      header = T,
      stringsAsFactors = F
    )
  
  datatable(
    tmpData,
    filter = "none",
    rownames = F,
    extensions = c('Buttons'),
    options = list(
      searching = FALSE,
      dom = 'Bfrtip',
      pageLength = 10,
      autoWidth = TRUE,
      buttons = list(#'csv',
        list(
          extend = 'csv',
          fieldSeparator = '\t',
          fieldBoundary = ''
        ))
    )
  )
}, server = F)