output$brewercols <- renderPlot({
  display.brewer.all()
})


output$brewerOptionsQual <- renderUI({
  req(reacValues$proteinData)
  selectInput(
    "brewerOptionsQual",
    "Choose a color palette for qualitative data:",
    row.names(brewer.pal.info),
    selected = "Set2"
  )
})

output$brewerOptionsDiv <- renderUI({
  req(reacValues$proteinData)
  selectInput(
    "brewerOptionsDiv",
    "Choose a color palette for diverging data:",
    row.names(brewer.pal.info[brewer.pal.info$category != "qual", ]),
    selected = "RdYlBu"
  )
})

myColors <- reactive({
  color <-
    ifelse(is.null(input$brewerOptionsQual),
           "Accent",
           input$brewerOptionsQual)
  reverse <- ifelse(input$revQual == "yes", TRUE, FALSE)
  nc <- nrow(colData(reacValues$proteinData))
  
  mycols <- c()
  maxcol <- brewer.pal.info[color, "maxcolors"]
  if (maxcol < nc) {
    mycols <- colorRampPalette(brewer.pal(maxcol, color))(nc)
  } else {
    nc <- ifelse(nc < 3, 3, nc)
    mycols <- brewer.pal(nc, color)
  }
  if (reverse)
    mycols <- rev(mycols)
  mycols
})

mySelectionColors <- reactive({
  color <- ifelse(is.null(input$brewerOptionsQual), "Set2", input$brewerOptionsQual)
  reverse <- ifelse(input$revQual=="yes", TRUE, FALSE)
  #nc <- length(unique(reacValues$proteinData$groups))
  
  nc <- length(unique(colData(reacValues$proteinData)$groups))
  
  mycols <- c()
  maxcol <- brewer.pal.info[color, "maxcolors"]
  if (maxcol < nc) {
    mycols <- colorRampPalette(brewer.pal(maxcol, color))(nc)
  } else {
    nc <- ifelse(nc<3, 3, nc)
    mycols <- brewer.pal(nc, color)
  }
  
  names(mycols) <- unique(colData(reacValues$proteinData)$groups)
  
  if (reverse) mycols <- rev(mycols)
  mycols
})

output$myColorPanel <- renderUI({
  req(reacValues$proteinData)
  
  groups <- unique(colData(reacValues$proteinData)$groups)
  
  lapply(seq_along(groups), function(i) {
    colourInput(paste("col", i, sep="_"), paste0(groups[i]), mySelectionColors()[i])
  })
})

output$groupFactors <- renderUI({
  req(reacValues$analysisSuccess)
  selectizeInput(
    "groupFactors",
    "Select the ordering of groups on the x-axis of plots (QC-plots, profile plots).
      Ordering will be done from first to last selected.",
    unique(reacValues$expDesign$groups),
    multiple = T,
    selected = NULL
  )
})

observeEvent(input$submitGroupFactors, {
  req(input$groupFactors)
  reacValues$groupFactors <- input$groupFactors
})

observeEvent(input$resetGroupFactors, {
  reacValues$groupFactors <- NULL
})

output$groupFactorsSummary <- renderText({
  req(reacValues$groupFactors)
  out <- "Group factors:\n\n"
  for (idx in seq_along(reacValues$groupFactors)) {
    out <- paste0(out, "\t", idx, ": ", reacValues$groupFactors[idx], "\n")
  }
  out
})

output$volcanoMAColors <- renderUI({
  pal <- brewer.pal(3, "Set2")
  
  lapply(seq_along(1:2), function(i) {
    colourInput(paste("volcanocol", i, sep="_"), paste0("Color ", i, ":"), pal[i])
  })
})

output$fcPlotColors <- renderUI({
  pal <- myScatterColors()
  
  lapply(seq_along(pal), function(i) {
    colourInput(paste("fcplotcol", i, sep="_"), paste0("Color ", i, ":"), pal[i])
  })
})

output$eulerColors <- renderUI({
  req(input$upset1Sample)
  pal <- myScatterColors()
  
  lapply(seq_along(input$upset1Sample), function(i) {
    colourInput(paste("eulercol", i, sep="_"), paste0("Color ", i, ":"), pal[i])
  })
})

myGroupColors <- reactive({
  groups <- unique(colData(reacValues$proteinData)$groups)
  
  colors <- lapply(seq_along(groups), function(i) {
    input[[paste("col", i, sep="_")]]
  })
  
  if (is.null(input$col_1)) colors <- mySelectionColors()
  
  names(colors) <- groups
  colors
})

myScatterColors <- reactive({
  color <- ifelse(is.null(input$brewerOptionsScatter), "Set2", input$brewerOptionsScatter)
  reverse <- ifelse(input$revScatter=="yes", TRUE, FALSE)
  
  pal <- brewer.pal(5, color)
  if (reverse) pal <- rev(pal)
  pal
})

heatColors <- reactive({
  color <- ifelse(is.null(input$brewerOptionsDiv), "RdYlBu", input$brewerOptionsDiv)
  maxcol <- brewer.pal.info[color, "maxcolors"]
  
  reverse <- ifelse(input$revDiv=="yes", TRUE, FALSE)
  
  if (reverse) {
    pal <- colorRampPalette(rev(brewer.pal(n = maxcol, name =
                                             color)))(100)
  } else {
    pal <- colorRampPalette(brewer.pal(n = maxcol, name =
                                         color))(100)
  }
  pal
})
