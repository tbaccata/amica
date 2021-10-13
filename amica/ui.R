ui <-
  fluidPage(
    tags$script(inactivity),
    useShinyjs(),
    navbarPage(
      "amica",
      id = "navbar",
      inverse = TRUE,
      theme = bs_theme(
        version = 4,
        #bg = "#21908dff",
        #fg = "black",
        primary = '#339999',
        secondary = '#669966',
        accent = "#339999",
        base_font = font_google("Fira Sans"),
        "font-size-base" = "0.9rem"
      ),
      source('R/ui/uiInput.R', local = TRUE)$value,
      source('R/ui/uiQC.R', local = TRUE)$value,
      source('R/ui/uiDiffAbundance.R', local = TRUE)$value,
      source('R/ui/uicompareAmica.R', local = TRUE)$value,
      source('R/ui/uiAbout.R', local = TRUE)$value
    )
  )