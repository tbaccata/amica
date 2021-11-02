tabPanel(
  title = 'Compare amica datasets',
  value = 'comparemicatab',
  br(),br(),br(),
  actionButton("showAmicasTutorial", "Tutorial", icon = icon("info") ),
  br(),br(),br(),
  fluidRow(
    column(
      width = 6,
      radioButtons(
        inputId = "amicaCompareSource",
        label = "Select the file input.",
        choices = c(#"Load in example" = "example",
          "Upload amica file" = "amica")
      ),
      conditionalPanel(
        condition = "input.amicaCompareSource == 'amica'",
        fileInput("amicaFile2", "Upload amica_proteinGroups.txt.",
                  width = "60%"),
        helpText(
          "Have you run amica before and want to compare it to the currently loaded in dataset?"
        )
      ),
      
      radioButtons(
        'mergeKey',
        'Key column',
        c('Gene.names', 'ProteinIDs'),
        selected = 'Gene.names'
      ),
      helpText(
        'This column determines the key that is used to merge the experiments.
                                  Use "ProteinIDs" when both experiments used the same search database.
                                           Use "Gene.names" when you want to compare experiments from different origins (e.g different organisms)'
      )
    ),
    
    column(
      width = 6,
      textInput("suffix1", "Suffix for original input", value = "exp1"),
      helpText(
        'Enter a suffix to better distinguish the column names of the experiments. For example you could enter "AP-MS" (without the quotation mark) if your first experiment was of that kind.'
      ),
      textInput("suffix2", "Suffix for uploaded input", value = "exp2"),
      helpText(
        'Enter a suffix to better distinguish the column names of the experiments. For example you could enter "TurboID" (without the quotation mark) if your first experiment was of that kind.'
      ),
      textInput("subPattern", "Pattern to substitute in ProteinID column", value = ";"),
      helpText(
        'Everything after the pattern will be removed. You might have the case where "ProteinA" is the stable id in your original data, while the uploaded file contains the id "ProteinA;ProteinB". The pattern ";" (enter without the quotation mark) would remove ";ProteinB".'
      ),
      actionButton("submitAmicaComparisons", strong("Submit"))
    )
  ),
  verbatimTextOutput("summaryMergedAmica"),
  
  shinyjs::hidden(
    div(
      id = "compareAmicaInput",
      fluidRow(column(width = 6,
                      uiOutput('download_merged_amica'),),
            #    column(
            #      width = 6,
            #      HTML(
            #        '
            # <p>
            # <h3>Usage:</h3>
            # On this tab you can correlate samples of both data sets. <br>
            # If you go back to the "Differential Abundance" tab you see a select
            # bar at the of the page:</p> <br>
            #      <img src="compare_amicas_tutorial/multi_amica_data.png"><br>
            #      <p>Select "multiple_amica_input" and press the Submit button.
            #      The merged data set is now selected and you can compare multiple
            #      group comparisons across experiment etc. Almost the full 
            #      functionality from the Differential Abundance tab is available 
            #      for the merged data sets.</p>
            #      '
            #      )
            #    )
               ), 
      br(),br(),br(),
      inline(uiOutput("assayNamesAmicas")),
      
      fluidRow(
        column(
          width = 6,
          h4('Scatter plot'),
          uiOutput("compareScatterPlotsAmica"),
          ###
          inline(actionButton("submitScatterAmicas", "Plot scatter plot", icon = icon("cog"))),
          inline(actionButton("scatteramicaParams", "", icon = icon("wrench"))),
          shinyjs::hidden(
            div(
              style = "display: grid;
            grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
            grid-gap: 30px;",
              
              id = 'toggle_scatteramica_params',
              numericInput(
                "scatteramica_width",
                "Width in pixel.",
                value = 768,
                min = 338,
                max = 1352,
                step = 10
              ),
              numericInput(
                "scatteramica_height",
                "Height in pixel.",
                value = 676,
                min = 338,
                max = 1352,
                step = 10
              ),
              numericInput(
                "scatteramica_base",
                "Base font size in pt.",
                value = 14,
                min = 4,
                max = 32,
                step = 1
              ),
              radioButtons(
                "scatteramica_showLine",
                "Show line in plot?",
                choices = c("straight line", "none", "linear regression"),
                selected = "straight line"
              ),
              radioButtons(
                "scatteramica_format",
                "Save image",
                choices = c("svg", "png"),
                selected = "svg"
              )
            )
          ),
          ###
          plotlyOutput("scatterPlotsAmica", height = 800)
        ),
        column(width = 6,
               h4("Correlation plot"),
               uiOutput("corrAmicasSamplesInput"),
               inline(actionButton("submitCorAmicas", "Plot Correlation", icon = icon("cog"))),
               inline(actionButton("corAmicasParams", "", icon = icon("wrench"))),
               shinyjs::hidden(
                 div(
                   style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                   
                   id = 'toggle_corAmicas_params',
                   numericInput(
                     "corAmicas_width",
                     "Width in pixel.",
                     value = 768,
                     min = 338,
                     max = 1352,
                     step = 10
                   ),
                   numericInput(
                     "corAmicas_height",
                     "Height in pixel.",
                     value = 676,
                     min = 338,
                     max = 1352,
                     step = 10
                   ),
                   radioButtons(
                     "corAmicas_format",
                     "Save image",
                     choices = c("svg", "png"),
                     selected = "svg"
                   )
                 )
               ),
               plotlyOutput("corrAmicasPlotly", height = 800)
        )
      )
    )
  ),
  footer()
)