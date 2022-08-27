tabPanel(
  title = 'QC',
  value = 'qctab',
  h1("Quality Control"),
  h2("Compare different intensities"),
  br(),
  inline(uiOutput("assayNames")),
  inline(actionButton('assaysHelp', label = '', icon = icon("info"))),
  uiOutput("assaysHelpBox"),
  helpText("ImputedIntensity are the intensities used for quantification"),
  verbatimTextOutput("quantSummary"),
  downloadButton('qcReport', 'Generate QC Report', icon = icon("file")),
  tabsetPanel(
    type = "tabs",
    tabPanel(
      h3("PCA"),
      uiOutput("pcaSamplesInput"),
      actionButton("submitPCA", "Plot PCA", icon = icon("cog")),
      actionButton("pcaParams", "", icon = icon("wrench")),
      shinyjs::hidden(div(
        style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
        
        id = 'toggle_pca_params',
        numericInput(
          "pca_width",
          "Width in pixel.",
          value = 768,
          min = 338,
          max = 1352,
          step = 10
        ),
        numericInput(
          "pca_height",
          "Height in pixel.",
          value = 676,
          min = 338,
          max = 1352,
          step = 10
        ),
        radioButtons(
          "pca_format",
          "Save image",
          choices = c("svg", "png"),
          selected = "svg"
        ),
        numericInput(
          "pca_base",
          "Base font size in pt.",
          value = 14,
          min = 4,
          max = 32,
          step = 1
        ),
        numericInput(
          "pca_legend",
          "Legend font size in pt.",
          value = 10,
          min = 4,
          max = 16,
          step = 1
        ),
        numericInput(
          "pca_pointsize",
          "Point size.",
          value = 3,
          min = 1,
          max = 8,
          step = 1
        )
      )),
      fluidRow(
        column(width = 2),
        column(width = 8,
               plotlyOutput("pca", height = 800)),
        column(width = 2)
      )
      
    ),
    tabPanel(h3("Intensity statistics"),
             fluidRow(
               column(
                 width = 6,
                 h4("Boxplot"),
                 inline(actionButton("submitBoxplot", "Plot boxplot", icon = icon("cog"))),
                 inline(actionButton("boxplotHelp", icon = icon("info") , label = NULL)),
                 inline(uiOutput("boxplotHelpBox")),
                 ###
                 inline(actionButton("boxplotParams", "", icon = icon("wrench"))),
                 shinyjs::hidden(
                   div(
                     style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                     
                     id = 'toggle_boxplot_params',
                     numericInput(
                       "boxplot_width",
                       "Width in pixel.",
                       value = 768,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     numericInput(
                       "boxplot_height",
                       "Height in pixel.",
                       value = 676,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     radioButtons(
                       "boxplot_format",
                       "Save image",
                       choices = c("svg", "png"),
                       selected = "svg"
                     ),
                     numericInput(
                       "boxplot_base",
                       "Base font size in pt.",
                       value = 14,
                       min = 4,
                       max = 32,
                       step = 1
                     ),
                     numericInput(
                       "boxplot_legend",
                       "Legend font size in pt.",
                       value = 10,
                       min = 4,
                       max = 16,
                       step = 1
                     )
                   )
                 ),
                 ###
                 
                 plotlyOutput("boxPlot", height = 800)
               ),
               column(
                 width = 6,
                 h4("Density plot"),
                 inline(actionButton("submitDensity", "Plot density plot", icon = icon("cog"))),
                 inline(actionButton("densityHelp", icon = icon("info") , label = NULL)),
                 inline(uiOutput("densityHelpBox")),
                 ###
                 inline(actionButton("densityParams", "", icon = icon("wrench"))),
                 shinyjs::hidden(
                   div(
                     style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                     
                     id = 'toggle_density_params',
                     numericInput(
                       "density_width",
                       "Width in pixel.",
                       value = 768,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     numericInput(
                       "density_height",
                       "Height in pixel.",
                       value = 676,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     radioButtons(
                       "density_format",
                       "Save image",
                       choices = c("svg", "png"),
                       selected = "svg"
                     ),
                     numericInput(
                       "density_base",
                       "Base font size in pt.",
                       value = 14,
                       min = 4,
                       max = 32,
                       step = 1
                     ), numericInput(
                       "density_legend",
                       "Legend font size in pt.",
                       value = 10,
                       min = 4,
                       max = 16,
                       step = 1
                     )
                   )
                 ),
                 ###
                 plotlyOutput("densityPlot", height = 800)
               )
             )),
    
    tabPanel(
      h3("Correlation and CVs"),
      fluidRow(
        column(width = 6,
               h4("Correlation"),
               uiOutput("corrSamplesInput"),
               inline(actionButton("submitCor", "Plot Correlation", icon = icon("cog"))),
               inline(actionButton("corHelp", icon = icon("info") , label = NULL)),
               inline(uiOutput("corHelpBox")),
               ###
               inline(actionButton("corParams", "", icon = icon("wrench"))),
               shinyjs::hidden(
                 div(
                   style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                   
                   id = 'toggle_cor_params',
                   numericInput(
                     "cor_width",
                     "Width in pixel.",
                     value = 768,
                     min = 338,
                     max = 1352,
                     step = 10
                   ),
                   numericInput(
                     "cor_height",
                     "Height in pixel.",
                     value = 676,
                     min = 338,
                     max = 1352,
                     step = 10
                   ),
                   radioButtons(
                     "cor_format",
                     "Save image",
                     choices = c("svg", "png"),
                     selected = "svg"
                   ),
                   checkboxInput("cor_annot",
                                 "Annotate samples in dendrogram?",
                                 value = TRUE),
                 )
               ),
               
               plotlyOutput("corrPlotly", height = 800)
        ),
        column(width = 6,
               h4("Coefficient of variation"),
               inline(actionButton("submitCVs", "Plot CVs", icon = icon("cog"))),
               inline(actionButton("cvHelp", icon = icon("info") , label = NULL)),
               inline(uiOutput("cvHelpBox")),
               ###
               inline(actionButton("cvParams", "", icon = icon("wrench"))),
               shinyjs::hidden(
                 div(
                   style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                   
                   id = 'toggle_cv_params',
                   numericInput(
                     "cv_width",
                     "Width in pixel.",
                     value = 768,
                     min = 338,
                     max = 1352,
                     step = 10
                   ),
                   numericInput(
                     "cv_height",
                     "Height in pixel.",
                     value = 676,
                     min = 338,
                     max = 1352,
                     step = 10
                   ),
                   radioButtons(
                     "cv_format",
                     "Save image",
                     choices = c("svg", "png"),
                     selected = "svg"
                   ),
                   numericInput(
                     "cv_base",
                     "Base font size in pt.",
                     value = 14,
                     min = 4,
                     max = 32,
                     step = 1
                   ),
                   numericInput(
                     "cv_legend",
                     "Legend font size in pt.",
                     value = 10,
                     min = 4,
                     max = 16,
                     step = 1
                   )
                 )
               ),
               ###
               plotlyOutput("boxplotCV", height = 800)
        )
      )
    )
    
  ),
  
  
  tabsetPanel(
    type = "tabs",
    tabPanel(h3("Contaminants"),
             fluidRow(
               column(
                 h4("Contaminants"),
                 width = 6,
                 inline(actionButton("submitConts", "Plot Contamiants", icon = icon("cog"))),
                 inline(actionButton("contaminantHelp", icon = icon("info") , label = NULL)),
                 inline(uiOutput("contaminantHelpBox")),
                 ###
                 inline(actionButton("contaminantsParams", "", icon = icon("wrench"))),
                 shinyjs::hidden(
                   div(
                     style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                     
                     id = 'toggle_contaminants_params',
                     numericInput(
                       "contaminants_width",
                       "Width in pixel.",
                       value = 768,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     numericInput(
                       "contaminants_height",
                       "Height in pixel.",
                       value = 676,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     radioButtons(
                       "contaminants_format",
                       "Save image",
                       choices = c("svg", "png"),
                       selected = "svg"
                     ),
                     numericInput(
                       "contaminants_base",
                       "Base font size in pt.",
                       value = 14,
                       min = 4,
                       max = 32,
                       step = 1
                     ),
                     numericInput(
                       "contaminants_legend",
                       "Legend font size in pt.",
                       value = 10,
                       min = 4,
                       max = 16,
                       step = 1
                     )
                   )
                 ),
                 ###
                 plotlyOutput("contaminants", height = 600),
               ),
               column(
                 width = 6,
                 h4("Most abundant proteins"),
                 
                 
                 inline(uiOutput("samples")),
                 inline(actionButton("abundantHelp", icon = icon("info") , label = NULL)),
                 inline(uiOutput("abundantHelpBox")),
                 ###
                 inline(actionButton("abundantParams", "", icon = icon("wrench"))),
                 shinyjs::hidden(
                   div(
                     style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                     
                     id = 'toggle_abundant_params',
                     numericInput(
                       "abundant_width",
                       "Width in pixel.",
                       value = 768,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     numericInput(
                       "abundant_height",
                       "Height in pixel.",
                       value = 676,
                       min = 338,
                       max = 1352,
                       step = 10
                     ),
                     radioButtons(
                       "abundant_format",
                       "Save image",
                       choices = c("svg", "png"),
                       selected = "svg"
                     ),
                     numericInput(
                       "abundant_base",
                       "Base font size in pt.",
                       value = 14,
                       min = 4,
                       max = 32,
                       step = 1
                     )
                     # numericInput(
                     #   "abundant_legend",
                     #   "Legend font size in pt.",
                     #   value = 10,
                     #   min = 4,
                     #   max = 16,
                     #   step = 1
                     # )
                   )),
                 ###
                 plotlyOutput("mostAbundant", height = 600),
               )
             )),
    tabPanel(
      h3("Protein groups"),
      fluidRow(column(width = 6,
                      h4(
                        "Identified proteins"
                      ),),
               column(width = 6,
                      h4("Missing values"),)),
      fluidRow(
        column(
          width = 6,
          inline(actionButton("submitNumProts", "Plot barplot", icon = icon("cog"))),
          inline(actionButton("idHelp", icon = icon("info") , label = NULL)),
          inline(uiOutput("idHelpBox")),
          ###
          inline(actionButton("barplotIdParams", "", icon = icon("wrench"))),
          shinyjs::hidden(
            div(
              style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
              
              id = 'toggle_barplotId_params',
              numericInput(
                "barplotId_width",
                "Width in pixel.",
                value = 768,
                min = 338,
                max = 1352,
                step = 10
              ),
              numericInput(
                "barplotId_height",
                "Height in pixel.",
                value = 676,
                min = 338,
                max = 1352,
                step = 10
              ),
              radioButtons(
                "barplotId_format",
                "Save image",
                choices = c("svg", "png"),
                selected = "svg"
              ),
              numericInput(
                "barplotId_base",
                "Base font size in pt.",
                value = 14,
                min = 4,
                max = 32,
                step = 1
              ),
              numericInput(
                "barplotId_legend",
                "Legend font size in pt.",
                value = 10,
                min = 4,
                max = 16,
                step = 1
              )
            )
          ),
          ###
          plotlyOutput("barplotProteins", height = 600),
        ),
        column(
          width = 6,
          
          inline(actionButton("submitNumMVs", "Plot barplot", icon = icon("cog"))),
          inline(actionButton("mvHelp", icon = icon("info") , label = NULL)),
          inline(uiOutput("mvHelpBox")),
          ###
          inline(actionButton("barplotMvParams", "", icon = icon("wrench"))),
          shinyjs::hidden(
            div(
              style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
              
              id = 'toggle_barplotMv_params',
              numericInput(
                "barplotMv_width",
                "Width in pixel.",
                value = 768,
                min = 338,
                max = 1352,
                step = 10
              ),
              numericInput(
                "barplotMv_height",
                "Height in pixel.",
                value = 676,
                min = 338,
                max = 1352,
                step = 10
              ),
              radioButtons(
                "barplotMv_format",
                "Save image",
                choices = c("svg", "png"),
                selected = "svg"
              ),
              numericInput(
                "barplotMv_base",
                "Base font size in pt.",
                value = 14,
                min = 4,
                max = 32,
                step = 1
              ),
              numericInput(
                "barplotMv_legend",
                "Legend font size in pt.",
                value = 10,
                min = 4,
                max = 16,
                step = 1
              )
            )
          ),
          ###
          plotlyOutput("barplotMissingValues", height = 600),
        )
      )
    ),
    tabPanel(
    h3("Protein Overlap"),
    inline(actionButton("submitOverlapHeatmap", "Plot overlap", icon = icon("cog"))),
    inline(actionButton("overlapHeatmapHelp", icon = icon("info") , label = NULL)),
    inline(uiOutput("overlapHeatmapHelpBox")),
    ###
    inline(actionButton("overlapHeatmapParams", "", icon = icon("wrench"))),
    shinyjs::hidden(
      div(
        style = "display: grid;
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",

        id = 'toggle_overlapHeatmap_params',
        numericInput(
          "overlapHeatmap_width",
          "Width in pixel.",
          value = 768,
          min = 338,
          max = 1352,
          step = 10
        ),
        numericInput(
          "overlapHeatmap_height",
          "Height in pixel.",
          value = 676,
          min = 338,
          max = 1352,
          step = 10
        ),
        radioButtons(
          "overlapHeatmap_format",
          "Save image",
          choices = c("svg", "png"),
          selected = "svg"
        ),
        radioButtons(
          "overlapHeatmap_metric",
          "Which metric to plot?",
          choices = c("overlap_coefficient",
                      "jaccard_index",
                      "num_shared")
        ),
        checkboxInput("overlapHeatmap_annot",
                      "Annotate samples in dendrogram?",
                      value = FALSE)
      )
     ),
    ###
    fluidRow(
      column(width = 6,
             h4("Overlap Heatmap"),
             plotlyOutput("overlapHeatmapPlotly", height = 600)
      ),
      column(
        width = 6,
        h4("Overlap table"),
        DTOutput("overlapSummaryDT")
      )
    )
    )
  ),
  
  br(),
  br(),
  br(),
  fluidRow(column(width = 4,
                  h3("Scatter plot"),)),
  
  fluidRow(
    column(
      width = 2,
      uiOutput("assayName1"),
      uiOutput("compareScatterPlot1"),
      helpText("Select the sample from the intensity to plot on the x-axis.
                                                ")
      
    ),
    column(
      width = 2,
      uiOutput("assayName2"),
      uiOutput("compareScatterPlot2"),
      helpText("Select the sample from the intensity to plot on the y-axis.
                                                "),
      actionButton("submitScatterPlot", "Plot Scatterplot"),
      ###
      actionButton("scatterParams", "", icon = icon("wrench")),
      shinyjs::hidden(
        div(
          style = "display: grid; 
          grid-template-columns: 40% repeat(2, 40%); ## same as repeat(4, 20%)
          grid-gap: 20px;",
          
          id = 'toggle_scatter_params',
          numericInput(
            "scatter_width",
            "Width in pixel.",
            value = 768,
            min = 338,
            max = 1352,
            step = 10
          ),
          numericInput(
            "scatter_height",
            "Height in pixel.",
            value = 676,
            min = 338,
            max = 1352,
            step = 10
          ),
          radioButtons(
            "scatter_showLine",
            "Show line in plot?",
            choices = c("straight line", "none", "linear regression"),
            selected = "straight line"
          ),
          radioButtons(
            "scatter_format",
            "Save image",
            choices = c("svg", "png"),
            selected = "svg"
          ),
          numericInput(
            "scatter_base",
            "Base font size in pt.",
            value = 14,
            min = 4,
            max = 32,
            step = 1
          ),
          numericInput(
            "scatter_legend",
            "Legend font size in pt.",
            value = 10,
            min = 4,
            max = 16,
            step = 1
          )
        )
      )
      ###
    ),
    column(
      width = 7,
      plotlyOutput("scatterPlot", height = 800)
    ), column(width = 1)
  ),
  footer()
)