inline = function (x) {
  tags$div(style="display:inline-block;", x)
}

footer = function(x) {
  tags$div(
    style="footer{position: absolute; bottom:5%; left: 33%; padding:5px;}",
    HTML("<h5>How to cite us</h5>
         <p>Please cite amica in your publications, 
         we will provide a PMID pointing to this site shortly.\n
         </p>
         <p>Many thanks to the MS facility and the Max Perutz Labs IT team.
         If you have questions, feedback, suggestions, etc. please mail
         them to sebastian.didusch@univie.ac.at.
         </p>")
  )
}


# Define UI for application
ui <- #secure_app(head_auth = tags$script(inactivity),
  fluidPage(
    tags$script(inactivity),
    useShinyjs(),
    navbarPage(
      "amica v21.09.01",
      id = "navbar",
      inverse = TRUE,
      theme = bs_theme(
        version = 4,
        #bg = "#21908dff",
        #fg = "black",
        primary = '#339999', #primary = "#5ac865ff",
        secondary = '#669966', #secondary = "#21908dff",
        accent = "#339999",
        base_font = font_google("Fira Sans"),
        "font-size-base" = "0.9rem"
      ),
      tabPanel('Input',
               sidebarLayout(
                 sidebarPanel(
                   radioButtons(
                     inputId = "source",
                     label = "Select the file input.",
                     choices = c(
                       
                       "Upload output from MaxQuant or FragPipe" = "maxquant",
                       "Upload amica format" = "amica",
                       "Upload custom format" = "custom",
                       "Load in example" = "example"
                     )
                   ),
                   helpText(
                     "If 'amica' is selected previous output gets loaded in. When 'MaxQuant or FragPipe', 
                     or 'custom' are selected differential abundance get's computed with limma or DEqMS, 
                     for this you have to upload an additional 'contrast matrix' to your experimental 'design'."
                   ),
                   conditionalPanel(
                     condition = "input.source == 'example'",
                     
                     actionButton('exampleHelp', label = 'Example', icon = icon("info")),
                     uiOutput("exampleHelpBox"),
                     br(),br()
                   ),
                   # Wrap the file input in a conditional panel
                   

                   conditionalPanel(
                     # The condition should be that the user selects
                     # "file" from the radio buttons
                     condition = "input.source == 'amica'",
                     fileInput("amicaFile", "Upload amica_proteinGroups.txt."),
                     helpText(
                       "Have you run amica before? Input the amica output from a previous session."
                     )
                   ),
                   conditionalPanel(
                     # The condition should be that the user selects
                     # "file" from the radio buttons
                     condition = "input.source == 'maxquant'",
                     fileInput(
                       "maxquantFile",
                       "Upload 'proteinGroups.txt' file from MaxQuant or 'combined_protein.tsv' from FragPipe."
                     )
                   ),
                   conditionalPanel(
                     # The condition should be that the user selects
                     # "file" from the radio buttons
                     condition = "input.source == 'custom'",
                     fileInput("customFile", "Upload custom tab delimeted file."),
                     helpText(
                       "File needs to be tab-delimited and it needs to contain a proteinId, Gene.name, intensities and peptide counts."
                     )
                   ),
                   conditionalPanel(
                     # The condition should be that the user selects
                     # "file" from the radio buttons
                     condition = "input.source == 'custom'",
                     fileInput("specFile", "Upload a specification file in tab delimited format."),
                     helpText(
                       "File needs to have two columns: the column 'Variable' which defines the variable and 
                       'Pattern' which is the column name or a pattern."
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.source != 'example'",
                     fileInput(
                       "groupSpecification",
                       "Experimental design",
                       c('text/csv',
                         'text/comma-separated-values,text/plain'),
                       multiple = F
                     ),
                     helpText(
                       "A tab-separated file containing the sample name ordered in appearance in the uploaded file and its corresponding group."
                     ),
                   ),
                   
                   conditionalPanel(
                     condition = "input.source == 'custom' || input.source == 'maxquant'",
                     fileInput(
                       "contrastMatrix",
                       "Contrast matrix (group comparisons)",
                       c('text/csv',
                         'text/comma-separated-values,text/plain'),
                       multiple = F
                     ),
                     helpText(
                       "A tab-separated file with the sample-of-interest in the first column and the sample it is getting compared to in the second column."
                     ),
                   ),
                   
                   ######## PARAMETERS BEGIN
                   actionButton("submitAnalysis", "Upload"),
                   actionButton("resetAnalysis", "Reset"),
                   
                   conditionalPanel(
                     condition = "input.source != 'example' && input.source != 'amica' && input.submitAnalysis!=0",
                     h3("Analysis Parameters"),
                     bsCollapse(
                       id = "collapseExample",
                       open = "Panel 2",
                       bsCollapsePanel(
                         "Advanced",
                         #style="info",
                         h4("Filter on count values"),
                         numericInput(
                           "minRazor",
                           "min. razor peptides",
                           min = 1,
                           value = 2
                         ),
                         helpText(
                           "Default is 2. It is difficult to distinguish
                     quantitative effects on the protein level from peptide
                     specific effects for 'one-hit wonders'."
                         ),
                         numericInput("minMSMS", "min. MSMS count", min = 1, value = 3),
                         
                         h4("Filter on valid values per group"),
                         uiOutput("filterValuesInput"),
                         radioButtons(
                           "validValuesGroup",
                           "Require min valid values in at least one group or in all groups?",
                           choices = c("in_one_group", "in_each_group"),
                           selected = "in_one_group"
                         ),
                         numericInput(
                           "minValidValue",
                           "Min. valid values in group",
                           min = 1,
                           value = 3
                         ),
                         
                         h4("Normalization and statistical testing"),
                         radioButtons(
                           "renormalizationMethod",
                           "Normalization method",
                           choices = c("None", "Quantile", "VSN", "Median.centering"),
                           selected = "None"
                         ),
                         helpText("The renormalization utilzes LFQ intensities."),
                         checkboxInput(
                           "limmaTrend",
                           "use DEqMS? If not selected differential abundant proteins are detected by limma.",
                           value = FALSE
                         ),
                         helpText(
                           "DEqMS takes the number of PSMs of a protein into account for quantitation."
                         ),
                         
                         h4("Imputation"),
                         radioButtons(
                           "impMethod",
                           "Imputation method",
                           choices = c("normal", "min", "global"),
                           selected = "normal"
                         ),
                         
                         
                         helpText(
                           "'normal' samples from the
                     normal distribution of all intensities
                     from the same sample but dowshifted.
                     Min takes the minimum observed value
                     which is useful for analyzing pilots"
                         ),
                         numericInput(
                           "downshift",
                           "Imputation downshift",
                           value = 1.8,
                           min = 0,
                           step = 0.01
                         ),
                         numericInput(
                           "width",
                           "Imputation width",
                           value = 0.3,
                           min = 0,
                           step = 0.01
                         ),
                         style = "info"
                       )
                       
                     ),
                     actionButton("runAnalysis", "Analyze"),
                   )
                   
                   ######## PARAMETERS END
                 )
                 
                 ,
                 mainPanel(
                   #h1("amica"),
                   HTML('<center><img src="ga_amica.png" width="50%"></center>'),
                   #img(src = 'ga_amica.svg',  align = "center"),
                   br(),
                   br(),
                   br(),
                   p(
                     "That's latin and roughly translates to 'user-friendly software for MS-based proteomics data analysis'.
                                    amica's capabilities can deal with analyzing and visualizing common subjects in proteomics such as QC,
                                    quantitative proteomics and protein-protein interaction networks, but one of amica's
                                    main advantage is its comparability; you can compare different biological conditions
                                    from a single experiment or even multiple similar or distinct conditions from different experiments."
                   ),
                   p(
                     "Upload the required input files (explained on the sidebar) and the full functionality will be revealed."
                   ),
                   br(),
                   br(),
                   uiOutput("uploadSuccessMsg"),
                   uiOutput("analysisSuccessMsg"),
                   verbatimTextOutput("summaryText", placeholder = F),
                   
                   shinyjs::hidden(div(
                     id = 'hide_before_input',
                     bsCollapse(
                       id = "colorCollapse",
                       open = "Panel beneath",
                       
                       bsCollapsePanel(
                         style = "info",
                         "Choose colors v",
                         p(
                           "All visualizations are interactive and created with plotly.
                                                             When you hover over a plot with your mouse a toolbar is visable in the top right corner.
                                                             For some plots (volcano -, MA - and fold change plots) you can select data points with the select box or lasso functions which annotates the data points.
                                                             Plots can be downloaded in svg format when clicking on the camera icon."
                         ),
                         tabsetPanel(
                           type = "tabs",
                           tabPanel(h4("Available color palettes"),
                                    plotOutput("brewercols", height = "800px")),
                           tabPanel(h4("Qualitative (groups)"),
                                        uiOutput("brewerOptionsQual"),
                                        radioButtons(
                                          "revQual",
                                          "Reverse colors?",
                                          choices = c("yes", "no"),
                                          selected = "no"
                                        ),
                                        plotlyOutput("mtcarsBar", height = 800),
                                        helpText(
                                          "The Motor Trend Car Road Tests was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles (1973–74 models)."
                                        )
                                    ),
                           tabPanel(h4("Qualitative (scatter)"),
                                    selectInput(
                                      "brewerOptionsScatter",
                                      "Scatter plot colors",
                                      choices = c(
                                        "RdYlBu",
                                        "RdYlGn",
                                        "Spectral",
                                        "Accent",
                                        "Dark2",
                                        "Paired",
                                        "Pastel1",
                                        "Pastel2",
                                        "Set1",
                                        "Set2",
                                        "Set3"
                                      ),
                                      multiple = F,
                                      selected = "Set2"
                                    ),
                                    radioButtons(
                                      "revScatter",
                                      "Reverse colors?",
                                      choices = c("yes", "no"),
                                      selected = "no"
                                    ),
                                    plotlyOutput("irisScatter", height = 800),
                                    helpText("This famous (Fisher's or Anderson's) iris data set gives
                                             the measurements in centimeters of the variables sepal 
                                             length and width and petal length and width, respectively, 
                                             for 50 flowers from each of 3 species of iris. 
                                             The species are Iris setosa, versicolor, and virginica.")
                                    ),
                           tabPanel(
                             h4("Diverging"),
                             uiOutput("brewerOptionsDiv"),
                             radioButtons(
                               "revDiv",
                               "Reverse colors?",
                               choices = c("yes", "no"),
                               selected = "yes"
                             ),
                             plotlyOutput("mtcarsHeatmap", height = 800),
                             helpText(
                               "The Motor Trend Car Road Tests was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles (1973–74 models)."
                             )
                           )
                         )
                       )
                     )
                   )),
                   br(), br(), br(),
                   #actionButton("showMemory", "Show memory"),
                   # verbatimTextOutput("printMemory", placeholder = F),
                   uiOutput("download_amica"),
                   #p("Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam erat, sed diam voluptua. At vero eos et accusam et justo duo dolores et ea rebum. Stet clita kasd gubergren, no sea takimata sanctus est Lorem ipsum dolor sit amet. Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam erat, sed diam voluptua. At vero eos et accusam et justo duo dolores et ea rebum. Stet clita kasd gubergren, no sea takimata sanctus est Lorem ipsum dolor sit amet."),
                   #verbatimTextOutput("uploadSummary", placeholder = F),
                   
                   br(),br(),br(),
                   h3(textOutput("designTitle")),
                   DTOutput("expDesignDT"),
                   footer()#p("This is footer"))
                 )
               )),
      tabPanel(
        title = 'QC',
        value = 'qctab',
        h1("Quality Control"),
        h2("Compare different intensities"),
        uiOutput("assayNames"),
        helpText("ImputedIntensity are the intensities used for quantification"),
        verbatimTextOutput("quantSummary"),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            h3("PCA"),
            uiOutput("pcaSamplesInput"),
            actionButton("submitPCA", "Plot PCA", icon = icon("cog")),
            actionButton("pcaParams", "", icon = icon("wrench")),
            shinyjs::hidden(div(
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
              ),
              checkboxInput(
                "pca_show_label",
                "Show labels of samples?",
                value = FALSE
              ),
              
            )),
            fluidRow(
              column(width = 2),
              column(width = 8,
                     plotlyOutput("pca", height = 600)),
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
                       
                       plotlyOutput("boxPlot", height = 600)
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
                       plotlyOutput("densityPlot", height = 600)
                     )
                   )),
          tabPanel(
            h3("Coefficient of variation"),
            inline(actionButton("submitCVs", "Plot CVs", icon = icon("cog"))),
            inline(actionButton("cvHelp", icon = icon("info") , label = NULL)),
            inline(uiOutput("cvHelpBox")),
            ###
            inline(actionButton("cvParams", "", icon = icon("wrench"))),
            shinyjs::hidden(
              div(
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
            plotlyOutput("boxplotCV", height = 600)
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
                       shinyjs::hidden(div(
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
                         numericInput(
                           "abundant_base",
                           "Base font size in pt.",
                           value = 14,
                           min = 4,
                           max = 32,
                           step = 1
                         ),
                         numericInput(
                           "abundant_legend",
                           "Legend font size in pt.",
                           value = 10,
                           min = 4,
                           max = 16,
                           step = 1
                         )
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
                
                inline( actionButton("submitNumMVs", "Plot barplot", icon = icon("cog"))),
                inline(actionButton("mvHelp", icon = icon("info") , label = NULL)),
                inline(uiOutput("mvHelpBox")),
                ###
                inline(actionButton("barplotMvParams", "", icon = icon("wrench"))),
                shinyjs::hidden(
                  div(
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
            width = 8,
            plotlyOutput("scatterPlot", height = 600)
          )
        ),
        footer()
      ),
      
      tabPanel(
        title = 'Differential abundance',
        value = 'quanttab',
        h2("Analyze differentially abundant proteins"),
        p(
          "Get an overview, compare multiple conditions, find proteins specific to certain conditions and find functionally related biological terms."
        ),
        wellPanel(
          h3("Global parameters"),
          p(
            "Select the fold change threshold and whether to use multiple-testing correction or rely on your gut feeling for all subsequent analysis in this tab."
          ),
          fluidRow(
            column(
              width = 4,
              numericInput(
                "fcCutoff",
                "Fold change threshold",
                value = 1.5,
                min = 0,
                step = 0.01,
                width = '30%'
              ),
              helpText(
                "While the choice is arbitrary, lower thresholds might result in more false positives. Usual choices are between 1 and 2. You can lower the threshold and try to make sense of the resulting list of proteins (e.g if enriched terms make sense in the over-representation analysis below)."
              )
            ),
            column(
              width = 4,
              radioButtons(
                "sigCutoffValue",
                "Significance cutoff (which value to use)",
                choices = c("adj.p-value", "p-value", "none"),
                selected = "adj.p-value"
              ),
              helpText(
                "Adjusted p-value is (highly) recommended. Raw p-value cutoffs yield many false positives and give a general trend, not any statistical significance."
              )
            ),
            column(
              width = 4,
              radioButtons(
                "enrichmentChoice",
                "Which proteins to chose?",
                choices = c("enriched", "absolute", "reduced"),
                selected = "absolute"
              ),
              helpText(
                "If 'enriched' is selected only proteins above the FC threshold are selected. 'absolute' selects proteins above and below that negative threshold (e.g if FC-threshold equals 2 'absolute' selects proteins in the ranges [-inf - (-2)] and [2 - inf]. 'reduced' return all significant proteins below the negative FC cutoff."
              )
            )
          ),
          
        ),
        h3("Summary"),
        p(
          "Compare intersections of your data, visualize your differentially abundant proteins, 
          or select the gene name of a single protein to inspect its intensities in your different 
          conditions."
        ),
        tabsetPanel(
          type = "pills",
          
          tabPanel(
            h4("Analyze single comparison"),
            
            fluidRow(
              column(
                width = 3,
                uiOutput("maVolcanoSample"),
                inline(actionButton("maVolcanoSubmit", label = "Submit Analysis")),
                inline(actionButton("volcanoHelp", icon = icon("info") , label = NULL)),
                inline(uiOutput("VolcanoHelpBox")),
                inline(actionButton("volcanoParams", "", icon = icon("wrench"))),
                ###
                shinyjs::hidden(
                  div(
                    id = 'toggle_volcano_params',
                    numericInput(
                      "volcano_width",
                      "Width in pixel.",
                      value = 768,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "volcano_height",
                      "Height in pixel.",
                      value = 676,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "volcano_base",
                      "Base font size in pt.",
                      value = 12,
                      min = 4,
                      max = 32,
                      step = 1
                    ),
                    numericInput(
                      "volcano_legend",
                      "Legend font size in pt.",
                      value = 10,
                      min = 4,
                      max = 16,
                      step = 1
                    ),
                    numericInput(
                      "volcano_pointsize",
                      "Scatter point size.",
                      value = 2,
                      min = 1,
                      max = 4,
                      step = 1
                    ),
                    radioButtons(
                      'volcano_padj_y',
                      'Which p-values to plot on y-axis?',
                      choices = c('p-values',
                                  'adj. p-values'),
                      selected = 'p-values'
                    )
                  )
                )
              ),
              column(width = 9, )
            ),
            
            
            shinyjs::hidden(div(
              id = 'hide_before_single_submit',
              #inline(actionButton("maVolcanoPlot", label = "Update plots")),
              
              ###
              fluidRow(
                column(
                  width = 6,
                  h2("Volcano plot"),
                  plotlyOutput("volcanoPlot", height = 600)
                ),
                column(
                  width = 6,
                  h2("MA plot"),
                  plotlyOutput("maplot", height = 600)
                )
              ),
              br(), br(), br()
            ))
          ),
          tabPanel(
            h4("Analyze multiple comparisons"),
            
            fluidRow(
              column(
                width = 6,
                uiOutput("upset1Sample"),
                helpText("Select at least 2 comparisons."),
                inline(actionButton("submitMultiComp", label = "Submit Analysis")),
                inline(actionButton("upsetHelp", icon = icon("info"), label = NULL)),
                inline(uiOutput("UpsetHelpBox")),
                inline(actionButton("upsetParams", "", icon = icon("wrench"))),
                ###
                shinyjs::hidden(
                  div(
                    id = 'toggle_upset_params',
                    numericInput(
                      "upset_width",
                      "Width in inch.",
                      value = 12,
                      min = 3,
                      max = 18,
                      step = 1
                    ),
                    numericInput(
                      "upset_height",
                      "Height in inch.",
                      value = 6,
                      min = 3,
                      max = 18,
                      step = 1
                    ),
                    numericInput(
                      "upset_pointsize",
                      "Size of points in matrix plot.",
                      value = 4,
                      min = 1,
                      max = 8,
                      step = 1
                    ),
                    numericInput(
                      "upset_ratio",
                      "Ratio between main bar plot to matrix plot.",
                      value = 0.6,
                      min = 0.3,
                      max = 0.7,
                      step = 0.1
                    ),
                    radioButtons(
                      'upset_sorted',
                      'How to sort bars?',
                      choices = c('Frequency', 'Degree'),
                      selected = 'Frequency'
                    )
                  )
                )
              )
            ),
            shinyjs::hidden(
              div(
                id = 'hide_before_multi_submit',
                #inline(actionButton("plotMultiComp", label = "Update plot")),
                
                
                plotOutput('upsetPlot', height = 600),
                uiOutput("download_button_upset"),
                br(),
                br(),
                br()
              )
            )
          )
        ),
        br(),
        #################
        h3(textOutput("compSummary")),
        br(),
        verbatimTextOutput("parameterSummary"),
        br(),
        br(),
        br(),
        shinyjs::hidden(
          div(
            id = 'hide_before_comparisons',
        p(
          "The output table can be further filtered.
          You can for example use regular expressions (ProteinA|ProteinB|ProteinC)
          to select only proteins of interest or apply additional column filters.
          Subsequent visualizations (heatmap, fold change plot, PPI network) 
          and ORA will be computed on the proteins in the output table."
        ),
        
        div(style = 'overflow-x: scroll; max-width: 100%',
            DTOutput("groupComparisonsDT")),
        br(),
        
        
            verbatimTextOutput("filterDTSummary"),
            tabsetPanel(
              type = "tabs",
              tabPanel(
                h3("Heatmap"),
                p(
                  "Press 'Submit' to visualize your selection as a heatmap. Imputed intensities are used."
                ),
                div(
                  style = "display:inline-block",
                  uiOutput("heatmapSamplesInput"),
                  helpText(
                    "Select at least 2 groups If none is selected all comparisons are considered."
                  ),
                  actionButton("heatmapParams", "", icon = icon("wrench")),
                  shinyjs::hidden(
                    div(
                      id = "toggle_heatmap_params",
                      checkboxInput("clusterRows",
                                    "Cluster rows (proteins) of heatmap?",
                                    value = TRUE),
                      checkboxInput("clusterCols",
                                    "Cluster columns (samples) of heatmap?",
                                    value = TRUE),
                      checkboxInput("heatmap_labels",
                                    "Show gene names of proteins?",
                                    value = TRUE),
                      radioButtons(
                        "scaleHeatmap",
                        "Scale rows (proteins) or columns (samples) on z-scores? If none is selected intensity values are plotted.",
                        choices = c("row", "col", "none"),
                        selected = "row"
                      ),
                      numericInput(
                        "heatmap_width",
                        "Width in pixel.",
                        value = 768,
                        min = 338,
                        max = 1352,
                        step = 10
                      ),
                      numericInput(
                        "heatmap_height",
                        "Height in pixel.",
                        value = 676,
                        min = 338,
                        max = 1352,
                        step = 10
                      ),
                      numericInput(
                        "heatmap_base",
                        "Base font size in pt.",
                        value = 14,
                        min = 4,
                        max = 32,
                        step = 1
                      )
                    )
                  ),
                  actionButton("submitHeatmap", strong("Submit"))
                ),
                shinyjs::hidden(div(
                  id = 'hide_heatmap_before_submit',
                  plotlyOutput("compareHeatmap", height = 800), #height = "auto"),
                ))
              ),
              
              tabPanel(
                h3("Fold change plot"),
                inline(uiOutput("foldChangeSelection")),
                inline(actionButton("sumbitFoldChangePlot", label = "Compare fold changes")),
                inline(actionButton(
                  "FoldChangePlotHelp",
                  icon = icon("info"),
                  label = NULL
                )),
                inline(uiOutput("FoldChangePlotHelpBox")),
                ###
                inline(actionButton("fcParams", "", icon = icon("wrench"))),
                shinyjs::hidden(
                  div(
                    id = 'toggle_fc_params',
                    numericInput(
                      "fc_width",
                      "Width in pixel.",
                      value = 768,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "fc_height",
                      "Height in pixel.",
                      value = 676,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "fc_base",
                      "Base font size in pt.",
                      value = 14,
                      min = 4,
                      max = 32,
                      step = 1
                    ),
                    numericInput(
                      "fc_pointsize",
                      "Scatter point size.",
                      value = 2,
                      min = 1,
                      max = 4,
                      step = 1
                    )
                  )
                ),
                ###
                fluidRow(
                  column(width = 1),
                  column(width = 10,
                         plotlyOutput("foldChangePlot", height = 600)),
                  column(width = 1)
                )
                
              ),
              
              tabPanel(
                h3("Profile plot"),
                inline(selectizeInput(
                  inputId = "selectProfilePlotGene",
                  label = ("Select gene for profile plot"),
                  choices = NULL,
                  multiple = F,
                  options = list(maxItems = 1,
                                 placeholder =
                                   "Please choose...")
                )),
                
                helpText(
                  "The profile plot shows the imputed intensities with their mean and standard error."
                ),
                inline(actionButton("profileParams", "", icon = icon("wrench"))),
                shinyjs::hidden(
                  div(
                    id = 'toggle_profile_params',
                    numericInput(
                      "profile_width",
                      "Width in pixel.",
                      value = 768,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "profile_height",
                      "Height in pixel.",
                      value = 676,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "profile_base",
                      "Base font size in pt.",
                      value = 14,
                      min = 4,
                      max = 32,
                      step = 1
                    ) ,
                    numericInput(
                      "profile_pointsize",
                      "Point size.",
                      value = 3,
                      min = 1,
                      max = 8,
                      step = 1
                    )
                  )
                ),
                ###
                plotlyOutput("profilePlot", height = 600),
                div(style = 'overflow-x: scroll; max-width: 100%',
                    DTOutput("geneSummary"))
              ),
              tabPanel(
                h3("PPI Network"),
                inline(actionButton("networkHelp", icon = icon("info") , label = NULL)),
                inline(uiOutput("NetworkHelpBox")),
                inline(
                  numericInput(
                    "edgeWeightThresh",
                    "Min edge weight",
                    value = 0,
                    min = 0,
                    max = 1,
                    step = 0.1,
                    width = '50%'
                  )
                ),
                inline(uiOutput("download_button_spec_network")),
                visNetworkOutput("network", width = "100%", height = "600px")
              )
            ),
            
            h3("Over-Representation Analysis (ORA)"),
            wellPanel(fluidRow(
              column(
                width = 6,
                checkboxInput("showGenes", "Show genes in functional enrichment?"),
                helpText("Only select this feature if your gene set isn't too large."),
                checkboxInput("significantORA", "Only show significant terms?", value = TRUE),
                helpText(
                  "Only deselect this box if you are certain. The running time can increase dramatically if your gene list is too long."
                ),
              ),
              column(width = 6,
                     radioButtons(
                       "species",
                       "Organism",
                       choices = c(
                         "hsapiens",
                         "mmusculus",
                         "scerevisiae",
                         "athaliana",
                         "celegans",
                         "dmelanogaster"
                       ),
                       selected = "hsapiens"
                     ),
                     actionButton("submitORA", strong("Submit"))),
            )),
            shinyjs::hidden(
              div(
                id = 'hide_ora_before_submit',
                
                tabsetPanel(
                  type = "tabs",
                  tabPanel(
                    h3("Manhattan plot"),
                    fluidRow(
                      column(width = 1),
                      column(width = 10,
                             plotlyOutput("gostplot", height = 600)),
                      column(width = 1)
                    )
                  ),
                  tabPanel(
                    h3("Bar plot"),
                    radioButtons(
                      "orasource",
                      "Source",
                      choices = c(
                        "GO:MF",
                        "GO:CC",
                        "GO:BP",
                        "KEGG",
                        "REAC",
                        "CORUM",
                        "WP"
                      ),
                      selected = "GO:MF"
                    ),
                    colourInput("oraBar_color", "Select bar color", "#66c2a5"),
                    actionButton("submitORABar", strong("Submit")),
                    plotlyOutput("oraBarplot", height = 800)
                  )
                ),
                br(),
                uiOutput("download_button_ora"),
                #downloadLink("download1","Download table"),
                div(style = 'overflow-x: scroll; max-width: 100%',
                    DTOutput("gprofilerDT"))
              )
            )
          )
        ),
        footer()
      ),
      tabPanel(
        title = 'Compare amica datasets',
        value = 'comparemicatab',
        
        bsCollapse(
          id = "collapseUpload",
          open = "Upload amica",
          bsCollapsePanel("Toggle Upload",
                          style = "info",
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
                                fileInput("amicaFile2", "Upload amica_proteinGroups.txt."),
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
                          ))
        ),
        
        
        verbatimTextOutput("summaryMergedAmica"),
        
        
        shinyjs::hidden(
          div(
            id = "compareAmicaInput",
            
            uiOutput('download_merged_amica'),
            br(),
            br(),
            br(),
            
            fluidRow(
              column(
                width = 6,
                numericInput(
                  "fcThresh3",
                  "Fold change threshold",
                  value = 1.5,
                  min = 0,
                  step = 0.01
                ),
                helpText(
                  "While the choice is arbitrary, lower thresholds might result in more false positives. Usual choices are between 1 and 2. You can lower the threshold and try to make sense of the resulting list of proteins (e.g if enriched terms make sense in the over-representation analysis below)."
                ),
                
                radioButtons(
                  "enrichmentChoice3",
                  "Which proteins to chose?",
                  choices = c("enriched", "absolute", "reduced"),
                  selected = "absolute"
                ),
                helpText(
                  "If 'enriched' is selected only proteins above the FC threshold are selected. 'absolute' selects proteins above and below that negative threshold (e.g if FC-threshold equals 2 'absolute' selects proteins in the ranges [-inf - (-2)] and [2 - inf]. 'reduced' return all significant proteins below the negative FC cutoff."
                )
                
              ),
              column(
                width = 6,
                radioButtons(
                  "sigCutoffValue3",
                  "Significance cutoff (which value to use)",
                  choices = c("p-value", "adj.p-value"),
                  selected = "adj.p-value"
                ),
                helpText(
                  "Adjusted p-value is recommended. Raw p-value cutoffs yield many false positives and give a general trend, not any statistical significance."
                )
              )
            ),
            
            # inline(actionButton("plotMultiComp", label = "Update plot")),
            # inline(actionButton("upsetHelp", icon = icon("info"), label = NULL)),
            # inline(uiOutput("UpsetHelpBox")),
            
            uiOutput('upsetSelection2'),
            inline(actionButton("submitMultiComp2", label = "Submit Analysis")),
            inline(actionButton("upsetHelp2", icon = icon("info"), label = NULL)),
            inline(uiOutput("UpsetHelpBox2")),
            inline(actionButton("upsetParams2", "", icon = icon("wrench"))),
            shinyjs::hidden(
              div(
                id = 'toggle_upset_amicas_params',
                numericInput(
                  "upset_amicas_width",
                  "Width in inch.",
                  value = 12,
                  min = 3,
                  max = 18,
                  step = 1
                ),
                numericInput(
                  "upset_amicas_height",
                  "Height in inch.",
                  value = 6,
                  min = 3,
                  max = 18,
                  step = 1
                ),
                numericInput(
                  "upset_amicas_pointsize",
                  "Size of points in matrix plot.",
                  value = 4,
                  min = 1,
                  max = 8,
                  step = 1
                ),
                numericInput(
                  "upset_amicas_ratio",
                  "Ratio between main bar plot to matrix plot.",
                  value = 0.6,
                  min = 0.3,
                  max = 0.7,
                  step = 0.1
                ),
                radioButtons(
                  'upset_amicas_sorted',
                  'How to sort bars?',
                  choices = c('Frequency', 'Degree'),
                  selected = 'Frequency'
                )
              )
            ),
            plotOutput('upsetPlotAmicas', height = 600),
            uiOutput("download_button_upset_amicas"),
            br(),
            br(),
            br(),
            helpText(
              "Set comparison of differentially abundant proteins from selected comparisons under selected thresholds. The dots show which sets are getting compared. A dot not connected to another dot shows the number of proteins specific to that comparisons. The top barplot depicts the number of intersecting proteins, and the barplot on the side shows how many proteins are differentially abundant in the comparison. Change the selected comparisons to your needs."
            ),
            h3(textOutput("compSummary2")),
            br(),
            verbatimTextOutput("parameterSummary2"),
            br(),
            div(style = 'overflow-x: scroll; max-width: 100%',
                DTOutput("groupComparisonsDTamica")),
            br(),
            br(),
            br(),
            
            tabsetPanel(
              type = "tabs",
              tabPanel(
                h3("Fold change plot"),
                uiOutput("foldChangeSelection2"),
                radioButtons(
                  'imputeScatterBoolean',
                  'Replace missing fold changes with constant value?',
                  choices = c('yes', 'no'),
                  selected = 'no'
                ),
                helpText(
                  'If enabled proteins not quantified in one condition can be visualized in the fold change plot.'
                ),
                conditionalPanel(
                  condition = "input.imputeScatterBoolean == 'yes'",
                  numericInput(
                    "imputeScatterPlot",
                    "Imputation value",
                    value = 0,
                    min = -15,
                    max = 15,
                    step = 0.1
                  )
                ),
                actionButton("sumbitFoldChangePlot2", label = "Compare fold changes"),
                actionButton(
                  "FoldChangePlotHelp2",
                  icon = icon("info"),
                  label = NULL
                ),
                uiOutput("FoldChangePlotHelpBox2"),
                ###
                actionButton("fcamicaParams", "", icon = icon("wrench")),
                shinyjs::hidden(
                  div(
                    id = 'toggle_fcamica_params',
                    numericInput(
                      "fcamica_width",
                      "Width in pixel.",
                      value = 768,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "fcamica_height",
                      "Height in pixel.",
                      value = 676,
                      min = 338,
                      max = 1352,
                      step = 10
                    ),
                    numericInput(
                      "fcamica_base",
                      "Base font size in pt.",
                      value = 14,
                      min = 4,
                      max = 32,
                      step = 1
                    )
                  )
                ),
                ###
                
                fluidRow(
                  column(width = 1),
                  column(width = 10,
                         plotlyOutput("foldChangePlot2", height = 600)),
                  column(width = 1)
                )
                
              ),
              tabPanel(
                h3('Scatter plot'),
                uiOutput("compareScatterPlotsAmica"),
                ###
                actionButton("scatteramicaParams", "", icon = icon("wrench")),
                shinyjs::hidden(
                  div(
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
                    )
                  )
                ),
                ###
                plotlyOutput("scatterPlotsAmica", height = 600)
              )
            )
          )
        ),
        footer()
      )
    )
  )