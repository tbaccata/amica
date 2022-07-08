tabPanel(
  title = 'Differential abundance',
  value = 'quanttab',
  h2("Analyze differentially abundant proteins"),
  p(
    "Get an overview, compare multiple conditions, find proteins specific 
          to certain conditions and find functionally related biological terms."
  ),
  br(),br(),br(),
  actionButton("showTutorial", "Tutorial", icon = icon("info") ),
  br(),br(),br(),
  conditionalPanel(
    "output.multiAmicasInput",
    selectInput("selectedDataset", 
                "Which data set to analyze?",
                choices = c("original_data",
                            "multiple_amica_data"),
                selected = "original_data",
                multiple = F),
    actionButton("submitDatasetSelection", "Submit"),
    verbatimTextOutput("summaryLoadedDataSet")
  ),
  br(),br(),br(),
  wellPanel(
    h3("Global parameters"),
    
    p(
      "In this section you can select parameters to subset your data to show only protein groups of your interest. 
            The selected parameters are applied to all selected group comparisons and can be further filtered in the output table.
            Following filters can be set:"
    ),
    HTML("<p><b>Fold change threshold</b><br>
               <b>Significance cutoff</b><br>
               <b>Enriched, reduced or absolute fold change.</b><br>
               </p>"),
    
    fluidRow(
      column(
        width = 4,
        numericInput(
          "fcCutoff",
          "Log2 fold change threshold",
          value = 1.5,
          min = 0,
          step = 0.01,
          width = '30%'
        ),
        helpText(
          "While the choice is arbitrary, lower thresholds might result in more false positives. 
                Usual choices are between 1 and 2. You can lower the threshold and try to make sense 
                of the resulting list of proteins (e.g if enriched terms make sense in the over-representation analysis below)."
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
          "Adjusted p-value is (highly) recommended. Raw p-value cutoffs yield 
                many false positives and give a general trend, not any statistical significance. 
                Default threshold is 0.05, but this can still be lowered in the output table.
                In case you want to visualize proteins only based on fold changes, 
                you can also set no p-value significance cutoff."
        ),
        conditionalPanel(
          condition = "input.sigCutoffValue != 'none'",
          numericInput(
            "pvalCutoff",
            "p-value threshold",
            value = 0.05,
            min = 0.000001,
            step = 0.001,
            max = 0.0501,
            width = '30%'
          )
        ), 
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
          "If 'enriched' is selected only proteins above the FC threshold are selected.
                This is useful for interaction proteomics experiments, when you have group comparisons like Bait - Control 
                and only want to see enriched interaction partners vs a control.
                
                'absolute' selects proteins above and below that negative threshold 
                (e.g if FC-threshold equals 2 'absolute' selects proteins in the ranges [-inf - (-2)] and [2 - inf]. 
                
                'reduced' return all significant proteins below the negative fold change threshold."
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
              style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
              
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
              radioButtons(
                "volcano_format",
                "Save image",
                choices = c("svg", "png"),
                selected = "svg"
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
              ),
              uiOutput('volcanoMAColors')
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
      
      width = 6,
      uiOutput("upset1Sample"),
      
      helpText("Select at least 2 comparisons."),
      inline(actionButton("submitMultiComp", label = "Submit Analysis")),
      shinyjs::hidden(
        div(
          id = 'hide_before_multi_submit',
          #inline(actionButton("plotMultiComp", label = "Update plot")),
          fluidRow(
            
            column(width = 6,
                   h4("UpSet plot"),
                   br(),
                   inline(actionButton("upsetHelp", icon = icon("info"), label = NULL)),
                   inline(uiOutput("UpsetHelpBox")),
                   inline(actionButton("upsetParams", "", icon = icon("wrench"))),
                   shinyjs::hidden(
                     div(
                       style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                       
                       id = 'toggle_upset_params',
                       numericInput(
                         "upset_width",
                         "Width in inch.",
                         value = 8,
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
                         "upset_scale",
                         "Scale text sizes",
                         value = 2,
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
                   ),
                   inline(
                     actionButton("multiNameChange", "Change labels?", icon = icon("wrench"))
                   ),
                   shinyjs::hidden(
                     div(
                       id = 'toggle_multi_name_change',
                       uiOutput("multiCompLabelsInput"),
                       actionButton("changeMultiCompNames", "Change labels")
                     )
                   ),
                   plotOutput('upsetPlot', height = 600),
                   uiOutput("download_button_upset")
            ),
            column(
              width = 6,
              h4("Euler diagram"),
              br(),
              #inline(actionButton("upsetHelp", icon = icon("info"), label = NULL)),
              # inline(uiOutput("UpsetHelpBox")),
              inline(actionButton("eulerParams", "", icon = icon("wrench"))),
              ###
              shinyjs::hidden(
                div(
                  style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                  
                  id = 'toggle_euler_params',
                  numericInput(
                    "euler_width",
                    "Width in inch.",
                    value = 6,
                    min = 3,
                    max = 18,
                    step = 0.5
                  ),
                  numericInput(
                    "euler_height",
                    "Height in inch.",
                    value = 6,
                    min = 3,
                    max = 18,
                    step = 0.5
                  ),
                  checkboxInput(
                    "euler_quant",
                    "Show quantities in circles?",
                    value = TRUE
                  ),
                  checkboxInput(
                    "euler_legend",
                    "Show legend next to plot?",
                    value = TRUE
                  ),
                  checkboxInput(
                    "euler_line",
                    "Plot circle outline?",
                    value = FALSE
                  ),
                  uiOutput('eulerColors'),
                )
              ),
              plotOutput("eulerrPlot", height = 600),
              uiOutput("download_button_eulerr")
            )
          ),
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
      
      fluidRow(
        column(width = 3),
        column(width = 6,
               HTML(
                 "<p><b>The output table can be further filtered</b>.
          You can for example use regular expressions (e.g type in a query like 
          ProteinA|ProteinB|ProteinC in the Gene.names column)
          to select only proteins of interest or apply additional column filters.
          When multiple group comparisons are selected you can subset the proteins that 
          are only specific to one group comparison (e.g by writing 'yes' in one significant column 
          and 'no' in the other column.)
          </p>
          <p>Subsequent visualizations (heatmap, fold change plot, PPI network) 
          and ORA will be computed on the proteins in the output table. </p>"
               ),
               actionButton("showQueryTutorial", "Advanced Queries", icon = icon("info") ),
        ),
        column(width = 3)
      ),
      
      div(style = 'overflow-x: scroll; max-width: 100%',
          DTOutput("groupComparisonsDT")),
      br(),br(),br(),
      h4(textOutput("filterDTSummary")),
      br(),br(),br(),
      # verbatimTextOutput("filterDTSummary"),
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
                style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                
                id = "toggle_heatmap_params",
                checkboxInput("clusterRows",
                              "Cluster rows (proteins) of heatmap?",
                              value = TRUE),
                checkboxInput("clusterCols",
                              "Cluster columns (samples) of heatmap?",
                              value = TRUE),
                checkboxInput("heatmap_row_labels",
                              "Show gene names of proteins?",
                              value = TRUE),
                checkboxInput("heatmap_col_labels",
                              "Show samples names at the bottom?",
                              value = TRUE),
                checkboxInput("heatmap_annot",
                              "Annotate samples in dendrogram?",
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
                radioButtons(
                  "heatmap_format",
                  "Save image",
                  choices = c("svg", "png"),
                  selected = "svg"
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
            actionButton("submitHeatmap", strong("Submit")),
            inline(uiOutput("download_button_heatmap"))
          ),
          shinyjs::hidden(div(id = 'hide_heatmap_before_submit',
                              fluidRow(
                                column(width = 2),
                                column(
                                  width = 8,
                                  div(style='scroll;height:800px;overflow-y: scroll;',
                                      plotlyOutput("compareHeatmap")
                                  )
                                  #height = "auto"),
                                ),
                                column(width = 2)
                              )))
        ),
        tabPanel(
          h3("Dotplot"),
          helpText("Dot plots integrate quantitive information together with 
                   log2 fold changes and (adj.) p-values into one visualization. 
                   Proteins are displayed as dots, with their circle size 
                   corresponding to relative abundance (average intensities or 
                   fold changes can be selected). Fold changes are mapped as 
                   color gradients on the dots. When proteins are clustered 
                   based on fold changes, users have the option to display 
                   proteins with a positive fold change only."),
          br(),br(),
          fluidRow(
            column(width = 4,
                   wellPanel(
                     uiOutput("dotplotGroupComparisons"),
                     uiOutput("dotplotGroups"),
                     inline(actionButton("dotplotSelection", label = "Select Groups")),
                     uiOutput("boolUniqueGroups")
                   )
                   ),
            column(
              width = 4,
              uiOutput("submitDotplot"),
              div(style = "position:relative",
                  uiOutput('plot.ui'),
                  uiOutput("hover_info"))
            ), 
            column(
              width = 4,
              shinyjs::hidden(
                div(
                  id = 'hide_dotplot_before_submit',
                  uiOutput("dotplot_color_gradient"),
                  helpText(
                    "Specify a min. and a max. fold change.
                       All fold changes out side the specified
                       range will be capped."
                  ),
                  uiOutput("dotplot_size_gradient"),
                  helpText(
                    "The size of dotplots are scaled to the maximum
                       avg. intensity or maximum fold change of the proteins
                       in your selection."
                  ),
                  h3("Clustering"),
                  uiOutput("dotplot_clustering_option"),
                  helpText(
                    "Which quantitive information should be displayed as
                       circle size? This value will also be used to cluster the
                       Dot plot."
                  ),
                  uiOutput("dotplot_ctrl_substraction"),
                  inline(
                    selectInput(
                      "dotplot_distance_metric",
                      "Distance method",
                      choices = c(
                        "canberra",
                        "euclidean",
                        "maximum",
                        "manhattan",
                        "binary",
                        "minkowski"
                      ),
                      selected = "canberra"
                    )
                  ),
                  inline(
                    selectInput(
                      "dotplot_clustering_method",
                      "Clustering method",
                      choices = c("complete", "average", "single", "ward.D", "median", "centroid"),
                      selected = "complete"
                    )
                  ),
                  uiOutput("dotplot_cluster_columns"),
                  helpText(
                    "If column clustering is disabled, columns are displayed
                       in the order they were inputted."
                  ),
                  h3("Colors"),
                  div(
                    id = 'dotplot_colors',
                    selectInput(
                      "dotplot_palette",
                      "Color palette",
                      choices = list(
                        viridis = c("viridis",
                                    "inferno",
                                    "plasma",
                                    "magma",
                                    "cividis"),
                        diverging = c(
                          "BrBG",
                          "PiYG",
                          "PRGn",
                          "PuOr",
                          "RdBu",
                          "RdGy",
                          "RdYlBu",
                          "RdYlGn",
                          "Spectral"
                        ),
                        sequential = c(
                          "Blues",
                          "BuGn",
                          "BuPu",
                          "GnBu",
                          "Greens",
                          "Greys",
                          "Oranges",
                          "OrRd",
                          "PuBu",
                          "PuBuGn",
                          "PuRd",
                          "Purples",
                          "RdPu",
                          "Reds",
                          "YlGn",
                          "YlGnBu",
                          "YlOrBr",
                          "YlOrRd"
                        )
                      ),
                      selected = "viridis"
                    ),
                    checkboxInput("dotplot_rev_colors",
                                  "Reverse color palette?",
                                  value = FALSE)
                  ),
                  br(),
                  br(),
                  downloadButton('downloadDotPlot', 'Download Plot'),
                  inline(downloadButton('downloadDotPlotData', 'Download data'))
                )
              )
            )
          )
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
              style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
              
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
              radioButtons(
                "fc_format",
                "Save image",
                choices = c("svg", "png"),
                selected = "svg"
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
                "fc_legend",
                "Legend font size in pt.",
                value = 10,
                min = 4,
                max = 16,
                step = 1
              ),
              numericInput(
                "fc_pointsize",
                "Scatter point size.",
                value = 2,
                min = 1,
                max = 4,
                step = 1
              ),
              radioButtons(
                "fc_showLine",
                "Show line in plot?",
                choices = c("straight line", "none", "linear regression"),
                selected = "straight line"
              ),
              uiOutput('fcPlotColors'),
              uiOutput("fcPlotLabelsInput")
            )
          ),
          ###
          fluidRow(
            column(width = 1),
            column(width = 10,
                   plotlyOutput("foldChangePlot", height = 800)),
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
              style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
              
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
              radioButtons("profile_plot_type",
                           "Plot type",
                           choices = c("error_bars",
                                       "violin",
                                       "data_points"),
                           selected = "error_bars"),
              checkboxInput(
                "profile_add_points",
                "Add individual data points to plot?",
                value = T
              ),
              numericInput(
                "profile_jittersize",
                "Jitter Point size.",
                value = 2,
                min = 1,
                max = 8,
                step = 1
              ),
              radioButtons(
                "profile_format",
                "Save image",
                choices = c("svg", "png"),
                selected = "svg"
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
                "Error bar mean point size.",
                value = 0,
                min = 0,
                max = 8,
                step = 1
              )
            )
          ),
          ###
          plotlyOutput("profilePlot", height = 800),
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
          visNetworkOutput("network", width = "100%", height = "800px"),
          inline(uiOutput("download_button_spec_network")),
          inline(uiOutput('download_network_button')),
          inline(actionButton("showNodeTable", strong("Show Node table"))),
          br(),br(),br(),
          shinyjs::hidden(div(id = 'toggle_node_table',
                              h4("Node table"),
                              DTOutput("networkDT")))
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
            "Only deselect this box if you are certain. 
                  The running time can increase dramatically if your gene list is too long."
          ),
          checkboxInput("oraCustom", "Use quantified proteins as custom background?", value = FALSE),
          helpText(
            "A defined background gene set is required in order to determine 
            over-represented functional terms. By default, all protein-coding 
            genes are taken as background genes. Setting only quantified proteins 
            as custom backgorund can be important in some cases, e.g when only a 
            subset of proteins can be measured in an experiment. Examples include 
            tissue-specific experiments (e.g liver-proteomics)."
          )
        ),
        column(width = 6,
               inline(selectizeInput(
                 inputId = "gprofilerOrganism",
                 label = ("Select Organism"),
                 choices = NULL,
                 multiple = F,
                 options = list(maxItems = 1,
                                placeholder =
                                  "Please choose..."),
                 width = '250px'
               )),
               helpText("Please enter the scientific name by concatenating
                              the first letter of the name and the family name. 
                              Example: human - 'hsapiens', mouse - 'mmusculus'."),
               verbatimTextOutput("organismSources", placeholder = F),
               checkboxGroupInput(
                 "oraSources",
                 "Selected sources",
                 choices = c(
                   "GO:MF",
                   "GO:CC",
                   "GO:BP",
                   "REAC",
                   "KEGG",
                   "CORUM",
                   "WP",
                   "TF",
                   "MIRNA",
                   "HP",
                   "HPA"
                 ),
                 selected = c(
                   "GO:MF",
                   "GO:CC",
                   "GO:BP",
                   "REAC",
                   "KEGG",
                   "CORUM",
                   "WP"
                 ),
                 inline = T
               ),
               inline(actionButton("submitORA", strong("Submit"))),
               inline(actionButton('oraHelp', label = '', icon = icon("info"))),
               uiOutput("oraHelpBox")
        )
      )),
      shinyjs::hidden(
        div(
          id = 'hide_ora_before_submit',
          
          tabsetPanel(
            type = "tabs",
            tabPanel(
              h3("Manhattan plot"),
              fluidRow(column(width = 2),
                       column(width = 8,
                              plotlyOutput("gostplot", height = 800)),
                       column(width = 2))
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
                  "REAC",
                  "KEGG",
                  "CORUM",
                  "WP",
                  "TF",
                  "MIRNA",
                  "HP",
                  "HPA"
                ),
                selected = "GO:MF",
                inline = T
              ),
              inline(actionButton("oraBarParams", "", icon = icon("wrench"))),
              shinyjs::hidden(
                div(
                  style = "display: grid; 
          grid-template-columns: 30% repeat(2, 30%); ## same as repeat(4, 20%)
          grid-gap: 30px;",
                  
                  id = 'toggle_oraBar_params',
                  
                  numericInput(
                    "oraBar_width",
                    "Width in pixel.",
                    value = 768,
                    min = 338,
                    max = 1352,
                    step = 10
                  ),
                  numericInput(
                    "oraBar_height",
                    "Height in pixel.",
                    value = 676,
                    min = 338,
                    max = 1352,
                    step = 10
                  ),
                  numericInput(
                    "oraBar_base",
                    "Base font size in pt.",
                    value = 12,
                    min = 4,
                    max = 32,
                    step = 1
                  ),
                  radioButtons(
                    "oraBar_format",
                    "Save image",
                    choices = c("svg", "png"),
                    selected = "svg"
                  ),
                  
                  sliderInput(
                    "oraBar_maxTerms",
                    "Max. number of terms to plot",
                    min = 0,
                    max = 50,
                    step = 1,
                    value = 0
                  ), helpText("If '0' is selected all terms are plotted."),
                  colourInput("oraBar_color", "Select bar color", "#66c2a5")
                )
              ),
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
)