# splash screen that will optionally be shown when the user opens the shiny app
if (amicaGlobalVars$splashScreenEnabled == TRUE) {
  
  htmlContent <- amicaUtil$replacePlaceholders(
    paste(readLines(amicaGlobalVars$splashScreenHtmlContentSource)), 
    amicaGlobalVars$placeholder_replacement_map
    )
  
  showModal(modalDialog(
  title = amicaGlobalVars$splashScreenTitle,
  HTML(htmlContent),
  easyClose = TRUE,
  footer = NULL
  ))
}

observeEvent(input$showFileInput, {
  showModal(modalDialog(
    title = "Accepted File input",
    
    HTML("
      <h4>Inspecting previously analyzed output:</h4>
      <p>
      If you just want to try out the software select the 'Load in example' option 
      and press 'Upload'. After the upload the QC -, Diff. abundance - and Compare 
      amica data sets - tabs open in the main tab bar.
      </p>
      <p>
      If you want to (re-)inspect previously analyzed amica output you can just upload 
      the amica_protein_groups.tsv file together with the experimental design and
      after pressing the 'Upload' button you can also inspect the aforementioned tabs.
      </p>
       <p>
           One property that can be changed in the 'Input' tab are the colors for 
           all upcoming plots. Press on 'Click to choose colors and ordering of 
           groups in plots' to toggle the color selection that is separated into
           different categories, depending on which plots are outputted. The 
           selected colors propagate to all subsequent visualizations (e.g 
           groups from the exp. design will always have the same color in a plot 
           legend.)
           <center><img src='input_tutorial/colorsOrderInput.png' width='100%'></center>
           </p>
      <p>
      <h4>Analyzing a data set</h4>
      Automatically recognized input files are MaxQuant's 'proteinGroups.txt' file and 
FragPipe's 'combined_proteins.txt' file, however you can also upload a generic 
tab-separated format which can easily be processed by the addition of a file 
that maps relevant search engine-specific columns to a standard format.
           </p>
           
           <p>
           After uploading all required files and pressing the 'Upload' button 
           a new field becomes accessible which allows for changing analysis parameters:
           <center><img src='input_tutorial/advanced.png' width='100%'></center>
           </p>
          
           "),
    size = "l",
    easyClose = TRUE,
    footer = NULL
  ))
})


observeEvent(input$showAmicaInput, {
  showModal(modalDialog(
    title = "amica file input",
    
    HTML("<p>
           amica’s tab-separated protein groups file has the following columns:
           </p><br>"), 
    DTOutput("exampleAmicaFile"),
    
    HTML('
      <p>
        <ul style="list-style-type:square">
        <li>IntensityPrefix, ImputedIntensityPrefix and abundancePrefix columns 
        are log2 transformed, all 0s need to be converted to NANs. 
        No INF values allowed. amica searches for all Intensity prefixes in 
        the column names, if you want to provide more than the default intensities.
        However, all intensity prefixes must have the same number of samples in 
        order to get processed.</li>
        <li>ImputedIntensityPrefix should only contain filtered, 
        imputed and normalized values.</li>
        <li>quantCol: All proteins passing spectraCount and 
        razorUniqueCount thresholds that have been quantified are set to "+" in 
        this column. Otherwise no value ("") is written in the
column</li>
        <li>
        comparisonInfix: The infix is important to retrieve the group ids 
        from a group comparison (e.g for downstream visualizations like heatmaps). 
        The groups before and after the "__vs__" infix should match with groups 
        defined in the uploaded experimental design.
        </li>
        <li>
        razorUniqueCount is a column, razorUniquePrefix is the prefix to the 
        count per sample, but they may very well have the same value 
        (just like in MaxQuant’s proteinGroups.txt)
        </li>
        </ul>
      </p><br>
           <p>
           Proteins inferred from reverse hits and peptides ”only identified by 
           site modifications” are not to be written into amica’s output. 
           Additional columns can be added in the future but are at the
moment not considered when uploaded.
           </p>'),
    size = "l",
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$showDesign, {
  showModal(modalDialog(
    title = "Example tab-separated experimental design",
    HTML("<p>
           The tab-separated design file has two columns: <b>samples</b> and <b>groups</b>. The sample names in the samples column need to match 
           the column names of the input file in the order of the input file.</p>"),
    DTOutput("exampleDesign"),
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$showContrasts, {
  showModal(modalDialog(
    title = "Example tab-separated contrast matrix",
    
    HTML("<p>
           The tab-separated contrast matrix tells amica which group comparisons to perform. The column names of
this file can be freely chosen, but column names must be provided. For each row in this file the
comparison group1-group2 is performed. If one wants to change the sign of the fold changes the
position of the groups needs to be switched in the file (e.g group2-group1 instead of group1-group2)</p>"),
    
    
    DTOutput("exampleContrasts"),
    
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$showSpecifications, {
  showModal(modalDialog(
    title = "Specification file",
    
    HTML("<p>
           Following variables can be parsed:
           </p>"),
    
    DTOutput("specificationsExplanation"),
    HTML("<br>The proteinId column must only contain unique entries.
      If razorUnique count is missing some functionality will be lost (DEqMS).<br>
      It is important that the provided intensities are not log2-transformed.
           <br>An example format is provided in the examples.zip file<br>"),
    DTOutput("exampleSpecifications"),
    
    HTML("<p>
           The specification file needs to be uploaded if a custom tab-delimited 
           file is analyzed. The file has two columns, Variable and Pattern, 
           these are used to change the prefixes (or post-fixes) to identify the
           relevant columns in your data.</p>"),
    
    easyClose = TRUE,
    footer = NULL
  ))
})


observeEvent(input$showTutorial, {
  showModal(modalDialog(
    title = "Example use cases",
    
    HTML('
      <h3>Use case 1: Single group comparison</h3>
      <p>
      These small examples have been produced with the provided example data set. 
      As the example data contains AP-MS data we set the selection parameter to 
      "enriched" to retrieve differentially abundant proteins against the control.
      <center><img src="da_tutorial/global_param.png" width="100%"></center>
      Then we select a single group comparison (PRGMC1 bait vs MIA PaCa-2 cell 
      background) and press "Submit Analysis".
      <center><img src="da_tutorial/single_comp.png" width="100%"></center>
      
      When you hover over the volcano - or MA - plot you can see features to 
      manipulate the plot. When we utilize the select box or lasso tool we can 
      annotate the highest enriched proteins, as seen in the figure above.
      </p>
      
      <p>
      Further information on most plots can be acquired when you press the "info"
      icon. Plot parameters can be changed when pressing the "wrench" icon and 
      can be saved upon hovering over the plot and clicking on the "camera" icon.
      </p>
      
      <h3>Use case 2: Multiple group comparisons</h3>
      <p>
      
      When we select the "Analyze multiple comparisons" tab pill we can compare 
      how the prey proteins change with and without the small molecule inhibitor 
      AG-205, just select the two comparisons of the bait versus the controls: 
      
      <center><img src="da_tutorial/multi_comp.png" width="100%"></center>
      </p>
      
      <p>
      Scrolling down we can figure out the quantitive changes of prey proteins 
      with and without AG-205 treatment in a fold change plot:
      <center><img src="da_tutorial/fcplot.png" width="100%"></center>
      </p>
           
           '),
    size = "l",
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$showAmicasTutorial, {
  showModal(modalDialog(
    title = "How to compare a previously analyzed dataset with the current input data",
    
    HTML('
      <p>
      This demonstration example was produced using the example dataset from the 
      Input tab, as well as the "amica_protein_groups.tsv" file from the FragPipe 
      directory in the examples.zip archive. Parameters for parsing and merging 
      both amica files can also be set here:
      </p>
      <p><br><br>
      <center><img src="compare_amicas_tutorial/multi_amicas.png" width="100%"></center>
      </p>
      <br><br><br>
      
      <p>
      After you have pressed the "Submit" button a small summary becomes visible, 
      together with two QC plot options for the combined dataset (scatter plot 
      and correlation plot). If you visit the "Differential abundance" tab a new 
      selection input appears at the top of the page:
      </p>
      <br><br>
      <p>
      <center><img src="compare_amicas_tutorial/multi_amica_data.png" width="100%"></center>
      </p>
      <br><br><br>
      <p>
      You can now use the "Differential abundance" tab select "multiple_amica_input" 
      and press the Submit button. Now you can inspect proteins across different 
      groups and conditions for the combined dataset. Almost the full functionality 
      from the Differential Abundance tab is available for the merged data sets.
      </p>
           '),
    size = "l",
    easyClose = TRUE,
    footer = NULL
  ))
})


observeEvent(input$showQueryTutorial, {
  showModal(modalDialog(
    title = "Advanced queries",
    
    HTML('
      <h3>Use case: Visualize proteins from over-represented functional term</h3>
      <p>
      This small example has been produced with the provided example data set 
      (using the same global parameters as in the above tutorial). 
      An over-representation analysis was conducted utilizing the 46 enriched 
      proteins from the comparison PGRMC1__vs__MIAPACA. The "Show genes in 
      functional enrichment?" button was selected.

      <center><img src="query_tutorial/ora.png" width="100%"></center>
      </p>
      
      <p>
      Sorting the output table from lowest p-value to highest p-value 
      we find the term "actin binding" on top of the list. 15 proteins from the 
      enriched proteins are annotated with this term.
      <center><img src="query_tutorial/oratable.png" width="100%"></center>
      </p>
      
      <p>
      
      All visualizations (heatmap, fold change plot and PPI network) work 
      only on the proteins selected in the above output table and we can filter that 
      table to only show proteins annotated with our term of interest. The output 
      table can parse "regular expressions", so all we need to do is to copy paste 
      the comma-delimited gene names into a text editor (or text processing tool 
      like MS Word) and replace all commas with viertical line symbols ("|" which 
      is the logical "or" operator) with the "Find and Replace" tool:
      
      <center><img src="query_tutorial/edit.png" width="100%"></center>
      </p>
      
      <p>
      We can now paste the the vertical line delimited gene names into the 
      "Gene.names" search bar in the output table and we could successfully subset 
      the data table to only show proteins of our interest:
      <center><img src="query_tutorial/filtereddt.png" width="100%"></center>
      
      Below the table there is now a text message telling us that the original 
      table has been filtered and that only the remaining proteins are used in 
      subsequent visualizations. As an example you can now observe how the selected 
      proteins compare across different group comparisons in a heatmap, fold 
      change plot or specificty/PPI - network (which only works for multiple 
      group comparisons).
      </p>
           '),
    size = "l",
    easyClose = TRUE,
    footer = NULL
  ))
})


### HELP BOXES
output$VolcanoHelpBox <- renderUI({
  if (input$volcanoHelp %% 2) {
    helpText(
      "You can highlight proteins by drawing a box or chosing the lasso option.
               This selection stays highlighted when you change to another volcano plot,
               by selecting a different comparison and pressing the 'Submit Analysis' button."
    )
  } else {
    return()
  }
})

output$assaysHelpBox <- renderUI({
  if (input$assaysHelp %% 2) {
    HTML(
      "
      <p><b>LFQIntensity</b> or <b>Intensity</b> are intensities
      that still contain missing values and potential contaminants.
      These intensities are used to calculate the fraction of missing data
       or the number of identified protein groups in a sample. If the data was quantified
      with MaxQuant or FragPipe these intensities are already normalized.</p>

      <p><b>ImputedIntensity</b> are normalized and imputed intensities
      that are also used to calculate differential abundance. If the re-normalization
      option was selected in the input tab the <b>LFQIntensities</b> were normalized
      after removing potential contaminants, reverse hits, proteins only identified by site
      and protein groups that had too few valid values per group.</p>

           <p><b>iBAQ</b> (intensity-based absolute quantification) values are obtained
           by dividing protein intensities by the number of theoretically observable tryptic peptides.
            This measure correlates well with protein abundance and is for example used to calculate
            the percentage of contamination in a sample (if available).
           </p>

           <p><b>RawIntensity</b> are non-normalized, summed peptide intensities per protein group.</p>
           "
    )
  }
})

output$UpsetHelpBox <- renderUI({
  if (input$upsetHelp %% 2) {
    HTML(
      "<p>Set comparison of differentially abundant proteins from selected comparisons
      under selected thresholds. The dots show which sets are getting compared.
      A dot not connected to another dot shows the number of proteins specific to that comparisons.
      The top barplot depicts the number of intersecting proteins, and the barplot on
      the side shows how many proteins are differentially abundant in the comparison.
               <a href='https://jku-vds-lab.at/tools/upset/' target='_blank'>
      <img src='https://jku-vds-lab.at/assets/images/projects/upset/matrix.png' alt='UpSet plot explained'>
      </a>"
    )
  } else {
    return()
  }
})


output$exampleHelpBox <- renderUI({
  if (!(input$exampleHelp %% 2)) {
    HTML(
      "
        <p>Teakel, S. L., Ludescher, M., Thejer, B. M., Poschmann, G., Forwood, J. K., Neubauer, H., & Cahill, M. A.
        (2020). Protein complexes including PGRMC1 and actin-associated proteins are disrupted by AG-205.
        Biochemical and biophysical research communications, 524(1), 64-69.
        </p>
        <p>PGRMC1 is a protein from the MAPR family with a range of cellular functions.
                       PGRMC1 has been described to play a role in regulating membrane trafficking and
                       in cancer progression and response to therapies. To further understand the functions
                       of PGRMC1 and the mechanism of the small molecule inhibitor of PGRMC1, AG-205,
                       proteins differentially bound to PGRMC1 were identified following AG-205
                       treatment of MIA PaCa-2 cells.</p>
                       Data have been re-analyzed from PRIDE identifier PXD016455.
                       <a href='https://pubmed.ncbi.nlm.nih.gov/31980178'target='_blank'>Further information</a>.</p>
                       "
    )
  }
})

output$NetworkHelpBox <- renderUI({
  if (input$networkHelp %% 2) {
    HTML(
      "<p>PPI (protein-protein interaction) Network from IntAct (2022-07-13).
      All interactions are derived from literature curation or direct user
      submissions. Edge weights can be further filtered, more information can
      be found <a href='https://www.ebi.ac.uk/intact/' target='_blank'>here</a>.
      When multiple group comparisons are selected two types of edges are
      created: edges from the group comparison to the proteins and PPI edges from
      IntAct between the proteins. Networks can be downloaded in GML format
      enabling the visualization in a network tool.</p>
      <p>For data from a single group comparison fold changes are color coded
      onto the nodes, for single and multiple group comparisons sub-cellular
      locations are retrieved from <a href='https://cell-map.org/' target='_blank'>
      Human Cell Map</a> (Go, Christopher D., et al. A proximity-dependent
      biotinylation map of a human cell. Nature (2021): 1-5.). Human Cell Map
      provides two subcellular localization predictions, one resulting from
      Spatial analysis of functional enrichment (SAFE) and a different one
      resulting from non-negative matrix factorization (NMF). A summarized subcellular
      localization from NMF predictions can be mapped onto the nodes, predictions
      from SAFE and NMF can be seen in the node table below and are also
      included in the gml download.
      </p>"
    )
  } else {
    return()
  }
})

output$FoldChangePlotHelpBox <- renderUI({
  if (input$FoldChangePlotHelp %% 2) {
    helpText(
      "The fold change plot gets plotted for the proteins in your selection.
      Points are colored on their specificity to the comparison. In some cases fold changes
      seem to be high in both comparisons, but are only significantly differentially abundant in one comparison.
               "
    )
  }
})

### QC HELP
output$boxplotHelpBox <- renderUI({
  if (input$boxplotHelp %% 2) {
    helpText(
      "Box plots show the distribution of selected intensities which
      gives an overview about their similarities.
               "
    )
  }
})

output$densityHelpBox <- renderUI({
  if (input$densityHelp %% 2) {
    helpText(
      "The density plot shows a smoothed version of a histogram.
      It is especially useful to compare density plots of the intensities before and after imputation.
               "
    )
  }
})

output$corHelpBox <- renderUI({
  if (input$corHelp %% 2) {
    helpText(
      "
      The Pearson correlation plot visualizes how well samples (e.g replicates)
      correlate.
               "
    )
  }
})

output$cvHelpBox <- renderUI({
  if (input$cvHelp %% 2) {
    helpText(
      "The Coefficient of Variation (CV) gets calculated by the standard deviation of replicates divided by their mean per protein,
                which gives an estimate on the reproducibility of the experiment.
               "
    )
  }
})

output$contaminantHelpBox <- renderUI({
  if (input$contaminantHelp %% 2) {
    helpText(
      "The percentage of contaminants is calculated with the iBAQ (Intensity based Absolute Quantification) values and shows
                the percentage of total signal from potential contaminants.
               "
    )
  }
})

output$abundantHelpBox <- renderUI({
  if (input$abundantHelp %% 2) {
    helpText(
      "The percentage of most abundant proteins from a sample is calculated with the iBAQ (Intensity based Absolute Quantification) values and shows
                the percentage of the 15 most abundant signals from the total signal.
               "
    )
  }
})

output$idHelpBox <- renderUI({
  if (input$idHelp %% 2) {
    helpText(
      "The number of identified proteins is calculated from the LFQ intensities before imputation.
               "
    )
  }
})

output$mvHelpBox <- renderUI({
  if (input$mvHelp %% 2) {
    helpText(
      "The percentage of missing values is calculated from the LFQ intensities before imputation.
               "
    )
  }
})

output$overlapHeatmapHelpBox <- renderUI({
  if (input$overlapHeatmapHelp %% 2) {
    HTML(
      "
      <strong>Jaccard index</strong>: The Jaccard index for two sets is defined 
      as the the size of the intersection divided by the size of the union: 
      J(X,Y) = |X∩Y| / |X∪Y|.<br>
      <strong>Overlap coefficient</strong>: Related to the Jaccard index, 
      defined as the size of the intersection divided by he smaller of the size 
      of the two sets: J(X,Y) = |X∩Y| / min(|X|,|Y|)<br>
      <strong>Number of shared proteins</strong>: The size of the intersection.<br>
      "
    )
  }
})

output$oraHelpBox <- renderUI({
  
  text <- "
        Information from <a href='https://biit.cs.ut.ee/gprofiler/gost'target='_blank'>
        gprofiler's web site</a>:
        <p>
        gprofiler performs functional enrichment analysis, also known as
        over-representation analysis (ORA) or gene set enrichment analysis, on
        input gene list. It maps genes to known functional information sources
        and detects statistically significantly enriched terms. gprofiler
        regularly retrieves data from Ensembl database and fungi, plants or
        metazoa specific versions of Ensembl Genomes, and parasite specific data
        from WormBase ParaSite. In addition to Gene Ontology, gprofiler includes
        pathways from KEGG Reactome and WikiPathways; miRNA targets from
        miRTarBase and regulatory motif matches from TRANSFAC; tissue specificity
        from Human Protein Atlas; protein complexes from CORUM and human disease
        phenotypes from Human Phenotype Ontology. g:GOSt supports close to 500
        organisms and accepts hundreds of identifier types.
        </p><h4>Source information:</h4>"
  
  out <- gprofiler2::get_version_info()
  srcNames <- names(out$sources)
  
  outText <- "<p>"
  for (src in srcNames) {
    src_info <- out$sources[[src]]
    outText <- paste0(outText, src, ':&emsp; ', src_info[[2]], '<br>')
  }
  
  outText <- paste0(outText, "</p>")
  
  if (input$oraHelp %% 2) {
    HTML(
      paste0(
        text, outText
      )
    )
  }
})




observeEvent(input$showReportingStandards, {
  out <- '
  <h3>What to report</h3>
  <p>Please report all analysis parameters in the "Input tab" used for data 
analysis, that are available for download after data upload and which are 
also depicted in the QC report:</p>

<ul>
	<li>Filtering proteins on MS/MS counts and min. razor peptides</li>
	<li>Filtering proteins on min. valid values per group</li>
	<li>Normalization settings</li>
	<li>Imputation settings</li>
</ul>



<p>Please also make sure to save the selected global parameters used for 
quantification in the "Differential abundance tab", which are also 
available in the Differential Abundance report:</p>

<ul>
	<li>Fold change threshold</li>
	<li>Significance cutoff (adj. p-value, p-value or none)</li>
	<li>(adj.) p-value threshold</li>
	<li>How to apply fold change threshold (absolute, show only 
	enriched, or show only reduced)</li>
</ul>

  <h3>Used libraries and ressources</h3>
         <p>Please cite following publications (if you have used them):</p>
         <hr>
         <h4>Upload Analysis options</h4>
         <ul>
         <li>(Differential expression analysis) <b>limma:</b> Ritchie, Matthew E., et al. "limma powers
         differential expression analyses for RNA-sequencing and microarray
         studies." Nucleic acids research 43.7 (2015): e47-e47.</li>
         <li>(Differential expression analysis) <b>DEqMS: </b> Zhu, Yafeng, et al. "DEqMS: a method for accurate
         variance estimation in differential protein expression analysis."
         Molecular & Cellular Proteomics 19.6 (2020): 1047-1057.</li>
         </ul>
         <hr>
         <h4>Data analysis and visualizations in multiple tabs</h4>
         <ul>
         <li>(Heatmaply) <b>heatmaply: </b>Galili, Tal, et al. "heatmaply: an R
         package for creating interactive cluster heatmaps for online
         publishing." Bioinformatics 34.9 (2018): 1600-1602.</li>
         <ul>
         Following visualizations are created by the <b>heatmaply</b> library:
         <li>QC tab: interactive correlation plot</li>
         <li>QC tab: interactive protein overlap plot</li>
         <li>Differential Abundance tab: interactive heatmap</li>
         </ul>
         </ul>
         <hr>
         <h4>Data analysis and visualizations in Differential Abundance tab</h4>
         <ul>
         <li>(UpSet plot) <b>UpSet: </b> Lex, Alexander, et al. "UpSet:
         visualization of intersecting sets." IEEE transactions on visualization
         and computer graphics 20.12 (2014): 1983-1992.</li>
         The UpSet plot is a powerful visualization for multiple set comparisons.
         <li>(PPI Networks) <b>IntAct: </b>Orchard, Sandra, et al. "The MIntAct
         project—IntAct as a common curation platform for 11 molecular
         interaction databases." Nucleic acids research 42.D1 (2014): D358-D363.</li>
         <li>(Subcell. localization) <b>Human CellMap: </b>Go, Christopher D., et al.
         "A proximity-dependent biotinylation map of a human cell." Nature (2021): 1-5.</li>
         <li>(ORA) <strong>gprofiler2: </strong>Raudvere, Uku, et al. "g: Profiler: a web
         server for functional enrichment analysis and conversions of gene
         lists (2019 update)." Nucleic acids research 47.W1 (2019): W191-W198.</li>
         </ul>
         </p>
         <br>
  '
  showModal(
    modalDialog(
      title = "Reporting standards",
      HTML(out),
      size = "l",
      easyClose = TRUE,
      footer = NULL
    )
  )
})

