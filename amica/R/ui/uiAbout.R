tabPanel(
  title = 'About',
  value = 'abouttab',
  HTML('<p>
    <h3>News</h3>
         <hr>
         <ul>
         <li>2023-05-30 - New version (3.0.0).</li>
         <li>2023-05-30 - Additional file upload option for DIA (Spectronaut and DIA-NN).</li>
         <li>2023-05-30 - Additional file upload option for TMT (FragPipe).</li>
         <li>2023-05-30 - New feature added - option to highlight proteins in volcano - and MA plots.</li>
         <li>2023-05-30 - Updated IntAct version (2022-07-13).</li>
         <li>2023-05-30 - Changed default plot colors.</li>
	 <li>22.12.09 - amica published in <a href="https://doi.org/10.1186/s12864-022-09058-7"
	 target="_blank">BMC Genomics</a>.</li>
         <li>2022-08-30 - New feature added - option to not impute missing values.</li>
         <li>2022-08-30 - New feature added - QC - and diff. abundance reports.</li>
         <li>2022-08-30 - Fixed error handling in upload.</li>
         <li>2022-02-28 - New feature added - Dot plot in Diff. abundance tab.</li>
         <li>2021-11-23 - amica posted on <a href="https://www.biorxiv.org/content/10.1101/2021.11.23.466958" 
         target="_blank">bioRxiv</a>.</li>
         <li>2021-10-14 - github went public.</li>
         </ul>
         </p>
         <p>
         <h3>Used libraries and ressources</h3>
         <p>Please cite following publications (if you have used them):</p>
         <hr>
         <ul>
         <li>(Differential expression analysis) <b>limma:</b> Ritchie, Matthew E., et al. "limma powers 
         differential expression analyses for RNA-sequencing and microarray 
         studies." Nucleic acids research 43.7 (2015): e47-e47.</li>
         <li>(Differential expression analysis) <b>DEqMS: </b> Zhu, Yafeng, et al. "DEqMS: a method for accurate 
         variance estimation in differential protein expression analysis." 
         Molecular & Cellular Proteomics 19.6 (2020): 1047-1057.</li>
         <li>(ORA) <strong>gprofiler2: </strong>Raudvere, Uku, et al. "g: Profiler: a web 
         server for functional enrichment analysis and conversions of gene 
         lists (2019 update)." Nucleic acids research 47.W1 (2019): W191-W198.</li>
         <li>(PPI Networks) <b>IntAct: </b>Orchard, Sandra, et al. "The MIntAct 
         project—IntAct as a common curation platform for 11 molecular 
         interaction databases." Nucleic acids research 42.D1 (2014): D358-D363.</li>
         <li>(Subcell. localization) <b>Human CellMap: </b>Go, Christopher D., et al. 
         "A proximity-dependent biotinylation map of a human cell." Nature (2021): 1-5.</li>
         <li>(Heatmaply) <b>heatmaply: </b>Galili, Tal, et al. "heatmaply: an R 
         package for creating interactive cluster heatmaps for online 
         publishing." Bioinformatics 34.9 (2018): 1600-1602.</li>
         <li>(UpSet plot) <b>UpSet: </b> Lex, Alexander, et al. "UpSet: 
         visualization of intersecting sets." IEEE transactions on visualization 
         and computer graphics 20.12 (2014): 1983-1992.</li>
         </ul>
         </p>
         <br>
             '),
  footer()
)
