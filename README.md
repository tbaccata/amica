[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![R](https://img.shields.io/badge/R-4-blue)

# amica
amica: an interactive and user-friendly web-platform for the analysis of proteomics data



![amica_logo](/amica/www/ga_amica.png)




## Functionality

- Faciliting interactive analyses and visualizations with just a couple of clicks


### Input

- MaxQuant's **proteinGroups.txt**, FragPipe's **combined_protein.tsv** or any tab-separated file
- Processed data can be downloaded in a developed **amica format** which can also be used as input
- Experimental design mapping samples to conditions
- Contrast matrix file for group comparisons in case of MaxQuant, FragPipe or custom upload
- specification file for mapping relevant columns in case of custom file upload


### Outputs

- Analyzed data downloadable as **amica format**
- Almost all plots prduced by plotly (hover over plot and download plot as svg or png with the camera icon)
- All plots have customizable plot parameters (width, height, file format) 
- Data tables

### Analysis options

- Remove decoys and proteins only identified by site (MaxQuant)
- Filter on minimum peptide count and spectral count values
- Filter on minimum valid values per group
- Select intensities to
- (Re-)normalize intensities (VSN, Quantile, Median)
- Imputate missing values from normal distribution or replace them by constant value (useful for pilots)


### QC-plots

#### For different intensities (Raw intensities, LFQ intensities, imputed intensities)
- PCA
- Box plots
- Density plots
- Correlation plots (Pearson correlation)
- Bar plots (identified proteins, % contaminants, most abundant proteins) per sample
- Scatter plots


### Differential abundance analysis

- Primary filter options (log2FC thresholds, multiple-testing correction, select enriched or reduced proteins)
- Analyze single or multiple group comparison
- Volcano - and MA - plots
- Set comparisons (UpSet plot)
- Customizable output data table (can be further filtered)
- Heatmap
- Fold change plot
- Profile plot
- Protein-protein interaction (PPI) network
- Over-Representation Analysis (ORA)


### Compare multiple amica files

- Upload a second amica file from another experiment/analysis to combine data sets
- Download combined data set
- Correlate intensities from combined data set
- Differential abundance analysis for combined amica data set


## Dependencies

### Session info

```

> sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 21.04

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=de_AT.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] tools     stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] colourpicker_1.1.0 RColorBrewer_1.1-2 dplyr_1.0.7        data.table_1.14.0  Rmisc_1.5          plyr_1.8.6         lattice_0.20-44    pheatmap_1.0.12   
 [9] colourvalues_0.3.7 UpSetR_1.4.0       visNetwork_2.0.9   igraph_1.2.6       reshape2_1.4.4     bslib_0.2.5.1      gprofiler2_0.2.0   DEqMS_1.10.0      
[17] limma_3.48.1       DT_0.18            heatmaply_1.2.1    viridis_0.6.1      viridisLite_0.4.0  plotly_4.9.4.1     ggfortify_0.4.12   ggplot2_3.3.5     
[25] shinyBS_0.61       shinyjs_2.0.0      shiny_1.6.0       

loaded via a namespace (and not attached):
 [1] httr_1.4.2        sass_0.4.0        tidyr_1.1.3       jsonlite_1.7.2    foreach_1.5.1     assertthat_0.2.1  yaml_2.2.1        pillar_1.6.1     
 [9] glue_1.4.2        digest_0.6.27     promises_1.2.0.1  colorspace_2.0-2  htmltools_0.5.1.1 httpuv_1.6.1      pkgconfig_2.0.3   purrr_0.3.4      
[17] xtable_1.8-4      scales_1.1.1      webshot_0.5.2     later_1.2.0       tibble_3.1.2      farver_2.1.0      generics_0.1.0    ellipsis_0.3.2   
[25] cachem_1.0.5      withr_2.4.2       lazyeval_0.2.2    magrittr_2.0.1    crayon_1.4.1      mime_0.11         fs_1.5.0          fansi_0.5.0      
[33] registry_0.5-1    lifecycle_1.0.0   stringr_1.4.0     munsell_0.5.0     compiler_4.1.1    jquerylib_0.1.4   rlang_0.4.11      grid_4.1.1       
[41] iterators_1.0.13  htmlwidgets_1.5.3 crosstalk_1.1.1   miniUI_0.1.1.1    labeling_0.4.2    gtable_0.3.0      codetools_0.2-18  DBI_1.1.1        
[49] TSP_1.1-10        R6_2.5.0          seriation_1.3.0   gridExtra_2.3     fastmap_1.1.0     utf8_1.2.1        dendextend_1.15.1 stringi_1.7.3    
[57] Rcpp_1.0.7        vctrs_0.3.8       tidyselect_1.1.1 

```

## Local installation

- Using git and Rstudio
```
## Clone the repository
git clone https://github.com/tbaccta/amica.git

## Move to the folder
cd amica

## Inside R console or R studio
> library("shiny")

> runApp()

```

- Using Docker

Have docker installed and running (www.docker.com/get-started)

```
## Clone the repository
git clone https://github.com/tbaccata/amica.git

## Move to the folder
cd amica

## Build amica, the -t flag is the name of the docker image
docker build -t amica .

## Start amica from terminal

docker run -p 3838:3838 amica

## Open local interface

https://localhost:3838/amica


```

## Deploy amica with ShinyProxy 

When deploying a Shiny application with ShinyProxy, the application is simply bundled as an R package and installed into a Docker image. Every time a user runs an application, a container spins up and serves the application.

Detailed documentation is provided here (https://www.shinyproxy.io/documentation/).

A minimum working example based on documentation (https://www.shinyproxy.io/documentation/deployment):

```

## install docker image for amica
git clone https://github.com/tbaccata/amica.git
cd amica
docker build -t amica .
 
## download latest version and install it
wget https://www.shinyproxy.io/downloads/shinyproxy_2.5.0_amd64.deb
sudo dpkg -i shinyproxy_2.5.0_amd64.deb

## enable system process
sudo systemctl enable shinyproxy

## Add amica into specs part of the server /etc/shinyproxy/application.yml:
## In this file you can also specify the port for shinyproxy.

specs:
  - id: amica
    display-name: amica Shiny App
    description: Analysis and visualization tool for quantitative MS
    container-cmd: ["R", "-e", "shiny::runApp('/root/amica')"]
    container-image: amica
    access-groups: [scientists, mathematicians]

```
