[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub](https://img.shields.io/github/license/tbaccata/amica?color=brightgreen)
![R](https://img.shields.io/badge/R-4-blue)

# amica
amica: an interactive and user-friendly web-platform for the analysis of proteomics data

![amica_logo](.www/ga_amica.png)


## Raison d'etre

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


### Local installation

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

## Build LFQ-Analyst (Any name after -t)
> docker build -t amica .

## Start amica from terminal

> docker run -p 3838:3838 amica

## Open local interface

https://localhost:3838/amica


```







