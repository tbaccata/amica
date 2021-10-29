install.packages(
  c(
    'shiny',
    'rmarkdown',
    'shinyjs',
    'shinymanager',
    'shinyBS',
    'DT',
    'shinycssloaders',
    'bslib',
    'colourpicker',
    'reshape2',
    'igraph',
    'visNetwork',
    'UpSetR',
    'dplyr',
    'pheatmap',
    'DT',
    'data.table',
    'RColorBrewer',
    'Rmisc',
    'eulerr',
    'heatmaply',
    'ggfortify',
    'colourvalues'
  ),
  repos = 'https://cloud.r-project.org/'
)

install.packages('BiocManager')
library('BiocManager')
BiocManager::install(c('limma', 'DEqMS', 'gprofiler2', 'vsn'))
