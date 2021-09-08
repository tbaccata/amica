FROM rocker/r-base

LABEL maintainer "Sebastian Didusch <sebastian.didusch@univie.ac.at>"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libxslt-dev

# basic shiny functionality
RUN R -e "install.packages(c('shiny', 'rmarkdown', 'shinyjs', 'shinymanager', 'shinyBS', 'DT', 'shinycssloaders', 'bslib', 'profvis', 'colourpicker'), repos='https://cloud.r-project.org/')"

RUN R -e "install.packages(c('reshape2', 'igraph', 'visNetwork', 'UpSetR', 'dplyr', 'pheatmap', 'DT', 'data.table', 'RColorBrewer', 'Rmisc'), repos='https://cloud.r-project.org/')"




# install dependencies of the euler app
RUN R -e "install.packages('BiocManager'); library('BiocManager'); BiocManager::install(c('limma', 'DEqMS', 'SummarizedExperiment', 'gprofiler2', 'vsn'))"

RUN R -e "install.packages(c('heatmaply', 'ggfortify', 'colourvalues', 'pryr'),  repos='https://cloud.r-project.org/')"


# copy the app to the image
RUN mkdir /root/amica
COPY amica /root/amica

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/amica')"]

