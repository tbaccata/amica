ARG AMICA_SOURCE_HYPERLINK="https://www.github.com/tbaccata/amica"

ARG AMICA_SHINY_IDLE_TIMEOUT_SECONDS=3600
ARG AMICA_SHINY_PORT_DOCKER_INTERNAL=3838
ARG AMICA_SHINY_HOST_DOCKER_INTERNAL="0.0.0.0"

ARG AMICA_VERSION_OVERRIDE="3.0.1"

############################################################################
FROM rocker/r-base:4.4.1 AS build_deps

LABEL maintainer="Sebastian Didusch <sebastian.didusch@univie.ac.at>"
LABEL maintainer="Juraj Ahel <juraj.ahel@vbcf.ac.at>"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libxslt-dev

WORKDIR /src
COPY install_dependencies.R .
RUN Rscript install_dependencies.R

# # basic shiny functionality
# RUN R -e "install.packages(c('shiny', 'rmarkdown', 'shinyjs', 'shinymanager', 'shinyBS', 'DT', 'shinycssloaders', 'bslib', 'profvis', 'colourpicker', 'shinyalert'), repos='https://cloud.r-project.org/')"
# RUN R -e "install.packages(c('reshape2', 'igraph', 'visNetwork', 'UpSetR', 'dplyr', 'pheatmap', 'DT', 'data.table', 'RColorBrewer', 'Rmisc', 'eulerr'), repos='https://cloud.r-project.org/')"
# RUN R -e "install.packages(c('heatmaply', 'ggfortify', 'colourvalues', 'pryr'),  repos='https://cloud.r-project.org/')"
# RUN R -e "install.packages(c('tidyr', 'Cairo', 'cowplot'),  repos='https://cloud.r-project.org/')"

# # install dependencies from bioconductor
# RUN R -e "install.packages('BiocManager'); library('BiocManager'); BiocManager::install(c('limma', 'DEqMS', 'gprofiler2', 'vsn'))"

############################################################################
FROM build_deps AS final

WORKDIR /

# copy the app to the image
RUN mkdir /root/amica
COPY amica /root/amica

COPY Rprofile.site /usr/lib/R/etc/

ARG AMICA_SHINY_PORT_DOCKER_INTERNAL
ENV AMICA_SHINY_PORT=${AMICA_SHINY_PORT_DOCKER_INTERNAL}
EXPOSE ${AMICA_SHINY_PORT_DOCKER_INTERNAL}

ARG AMICA_SHINY_HOST_DOCKER_INTERNAL
ENV AMICA_SHINY_HOST=${AMICA_SHINY_HOST_DOCKER_INTERNAL}

ARG AMICA_SHINY_IDLE_TIMEOUT_SECONDS
ENV AMICA_SHINY_IDLE_TIMEOUT_SECONDS=${AMICA_SHINY_IDLE_TIMEOUT_SECONDS}

ARG AMICA_VERSION_OVERRIDE
ENV AMICA_VERSION_OVERRIDE=${AMICA_VERSION_OVERRIDE}

ARG AMICA_SOURCE_HYPERLINK
ENV AMICA_SOURCE_HYPERLINK=${AMICA_SOURCE_HYPERLINK}

CMD ["R", "-e", "shiny::runApp('/root/amica')"]

