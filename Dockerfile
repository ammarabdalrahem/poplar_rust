FROM rocker/r-ver:4.4.1
LABEL org.opencontainers.image.source=https://github.com/YOUR_USERNAME/YOUR_REPO_NAME

ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /project

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran pkg-config git curl \
    cmake \
    libglpk-dev libgsl-dev libxml2-dev \
    libcurl4-openssl-dev libssl-dev libgit2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
    libblas-dev liblapack-dev libudunits2-dev \
    libgdal-dev libgeos-dev libproj-dev \
    && rm -rf /var/lib/apt/lists/*


RUN install2.r --error --skipinstalled -n 4 \
    tidyverse lme4 knitr readxl ggplot2 dplyr \
    ggrepel cowplot ade4 viridis scales sf \
    reshape2 svglite ggpubr adegenet ape poppr \
    pegas factoextra RColorBrewer genepop hierfstat \
    mapdata mapplots ggsci ggforce rnaturalearth \
    rnaturalearthdata BiocManager remotes devtools


RUN Rscript -e "BiocManager::install(c('ggtree', 'ggtreeExtra'), ask = FALSE, update = FALSE)"

RUN Rscript -e "remotes::install_github('dbailleul/RClone', dependencies = TRUE, upgrade = 'never')"

RUN Rscript -e "remotes::install_github('ropensci/rnaturalearthhires', dependencies = TRUE, upgrade = 'never')"

COPY Table_data.tsv /project/Table_data.tsv
COPY data_analysis_mlp_new.R /project/data_analysis_mlp_new.R

CMD ["Rscript", "/project/data_analysis_mlp_new.R"]