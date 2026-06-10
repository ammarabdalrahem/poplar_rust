# syntax=docker/dockerfile:1
#
# Reproducible environment for the population-genetics workflow of
# Melampsora larici-populina (microsatellite markers).
#
#   Build:  docker build -t poplar_rust .
#   Run:    docker run --rm -v "$(pwd)/output:/project/output" poplar_rust
#
# A pre-built image is published automatically to the GitHub Container Registry:
#   docker pull ghcr.io/ammarabdalrahem/poplar_rust:latest
#
# ---------------------------------------------------------------------------
# Base image
# ---------------------------------------------------------------------------
# rocker/geospatial:4.4.1 pins R 4.4.1 together with the full geospatial stack
# (GDAL / GEOS / PROJ, sf, terra, stars) and the tidyverse, all precompiled.
# This is the environment the PCI data editor used to verify the analysis, so
# the container reproduces that review setup exactly.
FROM rocker/geospatial:4.4.1

LABEL org.opencontainers.image.source="https://github.com/ammarabdalrahem/poplar_rust"
LABEL org.opencontainers.image.description="Reproducible population-genetics workflow for Melampsora larici-populina (microsatellite markers)."
LABEL org.opencontainers.image.licenses="CC-BY-4.0"

ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /project

# ---------------------------------------------------------------------------
# Reproducibility: freeze CRAN to a dated Posit Package Manager snapshot
# ---------------------------------------------------------------------------
# Every build then resolves to the SAME package versions, and pulls precompiled
# binaries for the image's Ubuntu release (fast, no source compilation).
# Bump PKG_SNAPSHOT only when you deliberately want newer packages.
ARG PKG_SNAPSHOT=2025-04-15
RUN . /etc/os-release && \
    echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/${PKG_SNAPSHOT}'))" \
      >> /usr/local/lib/R/etc/Rprofile.site

# ---------------------------------------------------------------------------
# System libraries not already present in the geospatial base image
# ---------------------------------------------------------------------------
#   libglpk-dev : graph optimisation used by igraph / poppr
#   libgsl-dev  : GNU Scientific Library used by some genetics dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
      libglpk-dev libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# CRAN packages  (sf, terra and the tidyverse are already in the base image)
# ---------------------------------------------------------------------------
RUN install2.r --error --skipinstalled -n -1 \
      knitr readxl cowplot viridis ggrepel ggsci scales \
      factoextra svglite reshape2 ggpubr RColorBrewer ggforce \
      rnaturalearth rnaturalearthdata mapdata mapplots \
      lme4 genepop hierfstat adegenet poppr pegas ape ade4 \
      BiocManager remotes

# ---------------------------------------------------------------------------
# Bioconductor  (pinned to release 3.20, the release paired with R 4.4)
# ---------------------------------------------------------------------------
RUN Rscript -e "BiocManager::install(version = '3.20', ask = FALSE, update = FALSE)" && \
    Rscript -e "BiocManager::install(c('ggtree', 'ggtreeExtra'), ask = FALSE, update = FALSE)"

# ---------------------------------------------------------------------------
# GitHub packages
# ---------------------------------------------------------------------------
# rnaturalearthhires is pinned to an exact commit for full reproducibility.
RUN Rscript -e "remotes::install_github('ropensci/rnaturalearthhires@e4736f636baa1c013d77d2ba028dd5bc334defee', upgrade = 'never')"

# RClone is a stable, rarely-updated package. For strict reproducibility, append
# a commit SHA to the ref below, e.g. 'dbailleul/RClone@<commit-sha>'.
RUN Rscript -e "remotes::install_github('dbailleul/RClone', dependencies = FALSE, upgrade = 'never')"

# ---------------------------------------------------------------------------
# Project files
# ---------------------------------------------------------------------------
COPY Table_S1.csv             /project/Table_S1.csv
COPY data_analysis_mlp_new.R  /project/data_analysis_mlp_new.R

# Results are written under /project/output. Mount a host directory there to
# collect figures and tables on the host, e.g.:
#   docker run --rm -v "$(pwd)/output:/project/output" ghcr.io/ammarabdalrahem/poplar_rust:latest
VOLUME ["/project/output"]

CMD ["Rscript", "data_analysis_mlp_new.R"]
