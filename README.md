# Population Genetics of *Melampsora larici-populina*

Microsatellite-based workflow for population genetic analysis of *Melampsora larici-populina* isolates.

This repository contains an R workflow to:

- Load isolate metadata and microsatellite genotypes
- Define multilocus genotypes (MLGs) and multilocus lineages (MLLs)
- Assign isolates to genetic clusters
- Infer sexual versus asexual reproduction
- Calculate population genetic indices
- Build a bootstrapped phylogenetic tree
- Generate spatial maps
- Test the effects of geography and year using regression models

All outputs are written to two subdirectories created automatically on first run:
`output/figures/` for all plots and `output/tables/` for all data tables.

---

## Associated study

**Long-lasting coexistence of multiple asexual lineages alongside their sexual counterparts in a fungal plant pathogen**

**Authors:**
Ammar Abdalrahem, Axelle Andrieux, Ronan Becheler, Sébastien Duplessis, Pascal Frey, Benoit Marcais, Kadiatou Schiffer-Forsyth, Solenn Stoeckel, Fabien Halkett

---

## Fully reproducible with Docker

The entire analysis environment is packaged in a Docker image. No R installation, no package management, no dependency conflicts.

### Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (free, works on macOS / Windows / Linux)

### First option

The analysis using the `rocker/geospatial:4.4.1` image:

```bash
# 1. Clone the repository
git clone https://github.com/ammarabdalrahem/poplar_rust.git
cd poplar_rust

# 2. Run the analysis
docker run --rm \
  -v $(pwd):/project \
  -w /project \
  rocker/geospatial:4.4.1 \
  Rscript data_analysis_mlp_new.R
```

> **On Windows PowerShell**, replace `$(pwd)` with `${PWD}`.

### Second option: custom pre-built image (all packages pre-installed)

A ready-to-run image with every dependency baked in is published to the GitHub
Container Registry. No build step is required:

```bash
# Pull the latest published image
docker pull ghcr.io/ammarabdalrahem/poplar_rust:latest

# Run it, collecting results into ./output on your machine
docker run --rm \
  -v "$(pwd)/output:/project/output" \
  ghcr.io/ammarabdalrahem/poplar_rust:latest
```

> **On Windows PowerShell**, replace `$(pwd)` with `${PWD}`.

After the run, outputs appear in `output/figures/` and `output/tables/` inside your local project folder.

### How the image is built and pinned

The image is defined by the `Dockerfile` and rebuilt automatically by GitHub
Actions (`.github/workflows/docker-publish.yml`) on every push to `main` and on
each version tag, then pushed to `ghcr.io/ammarabdalrahem/poplar_rust`.
Pushing a tag such as `v1.0` publishes both `:1.0` and `:latest`.

For long-term reproducibility, package versions are frozen:

- **Base image:** `rocker/geospatial:4.4.1` (R 4.4.1 + the full geospatial stack)
- **CRAN:** pinned to a dated [Posit Package Manager](https://packagemanager.posit.co/) snapshot (`PKG_SNAPSHOT` build arg) so the same versions resolve on every build
- **Bioconductor:** pinned to release `3.20`
- **GitHub packages:** `rnaturalearthhires` pinned to an exact commit

---

## Repository contents

| File | Description |
|------|-------------|
| `data_analysis_mlp_new.R` | Main R script — use for terminal / Docker execution |
| `data_analysis_mlp.Rmd` | R Markdown version — use for interactive work in RStudio |
| `Table_S1.csv` | Input data: isolate metadata and microsatellite genotypes (Table S1) |
| `Dockerfile` | Docker image definition for full reproducibility |
| `.dockerignore` | Files excluded from the Docker build context |
| `.github/workflows/docker-publish.yml` | CI: builds the image and publishes it to GHCR |
| `README.md` | This file |

---

## Manual setup (without Docker)

### Requirements

- R 4.4.1 or later
- RStudio recommended for interactive work

### R packages

Packages are installed automatically in dependency order when the script is first run. To install manually:

```r
# Layer 1 — CRAN core
install.packages(c(
  "knitr", "ggplot2", "readxl", "tidyverse", "cowplot",
  "viridis", "ggrepel", "ggsci", "scales", "dplyr",
  "factoextra", "grid", "svglite", "reshape2"
))

# Layer 2 — spatial stack
install.packages(c(
  "sf", "rnaturalearth", "rnaturalearthdata", "mapdata", "mapplots"
))
install.packages(
  "rnaturalearthhires",
  repos = "https://ropensci.r-universe.dev",
  type  = "source"
)

# Layer 3 — population genetics
install.packages(c(
  "adegenet", "poppr", "hierfstat", "pegas", "ape",
  "lme4", "genepop", "ade4"
))

# Layer 4 — Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggtreeExtra"), ask = FALSE, update = FALSE)

# Layer 5 — GitHub
remotes::install_github("dbailleul/RClone", dependencies = TRUE)
```

### Run the script

```bash
Rscript data_analysis_mlp_new.R
```

On first run the script will:
1. Create `output/figures/` and `output/tables/` directories
2. Write all outputs into those subdirectories

---

## Input format

```text
Table_S1.csv
```

Requirements:

- Microsatellite locus columns must contain `Mlp` in the column name
- Missing genotype values must be coded as `999`

---

## Outputs

### Tables — written to `output/tables/`

| File | Description | Manuscript |
|------|-------------|------------|
| `Table1_a_genetic_indices_cluster_approach.csv` | Genetic indices per genetic cluster (clustering approach) | Table 1 |
| `Table1_b_genetic_indices_resampling_approach.csv` | Genetic indices per reproduction mode (resampling approach) | Table 1 |
| `Table1_c_genetic_indices_combination_approaches.csv` | Genetic indices per reproduction mode (combined approach) | Table 1 |
| `Table2_contingency_clustering_vs_resampling.csv` | Contingency table: cluster vs resampling assignments | Table 2 |
| `Table3_GLMM_binomial_regression.csv` | GLMM fixed-effect coefficients (binomial, Lat + Long + Year) | Table 3 |
| `Table4_top7_asexual_MLLs.csv` | Characteristics of the seven most abundant asexual MLLs | Table 4 |
| `Table_asexual_MLL_genetic_indices.csv` | Full genetic indices for all asexual MLLs | — |
| `genetic_indices_per_sexual_MLL.csv` | Genetic indices per sexual MLL | — |
| `new_genotype_data.csv` | Final per-isolate table with cluster and reproduction labels | — |
| `filtered_mll_years.csv` | MLLs recurring across multiple sampling years | — |
| `MST_data_mlp_pop_as_Reproduction_for_cloneEstimate.csv` | Minimum-spanning-tree export for ClonEstiMate | — |

> **Note on Table 1 — Pareto β column:** this statistic is computed by the external
> software [GenAPoPop](https://forge.inrae.fr/solenn.stoeckel/genapopop1.0) and cannot be produced by
> this R script. Run GenAPoPop separately on the same isolate data and insert the
> resulting Pareto β values into the exported CSV before final publication.

### Figures — written to `output/figures/`

| File | Description | Manuscript |
|------|-------------|------------|
| `FigS1_geographic_distribution.svg` | Sampling map (all isolates) | Fig. S1 |
| `FigS2_Silhouette_kmeans.png` | Silhouette plot for k-means cluster evaluation | Fig. S2 |
| `FigS3_cluster_assignments.png` | Cluster assignment probability scatter plot | Fig. S3 |
| `FigS4_effect_of_latitude.png` | GLM: proportion of sexual reproduction vs latitude | Fig. S4 |
| `FigS4_effect_of_latitude_glmm.png` | GLMM: proportion of sexual reproduction vs latitude | Fig. S4 |
| `Fig2_geographical_distribution_2009_2011.png` | Map of sexual vs asexual proportions, 2009 and 2011 | Fig. 2 |
| `Fig3a_tree_plot1.png` | Circular NJ phylogenetic tree (branches only) | Fig. 3 |
| `Fig3_tree_nj_plot.png` | Annotated NJ tree with MLL rings | Fig. 3 |
| `Fig4A_asex_mll_Locations.png` | Asexual lineage abundance across sampling locations | Fig. 4A |
| `Fig4B_asex_mll_Year.png` | Asexual lineage abundance across years | Fig. 4B |
| `Fig4_asex_mll_Year_Locations.png` | Combined lineage abundance figure | Fig. 4 |
| `Dapc_xval.png` | DAPC cross-validation plot | — |
| `DAPC_scatter.png` | DAPC scatter plot | — |
| `DAPC_compoplot.png` | DAPC composition plot | — |
| `map_all_years.png` | Map of sexual vs asexual proportions across all years | — |

---

## Notes

- `data_analysis_mlp_new.R` and `data_analysis_mlp.Rmd` produce **identical outputs** with identical file names
- The Docker image `rocker/geospatial:4.4.1` and the custom image `ghcr.io/ammarabdalrahem/poplar_rust:1.0` both contain R 4.4.1 with all required packages
- When run outside Docker, the script installs any missing packages automatically in dependency-ordered layers (CRAN core → spatial → genetics → Bioconductor → GitHub); a first run on a clean R installation may take several minutes
- Isolates with uncertain cluster assignment are excluded from downstream analyses
- The **Pareto β** column in Table 1 requires a separate GenAPoPop run (see note above)

---

## Citation

This repository accompanies the article below. If you use this workflow or the
Docker image, please cite **both** the article and the archived software.

**Article (preprint)**

> Abdalrahem, A., Andrieux, A., Becheler, R., Duplessis, S., Frey, P.,
> Marçais, B., Schiffer-Forsyth, K., Stoeckel, S., & Halkett, F. (2025).
> *Long-lasting coexistence of multiple asexual lineages alongside their sexual
> counterparts in a fungal plant pathogen.* bioRxiv 2025.03.28.645883.
> https://doi.org/10.1101/2025.03.28.645883

**Software / data archive (this repository)**

> Abdalrahem, A. (2025). *poplar_rust: population-genetics workflow for*
> Melampsora larici-populina [Software]. Zenodo.
> https://doi.org/10.5281/zenodo.15100450

---

## License

This project is distributed under the **CC BY 4.0** license.

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
