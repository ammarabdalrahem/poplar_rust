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


This repository contains a fully reproducible R workflow with Docker for cross-platform.


---

## Fully Reproducible with Docker

The entire analysis environment is packaged in a Docker image. No R installation, no package management, no dependency conflicts.

### Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (free, works on macOS / Windows / Linux)

### Quick start

```bash
# 1. Clone the repository (contains the input data)
git clone https://github.com/ammarabdalrahem/poplar_rust.git
cd poplar_rust

# 2. Run the analysis — one command
docker run --rm \
  -v $(pwd):/project \
  ghcr.io/ammarabdalrahem/poplar_rust:1.0
```

> **On Windows PowerShell**, replace `$(pwd)` with `${PWD}`.

All outputs (figures, tables, CSV files) will appear in your local `poplar_rust/` folder after the run completes.

---

## Associated study

**Long-lasting coexistence of multiple asexual lineages alongside their sexual counterparts in a fungal plant pathogen**

**Authors:**  
Ammar Abdalrahem, Axelle Andrieux, Ronan Becheler, Sébastien Duplessis, Pascal Frey, Benoit Marcais, Kadiatou Schiffer-Forsyth, Solenn Stoeckel, Fabien Halkett

---

## Repository contents

| File | Description |
|------|-------------|
| `data_analysis_mlp_new.R` | R script version for terminal / Docker execution |
| `Table_data.tsv` | Input file with isolate metadata and microsatellite genotypes |
| `Dockerfile` | Docker image definition for full reproducibility |

---

## Manual setup (without Docker)

If you prefer to run the analysis manually in R:

### Requirements

- R 4.4.1 or later recommended
- RStudio recommended for interactive work

### R packages

Install the required CRAN packages:

```r
install.packages(c(
  "lme4","knitr", "ggplot2", "readxl", "tidyverse",
  "genepop", "hierfstat", "mapdata", "mapplots",
  "adegenet", "poppr", "pegas", "ape",
  "cowplot", "ade4", "viridis", "ggrepel", "ggsci",
  "scales", "dplyr", "factoextra", "sf",
  "rnaturalearth", "rnaturalearthdata",
  "svglite", "BiocManager", "devtools", "reshape2",
  "ggpubr", "ggforce"
))
```

Install required Bioconductor packages (`ggtree` and `ggtreeExtra`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggtreeExtra"), ask = FALSE, update = FALSE)
```

Install RClone from GitHub:

```r
remotes::install_github("dbailleul/RClone", dependencies = TRUE)
```

### Run the script

```bash
Rscript data_analysis_mlp_new.R
```

---

## Input format

The input file must be named:

```text
Table_data.tsv
```

Requirements:

- Microsatellite locus columns must contain `Mlp` in the column name
- Missing genotype values must be coded as `999`

---

## Main outputs

| File | Description |
|------|-------------|
| `new_genotype_data.csv` | Final genotype table with cluster and reproduction labels |
| `filtered_mll_years.csv` | MLLs recurring across multiple years |
| `Fig1A_map.svg` | Sampling map |
| `Silhouette_kmeans.png` | Silhouette plot for cluster evaluation |
| `Dapc_xval.png` | Cross-validation plot for DAPC |
| `DAPC_scatter.png` | DAPC scatter plot |
| `DAPC_compoplot.png` | DAPC composition plot |
| `cluster_assignments.png` | Cluster assignment probability plot |
| `tree_plot1.png` | Circular phylogenetic tree |
| `tree_plot2.png` | Annotated phylogenetic tree |
| `map_all_years.png` | Map of sexual versus asexual proportions across all years |
| `map_2009_2011.png` | Map for 2009 and 2011 only |
| `effect_of_latitude_filtered.png` | GLM-based latitude effect plot |
| `effect_of_latitude_glmm.png` | GLMM-based latitude effect plot |
| `MST_data_mlp_pop_as_Reproduction_for_cloneEstimate.txt` | Export for downstream clone analysis |
| `asex_mll_Year.png` | Asexual lineage abundance across years |
| `asex_mll_Locations.png` | Asexual lineage abundance across locations |
| `asex_mll_Year_Locations.png` | Combined lineage abundance figure |

---

## Notes

- The `.R` script is the main reproducible workflow for generating the report
- The Docker image (`ghcr.io/ammarabdalrahem/poplar_rust:1.0`) contains R 4.4.1 and all required packages pre-installed
- Isolates with uncertain cluster assignment are excluded from downstream analyses
- Most output files are written to the working directory

---

## Citation

If you use this workflow or the associated Docker image, please cite:

> Abdalrahem, A., et al. (2026). Long-lasting coexistence of multiple asexual lineages alongside their sexual counterparts in a fungal plant pathogen.

---

## License

This project is distributed under the **CC BY 4.0** license.


[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
