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

## How to run

There are three routes. All three produce the same files in `output/figures/` and
`output/tables/` — choose whichever suits you.

| Route | R needed on your machine? | Best for |
|-------|---------------------------|----------|
| **1. Pre-built Docker image** | no | reproducing the published results with one command |
| **2. Docker, base image + this code** | no | verifying the analysis from the plain `rocker` base |
| **3. Locally in R** | yes | reading the code, inspecting objects, re-using parts |

> **Apple Silicon Macs (M1/M2/M3/M4).** The `rocker` base image is published for
> `linux/amd64` only, so every Docker command below passes `--platform=linux/amd64`
> and the `Dockerfile` pins the same platform. Without it Docker stops with
> *"no matching manifest for linux/arm64/v8"*. The container then runs under Rosetta
> emulation and is slower than on an Intel machine. The flag is harmless on Intel
> Macs, Windows and Linux, and it guarantees everyone builds the same architecture.

> **On Windows PowerShell**, replace `$(pwd)` with `${PWD}` in every command.

---

### Route 1 — pre-built image (fastest, nothing to install)

Every dependency is already baked into the image, so the analysis starts immediately.

```bash
docker pull --platform=linux/amd64 ghcr.io/ammarabdalrahem/poplar_rust:latest

docker run --rm --platform=linux/amd64 \
  -v "$(pwd)/output:/project/output" \
  ghcr.io/ammarabdalrahem/poplar_rust:latest
```

Results appear in `output/` in the folder you ran the command from.

To check what the image contains without running the analysis:

```bash
docker run --rm --platform=linux/amd64 \
  ghcr.io/ammarabdalrahem/poplar_rust:latest ls -la /project
```

### Route 2 — plain `rocker/geospatial` base image

This is the setup the PCI data editor used. It mounts this repository into a clean
`rocker/geospatial:4.4.1` container and runs the script there.

```bash
git clone https://github.com/ammarabdalrahem/poplar_rust.git
cd poplar_rust

docker run --rm --platform=linux/amd64 \
  -v "$(pwd)":/project \
  -w /project \
  rocker/geospatial:4.4.1 \
  Rscript data_analysis_mlp_new.R
```

> This route is considerably slower than Route 1: the base image does not contain
> the genetics packages, so R installs `adegenet`, `poppr`, `ggtree` and the rest
> inside the container on every run.

### Route 3 — locally in R

Requires **R 4.4.1 or later**. Missing packages are installed automatically the
first time the script runs, in dependency order (CRAN core → spatial → genetics →
Bioconductor → GitHub), which can take several minutes on a clean installation.

Clone the repository and make sure your working directory is the repository root,
so that `Table_S1.csv` and `tree_nj_boot1000.rds` are found:

```bash
git clone https://github.com/ammarabdalrahem/poplar_rust.git
cd poplar_rust
```

**3a. From the command line**

```bash
Rscript data_analysis_mlp_new.R
```

Run it top to bottom in a clean session. Sections share objects, so they cannot be
run out of order. To also see the printed tables and model summaries scroll past,
use `R -e 'source("data_analysis_mlp_new.R", echo = TRUE)'` instead.

**3b. In RStudio, with Knit**

1. Open **`data_analysis_mlp.Rmd`** in RStudio
2. Make sure the working directory is the repository root — the file does this
   itself when knitting, and for chunk-by-chunk work you can set it with
   **Session → Set Working Directory → To Source File Location**
3. Click **Knit** (or press **Cmd/Ctrl + Shift + K**)

Knitting runs the whole analysis and produces `data_analysis_mlp.html`: a single
report with a numbered, floating table of contents, every figure and table inline,
and foldable code blocks. The same files are written to `output/` as in the other
routes.

To work through the analysis step by step instead, run the chunks in order with the
green arrow on each chunk, or **Run → Run All**. Do not start from a chunk in the
middle — later sections depend on objects created earlier.

> `data_analysis_mlp.Rmd` and `data_analysis_mlp_new.R` contain identical code in
> identical order. The `.Rmd` is the reference version: edit it first, then
> regenerate the script with
> `knitr::purl("data_analysis_mlp.Rmd", output = "data_analysis_mlp_new.R", documentation = 2)`
> so the two cannot drift apart.

### How the image is built and pinned

The image is defined by the `Dockerfile` and rebuilt automatically by GitHub
Actions (`.github/workflows/docker-publish.yml`) on every push to `main` and on
each version tag, then pushed to `ghcr.io/ammarabdalrahem/poplar_rust`.
Pushing a tag such as `v1.0` publishes both `:1.0` and `:latest`.

For long-term reproducibility, package versions are frozen:

- **Base image:** `rocker/geospatial:4.4.1` (R 4.4.1 + the full geospatial stack, the environment verified by the data editor)
- **CRAN:** pinned to a dated [Posit Package Manager](https://packagemanager.posit.co/) snapshot (`PKG_SNAPSHOT` build arg) so the same versions resolve on every build
- **Bioconductor:** pinned to release `3.20`
- **GitHub packages:** `rnaturalearthhires` pinned to an exact commit

---

## Repository contents

| File | Description |
|------|-------------|
| `data_analysis_mlp_new.R` | Main R script — use for terminal / Docker execution |
| `data_analysis_mlp.Rmd` | R Markdown version — use for interactive work in RStudio |
| `Table_S1.csv` | **Input 1:** isolate metadata and microsatellite genotypes (Table S1) |
| `tree_nj_boot1000.rds` | **Input 2:** neighbour-joining tree, 1,000 bootstrap replicates, precomputed (Figure 3) |
| `Dockerfile` | Docker image definition for full reproducibility |
| `.dockerignore` | Files excluded from the Docker build context |
| `.github/workflows/docker-publish.yml` | CI: builds the image and publishes it to GHCR |
| `README.md` | This file |

---

## Installing the R packages by hand

Route 3 installs everything automatically on first run, so this section is only
needed if you would rather control the installation yourself, or if a layer failed
and you want to retry it on its own. The order matters: packages with heavy system
dependencies must be built before the packages that depend on them.

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

---

## Inputs

The workflow takes exactly two input files, both in the repository root and both
shipped inside the Docker image. Everything else is generated by the run, into
`output/figures/` and `output/tables/`, which are created automatically.

```text
Table_S1.csv            genotypes and metadata, 2,122 individuals
tree_nj_boot1000.rds    neighbour-joining tree, 1,000 bootstrap replicates
```

`Table_S1.csv` requirements:

- Microsatellite locus columns must contain `Mlp` in the column name
- Missing genotype values must be coded as `999`

`tree_nj_boot1000.rds` is read with `readRDS()` and used only for Figure 3. It is
precomputed because `poppr::aboot` with 1,000 replicates is too slow to run on a
laptop or inside the container; it was produced on a compute server. The `aboot()`
call that generated it is kept in the code, commented out immediately above the
`readRDS()` line, so the tree can be regenerated on a machine with enough
resources.

---

## Outputs

### Tables — written to `output/tables/`

| File | Description | Manuscript |
|------|-------------|------------|
| `Table1_a_genetic_indices_cluster_approach.csv` | Genetic indices per genetic cluster (clustering approach) | Table 1 |
| `Table1_b_genetic_indices_resampling_approach.csv` | Genetic indices per reproduction mode (resampling approach) | Table 1 |
| `Table1_c_genetic_indices_combination_approaches.csv` | Genetic indices per reproduction mode (combined approach) | Table 1 |
| `Table2_contingency_clustering_vs_resampling.csv` | Contingency table: cluster vs resampling assignments | Table 2 |
| `Table2_b_fisher_test_clustering_vs_resampling.csv` | Fisher's exact test: p-value, odds ratio, CI | Results 3.2 |
| `Table3_GLMM_binomial_regression.csv` | GLMM fixed-effect coefficients (binomial, Lat + Long, 2009 and 2011) | Table 3 |
| `Table3_b_GLMM_random_effect_year.csv` | GLMM random effect of Year (variance, SD) | Results 3.3 |
| `Table4_top7_asexual_MLLs.csv` | Characteristics of the seven most abundant asexual MLLs | Table 4 |
| `Table_S2_cavalaire_february_samplings.csv` | February samplings at Cavalaire-sur-Mer, by year and lineage | Table S2 |
| `Table_MLL_delineation_summary.csv` | Number of MLGs, number of MLLs, MLL distance thresholds | Results 3.1 |
| `Table_DAPC_crossvalidation_summary.csv` | Optimal PCs, mean correct reassignment, variance explained | Methods 2.3 |
| `Table_fisher_MLL_year_region.csv` | Fisher's exact tests: asexual MLL x year and x region | Results 3.5 |
| `Table_contingency_asexual_MLL_by_year.csv` | Contingency table behind the temporal Fisher test | — |
| `Table_contingency_asexual_MLL_by_region.csv` | Contingency table behind the spatial Fisher test | — |
| `Table_asexual_MLL_genetic_indices.csv` | Full genetic indices for all asexual MLLs | — |
| `genetic_indices_per_sexual_MLL.csv` | Genetic indices per sexual MLL | — |
| `new_genotype_data.csv` | Final per-isolate table with cluster and reproduction labels | — |
| `filtered_mll_years.csv` | MLLs recurring across multiple sampling years | — |
| `MST_data_mlp_pop_as_Reproduction_for_cloneEstimate.csv` | ClonEstiMate input, comma-separated (archive copy) | — |
| `MST_data_mlp_pop_as_Reproduction_for_cloneEstimate.txt` | ClonEstiMate input, tab-separated (load this one into the software) | — |
| `my_data.RData` | Snapshot of the working objects, for `load("output/my_data.RData")` when returning to the analysis | — |

> **Note on Table 1 — Pareto β column:** this statistic is computed by the external
> software [GenAPoPop](https://www6.inrae.fr/genapopop) and cannot be produced by
> this R script. Run GenAPoPop separately on the same isolate data and insert the
> resulting Pareto β values into the exported CSV before final publication.

### Figures — written to `output/figures/`

| File | Description | Manuscript |
|------|-------------|------------|
| `FigS1_geographic_distribution.svg` | Sampling map (all isolates) | Fig. S1 |
| `FigS2_Silhouette_kmeans.png` | Silhouette plot for k-means cluster evaluation | Fig. S2 |
| `FigS3_cluster_assignments.png` | Cluster assignment probability scatter plot | Fig. S3 |
| `FigS4_effect_of_latitude_glmm.png` | GLMM: proportion of sexual reproduction vs latitude | Fig. S4 |
| `Fig2_geographical_distribution_2009_2011.svg` | Map of sexual vs asexual proportions, 2009 and 2011 | Fig. 2 |
| `Fig3_tree_nj_plot.png` | Annotated NJ tree with MLL rings | Fig. 3 |
| `Fig4A_asex_mll_Locations.png` | Asexual lineage abundance across French regions (panel A) | Fig. 4A |
| `Fig4B_asex_mll_Year.png` | Asexual lineage abundance across years (panel B) | Fig. 4B |
| `Fig4_asex_mll_Year_Locations.png` | Combined two-panel lineage abundance figure | Fig. 4 |

> **Figure 1 of the manuscript** (the *M. larici-populina* life cycle) is a drawn
> schematic and is not produced by this code.

> **Note on Figure 3 — `tree_nj_boot1000.rds`:** the neighbour-joining tree with
> 1,000 bootstrap replicates (`poppr::aboot`) is too slow to run on a laptop, so
> it was computed separately on a compute server and the result saved as
> `tree_nj_boot1000.rds`, which ships with this repository. The script loads that
> file instead of recomputing the tree, and Figure 3 is drawn from it. The
> `aboot()` call that produced it is kept in the script, commented out directly
> above the `readRDS()` line, so the tree can be regenerated from scratch on a
> machine with enough resources:
>
> ```r
> tree2 <- aboot(final_GenInd, distance = diss.dist, tree = "nj",
>                missing = "asis", sample = 1000, quiet = TRUE, threads = 20)
> ```

---

## Citation

If you use this workflow or the associated Docker image, please cite:

> Abdalrahem, A., et al. (2026). Long-lasting coexistence of multiple asexual lineages alongside their sexual counterparts in a fungal plant pathogen.

---

## License

This project is distributed under the **CC BY 4.0** license.

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
