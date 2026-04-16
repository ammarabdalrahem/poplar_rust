# Population Genetics of *Melampsora larici-populina* — Microsatellite Analysis

An R Markdown workflow for population genetic analysis of poplar rust (*M. larici-populina*) isolates using microsatellite markers, including lineage definition, clustering, reproductive mode inference, spatial mapping, and regression.


This code workflow is part of the study titled:
## Long-lasting coexistence of multiple asexual lineages alongside their sexual counterparts in a fungal plant pathogen
Ammar Abdalrahem, Axelle Andrieux, Ronan Becheler, Sébastien Duplessis, Pascal Frey, Benoit Marcais, Kadiatou Schiffer-Forsyth, Solenn Stoeckel, Fabien Halkett

---

## What the Code Does

1. **Data loading** — Reads isolate metadata and microsatellite genotypes from `Table_data.tsv`
2. **MLG/MLL definition** — Defines multilocus genotypes and multilocus lineages using distance-based thresholding (`poppr`)
3. **Clustering** — Assigns isolates to genetic clusters via K-means + DAPC (`adegenet`); removes uncertain assignments (<80% posterior probability)
4. **Reproductive mode inference** — Classifies isolates as **asexual** or **sexual** using two approaches:
   - *Cluster approach*: negative Fis → asexual (indirect indication)
   - *Resampling approach*: MLL persisting across ≥ 2 years → asexual
   - Final label: asexual if flagged by either approach
5. **Population genetic indices** — Calculates N, G/N, Ar, Ho, Hs, Fis, rbarD per group (`hierfstat`, `poppr`)
6. **Phylogenetic tree** — Builds bootstrapped UPGMA/NJ tree annotated by lineage and reproductive mode (`ggtree`)
7. **Spatial analysis** — Maps sexual vs. asexual proportions across France (`mapdata`, `mapplots`)
8. **Regression** — GLM and GLMM testing effect of latitude, longitude, and year on reproduction mode (`lme4`)

---

## Dependencies

**R version ≥ 4.3** and **RStudio** (required for automatic path detection) are recommended.

Install all CRAN packages:

```r
install.packages(c(
  "knitr", "ggplot2", "readxl", "tidyverse", "genepop", "hierfstat",
  "mapdata", "mapplots", "data.table", "colorspace", "adegenet", "poppr",
  "pegas", "ape", "cowplot", "ade4", "remotes", "ggtreeExtra", "viridis",
  "factoextra", "openxlsx", "ggrepel", "RClone", "ggsci", "scales",
  "lme4", "dplyr", "reshape2", "rstudioapi"
))
```

Install the Bioconductor package `ggtree`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")
```

---

## Input Data

| File | Description |
|------|-------------|
| `Table_data.tsv` | Isolate metadata + microsatellite genotypes (must be in the same folder as the `.Rmd` file) |

- Microsatellite locus columns must contain `Mlp` in the name
- Missing genotype values must be coded as `999`

---

## How to Run

**In RStudio (recommended):**
1. Place `data_analysis_mlp.Rmd` and `Table_data.tsv` in the same folder
2. Open the `.Rmd` file in RStudio
3. Install dependencies (see above)
4. Click **Knit** → renders a full HTML report

**From the R console:**
```r
setwd("path/to/your/folder")
rmarkdown::render("data_analysis_mlp.Rmd")
```

---

## Outputs

| File | Description |
|------|-------------|
| `data_analysis_mlp.html` | Full rendered report |
| `new_genotype_data.csv` | Final genotype table with reproduction labels |
| `filtered_mll_years.csv` | MLLs recurring across multiple years |
| `cluster_assignments.png` | DAPC cluster probability plot |
| `tree_plot1/2.png` | Annotated phylogenetic tree |
| `map_all_years.png`, `map_2009_2011.png` | Spatial maps of reproduction mode |
| `effect_of_latitude_*.png` | Regression plots |

---
