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

## Associated study

**Long-lasting coexistence of multiple asexual lineages alongside their sexual counterparts in a fungal plant pathogen**

**Authors:**  
Ammar Abdalrahem, Axelle Andrieux, Ronan Becheler, Sébastien Duplessis, Pascal Frey, Benoit Marcais, Kadiatou Schiffer-Forsyth, Solenn Stoeckel, Fabien Halkett

## Repository contents

| File | Description |
|------|-------------|
| `data_analysis_mlp.Rmd` | Main R Markdown workflow |
| `data_analysis_mlp.R` | Optional R script version for terminal execution |
| `Table_data.tsv` | Input file with isolate metadata and microsatellite genotypes |

## Requirements

- R 4.3 or later recommended
- RStudio recommended for knitting the `.Rmd` file

Install the required CRAN packages:

```r
install.packages(c(
  "knitr", "ggplot2", "readxl", "tidyverse", "genepop", "hierfstat",
  "mapdata", "mapplots", "adegenet", "poppr", "pegas", "ape",
  "cowplot", "ade4", "ggtreeExtra", "viridis", "ggrepel", "RClone",
  "ggsci", "scales", "lme4", "dplyr", "factoextra", "sf",
  "rnaturalearth", "rnaturalearthdata", "grid", "rstudioapi",
  "reshape2", "svglite"
))
```

Install `ggtree` from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")
```

## Input format

The input file must be named:

```text
Table_data.tsv
```

Requirements:

- Microsatellite locus columns must contain `Mlp` in the column name
- Missing genotype values must be coded as `999`


## Quick start

### Option 1 — Run the `.Rmd` file in RStudio

1. Put `data_analysis_mlp.Rmd` and `Table_data.tsv` in the same folder
2. Open `data_analysis_mlp.Rmd` in RStudio
3. Install the required packages
4. Click **Knit**

This will generate the HTML report and output files.

### Option 2 — Run the `.R` script from the terminal

If you use the script version of the workflow (`data_analysis_mlp.R`), you can run it directly from the terminal:

```bash
Rscript data_analysis_mlp.R
```


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

## Notes

- The `.Rmd` file is the main reproducible workflow for generating the report
- The `.R` file is useful for command-line or terminal-based execution
- Isolates with uncertain cluster assignment are excluded from downstream analyses
- Most output files are written to the working directory


## License

This project is distributed under the **CC BY 4.0** license.


[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
