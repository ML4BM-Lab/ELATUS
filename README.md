# ARKAS: Ascertaining Relevant lncRNAs from Kallisto Applying ScRNA-seq

This workflow is designed to, through droplet-based scRNA-seq, uncover functionally relevant lncRNAs that are only detected by Kallisto and missed by CellRanger. 

## Installation
Installation be performed with devtools. Some bioconductor packages will be requiered, and devtools will not automatically install them. Therefore, you will have to install them manually, for example:
```{r}
BiocManager::install("scran")
```
Now you can install ARKAS with:
```{r}
install.packages("devtools")
devtools::install_github("kikegoni/ARKAS_repository")
```

## Getting started
```{r}
library("ARKAS")
```

## Example
```{r}
ARKAS(kallisto_path=system.file("extdata", "kallisto_example_raw_matrix", package = "ARKAS"), kallisto_name="cells_genes_NO_multimapping", cellRanger_path=system.file("extdata", "cellRanger_example_raw_matrix", package = "ARKAS"), organism = "Mouse", lower_emptydrops = 1000, EmptyDrops_FDR_thres = 0.01, cells_mito_threshold= 15, cells_max_threshold = 30000, cells_min_genes_detected_threshold = 500, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = 5, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)
```


