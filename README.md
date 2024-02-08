# ELATUS: Elucidating biologically relevant lncRNA annotated transcripts using scRNA-seq

This workflow is designed to, through droplet-based scRNA-seq, not only retain highly-expressed lncRNAs robustly detected by Cell Ranger and Kallisto, but also uncover functionally relevant lncRNAs that are only detected by Kallisto. 
![image](https://github.com/kikegoni/ELATUS/assets/35218702/64948f75-3cc0-4b25-9ad1-a619e48bb5c2)

## Installation
Installation be performed with devtools. Some bioconductor packages will be requiered, and devtools will not automatically install them. Therefore, you will have to install them manually, for example:
```{r}
BiocManager::install("scran")
```
Now you can install ELATUS with:
```{r}
install.packages("devtools")
devtools::install_github("kikegoni/ELATUS")
```

## Getting started
```{r}
library("ELATUS")
```

## Example
```{r}
ELATUS(kallisto_path=system.file("extdata", "kallisto_example_raw_matrix", package = "ELATUS"), kallisto_name="cells_genes_NO_multimapping", cellRanger_path=system.file("extdata", "cellRanger_example_raw_matrix", package = "ELATUS"), organism = "Mouse", lower_emptydrops = 1000, EmptyDrops_FDR_thres = 0.01, cells_mito_threshold= 15, cells_max_threshold = 30000, cells_min_genes_detected_threshold = 500, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = 5, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)
```


