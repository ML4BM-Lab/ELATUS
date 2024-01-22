#' Filter SCE objects to remove cells with a lot of mitochrondrial counts, an extremely high number of counts and also cells having counts in very few genes
#' @param SCE SCE object to be filtered
#' @param cells_mito_threshold Mitochondrial content (%) threshold. Keep cells with less mitocondrial content than the defined threshold. 
#' @param cells_max_threshold Keep cells with less counts than the defined threshold. 
#' @param cells_min_genes_detected_threshold Keep cells with counts in more genes than the defined threhold. 
#' @return Filtered SCE object

#' @export

Filtering <- function(sce, cells_mito_threshold, cells_max_threshold, cells_min_genes_detected_threshold)
{
    total_max_mito_content <- sce$subsets_Mito_percent > cells_mito_threshold
    total_max_expression <- sce$sum > cells_max_threshold
    total_min_genes_detected <- sce$detected < cells_min_genes_detected_threshold

    discard <-  total_max_expression | total_max_mito_content | total_min_genes_detected
    print(table(discard))
    print(S4Vectors::DataFrame(total_max_expression=sum(total_max_expression),total_max_mito_content=sum(total_max_mito_content),total_min_genes_detected=sum(total_min_genes_detected), Total=sum(discard)))
    sce$discard <- discard
    #filter and normalize data
    sce_filt <- sce[,discard==F]
    sce_filt
}