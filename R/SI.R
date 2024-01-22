#' Calculate the Specificity Index
#' This function calculates the Specificity Index from a set of genes, knowing each expression across different defined groups of cells (normally the clusters).
#' @param sce Kallisto SCE object
#' @param group_by With respect to which colData is the SI being calculated (i.e. louvain clusters)
#' @param average_by To define if the sumCountsAcrossCells is standardized with respect to the mean, median or nothing
#' @return Candidate lncRNAs highly expressed according to Kallisto and not by CellRanger 
#' @export
SI <- function(sce,group_by, average_by)
{
    out <- scuttle::sumCountsAcrossCells(SingleCellExperiment::counts(sce), sce[[group_by]], average=average_by)
    counts_cell = as.data.frame(SummarizedExperiment::assay(out))
    counts_cell <- counts_cell[Matrix::rowSums(counts_cell[])>0,]
    counts_cell_specificity_index=t(apply(counts_cell, 1, function(x) x/(sum(x))))
    colnames(counts_cell_specificity_index) = out$ids

    cell_type_specificity_score <- c()
    for (i in 1:nrow(counts_cell_specificity_index))
    {
        s=0
        for (j in 1:ncol(counts_cell_specificity_index))
        {
            if(counts_cell_specificity_index[i,j]!=0)
            {
                s=s+(log(counts_cell_specificity_index[i,j],ncol(counts_cell_specificity_index))*counts_cell_specificity_index[i,j]) # Same results checked as the Shannon Entropy Specificity (HS): https://apcamargo.github.io/tspex/metrics/#fnref:9
            }
        }
        cell_type_specificity_score <- c(cell_type_specificity_score,1+s)
    }
    names(cell_type_specificity_score) <- rownames(counts_cell_specificity_index)

    return(list("cell_type_specificity_score" = cell_type_specificity_score,"counts_cell_specificity_index" = counts_cell_specificity_index))
    
}
