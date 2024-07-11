#' This function gets the highly expressed genes in snRNA-seq from each SCE and outputs those lncRNAs which are only detected by Kallisto and kallisto multimappers adding more information about them.
#' @param kallisto_sce Kallisto SCE object
#' @param cellRanger_sce CellRanger SCE object
#' @param kallisto_multi_sce Kallisto multimappers SCE object
#' @param threshold_minumun_gene_counts Select genes with more than this total counts 
#' @param threshold_cells_detected Select genes present in at least a number of cells higher than this threshold 
#' @param lncrna_names LncRNAs gene names
#' @param gtf GTF annotation
#' @param exclusive to get the highly-expressed lncRNAs exclusive detected with Kallisto or to get thoso commonly detected by Cell Ranger and Kallisto
#' @return Candidate lncRNAs highly expressed  
#' @export
get_candidates_premRNA <- function(kallisto_sce,cellRanger_sce,kallisto_multi_sce,threshold_minumun_gene_counts, threshold_cells_detected, lncrna_names, gtf, exclusive = T)
{
    kallisto_top_genes <- top_genes(kallisto_sce,threshold_minumun_gene_counts,threshold_cells_detected )
    cellRanger_top_genes <- top_genes(cellRanger_sce,threshold_minumun_gene_counts,threshold_cells_detected )
    kallisto_multi_top_genes <- top_genes(kallisto_multi_sce,threshold_minumun_gene_counts,threshold_cells_detected )

    # Uniquely vs common
    input_list <- list(CellRanger = rownames(cellRanger_top_genes), Kallisto = rownames(kallisto_top_genes), Kallisto_multi = rownames(kallisto_multi_top_genes))

    if (exclusive == TRUE)
    {
        candidates_kallisto <- setdiff(intersect(input_list[[2]], input_list[[3]]), input_list[[1]])
    }
    else 
    {
        candidates_kallisto <- c(intersect(input_list[[2]], input_list[[1]]), intersect(input_list[[3]], input_list[[1]]))
    }
    candidates_kallisto_lncRNAs <- intersect(candidates_kallisto, lncrna_names)

    expression_k <- Matrix::rowSums(SingleCellExperiment::logcounts(kallisto_sce))[candidates_kallisto_lncRNAs]
    expression_CR <- Matrix::rowSums(SingleCellExperiment::logcounts(cellRanger_sce))[candidates_kallisto_lncRNAs]
    ratio_candidates <- (expression_k+1)/(expression_CR+1)
    candidates_kallisto_lncRNAs_df <- data.frame(candidates =candidates_kallisto_lncRNAs, ratio =  ratio_candidates)

    candidates_kallisto_lncRNAs_df$gene <- gtf$gene_id[match(candidates_kallisto_lncRNAs_df$candidates, gtf$gene_name)]
    candidates_kallisto_lncRNAs_df$kallisto_total_expression <- Matrix::rowSums(SingleCellExperiment::logcounts(kallisto_sce))[rownames(candidates_kallisto_lncRNAs_df)]
    candidates_kallisto_lncRNAs_df$cellRanger_total_expression <- Matrix::rowSums(SingleCellExperiment::logcounts(cellRanger_sce))[rownames(candidates_kallisto_lncRNAs_df)]
    
    candidates_kallisto_lncRNAs_df
}