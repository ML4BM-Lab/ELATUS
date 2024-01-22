#' Select highly expressed genes from a SCE object.
#' @param data SCE object
#' @param threshold_minumun_gene_counts Select genes with more than this total counts 
#' @param threshold_cells_detected Select genes present in at least a number of cells higher than this threshold 
#' @return SCE object with highly expressed genes

#' @export
top_genes <- function(data,threshold_minumun_gene_counts,threshold_cells_detected)
{
    out_data <- data[((Matrix::rowSums(SingleCellExperiment::counts(data)) > threshold_minumun_gene_counts) & (Matrix::rowSums(SingleCellExperiment::counts(data) != 0) > threshold_cells_detected)),]
    out_data
}