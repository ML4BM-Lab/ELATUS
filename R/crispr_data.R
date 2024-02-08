#' This function checks if highly-expressed exclusive lncRNAs have been functionally validated with external CRISPRi studies 
#' @param crispr_data Public CRISPRi data with lncRNAs that have been externally validated
#' @param candidates_lncRNAs Highly-expressed exclusive lncRNAS
#' @return Boolean vector on whether or not the highly-expressed exclusive lncRNAs has a hit in public CRISPRi studies
#' @export
crispr_info <- function(crispr_data, candidates_lncRNAs)
{
    candidates_lncRNAs$crispr_intersection <- gsub("\\..*","",candidates_lncRNAs$gene) %in% gsub("\\..*","",crispr_data$gene_id)
    return(candidates_lncRNAs)
}
