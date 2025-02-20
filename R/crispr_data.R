#' This function checks if highly-expressed exclusive lncRNAs have been functionally validated with external CRISPRi studies 
#' @param crispr_data Public CRISPRi data with lncRNAs that have been externally validated
#' @param candidates_lncRNAs Highly-expressed exclusive lncRNAS
#' @return Boolean vector on whether or not the highly-expressed exclusive lncRNAs has a hit in public CRISPRi studies
#' @export
crispr_info <- function(crispr_data, candidates_lncRNAs)
{
    v <- gsub("\\..*","",candidates_lncRNAs$gene)
    y <- gsub("\\..*","",crispr_data$gene_id)
    candidates_lncRNAs$crispr_intersection <- sapply(v, function(x) any(grepl(x, y)))
    return(candidates_lncRNAs)
}
