#' Extract the biologically relevant lncRNAs
#' This function extracts the most biologically relevant lncRNAs from a list of highly expressed lncRNAs only detected by Kallisto.
#' @param candidates List of highly expressed lncRNAs only detected by Kallisto
#' @param ratio_threshold Threshold to remove lncRNAs whose ratio of expression between Kallisto and CellRanger is smaller than this defined threshold.
#' @param CR_threshold Threshold to remove lncRNAs that have a CellRanger expression higher than this defined threshold.
#' @param SI_threshold Threshold to remove lncRNAs whose SI is smaller than this defined threshold.
#' @return most biologically relevant lncRNAs
#' @export
biologically_relevant_lncRNAs <- function(candidates,ratio_threshold,CR_threshold,SI_threshold)
{
    ratio_out <- candidates$ratio <= ratio_threshold
    CR_out <- candidates$cellRanger_total_expression >= CR_threshold
    SI_out <- candidates$SI <= SI_threshold
    discard <-  ratio_out | CR_out | SI_out
    filtered_candidates <- candidates[discard==F,]
    return(filtered_candidates)
}