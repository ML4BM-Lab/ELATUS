#' Calculate standard QC metrics
#' This function calculates standard QC metrics.
#' @param pipeline Kallisto or CellRanger SCE object 
#' @param mitochondrial_ens_ids Mitochondrial genes (ENS_IDs) 
#' @return SCE object with quality metrics 

#' @export
qc_metrics <- function(pipeline, mitochondrial_ens_ids )
{
  sce <- pipeline  
  sce <- scuttle::addPerCellQCMetrics(sce, subsets=list(Mito=which(rownames(sce)%in% mitochondrial_ens_ids)))
  sce
}