#' Perform emptydrops filtering using https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html
#' @param sce SCE object 
#' @param lower A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets
#' @param EmptyDrops_FDR_thres FDR threshold to classifly empty droplets
#' @return SCE object without empty droplets 

#' @export
emptydrops_filt <- function(sce, lower, EmptyDrops_FDR_thres)
{
    e.out <- DropletUtils::emptyDrops(SingleCellExperiment::counts(sce), lower = lower, niters = 10000)
    is.cell <- e.out$FDR <= EmptyDrops_FDR_thres
    print(paste("cells_detected: ",sum(is.cell, na.rm=TRUE)))
    non_empty_drops <- which(e.out$FDR < EmptyDrops_FDR_thres)
    sce_filt <- sce[,non_empty_drops]
    sce_filt
}
