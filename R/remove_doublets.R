#' Remove multiplets
#' This function import Kallisto raw count matrix and save it in a Seurat object.
#' @param data SCE object already filtered with EmptyDrops
#' @return SCE object without doublets

#' @export
remove_doublets <- function(data)
{
  data2 <- scDblFinder::scDblFinder(data, clusters=TRUE)
  data$isDoublet <- data2$scDblFinder.class == "doublet"
  data
}