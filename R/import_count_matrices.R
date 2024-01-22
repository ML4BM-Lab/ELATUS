#' Import Kallisto raw count matrix
#' This function import Kallisto raw count matrix and save it in a Seurat object.
#' @param path Path to the Kallisto raw count matrix 
#' @param name Name of the Kallisto raw count matrix 
#' @return A Seurat object 

#' @export
import_kallisto_sc <- function(path, name)
{
    kallisto_data <- BUSpaRse::read_count_output(path, name = name)
    kallisto <- Seurat::CreateSeuratObject(kallisto_data, project = "kallisto")
    kallisto
}

#' Import CellRanger raw count matrix
#' This function import CellRanger raw count matrix and save it in a Seurat object.
#' @param path Path to the CellRanger raw count matrix 
#' @param name Name of the CellRanger raw count matrix 
#' @return A Seurat object 
#' @export
import_CellRanger_sc <- function(path)
{
    CellRanger_data <- Seurat::Read10X(data.dir = path, gene.column = 1)
    CellRanger <- Seurat::CreateSeuratObject(CellRanger_data, project = "CellRanger")
    CellRanger
}


