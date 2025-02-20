#' ELATUS workflow ofor snRNA-seq using a pre-mRNA reference:This workflow loads the Kallisto and CellRanger raw count matrices as well as Kallisto allowing for multimappers, performs emptydrops, then doublet removal. Next it gets the highly expressed lncRNAs only detected by both Kallisto and Kallisto multimappers and the ratio of their expression between Kallisto/CellRanger and finally, after normalization and clustering it calculates the Specificity Index of each of them. The last step is generating the list of biologically relevant lncRNAs. In addition it identifies highly-expressed lncRNAs exclusively detected by Kallisto whose functionality has been proven by external biobliography (CRISPR screenings...) and those highly-expressed lncRNAs robustly detected by both Cell Ranger and Kallisto and by Cell Ranger and Kallisto multimappers.
#' @param kallisto_path Path to the Kallisto raw count matrix (As an example: kallisto_path=/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.Kallisto/output_bus_transcriptome_introns/bustools_results_no_multimappers)
#' @param kallisto_name Name of the Kallisto raw count matrix (As an example kallisto_name="cells_genes_NO_multimapping")
#' @param kallisto_multimappers_path Path to the Kallisto raw count matrix (As an example: kallisto_multimappers_path="/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.Kallisto/output_bus_transcriptome_introns/bustools_multimappers")
#' @param kallisto_multimappers_name Name of the Kallisto raw count matrix (As an example kallisto_multimappers_name="cells_genes_multimapping")
#' @param cellRanger_path Path to the CellRanger raw count matrix (As an example :cellRanger_path="/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.CellRanger/pbmc_granulocyte_sorted_3k_alldata_cellRanger_modified_GEX_test2/outs/raw_feature_bc_matrix")
#' @param organism Human or Mouse dataset (for the example "Human") 
#' @param lower_emptydrops A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets. For the example, 1000
#' @param EmptyDrops_FDR_thres FDR threshold to classifly empty droplets (For the example, 0.01)
#' @param threshold_minumun_gene_counts Select genes with more than this total counts  (For the paper = 250 counts)
#' @param threshold_cells_detected Select genes present in at least a number of cells higher than this threshold (For the paper = 25)
#' @param dimred_clustering dimensionality reduction (For the example "PCA")
#' @param k_neighbors the number of nearest neighbors used to construct the graph. Choose a smaller k to more but smaller clusters as lncRNAs tend to be expressed in small subpopulations. (in this example, k=5)
#' @param ratio_threshold Threshold to remove lncRNAs whose ratio of expression between Kallisto and CellRanger is smaller than this defined threshold (For the paper = 20)
#' @param CR_threshold Threshold to remove lncRNAs that have a CellRanger expression higher than this defined threshold (For the paper = 20)
#' @param SI_threshold Threshold to remove lncRNAs whose SI is smaller than this defined threshold (For the paper = 0.15)
#' @return A list with most biologically relevant lncRNAs

#' @export
ELATUS_premRNA <- function(kallisto_path, kallisto_name, cellRanger_path,kallisto_multimappers_path,kallisto_multimappers_name,  organism, lower_emptydrops, EmptyDrops_FDR_thres, threshold_minumun_gene_counts, threshold_cells_detected, dimred_clustering, k_neighbors, ratio_threshold, CR_threshold, SI_threshold)
{
    if (organism == "Human")
    {
        gencode_path <- system.file("extdata", "hg38_v37.rds", package = "ELATUS")
    } 
    if (organism == "Mouse")
    {
        gencode_path <- system.file("extdata", "mm10_vM27.rds", package = "ELATUS")
    }
    crispr_data1 <- readRDS(system.file("extdata", "hits_info_Liu_science_2015_ensids.rds", package = "ELATUS"))
    crispr_data2 <- readRDS(system.file("extdata", "essential_lncRNAs_hits_info_Wei_Santana_cell_2024.rds", package = "ELATUS"))
    crispr_data <- as.data.frame(c(crispr_data1$gene_id, crispr_data2$gene_id))
    colnames(crispr_data) <- "gene_id"
    gtf <- readRDS(gencode_path)
    gtf$gene_id <- gsub("_","-",gtf$gene_id)
    mitochondrial_ens_ids <- unique(gtf$gene_id[grep("^MT-",gtf$gene_name)])
    mitochondrial_names <- unique(gtf$gene_name[grep("^MT-",gtf$gene_name)])
    lncrna_ens_ids <- unique(c(gtf$gene_id[grep("lncRNA",gtf$gene_type)]))
    protein_coding_ens_ids <- unique(c(gtf$gene_id[gtf$gene_type=="protein_coding"]))
    lncrna_names <- unique(gtf$gene_name[gtf$gene_id %in% lncrna_ens_ids])
    protein_coding_names <-  unique(gtf$gene_name[gtf$gene_id %in% protein_coding_ens_ids])

    kallisto <- import_kallisto_sc(kallisto_path, kallisto_name)
    kallisto <- Seurat::as.SingleCellExperiment(kallisto)
    kallisto_sce <- qc_metrics(kallisto, mitochondrial_ens_ids)

    cellRanger <- import_CellRanger_sc(cellRanger_path)
    cellRanger <- Seurat::as.SingleCellExperiment(cellRanger)
    cellRanger_sce <- qc_metrics(cellRanger, mitochondrial_ens_ids)

    kallisto_multimappers <- import_kallisto_sc(kallisto_multimappers_path, kallisto_multimappers_name)
    kallisto_multimappers <- Seurat::as.SingleCellExperiment(kallisto_multimappers)
    kallisto_multimappers_sce <- qc_metrics(kallisto_multimappers, mitochondrial_ens_ids)

    #EmptyDrops filtering
    #for the example given (human PBMCs) the parameters are lower_emptydrops=1000 & EmptyDrops_FDR_thres = 0.01
    kallisto_filt_sce_ed <- emptydrops_filt(kallisto_sce, lower = lower_emptydrops, EmptyDrops_FDR_thres = EmptyDrops_FDR_thres)
    kallisto_multimappers_sce_ed <- emptydrops_filt(kallisto_multimappers_sce, lower = lower_emptydrops, EmptyDrops_FDR_thres = EmptyDrops_FDR_thres)
    cellRanger_filt_sce_ed <- emptydrops_filt(cellRanger_sce, lower = lower_emptydrops, EmptyDrops_FDR_thres = EmptyDrops_FDR_thres)

    # Remove doublets
    kallisto_filt_sce_ed_nodoubs <- remove_doublets(kallisto_filt_sce_ed)
    kallisto_filt_sce_ed_nodoubs <- kallisto_filt_sce_ed_nodoubs[,kallisto_filt_sce_ed_nodoubs$isDoublet == F]
    cellRanger_filt_sce_ed_nodoubs <- remove_doublets(cellRanger_filt_sce_ed)
    cellRanger_filt_sce_ed_nodoubs <- cellRanger_filt_sce_ed_nodoubs[,cellRanger_filt_sce_ed_nodoubs$isDoublet == F]
    kallisto_multimappers_sce_ed_nodoubs <- remove_doublets(kallisto_multimappers_sce_ed)
    kallisto_multimappers_sce_ed_nodoubs <- kallisto_multimappers_sce_ed_nodoubs[,kallisto_multimappers_sce_ed_nodoubs$isDoublet == F]

    # Now get the highly expressed lncRNAs only detected by Kallisto and the ratio of their expression between Kallisto/CellRanger. Also get the highly-expressedcommonly detected by Cell Ranger and Kallisto
    # uniquifyFeatures
    gene_name <- gtf$gene_name[match(rownames(kallisto_filt_sce_ed_nodoubs),gtf$gene_id)]
    rownames(kallisto_filt_sce_ed_nodoubs) <- scuttle::uniquifyFeatureNames(rownames(kallisto_filt_sce_ed_nodoubs), gene_name)
    gene_name <- gtf$gene_name[match(rownames(cellRanger_filt_sce_ed_nodoubs),gtf$gene_id)]
    rownames(cellRanger_filt_sce_ed_nodoubs) <- scuttle::uniquifyFeatureNames(rownames(cellRanger_filt_sce_ed_nodoubs), gene_name)
    gene_name <- gtf$gene_name[match(rownames(kallisto_multimappers_sce_ed_nodoubs),gtf$gene_id)]
    rownames(kallisto_multimappers_sce_ed_nodoubs) <- scuttle::uniquifyFeatureNames(rownames(kallisto_multimappers_sce_ed_nodoubs), gene_name)

    # We considered highly expressed lncRNAs as those with at least 250 counts in at least 25 cells
    candidate_lncRNAs_exclusive <- get_candidates_premRNA(kallisto_filt_sce_ed_nodoubs, cellRanger_filt_sce_ed_nodoubs, kallisto_multimappers_sce_ed_nodoubs,threshold_minumun_gene_counts = threshold_minumun_gene_counts, threshold_cells_detected = threshold_cells_detected,lncrna_names = lncrna_names,gtf=gtf)
    candidate_lncRNAs_common <- get_candidates_premRNA(kallisto_filt_sce_ed_nodoubs, cellRanger_filt_sce_ed_nodoubs, kallisto_multimappers_sce_ed_nodoubs, threshold_minumun_gene_counts = threshold_minumun_gene_counts, threshold_cells_detected = threshold_cells_detected,lncrna_names = lncrna_names,gtf=gtf, exclusive = F)

    # clustering
    set.seed(100100100)
    kallisto_filt_sce_ed_nodoubs <- scater::runPCA(kallisto_filt_sce_ed_nodoubs) 
    g <- scran::buildSNNGraph(kallisto_filt_sce_ed_nodoubs, use.dimred = dimred_clustering, k = k_neighbors ) # k is the number of nearest neighbors used to construct the graph. Choose a smaller k to more but smaller clusters as lncRNAs tend to be expressed in small subpopulations. (in this example, k=5). dimred_clustering is the dimensionality reduction (PCA here, but could be the corrected space after integrating samples)
    clust <- igraph::cluster_louvain(g)$membership
    print(table(clust))
    kallisto_filt_sce_ed_nodoubs$louvain_clusters <- factor(clust)

    kallisto_multimappers_sce_ed_nodoubs <- scater::runPCA(kallisto_multimappers_sce_ed_nodoubs) 
    g <- scran::buildSNNGraph(kallisto_multimappers_sce_ed_nodoubs, use.dimred = dimred_clustering, k = k_neighbors ) # k is the number of nearest neighbors used to construct the graph. Choose a smaller k to more but smaller clusters as lncRNAs tend to be expressed in small subpopulations. (in this example, k=5). dimred_clustering is the dimensionality reduction (PCA here, but could be the corrected space after integrating samples)
    clust <- igraph::cluster_louvain(g)$membership
    print(table(clust))
    kallisto_multimappers_sce_ed_nodoubs$louvain_clusters <- factor(clust)

    # Calculate the Specificity Index for each gene
    SI <- SI(kallisto_filt_sce_ed_nodoubs,group_by="louvain_clusters", average_by="mean")
    SI_kallisto_multi <- SI(kallisto_multimappers_sce_ed_nodoubs,group_by="louvain_clusters", average_by="mean") # for the lncRNAs that kallisto by default doesn't detect due to conflicts with the premRNA annnotation
    cell_type_specific_score <- SI[["cell_type_specificity_score"]]
    counts_cell_specificity_index <- SI[["counts_cell_specificity_index"]]
    cell_type_specific_score_multi <- SI_kallisto_multi[["cell_type_specificity_score"]]
    counts_cell_specificity_index_multi <- SI_kallisto_multi[["counts_cell_specificity_index"]]

    candidate_lncRNAs_exclusive$SI <- cell_type_specific_score[rownames(candidate_lncRNAs_exclusive)]
    to_subset <- setdiff(rownames(candidate_lncRNAs_common),names(cell_type_specific_score))
    cell_type_specific_score <- c(cell_type_specific_score, cell_type_specific_score_multi[to_subset])
    candidate_lncRNAs_common$SI <- cell_type_specific_score[rownames(candidate_lncRNAs_common)]

    # To know in which cluster the SI is the highest : Functionality to be added
    #candidate_lncRNAs_exclusive$cell_type_SI <- colnames(counts_cell_specificity_index[rownames(candidate_lncRNAs_exclusive),])[apply(counts_cell_specificity_index[rownames(candidate_lncRNAs_exclusive),],1,which.max)]
    candidate_lncRNAs_exclusive <- crispr_info(crispr_data, candidate_lncRNAs_exclusive)
    #candidate_lncRNAs_common$cell_type_SI <- colnames(counts_cell_specificity_index[intersect(rownames(candidate_lncRNAs_common), rownames(counts_cell_specificity_index)),])[apply(counts_cell_specificity_index[intersect(rownames(candidate_lncRNAs_common), rownames(counts_cell_specificity_index)),],1,which.max)]

    # Get the biologically relevant lncRNAs from these candidates (In the paper we used the following parameters: ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)
    #ratio_threshold = 20
    #CR_threshold = 20
    #SI_threshold = 0.15
    exclusive_lncRNAs_CRISPRi <- candidate_lncRNAs_exclusive[candidate_lncRNAs_exclusive$crispr_intersection == T,]
    exclusive_biologically_relevant_lncRNAs <- biologically_relevant_lncRNAs(candidate_lncRNAs_exclusive, ratio_threshold,CR_threshold,SI_threshold)
    candidate_lncRNAs_common$crispr_intersection = "NA"
    if (nrow(exclusive_lncRNAs_CRISPRi)>0)
    {
        exclusive_lncRNAs_CRISPRi$category = "Exclusive_lncRNA_CRISPRi"
    }
    exclusive_biologically_relevant_lncRNAs$category = "Exclusive_lncRNA"
    candidate_lncRNAs_common$category = "Common_lncRNA"

    biologically_relevant_lncRNAs <- rbind(exclusive_biologically_relevant_lncRNAs, candidate_lncRNAs_common,exclusive_lncRNAs_CRISPRi)
    biologically_relevant_lncRNAs


}