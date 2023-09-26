runGRaNIE <- function(dir_output = "output_GRaNIE", 
                      datasetName = "undescribed",
                      file_peaks, file_rna, file_metadata,
                      TFBS_source = "custom",
                      TFBS_folder = NULL,
                      TFBS_JASPAR_useSpecificTaxGroup = NULL,
                      genomeAssembly = "hg38",
                      normalization_peaks = "DESeq2_sizeFactors", 
                      idColumn_peaks = "peakID",
                      normalization_rna = "limma_quantile", 
                      idColumn_RNA =  "ENSEMBL",
                      includeSexChr = FALSE,
                      minCV = 0,
                      minNormalizedMean_peaks = 5,
                      minNormalizedMean_RNA = 1,
                      minSizePeaks = 5,
                      corMethod = "pearson",
                      promoterRange = 250000, 
                      useGCCorrection = FALSE,
                      TF_peak.fdr.threshold = 0.2,
                      peak_gene.fdr.threshold = 0.1,
                      runTFClassification = FALSE,
                      runNetworkAnalyses = FALSE,
                      nCores = 4,
                      forceRerun = TRUE
) {
  
  checkmate::assertSubset(TFBS_source, c("custom", "JASPAR"))
  
  
  file_GRN = paste0(dir_output, "/GRN.qs")
  
  startFresh = TRUE
  
  if (file.exists(file_GRN) & !forceRerun) {
    GRN = qs::qread(file_GRN)
    startFresh = FALSE
    
  }
  if (startFresh) {
    
    countsATAC   = readr::read_tsv(file_peaks)
    countsRNA    = readr::read_tsv(file_rna)
    metadata.all = readr::read_tsv(file_metadata) 
    
    
    # Arbitrary list with information and metadata that is stored within the GRN object
    metadata.l = list(name = datasetName,
                      file_peaks = file_peaks,
                      file_rna  = file_rna,
                      genomeAsembly = genomeAssembly,
                      file_metadata = file_metadata 
    )
    
    
    GRN = initializeGRN(objectMetadata = metadata.l, 
                        outputFolder = dir_output,
                        genomeAssembly = genomeAssembly)
    
    
    GRN = addData(GRN,
                  counts_peaks = countsATAC, normalization_peaks = normalization_peaks, idColumn_peaks = idColumn_peaks,
                  counts_rna = countsRNA, normalization_rna = normalization_rna, idColumn_RNA = idColumn_RNA,
                  sampleMetadata = metadata.all, allowOverlappingPeaks = TRUE, 
                  forceRerun = forceRerun)
    
    GRN = plotPCA_all(GRN, data = c("rna", "peaks"), topn = 500, type = "normalized", removeFiltered = FALSE, forceRerun = forceRerun)
    
    if (TFBS_source == "custom") {
      
      # Should the pipeline be run for only a subset of TFs or all? The special keyword "all" will use all TF that are found in the HOCOMOCO folder; however, if only a subset should be considered, specify the subset here with c() and the TF names, as shown below
      TFs = "all"
      
      if (is.null(TFBS_folder)) {
        # Base directory of the folder with the TFBS predictions. 
        # The TFBS predictions are expected as *.bed files as well as a translation table with the name translationTable.csv
        # We provide all files here: https://www.embl.de/download/zaugg/GRN/hg19_hg38_mm10_PWMScan.zip (7.5 GB)
        # Make sure they are in the same genome assembly as the peaks data
        hocomocoVersion = dplyr::if_else(genomeAssembly == "hg38", "v11", "v10")
        motifFolder = paste0("/g/zaugg/zaugg_shared/annotations/TFBS/", genomeAssembly, "/PWMScan_HOCOMOCO", hocomocoVersion)
        
      } else {
        motifFolder = TFBS_folder
      }
      
      
      GRN = addTFBS(GRN, source = TFBS_source, motifFolder = motifFolder, TFs = TFs, filesTFBSPattern = "_TFBS", fileEnding = ".bed", forceRerun = forceRerun)
      
    } else {
      
      GRN = addTFBS(GRN, source = TFBS_source, JASPAR_useSpecificTaxGroup = TFBS_JASPAR_useSpecificTaxGroup, JASPAR_removeAmbiguousTFs = TRUE, 
                    forceRerun = forceRerun)
    }
    
    ######################
    # PEAKS TFBS OVERLAP #
    ######################
    
    GRN = overlapPeaksAndTFBS(GRN, nCores = nCores, forceRerun = forceRerun)
    
    qs::qsave(GRN, file_GRN)
  }
  
  if (length(GRN@connections) == 0) {
    
    
    # Chromosomes to include for peaks and peak-gene associations. This should be a vector of chromosome names
    if (includeSexChr) {
      if (stringr::str_starts(genomeAssembly, "mm")) {
        chrToKeep_peaks = c(paste0("chr", 1:19), "chrX", "chrY")
      } else {
        chrToKeep_peaks = c(paste0("chr", 1:22), "chrX", "chrY")
      }
      
    } else {
      if (stringr::str_starts(genomeAssembly, "mm")) {
        chrToKeep_peaks = c(paste0("chr", 1:19))
      } else {
        chrToKeep_peaks = c(paste0("chr", 1:22))
      }
    }
    
    # If number of connections is 0 here, catch itand return an empty connctions dataframe
    GRN = tryCatch( {
      
      filterData(GRN, minNormalizedMean_peaks = minNormalizedMean_peaks, minNormalizedMeanRNA = minNormalizedMean_RNA, 
                 chrToKeep_peaks = chrToKeep_peaks, minSize_peaks = minSizePeaks,
                 minCV_peaks = minCV, minCV_genes = minCV, forceRerun = forceRerun)
    }, error = function(e) {
      warning("An error occured while executing the function filterData. This usually means there were no connections left. An empty file for the connections has been created.")
      file_connections = paste0(dir_output, "/connections_TFPeak", TF_peak.fdr.threshold, "_peakGene", peak_gene.fdr.threshold, ".tsv.gz")
      file.create(file_connections)
      qs::qsave(GRN, file_GRN)
      return(GRN)
    }
    )
    
    
    GRN = addConnections_TF_peak(GRN, connectionTypes = c("expression"), plotDiagnosticPlots = FALSE, plotDetails = FALSE, 
                                 corMethod = corMethod, maxFDRToStore = 0.3,
                                 useGCCorrection = useGCCorrection, percBackground_size = 75, percBackground_resample = TRUE,
                                 forceRerun = forceRerun)
    
    
    file_input_TADs = ""
    overlapTypeGene = "TSS"
    
    GRN = addConnections_peak_gene(GRN,
                                   overlapTypeGene = overlapTypeGene,
                                   corMethod = corMethod, shuffleRNACounts = TRUE,
                                   promoterRange = promoterRange, TADs = NULL,
                                   nCores = nCores, plotDiagnosticPlots = TRUE,
                                   forceRerun = forceRerun)
    
    
    GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = TF_peak.fdr.threshold, 
                                   TF_peak.connectionTypes = "expression" ,
                                   peak_gene.fdr.threshold = peak_gene.fdr.threshold,
                                   gene.types = c("protein_coding"),
                                   allowMissingTFs = FALSE, allowMissingGenes = FALSE,
                                   peak_gene.r_range = c(0,1))
    
    
    GRN = add_TF_gene_correlation(GRN, corMethod = corMethod, nCores = nCores)
    
    
    
    
    qs::qsave(GRN, file_GRN)
  } 
  
  file_connections = paste0(dir_output, "/connections_TFPeak", TF_peak.fdr.threshold, "_peakGene", peak_gene.fdr.threshold, ".tsv.gz")
  
  if (!file.exists(file_connections) | forceRerun) {
    
    connections.df = getGRNConnections(GRN, 
                                       include_TF_gene_correlations = TRUE, 
                                       include_peakMetadata = TRUE, 
                                       include_TFMetadata = TRUE, 
                                       include_geneMetadata = TRUE)
    readr::write_tsv(connections.df, file_connections)
  }
  
  
  
  if (runTFClassification) {
    
    if (forceRerun || is.null(GRN@data$TFs$classification)) {
      GRN = AR_classification_wrapper(GRN, corMethod = corMethod, forceRerun = forceRerun)
      qs::qsave(GRN, file_GRN)
    }
  }
  
  if (runNetworkAnalyses) {
    
    if (forceRerun || is.null(GRN@stats[["Enrichment"]]) || is.null(GRN@graph$TF_gene) | is.null(GRN@graph$TF_peak_gene) ) {
      GRN = performAllNetworkAnalyses(GRN, forceRerun = forceRerun)
      qs::qsave(GRN, file_GRN)
    } 
    
  }
  
  
  # GRN = visualizeGRN(GRN, plotAsPDF = FALSE)
  
  
  
  GRN
  
}


# algorithm: 1 = original Louvain algorithm
# algorithm: 2 = Louvain algorithm with multilevel refinement
# algorithm: 3 = SLM algorithm
# algorithm: 4 = Leiden algorithm). Leiden requires the leidenalg python.
# Algorithms discussed here: https://arxiv.org/pdf/1308.6604.pdf

#' @import checkmate
#' @import futile.logger
#' @import readr
#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
#' @importFrom ggplot2 `%+%`
#' @export
prepareSeuratData_GRaNIE <- function(seu.s, outputDir = "pseudobulk", saveSeuratObject = TRUE,
                                     file_RNA_features = "/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_hg38.tsv.gz", 
                                     assayName_RNA = "RNA", assayName_ATAC= "ATAC", 
                                     prepareData = TRUE, SCT_nDimensions = 50,  dimensionsToIgnore_LSI_ATAC = 1,
                                     pseudobulk_source = "cluster",
                                     countAggregation = "mean",
                                     clusteringAlgorithm = 3, clusterResolutions = c(0.1, seq(0.25, 1, 0.25), seq(2,20,2)), 
                                     minClusters = NULL,
                                     minCellsPerCluster = 25,
                                     subsample_percentage = 100, subsample_n = 1,
                                     doDimPlots = TRUE,
                                     forceRerun = FALSE
) {
    
    start = Sys.time()
    
    checkmate::assertClass(seu.s, "Seurat")
    checkmate::assertSubset(c(assayName_RNA, assayName_ATAC), names(seu.s@assays))
    checkmate::assertFlag(prepareData)
    checkmate::assertFileExists(file_RNA_features, access = "r")
    checkmate::assertSubset(pseudobulk_source, c("cluster", colnames(seu.s@meta.data)))
    
    checkmate::assertFlag(saveSeuratObject)
    checkmate::assertChoice(countAggregation, c("sum", "mean"))
    
    checkmate::assertSubset(clusteringAlgorithm, c(1:4), empty.ok = FALSE)
    checkmate::assertIntegerish(clusteringAlgorithm, len = 1)
    checkmate::assertNumeric(clusterResolutions, any.missing = FALSE, min.len = 1, lower = 0, null.ok = TRUE)
    
    checkmate::assertIntegerish(minClusters, lower = 1, null.ok = TRUE)
    checkmate::assertIntegerish(minCellsPerCluster, lower = 1)
    
    checkmate::assertIntegerish(SCT_nDimensions, lower = 5, upper = 500)
    checkmate::assertIntegerish(dimensionsToIgnore_LSI_ATAC, lower = 0, upper = SCT_nDimensions - 1)
    
    checkmate::assertNumber(subsample_percentage, lower = 0, upper = 100)
    checkmate::assertNumber(subsample_n, lower = 1)
    checkmate::assertFlag(doDimPlots)
    checkmate::assertFlag(forceRerun)
    
    if (!dir.exists(outputDir)) {
        dir.create(outputDir)
    }
    
    GRaNIE:::.startLogger(paste0(outputDir, "/prepareData_Seurat_GRaNIE.R.log") , "INFO",  removeOldLog = FALSE)
    
    if (saveSeuratObject & !GRaNIE:::is.installed("qs")) {
        stop("Package qs is not installed but needed due to saveSeuratObject = TRUE")
    }
    
    # Currently hard-coded options
    # See SCT transform
    returnOnlyVarGenes = FALSE
    
    # Needed to retrieve ENSEMBL IDs for gene names
    if (!is.null(file_RNA_features)) {
        futile.logger::flog.info(paste0("Reading and checking the RNA features file."))
        
        features = readr::read_tsv(file_RNA_features, col_names = FALSE, col_types = "ccffii")
        colnames(features) = c("ens_id", "gene", "view", "chr", "start", "end")
        
    }
    
    
    
    if (prepareData) {
        
        futile.logger::flog.info(paste0("Preparing RNA and ATAC data (RNA: SCTransform, PCA, UMPA; ATAC: TFID, FindTopFeatures, SVD, UMAP; combined: FindMultiModalNeighbors, UMAP). This may take a while."))
        
        Seurat::DefaultAssay(seu.s) <- assayName_RNA
        
        dimToUseRNA = seq_len(SCT_nDimensions)
        
        # Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. 
        # To save memory, we store these values only for variable genes, by setting the return.only.var.genes = TRUE by default in the SCTransform() function call.
        seu.s <- Seurat::SCTransform(seu.s, verbose = FALSE, return.only.var.genes = returnOnlyVarGenes) %>% 
            Seurat::RunPCA() %>% 
            Seurat::RunUMAP(dims = dimToUseRNA, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
        
        
        Seurat::DefaultAssay(seu.s) <- assayName_ATAC
        seu.s <- Signac::RunTFIDF(seu.s)
        seu.s <- Signac::FindTopFeatures(seu.s, min.cutoff = 'q0')
        seu.s <- Signac::RunSVD(seu.s)
        
        
        dimToUse_ATAC = setdiff(1:SCT_nDimensions, dimensionsToIgnore_LSI_ATAC)
        
        seu.s <- Seurat::RunUMAP(seu.s, reduction = 'lsi', dims = dimToUse_ATAC, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
        
        # We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
        seu.s <- Seurat::FindMultiModalNeighbors(seu.s, reduction.list = list("pca", "lsi"), dims.list = list(dimToUseRNA, dimToUse_ATAC),
                                                 weighted.nn.name = "weighted.nn",
                                                 knn.graph.name = "wknn",
                                                 snn.graph.name = "wsnn")
        
        seu.s <- Seurat::RunUMAP(seu.s, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        
        
    } # end if prepareData
    
    # Check whether the Seurat object is suffering from a bug introduced by the motifs object from scenic
    
    if ("motifs" %in% slotNames(seu.s@assays$ATAC) & !is.null(seu.s@assays$ATAC@motifs)) {
        names(seu.s@assays$ATAC@motifs@motif.names) <- seu.s@assays$ATAC@motifs@motif.names
        names(seu.s@assays$ATAC@motifs@pwm) <- seu.s@assays$ATAC@motifs@motif.names
        names(seu.s@assays$ATAC@motifs@positions) <- seu.s@assays$ATAC@motifs@motif.names
    }
    
    
    if (saveSeuratObject) {
        
        file_seurat = paste0(outputDir, "/seuratObject.qs")
        qs::qsave(seu.s,  file_seurat)
        
    } 
    
    
    ###############
    # Subsampling #
    ###############
    
    seu.orig.s = seu.s
    
    if (subsample_percentage == 100) subsample_n = 1
    
    for (subsampleCur in seq_len(subsample_n)) {
        
        
        subsample_suffix = dplyr::if_else(subsample_percentage != 100 | subsample_n > 1, 
                                          paste0("_subsampling", subsample_percentage, "_", subsampleCur), 
                                          "")
        
        
        # Currently the only option that works. Count aggregation
        sumCounts  = countAggregation == "sum"
        aggregationType = dplyr::if_else(sumCounts == TRUE, "_sum", "_mean")
        
        
        # Choose assayName according to aggregation type
        
        # Whether scale.data or counts is better remains to be determined
        slotName_aggregation_RNA  = dplyr::if_else(sumCounts, "counts", "scale.data")
        slotName_aggregation_RNA  = dplyr::if_else(sumCounts, "counts", "counts")
        
        # After reading Github issues, https://github.com/stuart-lab/signac/discussions/740 in particular, we use only raw counts for now
        # scale.data is not defined for ATAC, only counts and data after the TF-IDF normalization
        slotName_aggregation_ATAC  = dplyr::if_else(sumCounts, "counts", "counts")
        
        if (pseudobulk_source == "cluster") {
            
            searchClusterResolution = FALSE
            searchClusterResolutionStop = FALSE
            if (is.null(clusterResolutions)) {
                
                checkmate::assertIntegerish(minClusters, lower = 1)
                futile.logger::flog.info(paste0("Running clustering for the smallest resolution with at least: ", minClusters, ". This may take a while."))
                searchClusterResolution = TRUE
                clusterResolutions = c(0.1, seq(0.25, 1, 0.25), seq(2,20,2))
                
            } else {
                futile.logger::flog.info(paste0("Running clustering for the following resolutions: ", paste0(clusterResolutions, collapse = ", "), ". This may take a while."))
            }
            
            
            for (clusterResolution in unique(clusterResolutions)) {
                
                if (searchClusterResolutionStop) {
                    next
                }
                
                if (searchClusterResolution) {
                    clusterStr = paste0("minClusters", minClusters)
                } else {
                    clusterStr = paste0("res", clusterResolution)
                }
                
                file_in_rna   = paste0(outputDir, "/rna.pseudobulkFromClusters_", clusterStr, aggregationType, subsample_suffix, ".tsv.gz")
                file_in_atac  = paste0(outputDir, "/atac.pseudobulkFromClusters_", clusterStr, aggregationType, subsample_suffix, ".tsv.gz")
                file_metadata = paste0(outputDir, "/metadata_", clusterStr, aggregationType, subsample_suffix, ".tsv.gz")
                
                if (file.exists(file_in_rna) & file.exists(file_in_atac)  & file.exists(file_metadata) & !forceRerun) {
                    futile.logger::flog.info(paste0(" Files ", file_in_rna, ", ", file_in_atac, " and ", file_metadata, " already found, not overwriting due to forceRerun = TRUE. Returning Seurat object before clustering."))
                    next
                }
                futile.logger::flog.info(paste0(" Resolution ", clusterResolution))
                
                seu.s = seu.orig.s
                
                seu.s <- Seurat::FindClusters(seu.s, graph.name = "wsnn", 
                                              resolution = clusterResolution,
                                              algorithm = clusteringAlgorithm, verbose = FALSE)
                
                clusterColName = paste0("wsnn_res.", clusterResolution)
                Seurat::Idents(seu.s) = clusterColName
                seu.s@meta.data[[clusterColName]] = paste0("cluster", seu.s@meta.data[[clusterColName]])
                seu.s@meta.data$seurat_clusters = paste0("cluster", seu.s@meta.data$seurat_clusters)
                
                # Only needed once as NOT dependent on subsampling
                if (doDimPlots & subsampleCur == 1) {
                    pdfFile = paste0(outputDir, "/dimPlot_combined_clusterResolution", clusterResolution, ".pdf")
                    .doDimPlots(seu.s, groupBy = "seurat_clusters", pdfFile = pdfFile)
                }
                
                nClusters = length(which(table(Seurat::Idents(seu.s)) %>% as.vector() >= minCellsPerCluster))
                futile.logger::flog.info(paste0(" Number of clusters with at least ", minCellsPerCluster, " cells per cluster: ", nClusters))
                
                if (searchClusterResolution) {
                    if (nClusters < minClusters) {
                        next
                    } else {
                        futile.logger::flog.info(paste0("  Minimum resolution found with at least ", minClusters, " clusters"))
                        searchClusterResolutionStop = TRUE
                    }
                }
                
                # Subsampling
                seu.s = .subsampleSeuratObject(seu.s, subsample_percentage = subsample_percentage)
                if (subsample_percentage != 100) {
                    nClusters = length(which(table(Seurat::Idents(seu.s)) %>% as.vector() >= minCellsPerCluster))
                    futile.logger::flog.info(paste0(" Number of clusters with at least ", minCellsPerCluster, " cells per cluster after subsampling: ", nClusters))
                }
                
                # Metadata
                metadata = .writeMetadata(seu.s, file_metadata, minCellsPerCluster = minCellsPerCluster)
                
                ### Count aggregation
                futile.logger::flog.info(paste0(" Agggregate and prepare RNA counts for each cluster"))
                
                assayName_aggregation = dplyr::if_else(sumCounts, assayName_RNA, "SCT")
                rna.pseudobulk.clean = .aggregateCounts(seu.s, assayName = assayName_aggregation, slotName = slotName_aggregation_RNA,
                                                        groupBy = clusterColName, sumCounts = sumCounts, ID_column = "gene")
                
                # Delete clusters that contain too few cells
                rna.pseudobulk.clean = dplyr::select(rna.pseudobulk.clean, "gene", dplyr::one_of(metadata$sampleName))
                
                # rna.before = GetAssayData(object = seu.s, assay = "RNA", slot = "counts")
                # rna.before2 = GetAssayData(object = seu.s, assay = "RNA", slot = "data")
                # rna.before3 = GetAssayData(object = seu.s, assay = "RNA", slot = "scale.data")
                # 
                # # Sanity check
                # for (i in 1:10) {
                #   stopifnot(identical(sum(rna.before[i,]), sum(rna.pseudobulk$RNA[i,]))) 
                # }
                
                # Merge with the actual features and their correct mappings.
                rna.pseudobulk.clean = .addEnsemblIDs(rna.pseudobulk.clean, mapping = features)
                
                futile.logger::flog.info(paste0(" Writing RNA counts to file ", file_in_rna))
                readr::write_tsv(rna.pseudobulk.clean, file_in_rna)
                
                .printScarcity(rna.pseudobulk.clean, type = "RNA")
                
                ########
                # ATAC #
                ########
                
                futile.logger::flog.info(paste0(" Subsample, aggregate and prepare ATAC counts for each cluster"))
                atac.pseudobulk.clean = .aggregateCounts(seu.s, assayName = assayName_ATAC, slotName = slotName_aggregation_ATAC,
                                                         groupBy = clusterColName, sumCounts = sumCounts, ID_column = "peakID")
                
                # Delete clusters that contain too few cells
                atac.pseudobulk.clean = dplyr::select(atac.pseudobulk.clean, "peakID", dplyr::one_of(metadata$sampleName))
                
                
                # Replace the first hyphen with a colon
                atac.pseudobulk.clean$peakID = sub("-", ":", atac.pseudobulk.clean$peakID)
                
                futile.logger::flog.info(paste0(" Writing ATAC counts to file ", file_in_atac))
                readr::write_tsv(atac.pseudobulk.clean, file_in_atac)
                
                .printScarcity(atac.pseudobulk.clean, type = "ATAC")
                
                
            } # end for (clusterResolution in clusterResolutions)
            
            
        }  # end if (pseudobulk_source == "cluster")
        else {
            
            futile.logger::flog.info(paste0("Creating pseudobulk based on the following pre-existing cell identity column: ", pseudobulk_source))
            
            file_in_rna  = paste0(outputDir, "/rna.pseudobulkFromClusters_" , pseudobulk_source, aggregationType, subsample_suffix, ".tsv.gz")
            file_in_atac = paste0(outputDir, "/atac.pseudobulkFromClusters_", pseudobulk_source, aggregationType, subsample_suffix, ".tsv.gz")
            
            if (file.exists(file_in_rna) & file.exists(file_in_atac)) {
                futile.logger::flog.info(paste0(" Files ", file_in_rna, " and ", file_in_atac, " already found, not overwriting due to forceRerun = TRUE. Returning Seurat object before clustering."))
                
            } else {
                
                seu.s = seu.orig.s
                # Normalize levels and make the names compatible
                seu.s@meta.data[[pseudobulk_source]] = gsub(" ", "_", seu.s@meta.data[[pseudobulk_source]])
                seu.s@meta.data[[pseudobulk_source]] = gsub("-", "_", seu.s@meta.data[[pseudobulk_source]])
                
                Seurat::Idents(seu.s) =  pseudobulk_source
                
                if (doDimPlots & subsampleCur == 1) {
                    pdfFile = paste0(outputDir, "/dimPlot_combined_", pseudobulk_source, ".pdf")
                    .doDimPlots(seu.s, groupBy = pseudobulk_source, pdfFile = pdfFile)
                }
                
                nClusters = length(which(table(Seurat::Idents(seu.s)) %>% as.vector() >= minCellsPerCluster))
                futile.logger::flog.info(paste0(" Number of clusters with at least ", minCellsPerCluster, " cells per cluster: ", nClusters))
                
                # Subsampling
                seu.s = .subsampleSeuratObject(seu.s, subsample_percentage = subsample_percentage)
                if (subsample_percentage != 100) {
                    nClusters = length(which(table(Seurat::Idents(seu.s)) %>% as.vector() >= minCellsPerCluster))
                    futile.logger::flog.info(paste0(" Number of clusters with at least ", minCellsPerCluster, " cells per cluster after subsampling: ", nClusters))
                }
                
                # Metadata
                file_metadata = paste0(outputDir, "/metadata_", pseudobulk_source, subsample_suffix, ".tsv.gz")
                metadata = .writeMetadata(seu.s, file_metadata, minCellsPerCluster = minCellsPerCluster)
                
                
                # RNA #
                futile.logger::flog.info(paste0(" Aggregate and prepare RNA counts for each cluster"))
                
                assayName_aggregation = dplyr::if_else(sumCounts, assayName_RNA, "SCT")
                
                rna.pseudobulk.clean = .aggregateCounts(seu.s, assayName = assayName_aggregation, slotName = slotName_aggregation_RNA,
                                                        groupBy = pseudobulk_source, ID_column = "gene", 
                                                        sumCounts = sumCounts)
                
                # Delete clusters that contain too few cells
                rna.pseudobulk.clean = dplyr::select(rna.pseudobulk.clean, "gene", dplyr::one_of(metadata$sampleName))
                
                # Merge with the actual features and their correct mappings.
                rna.pseudobulk.clean2 = .addEnsemblIDs(rna.pseudobulk.clean, mapping = features)
                
                futile.logger::flog.info(paste0(" Writing RNA counts to file ", file_in_rna))
                readr::write_tsv(rna.pseudobulk.clean2, file_in_rna)
                
                .printScarcity(rna.pseudobulk.clean2, type = "RNA")
                
                
                # ATAC #
                futile.logger::flog.info(paste0(" Aggregate and prepare ATAC counts for each cluster"))
                
                atac.pseudobulk.clean = .aggregateCounts(seu.s,  assayName = assayName_ATAC, slotName = slotName_aggregation_ATAC,
                                                         groupBy = pseudobulk_source, ID_column = "peakID",
                                                         sumCounts = sumCounts)
                
                # Delete clusters that contain too few cells
                atac.pseudobulk.clean = dplyr::select(atac.pseudobulk.clean, "peakID", dplyr::one_of(metadata$sampleName))
                
                
                # Replace the first hyphen with a colon
                atac.pseudobulk.clean$peakID = sub("-", ":", atac.pseudobulk.clean$peakID)
                
                futile.logger::flog.info(paste0(" Writing ATAC counts to file ", file_in_atac))
                readr::write_tsv(atac.pseudobulk.clean, file_in_atac)
                
                .printScarcity(atac.pseudobulk.clean, type = "ATAC")
                
            }
            
        }
        
        
    } # end for all subsamples

    
    
    
    GRaNIE:::.printExecutionTime(start, prefix = "") 
    futile.logger::flog.info(paste0("Finished successfully, all output files have been saved in ", outputDir, ". Returning the Seurat object. Happy GRaNIE'ing!"))
    
    # Return object before it was potentially subsampled
    seu.orig.s
    
}

#' @import ggplot2
.doDimPlots <- function(seu.s, groupBy = "seurat_clusters", pdfFile) {
    
    
    if (!is.null(pdfFile)) {
        futile.logger::flog.info(paste0(" Writing DimPlot to file ", pdfFile))
        grDevices::pdf(pdfFile, width = 15, height = 10)
    }
    
    p1 <- Seurat::DimPlot(seu.s, reduction = "umap.rna", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA")
    p2 <- Seurat::DimPlot(seu.s, reduction = "umap.atac", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
    p3 <- Seurat::DimPlot(seu.s, reduction = "wnn.umap", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("WNN")
    plot(p1 + p2 + p3 & Seurat::NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
    
    if (!is.null(pdfFile)) {
        grDevices::dev.off()
    }
    
    
}

.removeClusters <- function(seu.s, minCellsPerCluster = 25) {
    
    # TODO
    stop("Not implemented yet")
    metadata = tibble::tibble(sampleName = levels(seu.s), 
                              nCells = table(Seurat::Idents(seu.s)) %>% as.vector())
    
    clustersToRemove = metadata %>% dplyr::filter(nCells < minCellsPerCluster) %>% dplyr::pull(sampleName)
    clustersToKeep   = metadata %>% dplyr::filter(nCells >= minCellsPerCluster) %>% dplyr::pull(sampleName)
    
    cellsRemove = names(Seurat::Idents(seu.s))[which(Seurat::Idents(seu.s) %in% clustersToRemove)]
    
    seu.s[,!colnames(seu.s) %in% toRemove]
    seu.s[,colnames(seu.s) %in% cellsRemove]
    
    seu2.s<- subset(seu.s, idents = clustersToKeep)
    
    
    
    seu.s
}

.writeMetadata <- function(seu.s, file_metadata, minCellsPerCluster = 25) {
    
    nClusters = levels(seu.s)
    futile.logger::flog.info(paste0(" Number of clusters found before filtering: ",  length(nClusters)))
    
    futile.logger::flog.info(paste0(" Cells per cluster before filtering: "), table(Seurat::Idents(seu.s)), capture = TRUE)
    
    metadata = tibble::tibble(sampleName = paste0("cluster", levels(seu.s)), 
                              nCells = table(Seurat::Idents(seu.s)) %>% as.vector()) %>%
        dplyr::filter(nCells >= minCellsPerCluster)
    
    
    futile.logger::flog.info(paste0(" Writing metadata after filtering  to file ", file_metadata))
    
    if (!dir.exists(dirname(file_metadata))) {
        dir.create(dirname(file_metadata), recursive = TRUE)
    }
    
    readr::write_tsv(metadata, file_metadata)
    
    metadata
}

.aggregateCounts <- function(seu.s, assayName, groupBy = "ident", sumCounts = TRUE, slotName = "counts", 
                             ID_column = "peakID") {
    
    checkmate::assertSubset(slotName, c("scale.data", "counts", "data"))
    
    
    
    
    if (sumCounts) {
        
        # If slot is set to 'data', this function assumes that the data has been log normalized and 
        # therefore feature values are exponentiated prior to aggregating so that sum is done in non-log space. 
        # Otherwise, if slot is set to either 'counts' or 'scale.data', no exponentiation is performed prior to 
        # aggregating
        
        
        pseudobulk = Seurat::AggregateExpression(
            seu.s,
            assays = assayName,
            features = NULL,
            return.seurat = FALSE,
            group.by = groupBy,
            add.ident = NULL,
            slot =  slotName,
        ) 
        
    } else {
        
        # When aggregating by mean, we should use the normalized counts as basis
        # and NOT the raw counts. See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9487674/ for details
        # The aggregation of the count values for pseudobulk methods can be performed using two approaches: 
        # cumulative summing of raw count values (sum) or averaging single-cell-normalized count values (mean). The sum aggregation is followed by bulk normalization, and it has achieved better performances in earlier studies than the mean aggregation [6].
        
        # seu.s[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA. 
        # The ‘corrected’ UMI counts are stored in seu.s[["SCT"]]@counts
        # We store log-normalized versions of these corrected counts in seu.s[["SCT"]]@data
        
        # assayRNA = "SCT"
        
        # The residuals for this model are normalized values, and can be positive or negative. 
        # Positive residuals for a given gene in a given cell indicate that we observed more UMIs than expected given the gene’s average expression 
        # in the population and cellular sequencing depth, while negative residuals indicate the converse.
        
        # Take the mean of the normalized values per cluster
        pseudobulk = Seurat::AverageExpression(
            seu.s,
            assays = assayName,
            features = NULL,
            return.seurat = FALSE,
            group.by = groupBy,
            add.ident = NULL,
            slot = slotName # If slot is set to 'data', this function assumes that the data has been log normalized and therefore feature values are exponentiated prior to aggregating so that sum is done in non-log space. Otherwise, if slot is set to either 'counts' or 'scale.data', no exponentiation is performed prior to aggregating
        ) 
        
    }
    
    # atac.before = GetAssayData(object = seu.s, assay = "ATAC", slot = "counts")
    # atac.before2 = GetAssayData(object = seu.s, assay = "ATAC", slot = "data")
    # atac.before3 = GetAssayData(object = seu.s, assay = "ATAC", slot = "scale.data")
    # 
    # # Sanity check
    # for (i in 1:10) {
    #   stopifnot(identical(sum(atac.before[i,]), sum(atac.pseudobulk$ATAC[i,]))) 
    # }
    
    # Prepare data for GRN format
    pseudobulk[[assayName]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(ID_column) %>%
        tibble::as_tibble() 
    
    
}

.subsampleSeuratObject <- function(seu.s, subsample_percentage) {
    
    if (subsample_percentage == 100) {
        return(seu.s)
    }
    
    cellBarcodes = colnames(seu.s)
    
    if (subsample_percentage < 100) {
        
        nColFinal = floor(ncol(seu.s) * (subsample_percentage / 100))
        stopifnot(nColFinal <= ncol(seu.s))
        futile.logger::flog.info(paste0(" Subsampling counts to ", subsample_percentage, " % before count aggregeation"))
        
        # This seems necessary but WHY?
        # https://stackoverflow.com/questions/52407517/random-number-generation-in-r-using-sample-always-returning-the-same-resul
        rm(.Random.seed, envir=globalenv())
        #rm(.Random.seed)
        cellBarcodes = sample(colnames(seu.s), size = nColFinal, replace = FALSE)
        
        
    } else if (subsample_percentage > 100) {
        nColFinal = ncol(seu.s) * (subsample_percentage / 100)
        stopifnot(nColFinal > ncol(seu.s))
        futile.logger::flog.info(paste0(" Subsampling counts to ", subsample_percentage, " % before count aggregeation"))
        cellBarcodes = sample(colnames(seu.s), size = nColFinal, replace = TRUE)
    } 
    
    seu.s[, cellBarcodes]
}



.addEnsemblIDs <- function(rna.pseudobulk.clean, mapping) {
    
    futile.logger::flog.info(paste0(" Mapping gene names to Ensembl IDs using the provided features file..."))
    rna.pseudobulk.clean = dplyr::left_join(rna.pseudobulk.clean, mapping, by = "gene")
    
    rna.pseudobulk.clean.filt =  dplyr::filter(rna.pseudobulk.clean, !is.na(.data$ens_id))
    nRows_missing = nrow(rna.pseudobulk.clean) - nrow(rna.pseudobulk.clean.filt)
    genesNA = rna.pseudobulk.clean$gene[which(!rna.pseudobulk.clean$gene %in% rna.pseudobulk.clean.filt$gene)]
    futile.logger::flog.info(paste0(" Deleted the following ", nRows_missing, " rows because no Ensembl ID could be found: ", paste0(genesNA, collapse = ", ")))
    
    # Check ambiguities for gene names, some gene names may be ambiguous
    geneNameFreq = table(rna.pseudobulk.clean.filt$gene)
    genesDelete = names(which(geneNameFreq > 1))
    
    futile.logger::flog.info(paste0(" Deleted the following ", length(genesDelete), " genes because the gene name was not unique and appeared multiple times: ", paste0(genesDelete, collapse = ",")))
    
    rna.pseudobulk.clean = dplyr::filter(rna.pseudobulk.clean.filt, !.data$gene %in% genesDelete) %>%
        dplyr::select("ens_id", tidyselect::everything(), -"view", -"chr", -"start", -"end", -"gene") %>%
        dplyr::rename(ENSEMBL = "ens_id")
    
    rna.pseudobulk.clean
    
}

.printScarcity <- function(pseudobulk, type = "RNA") {
    
    checkmate::assertSubset(type, c("RNA", "ATAC"))
    colname = dplyr::if_else(type == "RNA", "ENSEMBL", "peakID")
    
    pseudobulk.m = pseudobulk %>%
        tibble::column_to_rownames(colname) %>%
        as.matrix()
    
    spasity_fraction = (length(pseudobulk.m) - Matrix::nnzero(pseudobulk.m)) / length(pseudobulk.m)
    
    futile.logger::flog.info(paste0(" Sparsity ", type, ": ",  round(spasity_fraction * 100,1), " %"))
}
# 
# prepareSeuratDataTrajectory_GRaNIE <- function(seu.s, binsPerTrajectory) {
#   
#   # Check trajectory and pseudotime columns in metadata
#   checkmate::assertSubset(c("trajectory", "pseudotime"), colnames(seu.s@meta.data))
#   
#   # Add the bins
#   
#   for (traj in unique(seu.s@meta.data$trajectory)) {
#     
#     # Subset to cells only from the chosen trajectory / partition
#     subset.seu.s <- subset(seu.s, subset = cell_type %in% traj_celltypes[[traj]])
#     
#     
#     subset.seu.s@meta.data <- subset.seu.s@meta.data %>%
#       mutate(
#         bin_id = ggplot2::cut_number(!!!syms(pseudotime), n = binsPerTrajectory, labels=FALSE),
#         bin_interval = ggplot2::cut_number(!!!syms(pseudotime), n = binsPerTrajectory),
#         traj_bin = paste(bin_id, traj, sep = "__"),
#         pseudo_bin_id = paste(!!!syms(pseudo_source), bin_id, traj, sep = "__")
#       )
#     
#     bins <- subset.seu.s@meta.data %>% distinct(bin_id) %>% pull() %>% sort()
#     
#     # Prepare Seurat object for each bin and then run GRaNIE on it
#     for(i in 1:(length(bins)-1)) {
#       splitOverlap_Seurat(subset.seu.s, 
#                           bins[i], 
#                           bins[i+1], 
#                           traj, 
#                           pseudo_source, 
#                           proj_folder, 
#                           file_RNA_features = file_RNA_features_file,
#                           min_threshold = min_threshold,
#                           forceRerun = forceRerun)
#       
#       computeOverlap_GRaNIE(bins[i], 
#                             bins[i+1], 
#                             traj, 
#                             pseudo_source = pseudo_source, 
#                             dataset_name = paste(dataset_name, bins[i], bins[i+1], pseudo_source, sep = '__'), 
#                             proj_folder = proj_folder,
#                             forceRerun = forceRerun,
#                             TF_peak.fdr.threshold = TF_peak.fdr.threshold,
#                             peak_gene.fdr.threshold = peak_gene.fdr.threshold)
#     }
#   
#   print('Adding bin info to Seurat object')
#   
#   seu.s <- add_bininfo_Seurat(seu.s, trajectory_list, ident, 'pseudotime_order_RNA')
#   
#   # Run prepareSeuratData_GRaNIE for each trajectory + bin variant
#   # run GRaNIE_batch for each of them
# }


