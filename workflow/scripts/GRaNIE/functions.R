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
