source("../../scFunctions.R")

args_file <- fread("../../Data/WholeCohort.csv")

# atac requirements
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(all_versions = FALSE)
)
annotations <- readRDS("../../Data/annotations_ATAC.RDS")

# cell type annotation
BoroniReference <- readRDS("PATH/boroniCellTypeAnnotRef.RDS")

for(i in seq(1, args_file[-10])){
  mode <- args_file[[i, 2]]
  dir_path <- args_file[[i, 3]]
  folder_name <- args_file[[i, 1]]
  rna_path <- args_file[[i, 4]]
  file_path_atac <- args_file[[i, 5]]
  reads_cutoff <- args_file[[i, 6]]
  genes_cutoff <- args_file[[i, 7]]
  mito_cutoff <- args_file[[i, 8]]
  fragments_path <- args_file[[i, 9]]

  save_dir <- paste0(path, folder_name)
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }

  # processing RNA
  count_matrix <- read_10x_rna(rna_path)
  dupInd <- which(duplicated(rownames(count_matrix)))
  count_matrix <- count_matrix[-dupInd,]
  sce <- SingleCellExperiment(assays = list(counts = count_matrix))

  # identifying mitochondrial genes
  genes_list <- fread(paste0(rna_path, "/", "features.tsv.gz"), header = FALSE)
  is.mito <- which(grepl("^MT-", genes_list$V2))

  # calculating QC metrics
  sceQC <- qc_filter(sce, mito_list = is.mito)
  write_tsv(sceQC[[2]], paste0(save_dir, "/QCMetrics_byCell.tsv"))

  seurat <- CreateSeuratObject(counts = counts(sceQC[[1]]))
  seurat <- seurat %>%
    SCTransform() %>%
    RunPCA()
  
  # clustering
  DefaultAssay(seurat) <- "RNA"
  seurat <- seurat_rna_pipeline(seurat)

  f <- folder_name
  k <- dplyr::filter(clusteringParams, folder_name == f) %>%
    pull(k)
  res <- dplyr::filter(clusteringParams, folder_name == f) %>%
    pull(res)
  seurat <- FindNeighbors(seurat, dims = 1:15, k.param = k) %>%
    FindClusters(resolution = res, graph.name = "RNA_snn")
  
  # cell type annotation
  anchors <- FindTransferAnchors(reference = BoroniReference,
                                 query = seurat,
                                 # normalization.method = "SCT",
                                 dims = 1:30,
                                 reference.reduction = "rpca")
  predictions <- TransferData(anchorset = anchors, refdata = BoroniReference$cellType, dims = 1:30)
  seurat <- AddMetaData(seurat, predictions)

  saveRDS(seurat, paste0(path, folder_name, "/rna_ProcessedSeurat.RDS"))

  # processing ATAC
  df_atac <- read_10x_atac(file_path_atac) %>%
    as.matrix() %>%
    as.data.frame()

  chrom_assay <- make_chrom_assay_atac(df_atac, fragments_path, annotations)

  # convert barcodes
  bc <- convert_atac_indices(colnames(chrom_assay))
  chrom_assay <- RenameCells(chrom_assay, new.names = bc)

  # memory constraints > cannot use raw ATAC output
  # filtering for cells that are only present in both RNA and ATAC assays
  valid_cells <- intersect(colnames(chrom_assay), colnames(seurat))

  chrom_assay <- subset(chrom_assay, cells = valid_cells)
  seurat <- subset(seurat, cells = valid_cells)

  seurat_atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
  seurat_atac[["RNA"]] <- seurat[["RNA"]]
  seurat_atac <- AddMetaData(seurat_atac, seurat[[]])
  saveRDS(seurat_atac, paste0(path, folder_name, "/separate_ProcessedSeurat.RDS"))

  # run MACS2
  # calling cluster-specific peaks
  peaks <- CallPeaks(
    object = seurat_atac,
    assay = "ATAC",
    group.by = "seurat_clusters",
    macs2.path = "~/anaconda3/bin/macs2"
  )

  keep.peaks <- as.logical(seqnames(granges(peaks)) %in% main.chroms)
  peaks <- peaks[keep.peaks]

  peaks_assay <- FeatureMatrix(
    fragments = seurat_atac@assays[["ATAC"]]@fragments,
    features = peaks
  )
  peaks_assay_filtered <- peaks_assay[,which(colnames(peaks_assay) %in% colnames(seurat_atac))]

  seurat_atac[["PEAKS"]] <- CreateChromatinAssay(
    counts = peaks_assay_filtered,
    # sep = c(":", "-"),
    fragments = seurat_atac@assays[["ATAC"]]@fragments,
    annotation = annotations
  )
  DefaultAssay(seurat_atac) <- "PEAKS"
  seurat_atac <- seurat_atac_pipeline(seurat_atac)

  # adding motif info
  seurat_atac <- AddMotifs(
    object = seurat_atac,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm,
    assay = "PEAKS"
  ) %>%
    RunChromVAR(
      genome = BSgenome.Hsapiens.UCSC.hg38,
      assay = "PEAKS"
    )

  saveRDS(seurat_atac, paste0(path, folder_name, "/separate_ProcessedSeurat.RDS"))
}