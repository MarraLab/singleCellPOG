# libraries
library(tidyverse)
library(data.table)
library(Matrix)
library(Seurat)
library(Signac)
library(scater)
library(scRNAseq)
library(scran)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(biomaRt)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(scDblFinder)
library(infercnv)
library(harmony)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(SeuratWrappers)
library(SeuratExtend)
library(CellChat)
library(viridis)
library(slingshot)
library(glmmTMB)
library(emmeans)
library(rstatix)
library(httpgd)
library(qs2)
library(ggalluvial)
library(ggsankey)
library(circlize)
library(RColorBrewer)
library(paletteer)
library(glmmTMB)
library(org.Hs.eg.db)
library(annotatr)
library(aricode)  

path <- "PATH"
files <- list.files(path)
patient_ids <- paste0("POG", c("003", "130", "217", "415", "590", "643", "1128", "1329", "147", "196", "318", "326", "609", "650", "716", "732", "785"))

# converting between ATAC and RNA barcodes
atac_barcodes <- fread("PATH/cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz",
                       header = F) %>%
  mutate(V1 = paste(V1, "1", sep = "-"))
rna_barcodes <- fread("PATH/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz",
                      header = F) %>%
  mutate(V1 = paste(V1, "1", sep = "-"))

convert_rna_indices <- function(rna_bc){
  rna_indices <- match(rna_bc, rna_barcodes$V1)
  valid_atac_bc <- atac_barcodes$V1[rna_indices]
  return(valid_atac_bc)
}
convert_atac_indices <- function(atac_bc){
  atac_indices <- match(atac_bc, atac_barcodes$V1)
  valid_rna_bc <- rna_barcodes$V1[atac_indices]
  return(valid_rna_bc)
}

# create count matrices from 10X files
read_10x_rna <- function(file_path){
  all_files <- list.files(file_path)
  
  barcodes <- fread(paste0(file_path, "/", all_files[grepl("barcodes", all_files)]), header = FALSE)
  mat <- readMM(paste0(file_path, "/", all_files[grepl("mtx", all_files)]))
  colnames(mat) <- barcodes$V1
  
  features <- fread(paste0(file_path, "/", "features.tsv.gz"), header = FALSE)
  rownames(mat) <- features$V2
  
  return(mat)
}
read_10x_atac <- function(file_path){
  all_files <- list.files(file_path)
  
  barcodes <- fread(paste0(file_path, "/", all_files[grepl("barcodes", all_files)]), header = FALSE)
  mat <- readMM(paste0(file_path, "/", all_files[grepl("mtx", all_files)]))
  colnames(mat) <- barcodes$V1
  
  features <- fread(paste0(file_path, "/", "peaks.bed"), header = FALSE) %>%
    tidyr::unite(col = "V2", V1:V3, sep = ".", remove = T)
  rownames(mat) <- features$V2
  
  return(mat)
}

qc_filter <- function(sce, mito_list){
  metrics <- perCellQCMetrics(sce, subsets=list(Mito=mito_list))
  metrics_df <- as.data.frame(metrics) %>%
    rownames_to_column(var = "barcode")
  
  # filtering for cells
  lib <- isOutlier(metrics_df$sum, log=TRUE, type="lower")
  genes <- isOutlier(metrics_df$detected, log=TRUE, type="lower")
  # mito <- metrics_df$subsets_Mito_percent > 10
  mito <- metrics_df$subsets_Mito_percent > median(metrics_df$subsets_Mito_percent) + 3 * sd(metrics_df$subsets_Mito_percent)
  
  discard <- lib | genes | mito
  metrics_df$discard <- discard
  
  sce <- sce[,!discard]
  
  # removing doublets
  sce <- scDblFinder(sce)
  doublet_df <- data.frame("barcode" = colnames(sce),
                           "doublet" = sce@colData@listData[["scDblFinder.class"]])
  
  singletIndices <- which(sce$scDblFinder.class == "singlet")
  sce <- sce[,singletIndices]
  
  metrics_df <- metrics_df %>%
    left_join(doublet_df, by = "barcode")
  return(list(sce, metrics_df))
}

# seurat workflows
seurat_rna_pipeline <- function(seu){
  DefaultAssay(seu) <- "RNA"
  
  seu %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:10, reduction.name = "rna.umap")
}

seurat_atac_pipeline <- function(seu){
  DefaultAssay(seu) <- "PEAKS"
  
  seu %>%
    # NucleosomeSignal() %>%
    # TSSEnrichment(fast = FALSE) %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(reduction = "lsi", dims = 2:10, reduction.name = "peaks.umap") # edited for processing peak calls
}

make_chrom_assay_atac <- function(df_atac, fragments_path, annotations){
  # filtering for only peaks from standard chromosomes
  grange.counts <- StringToGRanges(rownames(df_atac), sep = c("\\.", "\\."))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- df_atac[as.vector(grange.use), ]
  
  # using fragment files for creating an assay object
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c("\\.", "\\."),
    genome = 'hg38',
    fragments = fragments_path,
    annotation = annotations
  )
  return(chrom_assay)
}

runInferCNV <- function(seu, sampleAnnot, save_dir, ref = NULL){
  counts_matrix <- seu[["RNA"]]$counts
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                      annotations_file = sampleAnnot,
                                      delim = "\t",
                                      gene_order_file = "PATH/Data/gencode_v21_gen_pos.complete.txt",
                                      ref_group_names = ref,
                                      chr_exclude = c("chrY", "chrM"))
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir = save_dir,  # dir is auto-created for storing outputs
                               # analysis_mode = 'subclusters',
                               cluster_by_groups = F,
                               cluster_references = F,
                               denoise = F,
                               HMM = F,
                               output_format = NA,
                               num_threads = 10,
                               no_prelim_plot = TRUE,
                               no_plot = TRUE)
}

# helper function for calculating SNV scores
snv_count <- function(cbsniffer_path){
  mut_calls <- fread(cbsniffer_path)
  mut_bc <- mut_calls %>%
    dplyr::filter(alt_count == 1) %>%
    group_by(barcode) %>%
    summarize(mut_count = n()) %>% 
    mutate(rna_bc = convert_atac_indices(barcode)) %>%
    dplyr::select(rna_bc, mut_count)
  
  return(mut_bc)
}