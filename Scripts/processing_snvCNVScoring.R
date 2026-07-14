source("scFunctions_forSubmission.R")

for(f in files){

  # running inferCNV
  seu <- readRDS(paste0(path, f, "/rna_ProcessedSeurat.RDS"))
  DefaultAssay(seu) <- "RNA"

  # annotations
  sample_annot <- data.frame(sample_id = colnames(seu)) %>%
    mutate(label = ifelse(seu$predicted.id %in% c("B_cells", "Myeloid", "Mast_cells", "NK/T_cells"),
                          "Reference", "Observation"))
  write_tsv(sample_annot,
            paste0(path, f, "/rna_SampleAnnot.tsv"),
            col_names = F)

  runInferCNV(seu, paste0(path, f, "/rna_SampleAnnot.tsv"), paste0(path, f, "/rna_inferCNV"), ref = "Reference")
}

# calculating snv scores for individual samples
args_file <- fread("PATH/WholeCohort.csv")
# This file is not included for privacy/security. The relevant columns are:
# "folder_name": Sample ID    
# "file_path": CellRanger outputs for snRNA-seq    
# "file_path_atac": CellRanger outputs for snATAC-seq 
# "fragments_path": snATAC-seq fragments file 
# "pog_path": SNVs from bulk WGS        
# "bam_path": Raw BAM file for snATAC-seq 

for(i in seq(1, nrow(args_file))){
  sample <- args_file$folder_name[i]
  
  atac_bam_path <- args_file$bam_path[i]
  file_path <- args_file$file_path[i]
  rna_bam_path <- gsub("filtered_feature_bc_matrix",
                       "possorted_genome_bam.bam",
                       file_path)

  # calculating atac-based score
  fragments_path <- args_file$fragments_path[i]

  save_dir <- paste0(path, sample)
  seu <- readRDS(paste0(save_dir, "/rna_ProcessedSeurat.RDS"))

  # run cbSniffer
  data.frame(barcodes = convert_rna_indices(colnames(seu))) %>%
    write_tsv(paste0(save_dir, "/atacBarcodes.tsv"))
  
  cb_string <- paste("python3 sander_mutationcalling_cb_sniffer.py",
                     bam_path,
                     paste0(save_dir, "/self_POGSNVcalls.tsv"), # SNVs from bulk
                     paste0(save_dir, "/atacBarcodes.tsv"),
                     paste0(save_dir, "/self_cbSniffer"),
                     "-bq 30") 
  system(cb_string)

  fragment coverage
  pogSNVs_path <- paste0(path, sample, "/self_POGSNVcalls.tsv")
  
  fragments <- fread(fragments_path)
  colnames(fragments) <- c("V1", "V2", "V3", "V4", "V5")
  fragments_GR <- GRanges(seqnames = fragments$V1,
                          ranges = IRanges(start = fragments$V2, end = fragments$V3),
                          atac_barcode = fragments$V4)
  pogSNVs <- fread(pogSNVs_path)
  pogSNVs_GR <- GRanges(seqnames = pogSNVs$chrm,
                        ranges = IRanges(start = pogSNVs$start, end = pogSNVs$stop),
                        SNV_ID = pogSNVs$gene_name)
  
  overlap <- findOverlaps(fragments_GR, pogSNVs_GR)
  intersectDF <- cbind(as.data.frame(fragments_GR[overlap@from]), as.data.frame(pogSNVs_GR[overlap@to]))
  colnames(intersectDF) <- c("chrFrag", "startFrag", "endFrag", "widthFrag", "strandFrag", "atac_barcode", "chrSNV", "startSNV", "endSNV", "widthSNV", "strandSNV", "SNV_ID")
  write_csv(intersectDF, paste0(path, sample, "/perCell_SNVFragmentCoverage.csv"))

  # calculate and plot SNV scores
  frags_perCell <- fread(paste0(path, sample, "/perCell_SNVFragmentCoverage.csv")) %>%
    mutate(rna_bc = convert_atac_indices(atac_barcode)) %>%
    group_by(rna_bc) %>%
    summarize(numFragsOverlapping = n()) %>%
    dplyr::filter(numFragsOverlapping != 0)

  cbsniffer_path <- paste0("PATH/",
                           sample,
                           "/cbSnifferOuts/allCells_cbATAC_counts_CB.tsv")
  SNVs_perCell <- snv_count(cbsniffer_path)
  SNVs_perCell2 <- left_join(data.frame(rna_bc = colnames(seu)),
                             SNVs_perCell)
  snvScore_byCell <- left_join(SNVs_perCell2, frags_perCell, by = c("rna_bc"))

  snvScore_byCell[is.na(snvScore_byCell)] <- 0
  snvScore_byCell <- snvScore_byCell %>%
    mutate(denominator = numFragsOverlapping,
           mutscore = mut_count / denominator)
  snvScore_byCell[is.na(snvScore_byCell)] <- 0
  
  snvScore_byCell <- fread(cbsniffer_path) %>%
    group_by(barcode) %>%
    summarize(snvScore = sum(alt_count) / sum(ref_count)) %>%
    dplyr::filter(!is.infinite(snvScore)) %>%
    mutate(barcode = convert_atac_indices(barcode))
  
  write_csv(snvScore_byCell,
            paste0("PATH/",
                   sample,
                   "/cbSnifferOuts/allCells_atacScore_fragments.tsv"))
  
  metadata <- snvScore_byCell %>%
    column_to_rownames(var = "barcode")
  seu <- AddMetaData(seu, metadata)
  seu$snvScore[is.na(seu$snvScore)] <- 0
  
  p1 <- DimPlot(seu, reduction = "rna.umap", group.by = "seurat_clusters") +
    ggtitle("RNA-based clusters")
  ggsave(paste0(path, sample, "/rnaClusters.png"), plot = p1, width = 6, height = 5)
  
  p2 <- FeaturePlot(seu, reduction = "rna.umap", features = "snvScore") +
    ggtitle("scATAC-based SNV score")
  ggsave(paste0(path, sample, "/atacSNVScore.png"), plot = p2, width = 6, height = 5)
}

# CNVs across multiple biopsies
# identify CNVs called by ploidetect in later biopsies that are not present in earlier ones
genes <- read_tsv("PATH/Data/gencode_v21_gen_pos.complete.txt", 
                  col_names = c("gene_symbol", "chr", "start", "end")) %>%
          mutate(chr = gsub("chr", "", chr))
gr_genes <- GRanges(seqnames = genes$chr,
                      ranges = IRanges(start = genes$start, end = genes$end),
                      gene_symbol = genes$gene_symbol)

compare_cnv_with_genes <- function(file1, file2) { 
  # Read CNV files
  cnv1 <- read_tsv(paste0(file1, "/cna.txt"), col_types = cols())
  cnv2 <- read_tsv(paste0(file2, "/cna.txt"), col_types = cols())
  
  # Create GRanges objects for CNVs
  gr1 <- GRanges(seqnames = cnv1$chr,
                 ranges = IRanges(start = cnv1$pos, end = cnv1$end), strand = "*")
  
  gr2 <- GRanges(seqnames = cnv2$chr,
                 ranges = IRanges(start = cnv2$pos, end = cnv2$end), strand = "*")
  
  # Find non-overlapping CNVs in sample2
  overlaps <- findOverlaps(gr2, gr1)
  # non_overlapping_indices <- setdiff(seq_along(gr2), queryHits(overlaps))
  # non_overlapping_cnvs <- cnv2[non_overlapping_indices, ]
  # non_overlapping_cnvs$size <- non_overlapping_cnvs$end - non_overlapping_cnvs$pos
  # filtered_cnvs <- non_overlapping_cnvs[non_overlapping_cnvs$size > threshold, ]
  

  # Compare CN values for overlapping regions
  cnv_diff <- data.frame(
    sample2_index = queryHits(overlaps),
    sample1_index = subjectHits(overlaps),
    chr = as.character(seqnames(gr2[queryHits(overlaps)])),
    start = start(gr2[queryHits(overlaps)]),
    end = end(gr2[queryHits(overlaps)]),
    CN_sample2 = cnv2$CN[queryHits(overlaps)],
    CN_sample1 = cnv1$CN[subjectHits(overlaps)]
  ) %>%
    mutate(size = end - start,
           CN_diff = CN_sample2 - CN_sample1,
           change = case_when(
             CN_diff > 0 ~ "gain",
             CN_diff < 0 ~ "loss",
             TRUE ~ "no_change"
           )) %>%
    dplyr::filter(change != "no_change")


  # Create GRanges for filtered CNVs
  # gr_filtered <- GRanges(seqnames = filtered_cnvs$chr,
  #                        ranges = IRanges(start = filtered_cnvs$pos, end = filtered_cnvs$end))
  
  # Create GRanges for filtered CNV differences
  gr_diff <- GRanges(seqnames = cnv_diff$chr,
                     ranges = IRanges(start = cnv_diff$start, end = cnv_diff$end))

  # Find genes fully contained in CNV regions
  gene_hits <- findOverlaps(gr_genes, gr_diff)
  matched_genes <- genes[queryHits(gene_hits), ]
  matched_genes$cnv_index <- subjectHits(gene_hits)
  

  # Combine CNV differences with gene annotations
  cnv_diff$cnv_index <- seq_len(nrow(cnv_diff))
  result <- left_join(cnv_diff, matched_genes, by = "cnv_index") %>%
    dplyr::select(-cnv_index)
  
  return(result)
}

cnvComp_list <- list()
for(i in seq(1, length(multiBiop_list))){
  mb <- names(multiBiop_list)[i]
  mbs <- files[grepl(mb, files)]

  comb <- combn(mbs, 2, simplify = FALSE)
  for(i in seq(1, length(comb))){
    file1 <- wholeCohort[[which(wholeCohort$folder_name == comb[[i]][1]), "pd_path"]]
    file2 <- wholeCohort[[which(wholeCohort$folder_name == comb[[i]][2]), "pd_path"]]

    index <- paste(comb[[i]][1], comb[[i]][2], sep = ".")

    cnvComp_list[[index]] <- compare_cnv_with_genes(file1, file2)
  }
}

# for CNVs also detected in sc second biop, determine what percentage of first biop malig cells also contain it
# Function to calculate directional CNV percentage
directional_cnv_percent <- function(matrix, gene, direction) {
  if (!(gene %in% colnames(matrix))) return(NA)
  values <- matrix[,gene]
  neutral <- median(matrix, na.rm = T)
  if (direction == "gain") {
    return(mean(values > neutral, na.rm = TRUE) * 100)
  } else if (direction == "loss") {
    return(mean(values < neutral, na.rm = TRUE) * 100)
  } else {
    return(NA)
  }
}

compare_cnv_in_cells_directional <- function(cnv_list, cnv_matrix) {
  # Filter CNV list for valid genes and directions
  cnv_list <- cnv_list %>%
    dplyr::filter(!is.na(gene_symbol), change %in% c("gain", "loss")) %>%
    dplyr::filter(gene_symbol %in% colnames(cnv_matrix),
                  !duplicated(gene_symbol))
  
  # Extract biopsy labels from cell names
  biopsy_labels <- sapply(rownames(cnv_matrix), function(x) strsplit(x, "_")[[1]][3])
  biopsy1_cells <- cnv_matrix[biopsy_labels == min(biopsy_labels), , drop = FALSE]
  biopsy2_cells <- cnv_matrix[biopsy_labels == max(biopsy_labels), , drop = FALSE]
  
  # Neutral baseline (median CN across all cells)
  neutral <- median(as.matrix(cnv_matrix), na.rm = TRUE)
  
  # Match CNV genes to matrix columns
  gene_idx <- match(cnv_list$gene_symbol, colnames(cnv_matrix))
  
  # Vectorized computation
  biopsy1_percent <- numeric(length(gene_idx))
  biopsy2_percent <- numeric(length(gene_idx))
  biopsy1_ploidy <- numeric(length(gene_idx))
  biopsy2_ploidy <- numeric(length(gene_idx))
  
  for (i in seq_along(gene_idx)) {
    idx <- gene_idx[i]
    if (is.na(idx)) {
      biopsy1_percent[i] <- NA
      biopsy2_percent[i] <- NA
      biopsy1_ploidy[i] <- NA
      biopsy2_ploidy[i] <- NA
      next
    }
    
    dir <- cnv_list$change[i]
    
    # Masks for gain/loss
    if (dir == "gain") {
      mask_b1 <- biopsy1_cells[, idx] > neutral
      mask_b2 <- biopsy2_cells[, idx] > neutral
    } else {
      mask_b1 <- biopsy1_cells[, idx] < neutral
      mask_b2 <- biopsy2_cells[, idx] < neutral
    }
    
    # Percent of cells
    biopsy1_percent[i] <- mean(mask_b1, na.rm = TRUE) * 100
    biopsy2_percent[i] <- mean(mask_b2, na.rm = TRUE) * 100
    
    # Average ploidy among cells with gain/loss
    biopsy1_ploidy[i] <- if (any(mask_b1, na.rm = TRUE)) mean(biopsy1_cells[mask_b1, idx], na.rm = TRUE) else NA
    biopsy2_ploidy[i] <- if (any(mask_b2, na.rm = TRUE)) mean(biopsy2_cells[mask_b2, idx], na.rm = TRUE) else NA
  }
  
  # Combine results
  result <- cnv_list %>%
    mutate(biopsy1_percent = biopsy1_percent,
           biopsy2_percent = biopsy2_percent,
           biopsy1_ploidy = biopsy1_ploidy,
           biopsy2_ploidy = biopsy2_ploidy)
  
  return(result)
}

cnvComp_propCells_list <- list()
for(i in seq(1, length(multiBiop_list))){
  mb <- names(multiBiop_list)[i]
  mbs <- files[grepl(mb, files)]

  seu <- multiBiop_list[[mb]]
  seu <- AddMetaData(seu, allSamples_merged[[]])
  seu_malig <- subset(seu, identity == "Malignant")

  maligCellNames <- paste(sapply(colnames(seu_malig), function(x) convert_rna_indices(str_split(x, "_")[[1]][1])),
                          seu_malig$sample,
                          sep = "_")


  cnvObj <- multiBiopCNV_list[[mb]]
  cnv_mat <- t(cnvObj@expr.data)

  comb <- combn(mbs, 2, simplify = FALSE)
  for(i in seq(1, length(comb))){
    index <- paste(comb[[i]][1], comb[[i]][2], sep = ".")
    cnv_list <- cnvComp_list[[index]]

    cells <- intersect(which(rownames(cnv_mat) %in% maligCellNames), 
                       which(grepl(paste(comb[[i]], collapse = "|"), rownames(cnv_mat))))
    cnv_matrix <- cnv_mat[cells,]
    cnvComp_propCells_list[[index]] <- compare_cnv_in_cells_directional(cnv_list, cnv_matrix) %>%
      mutate(pair = index)
  }
}
cnvComp_propCells <- Reduce(rbind, cnvComp_propCells_list)

# SNVs across multiple biopsies
# detection of SNVs from other samples
all_samples <- list.files(path)
multi_biop_samples <- unique(grep(paste(multi_biop,collapse="|"), all_samples, value=TRUE))

# data frame for determining self vs others
patientLabels <- data.frame(sample = all_samples) %>%
  separate(col = sample, into = c("patient", "biopsy"), remove = F)

get_unique_snvs <- function(a_snvs, b_snvs){
  # comparator; a_snvs is already a df
  b_snvs_df <- fread(b_snvs) %>%
    dplyr::filter(chrm %in% chr_noContigs) %>%
    mutate(SNV_ID = paste(chrm, start, var, sep = "_"))
  
  b_snvs_unique <- b_snvs_df %>%
    dplyr::filter(!(SNV_ID %in% a_snvs$SNV_ID)) %>%
    dplyr::select(-SNV_ID)
  
  return(list(
    nrow(b_snvs_unique), # used only for same patient, to establish baseline of snvs to input
    b_snvs_unique
  ))
}

for(mb in multi_biop_samples){
  sample_folder <- paste0(path, mb)
  bam_path <- args_file %>%
    dplyr::filter(folder_name == mb) %>%
    dplyr::select(bam_path) %>%
    pull()
  
  # get sample(s) with same patient
  p <- str_split(mb, "_")[[1]][1]
  b <- str_split(mb, "_")[[1]][2]
  
  same_pt <- patientLabels %>%
    dplyr::filter(patient == p, biopsy != b)
  
  a_snvs <- fread(paste0(sample_folder, "/self_POGSNVcalls.tsv")) %>%
    mutate(SNV_ID = paste(chrm, start, var, sep = "_"))
  
  # calculate how many SNVs to use as input
  if(nrow(same_pt) > 0){
    numSNVs_samePT <- c()
    for(i in seq(1, nrow(same_pt))){
      comp_folder <- paste0(path, same_pt$sample[i])
      unique_snvs <- get_unique_snvs(a_snvs, paste0(comp_folder, "/self_POGSNVcalls.tsv"))
      numSNVs_samePT <- c(numSNVs_samePT, unique_snvs[[1]])
    }
    numInput <- round(mean(numSNVs_samePT), digits = 0)
  } else{
    numInput <- 0
  }
  
  # iterate through all samples
  if(!dir.exists(paste0(sample_folder, "/SNVNullDist"))){
    system(paste0("mkdir ", sample_folder, "/SNVNullDist"))
  } 
  else{
    system(paste0("rm -r ", sample_folder, "/SNVNullDist"))
    system(paste0("mkdir ", sample_folder, "/SNVNullDist")) # !!!
  }

  s_folder <- paste0(sample_folder, "/SNVNullDist")
  
  if(file.exists(paste0(sample_folder, "/self_POGSNVcalls.tsv"))){
    sample_snvs <- fread(paste0(sample_folder, "/self_POGSNVcalls.tsv")) %>%
      mutate(SNV_ID = paste(chrm, start, var, sep = "_"))
    
    for(s in all_samples){
      comp_folder <- paste0(path, s)
      
      if(!file.exists(paste0(comp_folder, "/SNVNullDist_SNVCalls"))){
        system(paste0("mkdir ", comp_folder, "/SNVNullDist_SNVCalls"))
      } 
      else{
        system(paste0("rm -r ", comp_folder, "/SNVNullDist_SNVCalls"))
        system(paste0("mkdir ", comp_folder, "/SNVNullDist_SNVCalls"))
      }
      c_folder <- paste0(comp_folder, "/SNVNullDist_SNVCalls")
      
      if(file.exists(paste0(comp_folder, "/self_POGSNVcalls.tsv")) & mb != s){
        
        # identify SNVs unique to comparator
        comp_snvs_unique <- get_unique_snvs(sample_snvs, paste0(comp_folder, "/self_POGSNVcalls.tsv"))
        
        if(numInput == 0){
          numInput <- nrow(comp_snvs_unique[[2]])
        }

        # extract at most numInput
        if(comp_snvs_unique[[1]] > numInput){
          output <- slice_sample(comp_snvs_unique[[2]], n = numInput)
        } else{output <- comp_snvs_unique[[2]]}
        
        write_tsv(output,
                  paste0(c_folder, "/", mb, "_unique_POGSNVcalls.tsv"))
        
        # run cbSniffer
        cb_string <- paste("python3 sander_mutationcalling_cb_sniffer.py",
                           bam_path,
                           paste0(c_folder, "/", mb, "_unique_POGSNVcalls.tsv"),
                           paste0(sample_folder, "/"),
                           paste0(s_folder, "/", s, "_cbSniffer"),
                           # "-mq 10",
                           "-bq 30")
        system(cb_string)
      }else{next}
    }
  }else{next}
}

# scoring of SNVs detected from other samples
for(mb in multi_biop_samples){
  frag_file <- dplyr::filter(args_file, folder_name == mb) %>%
    pull(fragments_path)
  fragments <- fread(frag_file)
  colnames(fragments) <- c("V1", "V2", "V3", "V4", "V5")
  fragments_GR <- GRanges(seqnames = fragments$V1,
                          ranges = IRanges(start = fragments$V2, end = fragments$V3),
                          atac_barcode = fragments$V4)
  
  snvCalls_path <- paste0(path, mb, "/SNVNullDist/")
  clusterLabels <- fread(paste0(path,
                                mb,
                                "/barcodesClusters_labelled.csv"))
  
  for(i in seq(1, nrow(args_file))){
    f <- args_file$folder_name[i]
    snv_path <- args_file$pog_path[i]
    nullSNVs <- paste0(path, f, "/SNVNullDist_SNVCalls/", mb, "_unique_POGSNVcalls.tsv")
    
    pogSNVs <- fread(nullSNVs)
    pogSNVs_GR <- GRanges(seqnames = pogSNVs$chrm,
                          ranges = IRanges(start = pogSNVs$start, end = pogSNVs$stop),
                          SNV_ID = pogSNVs$gene_name)
    
    overlap <- findOverlaps(fragments_GR, pogSNVs_GR)
    intersectDF <- cbind(as.data.frame(fragments_GR[overlap@from]), as.data.frame(pogSNVs_GR[overlap@to]))
    colnames(intersectDF) <- c("chrFrag", "startFrag", "endFrag", "widthFrag", "strandFrag", "atac_barcode", "chrSNV", "startSNV", "endSNV", "widthSNV", "strandSNV", "SNV_ID")
    write_csv(intersectDF, paste0(nullDist, f, "_perCell_SNVFragmentCoverage.csv"))
    
    if(mb != f){
      intersectDF <- fread(paste0(snvCalls_path, f, "_perCell_SNVFragmentCoverage.csv"))
      cbsniffer_path <- paste0(snvCalls_path, f, "_cbSniffer_counts_CB.tsv")
    } else{
      intersectDF <- fread(paste0(path, f, "/perCell_SNVFragmentCoverage.csv"))
      cbsniffer_path <- paste0(path, f, "/self_cbSniffer_counts_CB.tsv")
    }

    frags_perCell <- intersectDF %>%
      mutate(rna_bc = convert_atac_indices(atac_barcode)) %>%
      group_by(rna_bc) %>%
      summarize(numFragsOverlapping = n())
    
    cb <- fread(cbsniffer_path) %>%
      mutate(mutation_id = paste(chrm, start, sep = "_")) %>%
      dplyr::filter(alt_count == 1) %>%
      mutate(rna_bc = convert_atac_indices(barcode)) %>%
      group_by(rna_bc) %>%
    
    SNVs_perCell <- snv_count(cbsniffer_path)
    
    snvScore_byCell <- left_join(clusterLabels, frags_perCell, by = c("barcode" = "rna_bc")) %>%
      left_join(SNVs_perCell, by = c("barcode" = "rna_bc")) %>%
      left_join(cb, by = c("barcode" = "rna_bc"))
    
    snvScore_byCell[is.na(snvScore_byCell)] <- 0
    snvScore_byCell <- snvScore_byCell %>%
      mutate(denominator = numFragsOverlapping,
             mutscore = mut_count / denominator) %>%
      dplyr::filter(mutscore <= 1)
    
    write_csv(snvScore_byCell, paste0(path, mb, "/snvScoreComp/", f, "_scores.csv"))
}}

# calculate difference in SNV scores between biopsies, minus background
patientLabels <- data.frame(sample = files) %>%
  separate(col = sample, into = c("patient", "biopsy"), remove = F)
snvScoresOut <- data.frame(barcode = character(0),
                           snvScore = numeric(0),
                           comparator = character(0),
                           category = character(0),
                           sample = character(0))

for(mb in multi_biop_samples){
  scoresPath <- paste0(path, mb, "/snvScoreComp/")
  
  # self
  selfSNVs <- fread(paste0(path, mb, "/snvScores_perCell.csv")) %>%
    right_join(clusterAnnot, by = c("rna_bc" = "barcode")) %>%
    mutate(sample = mb,
           comparator = mb,
           category = "self") %>%
    dplyr::select(rna_bc, mutscore, comparator, category, cluster, sample) %>%
    `colnames<-`(c("barcode", "snvScore", "comparator", "category", "cluster", "sample"))
  selfSNVs[is.na(selfSNVs)] <- 0
  snvScoresOut <- rbind(snvScoresOut, selfSNVs)
  
  # same pt
  p <- str_split(mb, "_")[[1]][1]
  b <- str_split(mb, "_")[[1]][2]
  
  same_pt <- patientLabels %>%
    dplyr::filter(patient == p, biopsy != b) %>%
    pull(sample)

  for(spt in same_pt){
    scores <- fread(paste0(path, mb, "/snvScoreComp/", spt, "_scores.csv")) %>%
      right_join(clusterAnnot, by = c("barcode")) %>%
      mutate(sample = mb,
             comparator = spt,
             category = "same patient") %>%
      dplyr::select(barcode, mutscore, comparator, category, cluster.x, sample) %>%
      `colnames<-`(c("barcode", "snvScore", "comparator", "category", "cluster", "sample"))
    scores[is.na(scores)] <- 0
    snvScoresOut <- rbind(snvScoresOut, scores)
  }
  
  # others
  others <- list.files(paste0(path, mb, "/snvScoreComp/"))
  same_pts <- paste(same_pt, collapse = "|")
  for(o in others){
    if(!grepl(same_pts, o)){
      scores <- fread(paste0(path, mb, "/snvScoreComp/", o)) %>%
        right_join(clusterAnnot, by = c("barcode")) %>%
        mutate(sample = mb,
               comparator = spt,
               category = "other") %>%
        dplyr::select(barcode, mutscore, comparator, category, cluster.x, sample) %>%
        `colnames<-`(c("barcode", "snvScore", "comparator", "category", "cluster", "sample"))
      scores[is.na(scores)] <- 0
     snvScoresOut <- rbind(snvScoresOut, scores)
   }
  }
}

snvScores_forSamePt_minusOther <- snvScoresOut %>%
  dplyr::select(sample, barcode, category, snvScore) %>%
  pivot_wider(names_from = category, values_from = snvScore, values_fn = mean)
snvScores_forSamePt_minusOther[is.na(snvScores_forSamePt_minusOther)] <- 0
snvScores_forSamePt_minusOther <- mutate(snvScores_forSamePt_minusOther,
                                         selfScore = self - other,
                                         samePtScore = same_pt - other)

# background SNVs in normal cells
# create a pool of SNVs from all samples
snvList <- list()
for(f in files){
    sample_folder <- paste0(path, f)
    snvList[[f]] <- fread(paste0(sample_folder, "/self_POGSNVcalls.tsv")) %>%
        mutate(SNV_ID = paste(chrm, start, var, sep = "_"),
               sample = f)
}
snvsDF <- Reduce(rbind, snvList)
write_tsv(snvsDF, "Data/snvs_allSamples.tsv")
snvsDF <- fread("/projects/marralab/cayan_prj/PrecisionMed/Data/snvs_allSamples.tsv")

# randomly select SNVs from other samples, and look for them in current sample
for(f in files){
    sample_folder <- paste0(path, f)
    bam_path <- args_file %>%
        dplyr::filter(folder_name == f) %>%
        dplyr::select(bam_path) %>%
        pull()
    
    a_snvs <- fread(paste0(sample_folder, "/self_POGSNVcalls.tsv")) %>%
        mutate(SNV_ID = paste(chrm, start, var, sep = "_"))
  
    # calculate how many SNVs to use as input
    numInput <- nrow(a_snvs)
  
    # randomly select from pool of all SNVs
    output <- slice_sample(dplyr::filter(snvsDF, !(gene_name %in% a_snvs$gene_name)), n = numInput)

    write_tsv(output,
              paste0(sample_folder, "/null_POGSNVcalls.tsv"))
        
    # run cbSniffer
    cb_string <- paste("python3 sander_mutationcalling_cb_sniffer.py",
                        bam_path,
                        paste0(sample_folder, "/null_POGSNVcalls.tsv"),
                        paste0(sample_folder, "/separate_barcodes.tsv"),
                        paste0(sample_folder, "/null_cbSniffer"),
                        # "-mq 10",
                        "-bq 30")
    system(cb_string)
}

# scoring of SNVs detected from other samples
for(f in files){
    frag_file <- dplyr::filter(args_file, folder_name == f) %>%
        pull(fragments_path)
    fragments <- fread(frag_file)
    colnames(fragments) <- c("V1", "V2", "V3", "V4", "V5")
    fragments_GR <- GRanges(seqnames = fragments$V1,
                            ranges = IRanges(start = fragments$V2, end = fragments$V3),
                            atac_barcode = fragments$V4)

    clusterLabels <- fread(paste0(path,
                                    f,
                                    "/barcodesClusters_labelled.csv"))
  
    nullSNVs <- paste0(path, f, "/null_POGSNVcalls.tsv")
    pogSNVs <- fread(nullSNVs)
    pogSNVs_GR <- GRanges(seqnames = pogSNVs$chrm,
                          ranges = IRanges(start = pogSNVs$start, end = pogSNVs$stop),
                          SNV_ID = pogSNVs$gene_name)
    
    overlap <- findOverlaps(fragments_GR, pogSNVs_GR)
    intersectDF <- cbind(as.data.frame(fragments_GR[overlap@from]), as.data.frame(pogSNVs_GR[overlap@to]))
    colnames(intersectDF) <- c("chrFrag", "startFrag", "endFrag", "widthFrag", "strandFrag", "atac_barcode", "chrSNV", "startSNV", "endSNV", "widthSNV", "strandSNV", "SNV_ID")
    write_csv(intersectDF, paste0(path, f, "/null_SNVFragmentCoverage.csv"))

    frags_perCell <- intersectDF %>%
      mutate(rna_bc = convert_atac_indices(atac_barcode)) %>%
      group_by(rna_bc) %>%
      summarize(numFragsOverlapping = n())
    
    cbsniffer_path <- paste0(path, f, "/null_cbSniffer_counts_CB.tsv")
    cb <- fread(cbsniffer_path) %>%
      mutate(mutation_id = paste(chrm, start, sep = "_")) %>%
      dplyr::filter(alt_count == 1) %>%
      mutate(rna_bc = convert_atac_indices(barcode)) %>%
      group_by(rna_bc)
    
    SNVs_perCell <- snv_count(cbsniffer_path)
    
    snvScore_byCell <- left_join(clusterLabels, frags_perCell, by = c("barcode" = "rna_bc")) %>%
      left_join(SNVs_perCell, by = c("barcode" = "rna_bc")) %>%
      left_join(cb, by = c("barcode" = "rna_bc"))
    
    snvScore_byCell[is.na(snvScore_byCell)] <- 0
    snvScore_byCell <- snvScore_byCell %>%
      mutate(denominator = numFragsOverlapping,
             mutscore = mut_count / denominator)
    
    write_csv(snvScore_byCell, paste0(path, f, "/nullScores.csv"))
}