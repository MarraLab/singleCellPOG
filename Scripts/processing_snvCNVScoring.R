source("../scFunctions.R")

# helper functions
runInferCNV <- function(seu, sampleAnnot, save_dir, ref = NULL){
  counts_matrix <- seu[["RNA"]]$counts
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                      annotations_file = sampleAnnot,
                                      delim = "\t",
                                      gene_order_file = "gencode_v21_gen_pos.complete.txt",
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

files <- files[grepl(paste(patient_ids, collapse="|"), files)]
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
args_file <- fread("../Data/WholeCohort.csv")
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

for(mb in multi_biop_samples){# rev(all_samples)){ # run the faster samples first
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
}

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
