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
library(httpgd)
library(qs2)
library(ggalluvial)
library(ggsankey)
library(circlize)
library(RColorBrewer)
library(paletteer)
library(glmmTMB)
library(org.Hs.eg.db)

# Figure 1
# A
patientChar_long <- make_long(patient_char_mat_clean, `Patient ID`, `Broad Cancer Type`, `Stage at Biopsy`)
ggplot(patientChar_long, aes(x = x, 
                             next_x = next_x, 
                             node = node, 
                             next_node = next_node,
                             fill = factor(node),
                             label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 2.5, color = 1, fill = "white") +
  theme_sankey(base_size = 16) +
  scale_fill_viridis_d(option = "inferno", alpha = 0.95) +
  theme(legend.position = "none")

# B
DimPlot(allSamples_merged, reduction = "rna.umap", group.by = c("Comparator Broad_oncotree_type")) +
  paletteer::scale_color_paletteer_d("ggsci::default_uchicago") +
  ggtitle("Cancer Type")

# C
DimPlot(allSamples_merged, reduction = "rna.umap", group.by = c("identity")) +
  ggtitle("Cell Type") +
  scale_colour_nejm()

# D
FeaturePlot(allSamples_merged, reduction = "rna.umap", features = "snvScore") +
  ggtitle("SNV Score")

# E
comp <- list(c("Malignant", "Non-malignant"))
data.frame(
  snvScore = allSamples_merged$snvScore,
  cellType = allSamples_merged$identity
) %>%
  mutate(cellType = ifelse(cellType == "Malignant", "Malignant", "Non-malignant")) %>%
  ggplot(aes(x = cellType, y = snvScore)) +
  geom_violin(scale = "width", width = 1.1) +
  geom_boxplot(width = 0.1) +
  stat_compare_means(comparisons = comp, method = "wilcox.test",
                     label = "p.format", size = 5) +
  theme_linedraw() +
  labs(x = "Cell Type", y = "SNV Score") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# F
cohort %>%
  ggscatter(x = "cellTypeProp", y = "bulkTC_seq", 
            shape = 21, fill = "Comparator Broad_oncotree_type", 
            size = 5, alpha = 0.8,
            add = "reg.line", conf.int = F, cor.coef = F,
            xlab = "Proportion of Malignant Cells from Single-cell",
            ylab = "Inferred Tumour Content from Bulk WGS") +
  paletteer::scale_fill_paletteer_d("ggsci::default_uchicago") +
  geom_smooth(aes(group = 1), method = "lm", color = "black") +
  stat_cor(aes(group = 1), label.y = 0.3, label.x.npc = "left", size = 6, method = "spearman") +
  labs(fill = "Cancer Type") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  theme(strip.text = element_text(size = 14, angle = 90),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

# Figure 2
# A
# mark CD8+ T cells
cd8TCells_indices <- which(allSamples_merged$identity == "NK/T" & (allSamples_merged[["RNA"]]$data["CD8A",] > 0 | allSamples_merged[["RNA"]]$data["CD8B",] > 0))

# exclude NK cells
cd3TCells_indices <- which(allSamples_merged$identity == "NK/T" & allSamples_merged[["RNA"]]$data[c("CD3E", "CD3D", "CD3G"),] %>% colSums() > 0)

cd8TCells_indices <- intersect(cd8TCells_indices, cd3TCells_indices)
cellType_cd8TCellsLabelled <- data.frame(
  sampleName = allSamples_merged$sample,
  barcode = colnames(allSamples_merged),
  cellType = allSamples_merged$identity
) %>%
  mutate(cellType = ifelse(cellType == "NK/T",
                           ifelse(barcode %in% colnames(allSamples_merged)[cd8TCells_indices],
                                  "CD8+ T cells",
                                  "Other NK/T cells"),
                           cellType))

cellTypeProp <- cellType_cd8TCellsLabelled %>%
  group_by(sampleName) %>%
  mutate(numCells = n()) %>%
  ungroup() %>%
  group_by(sampleName, cellType) %>%
  summarize(scProp = n() / numCells) %>%
  ungroup() %>%
  dplyr::filter(!duplicated(cbind(sampleName, cellType)))
cellTypeProp <- inner_join(cellTypeProp, patient_char, by = c("sampleName"))
cellTypeProp <- inner_join(cellTypeProp, idMapping, by = c("sampleName" = "sample"))

cellTypeProp %>%
  ggplot(aes(x = newID, y = scProp, fill = cellType)) +
  geom_bar(stat = "identity", position = "stack") + 
  facet_grid(~`Comparator Broad_oncotree_type`, scales = "free_x", space = "free_x") + 
  labs(x = "Sample", fill = "Cell Type", y = "Proportion of Cells") +
  theme_bw() +
  paletteer::scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme(strip.text = element_text(size = 14, angle = 90),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# B
cd8_comp <- full_join(dplyr::filter(cellTypeProp, cellType == "CD8+ T cells"), 
                      cd8_scores)
cd8_comp$scProp[is.na(cd8_comp$scProp)] <- 0

cd8_comp %>%
  ggscatter(x = "T_cell_d8_score", y = "scProp", 
            shape = 21, fill = "Comparator Broad_oncotree_type", 
            size = 5, alpha = 0.8,
            add = "reg.line", conf.int = F, cor.coef = F,
            xlab = "CIBERSORT score for CD8+ T-cells",
            ylab = "Proportion of CD8+ T-cells") +
  paletteer::scale_fill_paletteer_d("ggsci::default_uchicago") +
  geom_smooth(aes(group = 1), method = "lm", color = "black") +
  stat_cor(aes(group = 1), label.y = 0.3, label.x.npc = "left", size = 6, method = "spearman") +
  labs(fill = "Cancer Type") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  theme(strip.text = element_text(size = 14, angle = 90),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

# C
seu$CD8_pos <- seu[["RNA"]]$data["CD8A",] > 0 | seu[["RNA"]]$data["CD8B",] > 0
seu$CD4_pos <- seu[["RNA"]]$data["CD4",] > 0
seu$doublePos <- seu$CD8_pos & seu$CD4_pos & seu$identity == "NK/T"

print(sum(seu$doublePos) / sum(seu$identity == "NK/T"))

sample <- str_split(colnames(seu)[1], "_")[[1]][2]
seu$newTIdents <- ifelse(seu$identity == "NK/T",
                        case_when(
                        seu$doublePos ~ "CD8+CD4+",
                        # seu$CD8_pos ~ "CD8+",
                        # seu$CD4_pos ~ "CD4+",
                        .default = "NK/T"
                        ),
                        "Non-T")

doublePos_umapCoord <- as.data.frame(Embeddings(object = seu[["rna.umap"]]))[names(which(seu$doublePos)),] %>%
    mutate(newTIdents = "CD8+CD4+")

DimPlot(seu, reduction = "rna.umap", group.by = "newTIdents") +
    geom_point(data = doublePos_umapCoord, aes(x = rnaumap_1, y = rnaumap_2, colour = newTIdents), size = 0.5, alpha = 0.8) +
    scale_colour_manual(values = pal_t) +
    ggtitle(NULL) +
    theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 20),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20),
    legend.position = "bottom"
    )

# D
cd8Tcells_pog1329 <- intersect(colnames(pog1329), colnames(allSamples_merged)[cd8TCells_indices])
tCells_traj <- subset(pog1329, cells = cd8Tcells_pog1329)

# regressing out cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(tCells_traj) <- "RNA"
cd8_joined <- CellCycleScoring(tCells_traj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cd8_joined$CC.Difference <- cd8_joined$S.Score - cd8_joined$G2M.Score
cd8_joined <- ScaleData(cd8_joined, vars.to.regress = "CC.Difference", features = rownames(cd8_joined))
cd8_joined <- RunPCA(cd8_joined) %>%
  RunUMAP(dims = 1:15)

cd8_joined <- FindNeighbors(cd8_joined) %>%
  FindClusters(resolution = 0.5)

# annotating at either single- or double-positive CD8+ T cells
dpt <- subset(cd8_joined, doublePos) %>%
  seurat_rna_pipeline()
dpt <- FindNeighbors(dpt) %>%
  FindClusters(resolution = 0.2)

cd8_joined$dptID <- dpt$seurat_clusters
cd8_joined$dptID <- ifelse(is.na(cd8_joined$dptID), "CD4-", cd8_joined$dptID)

DimPlot(cd8_joined, reduction = "umap", group.by = "dptID") +
  ggtitle(NULL) +
  scale_colour_jama() +
  theme(legend.position = "bottom",
        axis.text       = element_text(size = 16),
        axis.title      = element_text(size = 20),
        legend.text     = element_text(size = 16))

# E
# Find differentially expressed genes
dptMarkers <- FindAllMarkers(dpt, logfc.threshold = 0.3, min.pct = 0.1, only.pos = T) %>%
  dplyr::filter(p_val_adj <= 0.05)

dptMarkers_toPlot <- dptMarkers %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  dplyr::slice(1:5)

colAnnot <- data.frame(cluster = dpt$seurat_clusters) %>%
              arrange(cluster)
mat <- dpt[["RNA"]]$data[dptMarkers_toPlot, rownames(colAnnot)] %>%
  as.matrix()

ann_colors <- list(
  cluster = setNames(
    pal_jama("default")(3),  
    c("0", "1", "2")          
  )
)

pheatmap(mat, cluster_rows = F, cluster_cols = F,
         annotation_col = colAnnot,
         annotation_colors = ann_colors,
         show_colnames = F, 
         filename = "UpdatedSet/pog1329_dptDEGs.jpg", width = 8, height = 4)

# Figure 3
# B
data.frame(
  patient = sapply(malig_mbSeu_earlyBiop$newID, function(x) str_split(x, "_")[[1]][1]),
  sample = malig_mbSeu_earlyBiop$newID,
  laterBiopScore = apply(malig_mbSeu_earlyBiop[[]][,metadataField], 1, max, na.rm = TRUE)
) %>%
  rbind(data.frame(
    patient = sapply(malig_mbSeu_lateBiop$newID, function(x) str_split(x, "_")[[1]][1]),
    sample = malig_mbSeu_lateBiop$newID,
    laterBiopScore = malig_mbSeu_lateBiop$snvScore
  )) %>%
  ggplot(aes(x = sample, y = laterBiopScore)) +
  geom_violin(scale = "width", width = 1) +
  geom_boxplot(width = 0.1) +
  facet_grid(cols = vars(patient), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title      = element_text(size = 16),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20),
    strip.text = element_text(size = 16)
  ) +
  labs(x = "Sample", y = "SNV Score for Latest Bx")

# C
# identify CNVs called by ploidetect in later biopsies that are not present in earlier ones
genes <- read_tsv("../../Data/gencode_v21_gen_pos.complete.txt", 
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

cnvComp_propCells %>%
  dplyr::filter(CN_sample1 == 2) %>% # filter for neutral CNVs
  mutate(diff = biopsy2_percent - biopsy1_percent) %>%
  ggplot(aes(x = diff)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 20),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20)
  ) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 50, colour = "red") +
  labs(x = "Difference in Proportion (Biopsy 2 - Biopsy 1)", y = "Count")

# D & E
# circos for regions of large expansion
# genomic locations
largeExp_coords <- cnvComp_propCells %>%
  dplyr::filter(CN_sample1 == 2) %>% # filter for neutral CNVs
  mutate(diff = biopsy2_percent - biopsy1_percent,
         group = ifelse(diff > 50, "Large Expansion", ifelse(diff > 0, "Expansion", "Contraction")),
         chr.y = paste0("chr", chr.y)) %>%
  dplyr::filter(group == "Large Expansion")


# Install and load circlize if needed
# install.packages("circlize")
library(circlize)

# cnv_status could be "gain", "loss", "neutral"
# Define colors for CNV status
cnv_colors <- c(
  gain = "#E5195D",      # pink/red for gain
  loss = "#197FE5"     # blue for loss
  # neutral = "#84E0A3"    # green for neutral
)


for(p in unique(largeExp_coords$pair)){
  largeExp_coords_filtered <- dplyr::filter(largeExp_coords, pair == p)
  dfPlot <- data.frame(
    chromosome = largeExp_coords_filtered$chr.y,
    start = largeExp_coords_filtered$start.y,
    end = largeExp_coords_filtered$end.y,
    cnv_status = largeExp_coords_filtered$change
  )

  # Convert to genomic format (list of regions and values)
  # genomic_data <- list(dfPlot[, 1:3], dfPlot[, "cnv_status", drop = FALSE])

  # Initialize circos plot
  # jpeg(paste0("UpdatedSet/", p, "_largeExpCircos.jpg"), width = 4, height = 4, units = "in", res = 300)
  circos.clear()
  circos.par("track.height" = 0.05)
  circos.initializeWithIdeogram(species = "hg38", chromosome.index = unique(dfPlot$chromosome)) # or your species


  # Add CNV regions as a track
  circos.genomicTrackPlotRegion(dfPlot, ylim = c(0, 1), track.height = 0.2,
    panel.fun = function(region, value, ...) {
      for (i in seq_len(nrow(region))) {
        status <- value$cnv_status[i]
        cnv_col <- cnv_colors[[status]]
        circos.genomicRect(region[i, , drop = FALSE], ytop = 1, ybottom = 0,
                           col = cnv_col, border = NA)
      }
    }
  )
  # dev.off()
}

# F - H
# hone in on specific driver CNVs
# function for generating CNV plots
plotCNV_multiBiopsy <- function(mb, cnv, xLims, yLims){
  seu <- multiBiop_list[[mb]]
  seu <- AddMetaData(seu, allSamples_merged[[]])
  seu_malig <- subset(seu, identity == "Malignant")

  cnvObj <- multiBiopCNV_list[[mb]]
  cnvOI <- cnvObj@expr.data[cnv,] %>%
    as.data.frame()
  colnames(cnvOI) <- c(cnv)

  # remapping samples
  seu_malig <- AddMetaData(seu_malig, inner_join(
    data.frame(sample = seu_malig$sample),
    idMapping,
    by = "sample"
  ))

  seu_malig <- RenameCells(seu_malig, new.names = paste(sapply(colnames(seu_malig), function(x) convert_rna_indices(str_split(x, "_")[[1]][1])),
                                                        seu_malig$sample,
                                                        sep = "_"))
  seu_malig <- AddMetaData(seu_malig, cnvOI)

  plot <- data.frame("UMAP_1" = seu_malig@reductions[["rna.umap"]]@cell.embeddings[,"rnaumap_1"],
                     "UMAP_2" = seu_malig@reductions[["rna.umap"]]@cell.embeddings[,"rnaumap_2"],
                     "Biopsy" = seu_malig$newID,
                     "CNV" = seu_malig[[]][[cnv]]) %>%
    drop_na() %>%
    arrange(CNV) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, fill = CNV), colour = "black") +
    geom_jitter(pch = 21, size = 1.5, alpha = 0.8) +
    xlim(xLims) +
    ylim(yLims) +
    # paletteer::scale_color_paletteer_d("rockthemes::secondlaw", direction = -1) +
    scale_colour_jama() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
    labs(fill = "InferCNV Score") +
    ggtitle(cnv) +
    theme_bw() +
    theme(
      axis.text       = element_text(size = 16),
      axis.title      = element_text(size = 20),
      legend.text     = element_text(size = 16, angle = 90),
      legend.title    = element_text(size = 20),
      title = element_text(size = 20),
      legend.position = "bottom"
    )
  plot2 <- DimPlot(seu_malig, reduction = "rna.umap", group.by = "newID") +
    ggtitle(NULL) +
    labs(x = "UMAP_1", y = "UMAP_2") +
    scale_colour_nejm() +
    theme_bw() +
    theme(
      axis.text       = element_text(size = 16),
      axis.title      = element_text(size = 20),
      legend.text     = element_text(size = 16),
      legend.title    = element_text(size = 20),
      title = element_text(size = 20),
      legend.position = "bottom"
    ) +
    guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size = 2.5)))
  p3 <- plot2 + plot
  print(p3)
  ggsave(paste0("UpdatedSet/", mb, "_", cnv, ".jpg"), plot = p3, width = 9, height = 5)
}

# get proportions of cells with gain
getCNVProp <- function(mb, cnv){
  seu <- multiBiop_list[[mb]]
  seu <- AddMetaData(seu, allSamples_merged[[]])
  seu_malig <- subset(seu, identity == "Malignant")

  cnvObj <- multiBiopCNV_list[[mb]]
  cnvOI <- cnvObj@expr.data[cnv,] %>%
    as.data.frame()
  colnames(cnvOI) <- c(cnv)
  med <- median(cnvObj@expr.data[cnv,])

  # remapping samples
  seu_malig <- AddMetaData(seu_malig, inner_join(
    data.frame(sample = seu_malig$sample),
    idMapping,
    by = "sample"
  ))

  seu_malig <- RenameCells(seu_malig, new.names = paste(sapply(colnames(seu_malig), function(x) convert_rna_indices(str_split(x, "_")[[1]][1])),
                                                        seu_malig$sample,
                                                        sep = "_"))
  seu_malig <- AddMetaData(seu_malig, cnvOI)

  data.frame(sample = seu_malig$sample,
             cnvValue = seu_malig[[]][,cnv]) %>%
    mutate(cnvStatus = ifelse(cnvValue > med, T, F)) %>%
    group_by(sample) %>%
    summarize(propGain = sum(cnvStatus) / n())
}

# Figure 4
# recoding treatment response
tx_history <- fread("/projects/POG/Clinical_Actions/Oasis_data/20231204_merge_20221014/drug_treatment.tab") %>%
  dplyr::filter(patient_id %in% multiBiop_ids,
                pog_informed != "N") %>% 
  mutate(patient_id = str_replace(patient_id, " ", ""),
         TreatmentResponse = case_when(
  grepl("PD", best_response) ~ "Progressive disease",
  grepl("PR", best_response) ~ "Partial response",
  grepl("SD", best_response) ~ "Stable disease",
  best_response == "" ~ "No information",
  .default = best_response
)) %>%
  inner_join(biopDate, by = c("patient_id" = "Participant_project_identifier")) %>%
  dplyr::filter(course_begin_on >= `Biopsy Defined_date`) %>%
  arrange(original_drug_list, desc(`Biopsy Defined_date`)) %>%
  dplyr::filter(!duplicated(cbind(course_begin_on, original_drug_list))) # drugs are associated with biopsy closest to course start

# calculating number of days from first biopsy
tx_history$course_begin_on <- as.Date(tx_history$course_begin_on)
tx_history$course_end_on <- as.Date(tx_history$course_end_on)
tx_history$`Biopsy Defined_date` <- as.Date(tx_history$`Biopsy Defined_date`)
tx_history$progression_on <- as.Date(tx_history$progression_on)

txHistory_toPlot <- tx_history %>%
  dplyr::filter(course_duration > 0) %>%
  arrange(patient_id, `Biopsy Defined_date`) %>%
  group_by(patient_id) %>%
  mutate(earliestBiopsyDate = min(`Biopsy Defined_date`)) %>%
  ungroup() %>%
  mutate(dateSinceEarliest = `Biopsy Defined_date` - earliestBiopsyDate) %>%
  group_by(patient_id) %>%
  mutate(biopsyID = dense_rank(dateSinceEarliest)) %>% 
  mutate(daysAfterBio_start = course_begin_on - earliestBiopsyDate,
         progressionDate = progression_on - earliestBiopsyDate)
