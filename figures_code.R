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
# A
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

txHistory_toPlot$patient_id <- factor(txHistory_toPlot$patient_id, 
                                      levels = rev(c("POG003", "POG217", "POG590", "POG130", "POG643")))
txHistory_toPlot %>%
  ggplot() +
  geom_segment(aes(x = daysAfterBio_start, xend = daysAfterBio_start + course_duration, 
                   y = patient_id, yend = patient_id, colour = TreatmentResponse),
               size = 6) +
  geom_point(aes(x = progressionDate, y = patient_id), shape = 21, fill = "grey", colour = "black", size = 3) +
  geom_point(aes(x = dateSinceEarliest, y = patient_id), shape = 23, fill = "red", colour = "black", size = 3) +
  facet_grid(~biopsyID, scales = "free", space = "free") +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Days after First Biopsy", y = "Patient ID", colour = "Response") +
  theme(legend.position = "top",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))

# B
seu <- multiBiop_list[["POG003"]]
seu <- AddMetaData(seu, allSamples_merged[[]])

# subset for malignant cells and re-cluster
DefaultAssay(seu) <- "RNA"
maligSeu <- subset(seu, identity == "Malignant") %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)

pog003_2biop <- subset(maligSeu, sample %in% c("POG003_3", "POG003_4"))
pog003_2biop <- FindNeighbors(pog003_2biop) %>%
  FindClusters(resolution = 0.5)
pog003_2biop <- AddMetaData(pog003_2biop, inner_join(
  data.frame(sample = pog003_2biop$sample),
  idMapping,
  by = "sample"
))
DimPlot(pog003_2biop, reduction = "rna.umap", group.by = c("newID")) +
  scale_colour_jama() +
  xlim(-8, 5) +
  ylim(-12,5) +
  ggtitle(NULL) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))

# C
DimPlot(pog003_2biop, reduction = "rna.umap", group.by = c("seurat_clusters"), alpha = 0.8) +
  paletteer::scale_colour_paletteer_d("calecopal::kelp1") +
  xlim(-8, 5) +
  ylim(-12,5) +
  ggtitle(NULL) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))

# D
library(SeuratExtend)
options(spe = "human")

pog003_2biop <- GeneSetAnalysis(pog003_2biop, genesets = hall50$human)
hm50 <- c("oxidative", "glycolysis", "oxygen", "fatty", "adipogenesis", "hypoxia", "apoptosis", "mtorc1", "tgf", "mapk")
hmIndices <- which(grepl(paste(hm50, collapse = "|"),
                          rownames(pog003_2biop@misc$AUCell$genesets),
                         ignore.case = T))

# sorted in order that makes sense
hmIndices <- c(26, 34, 35, 33, 12, 2, 36)
matr <- pog003_2biop@misc$AUCell$genesets[hmIndices,]
# stats <- CalcStats(matr, f = pog003_2biop$seurat_clusters)
Heatmap(CalcStats(matr, f = pog003_2biop$seurat_clusters, method = "mean")) +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)) + 
  labs(x = "Cluster", fill = "Mean") +
  paletteer::scale_fill_paletteer_c("grDevices::Inferno", limits = c(0, 0.1))

# E
pog217 <- multiBiop_list[["POG217"]]
pog217 <- AddMetaData(pog217, allSamples_merged[[]])

tSeu <- subset(pog217, identity == "NK/T") %>%
  seurat_rna_pipeline() %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)

# T cell umaps
tUMAP <- data.frame("UMAP_1" = tSeu@reductions[["rna.umap"]]@cell.embeddings[,"rnaumap_1"],
                    "UMAP_2" = tSeu@reductions[["rna.umap"]]@cell.embeddings[,"rnaumap_2"],
                    "Biopsy" = tSeu$sample,
                    "Cluster" = tSeu$seurat_clusters,
                    "Subtype" = tSeu$subtype)

DimPlot(tSeu, reduction = "rna.umap", group.by = c("sample")) +
  scale_colour_jama() +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(NULL) +
  theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 18),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 18),
    title = element_text(size = 20)
  )

# F
DimPlot(tSeu, reduction = "rna.umap", group.by = c("subtype")) +
  paletteer::scale_color_paletteer_d("MoMAColors::Rattner") +
  # theme(legend.position = "bottom") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(NULL) +
  theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 18),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 18),
    title = element_text(size = 20)
  )

# G
plotCluster_knownMarkers <- function(seu, geneList){
  meanExpr_sub <- data.frame(
    barcode = character(),
    cluster = numeric(),
    category = character(),
    meanExpr = numeric()
  )
  
  for(i in seq(1, length(geneList))){
    genes <- intersect(Features(seu), geneList[[i]])
    if(length(genes) <= 1){next}
    meanExpr <- seu[["RNA"]]$data[genes,] %>% colMeans(na.rm = T)
    out <- data.frame(
      barcode = colnames(seu),
      cluster = seu$seurat_clusters,
      category = names(geneList)[i],
      meanExpr = meanExpr
    )
    meanExpr_sub <- rbind(meanExpr_sub, out)
  }
  
  return(
    meanExpr_sub %>%
      ggplot(aes(x = cluster, y = meanExpr)) +
      geom_violin(width = 1.2, scale = "width") +
      geom_boxplot(width = 0.1) +
      facet_wrap(~category, ncol = 1) +
      theme_bw() +
      theme(legend.position = "bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(strip.text = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),) +
      labs(x = "Cluster", y = "Mean Expression of Marker Genes Per Cell")
    )
}

tCell_markers <- list(
  "CD8_Cytotoxic" = c("CD8A", "CD8B", "GZMB", "PRF1"),
  "CD4" = c("CD4", "IL7R"),
  "Tregs" = c("FOXP3", "IL2RA"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4"),
  "NK" = c("NKG7", "GNLY", "KLRB1"),
  "Proliferative" = c("MKI67", "TOP2A", "PCNA")
)

plotCluster_knownMarkers(tSeu, tCell_markers)  

tSeu$subtype <- case_when(
  tSeu$seurat_clusters == "0" ~ "CD8",
  tSeu$seurat_clusters == "1" ~ "CD4",
  tSeu$seurat_clusters == "2" ~ "CD8_ExhHi",
  tSeu$seurat_clusters == "3" ~ "Tregs",
  tSeu$seurat_clusters == "4" ~ "NK",
  tSeu$seurat_clusters == "5" ~ "Proliferating",
  tSeu$seurat_clusters == "6" ~ "Unknown",
  tSeu$seurat_clusters == "7" ~ "CD4"
)

# H
# running CellChat
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)

runCellChat <- function(seu, ident, savePath, savePath_cc){
  # Idents(seu) <- ident
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu, assay = "RNA")
  
  seu <- NormalizeData(seu)
  meta <- data.frame(
    samples = seu$orig.ident
  )
  
  cellChat <- createCellChat(object = seu, group.by = ident, assay = "RNA")
  cellChat <- addMeta(cellChat, meta = meta)
  
  cellChat@DB <- CellChatDB.use
  
  cellChat <- subsetData(cellChat)
  cellChat <- identifyOverExpressedGenes(cellChat)
  cellChat <- identifyOverExpressedInteractions(cellChat)
  
  cellChat <- smoothData(cellChat, adj = PPI.human)
  
  cellChat <- computeCommunProb(cellChat, raw.use = FALSE, type = "triMean")
  
  cellChat <- filterCommunication(cellChat, min.cells = 10)
  cellChat <- aggregateNet(cellChat)
  
  # sourceIndex <- which(grepl(paste(unique(endo_corrected$subtypes), collapse = "|"),
  #                            levels(cellChat@idents)))
  plot <- netVisual_bubble(cellChat, remove.isolate = FALSE)

  write_csv(as.data.frame(plot[["data"]]), savePath)
  saveRDS(cellChat, savePath_cc)
}

pog217_1_malig <- subset(maligSeu, sample == "POG217_1")  %>%
  seurat_rna_pipeline() %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.1)
pog217_1_malig$subtype <- paste("Malignant", pog217_1_malig$seurat_clusters, sep = "_")

pog217_1_cellchat <- subset(pog217, identity %in% c("NK/T", "Malignant") & sample == "POG217_1")
pog217_1_cellchat <- AddMetaData(pog217_1_cellchat, rbind(tSeu[[]], pog217_1_malig[[]]))
pog217_1_cellchat$subtype <- as.factor(pog217_1_cellchat$subtype)

runCellChat(pog217_1_cellchat, 
            "subtype", 
            paste0("../../Objects/POG217_1/tCellMalig_interactions.csv"), 
            paste0("../../Objects/POG217_1/tCellMalig_cellChat.RDS"))

pog217_2_malig <- subset(maligSeu, sample == "POG217_2")  %>%
  seurat_rna_pipeline() %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.1)
pog217_2_malig$subtype <- paste("Malignant", pog217_2_malig$seurat_clusters, sep = "_")

pog217_2_cellchat <- subset(pog217, identity %in% c("NK/T", "Malignant") & sample == "POG217_2")
pog217_2_cellchat <- AddMetaData(pog217_2_cellchat, rbind(tSeu[[]], pog217_2_malig[[]]))
pog217_2_cellchat$subtype <- as.factor(pog217_2_cellchat$subtype)

runCellChat(pog217_2_cellchat, 
            "subtype", 
            paste0("../../Objects/POG217_2/tCellMalig_interactions.csv"), 
            paste0("../../Objects/POG217_2/tCellMalig_cellChat.RDS"))

pog217_1_ccInter <- fread("../../Objects/POG217_1/tCellMalig_interactions.csv")
pog217_2_ccInter <- fread("../../Objects/POG217_2/tCellMalig_interactions.csv")

# interactions
pog217_inter <- rbind(
  mutate(pog217_1_ccInter, sample = "POG217_Bx1"),
  mutate(pog217_2_ccInter, sample = "POG217_Bx2")
) %>%
  dplyr::filter(grepl("Malignant", source), 
                target %in% tSeu$subtype,
                receptor %in% c(tCell_markers[["Exhaustion"]])) %>%
  group_by(target, ligand, receptor, sample) %>%
  summarize(avgProb = mean(prob)) %>%
  mutate(interaction = paste(ligand, receptor, sep = "->"))

pog217_inter_controls <- rbind(
  mutate(pog217_1_ccInter, sample = "POG217_Bx1"),
  mutate(pog217_2_ccInter, sample = "POG217_Bx2")
) %>%
  dplyr::filter(ligand %in% c("HLA-A", "HLA-B"),
                receptor == "CD8A") %>%
  group_by(target, ligand, receptor, sample) %>%
  summarize(avgProb = mean(prob)) %>%
  mutate(interaction = paste(ligand, receptor, sep = "->"))

ggplot(rbind(pog217_inter, pog217_inter_controls), aes(x = sample, y = interaction, colour = avgProb, size = avgProb)) +
  geom_point() +
  facet_grid(cols = vars(target)) +
  theme_bw() +
  theme(strip.text = element_text(angle = 90, size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)) + 
  paletteer::scale_color_paletteer_c("grDevices::ag_Sunset") + 
  labs(x = "Sample", y = "Interaction", colour = "Interaction\nProbability", size = "Interaction\nProbability")
ggsave("UpdatedSet/pog217_ccImmunosupp.jpg", width = 9, height = 10)

# I
pog590 <- multiBiop_list[["POG590"]]
DefaultAssay(pog590) <- 'chromvar'

pog590_differential.activity <- FindMarkers(
  object = pog590,
  ident.1 = '1', # second biopsy
  ident.2 = '3', # first biopsy
  # only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0
)

motifGenes <- fread("/projects/marralab/cayan_prj/PrecisionMed/Data/motifID_geneName.tsv")
pog590_differential.activity <- pog590_differential.activity %>%
  rownames_to_column(var = "jaspar_matrix") %>%
  inner_join(motifGenes)

pog590_differential.activity_favoursB2 <- pog590_differential.activity %>%
  dplyr::filter(avg_diff > 0) %>%
  arrange(p_val_adj) %>%
  rownames_to_column(var = "Rank") %>%
  mutate(Rank = as.numeric(Rank),
         logpval_adj = -log(p_val_adj))
pog590_differential.activity_favoursB1 <- pog590_differential.activity %>%
  dplyr::filter(avg_diff < 0) %>%
  arrange(desc(p_val_adj)) %>%
  rownames_to_column(var = "Rank") %>%
  mutate(Rank = as.numeric(Rank) + max(as.numeric(pog590_differential.activity_favoursB2$Rank)),
         logpval_adj = log(p_val_adj))

# plotting
pog590_differential.activity_fox <- rbind(pog590_differential.activity_favoursB2, pog590_differential.activity_favoursB1) %>%
  dplyr::filter(grepl("FOX", gene_name, ignore.case = T))

rbind(pog590_differential.activity_favoursB2, pog590_differential.activity_favoursB1) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "Yes", "No")) %>%
  ggplot(aes(x = Rank, y = logpval_adj, colour = Significant)) +
  geom_point(alpha = 0.25) +
  geom_point(data = pog590_differential.activity_fox,
             aes(x = Rank, y = logpval_adj), colour = "#0A74B2") + 
  geom_text_repel(data = dplyr::filter(pog590_differential.activity_fox, logpval_adj < 200),
                  aes(x = Rank, y = logpval_adj, label = gene_name),
                  colour = "black", max.overlaps = 50, point.padding = 0.1,
                  min.segment.length = 1, box.padding = 0.15, size = 6) +
  scale_color_manual(values = c("light grey", "black")) +
  labs(x = "", y = "Log(Adj. P-value)", colour = "Significant") +
  theme_linedraw() +
  theme(legend.position = "top") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 20),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20)
  )

ggplot(data = dplyr::filter(pog590_differential.activity_fox, logpval_adj > 200), 
       aes(x = Rank, y = logpval_adj)) +
  geom_point(colour = "#0A74B2") +
  geom_text_repel(aes(label = gene_name), colour = "black", max.overlaps = 50, point.padding = 0.1,
                  min.segment.length = 1, box.padding = 0.25, size = 4) +
  labs(x = "", y = "Log(Adj. P-value)", colour = "Significant") +
  theme_linedraw() +
  theme(legend.position = "top") +
  theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 20),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20)
  ) +
  scale_x_discrete(labels = NULL, breaks = NULL)
ggsave("UpdatedSet/pog590_motifEnrichment_foxOnly.jpg", width = 5, height = 5)

# J
pog590_1_maligCells <- subset(pog590, sample == "POG590_1" & identity == "Malignant")
DimPlot(pog590_1_maligCells, reduction = "rna.umap") # cluster 1 is the one that bridges into biopsy 2

# scoring for FOX family TFs
fox_ids <- rbind(pog590_differential.activity_favoursB2, pog590_differential.activity_favoursB1)  %>%
  dplyr::filter(grepl("^fox", gene_name, ignore.case = T)) %>%
  pull(jaspar_matrix)

fox_score <- colMeans(pog590[["chromvar"]]$data[fox_ids,], na.rm = T)
pog590$FOXMotifScore <- fox_score
pog590_malig <- subset(pog590, identity == "Malignant")
DimPlot(pog590_malig, reduction = "rna.umap")

pog590_1_maligCells$Label <- ifelse(Idents(pog590_1_maligCells) == 1,
                             ifelse(pog590_1_maligCells$sample == "POG590_1", "Biopsy 2-like", "Biopsy 2"),
                             "Biopsy 1")
DimPlot(pog590_malig, reduction = "rna.umap", group.by = "Label")
umap <- as.data.frame(pog590_malig@reductions[["rna.umap"]]@cell.embeddings) %>%
  mutate(Label = pog590_malig$Label)
biop2_like <- dplyr::filter(umap, Label == "Biopsy 2-like")
umap %>%
  dplyr::filter(Label != "Biopsy 2-like") %>%
  ggplot(aes(x = rnaumap_1, y = rnaumap_2, colour = Label)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = biop2_like, aes(x = rnaumap_1, y = rnaumap_2), colour = "#0A74B2", size = 1) +
  theme_classic() +
  theme(legend.position = "top") +
  theme(
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 20),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20)
  ) +
  xlim(-15, 0) +
  # scale_colour_continuous(high = "#250efe", low = "grey") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_colour_jama()

# K
# setting up comparisons
comp <- list(c("Biopsy 1", "Biopsy 2-like"), c("Biopsy 1", "Biopsy 2"), c("Biopsy 2-like", "Biopsy 2"))

data.frame(Label = pog590_malig$Label, FOXScore = pog590_malig$FOXMotifScore) %>%
  ggplot(aes(x = Label, y = FOXScore, fill = Label)) +
  geom_violin() +
  geom_boxplot(width = 0.075) +
  xlab("") +
  ylab("FOX TF Motif Score") +
  stat_compare_means(comparisons = comp, method = "wilcox.test", label = "p.format") +
  theme_minimal() +
  scale_x_discrete(limits = c('Biopsy 1', 'Biopsy 2-like', 'Biopsy 2')) +
  theme(
    axis.text       = element_text(size = 16),
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title      = element_text(size = 16),
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 20),
    title = element_text(size = 20)
  ) +
  scale_fill_manual(values = c("white", "white", "#0A74B2")) +
  theme(legend.position="none")

# Figure 5
# A
# from bulk RNA-seq
allSamples_BulkComp_zScore <- fread("../../Data/bulkRNA_zScores.csv") %>%
  dplyr::filter(hugo %in% relevant_genes,
                grepl("percentile|FC|tpmZscore", metric),
                !grepl("median", metric)) 

# from analyses_targetComparisons.R
geneExprComp_scRNA <- fread("../../Data/geneExprComp_scRNA_wscNormals_biopSite.csv") %>%
  rbind(geneExprComp_scRNA_origin) %>%
  dplyr::filter(!duplicated(cbind(sample, gene, cluster)))

pog_scComp <- allSamples_BulkComp_zScore %>%
  mutate(percentile = if_else(str_detect(metric, "percentile"),
                              value,
                              pnorm(z_score)*100)) %>%
  group_by(Sample, hugo) %>%
  dplyr::filter(percentile == max(percentile)) %>%
  ungroup() %>%
  dplyr::filter(!duplicated(cbind(Sample, hugo))) %>%
  full_join(geneExprComp_scRNA, by = c("hugo" = "gene", "Sample" = "sample")) %>%
  mutate(foundIn = ifelse(percentile >= 95,
                          ifelse(fdr_prop <= 0.05, "Both", "Bulk"),
                          ifelse(fdr_prop <= 0.05, "Single-cell", "Neither"))) %>%
  dplyr::filter(foundIn != "Neither")

pog_scComp %>%
  dplyr::select(Sample, hugo, cluster, prop_mal, mean_mal, foundIn) %>%
  dplyr::filter(!grepl("415", Sample)) %>%
  left_join(idMapping, by = c("Sample" = "sample")) %>%
  group_by(hugo) %>%
  mutate(numClusts = n()) %>%
  ggplot(aes(x = as.character(cluster), y = reorder(hugo, numClusts), shape = foundIn, colour = mean_mal, size = prop_mal)) +
  geom_point() +
  facet_grid(cols = vars(newID), scales = "free_x", space = "free_x") +
  scale_size_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Cluster", y = "Gene", colour = "Mean\nExpression", size = "Proportion\nof Cells", shape = "Significant in") +
  paletteer::scale_colour_paletteer_c("grDevices::Inferno") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme(strip.text = element_text(size = 20, angle = 90),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))
ggsave("UpdatedSet/allGraphKBGenes.jpg", width = 30, height = 18)

# bar plot of shared, sc, bulk across samples
# across genes
pog_scComp %>%
  group_by(hugo) %>%
  mutate(numClusts = n()) %>%
  ggplot(aes(x = reorder(hugo, numClusts), fill = foundIn)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)) +
  scale_fill_brewer(palette = "Dark2")

# across clusters
pog_scComp %>%
  dplyr::select(Sample, hugo, cluster, prop_mal, mean_mal, foundIn) %>%
  dplyr::filter(!grepl("415", Sample)) %>%
  ggplot(aes(x = cluster, fill = foundIn)) +
  geom_bar(position = "stack") +
  facet_grid(cols = vars(Sample), scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20)) +
  scale_fill_brewer(palette = "Dark2")

# B
pog_scComp %>%
  dplyr::filter(foundIn %in% c("Single-cell", "Both")) %>%
  group_by(Sample) %>%
  mutate(maxCluster = max(cluster) + 1) %>%
  ungroup() %>%
  group_by(Sample, hugo, maxCluster) %>%
  summarize(numClusts = n()) %>%
  ungroup() %>%
  mutate(propPresent = numClusts / maxCluster) %>%
  left_join(idMapping, by = c("Sample" = "sample")) %>%
  ggplot(aes(x = newID, y = propPresent)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)) +
  labs(y = "Proportion of Clusters", x = "Sample")

# C
pog1329 <- allSamples_list[["POG1329"]]
pog1329$maligClusts <- seu$seurat_clusters
pog1329$maligClusts <- as.character(pog1329$maligClusts)
pog1329$maligClusts[is.na(pog1329$maligClusts)] <- "Non-malignant"
DimPlot(pog1329, reduction = "rna.umap", group.by = "maligClusts") +
  ggtitle("POG1329 Malignant Clusters") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  paletteer::scale_color_paletteer_d("lisa::PabloPicasso_1") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.position = "bottom") +
  guides(colour = guide_legend(nrow=2, override.aes = list(size = 5)))

# D
FeaturePlot(pog1329, reduction = "rna.umap", features = "CD274") +
  ggtitle("POG1329 CD274 Expression") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.position = "bottom")
 # E
 pog_scComp %>% 
  dplyr::filter(!duplicated(cbind(Sample, hugo))) %>%
  group_by(hugo) %>%
  ggplot(aes(x = prop_mal, y = percentile, shape = foundIn, colour = foundIn)) +
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.position = "bottom") +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  labs(colour = "Significant in", shape = "Significant in", y = "Bulk Expression Percentile", x = "Proportion of Cells")

 # F
toPlot <- data.frame(RAF1 = allSamples_merged[["RNA"]]$data["RAF1",],
                    id = allSamples_merged$identity,
                    sample = allSamples_merged$sample) %>%
  mutate(category = ifelse(id == "Malignant", "Malignant", "Normal")) 
stat <- glmmTMB(RAF1 ~ category + sample, data = toPlot)
modelSumm <- coef(summary(stat))$cond

comp <- list(c("Malignant", "Normal"))
toPlot %>%
  dplyr::filter(sample %in% dplyr::filter(pog_scComp, hugo == "RAF1", prop_mal > prop_norm, foundIn %in% c("Single-cell", "Both"))$Sample) %>%
  group_by(sample, category) %>%
  summarize(prev = sum(RAF1 > 0) / n(),
            meanExpr = mean(RAF1)) %>%
  ungroup() %>%
  inner_join(idMapping, by = "sample") %>%
  # pivot_longer(cols = c(prev, meanExpr), names_to = "metric", values_to = "values") %>%
  # mutate(metric = ifelse(metric == "prev", "Proportion of Cells", "Mean Expression")) %>%
  ggplot(aes(x = category, y = prev)) +
  geom_boxplot() +
  geom_point(aes(colour = newID), size = 2, alpha = 1) +
  geom_line(aes(group = newID, color = newID)) +
  stat_compare_means(comparisons = comp, method = "wilcox.test",
                     label = "p.format", size = 5) +
  # geom_line(linewidth = 0.8, alpha = 0.5) +
  # facet_wrap(~metric, scales = "free") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)) +
  labs(x = "Cell Type", y = "Proportion of Cells Expressing RAF1", colour = "Sample")
