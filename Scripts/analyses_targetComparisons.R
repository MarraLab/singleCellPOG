source("scFunctions_forSubmission.R")

# relevant genes from gkb
graphkb <- fread("Data/GraphKBExpression.txt") 
relevance_reTreatment <- c("sensitivity", "response", "eligibility", "likely sensitivity", "targetable")

relevant_gkb <- relevant_genes <- graphkb %>%
  dplyr::filter(expression_type == "increased",
                relevance %in% relevance_reTreatment,
                !duplicated(gene))
relevant_genes <- relevant_gkb %>%
  pull(gene) %>%
  unique()

# load in all normal references
files <- list.files("PATH")
files <- files[grepl("_refFiltered.RDS", files)]

refList <- list()
for(f in files){
  seu <- readRDS(paste0("PATH/", f))
  seu$batch <- f
  refList[[f]] <- seu
}

refCombined <- Reduce(merge, refList) %>%
  JoinLayers() %>%
  seurat_rna_pipeline()
saveRDS(refCombined, "PATH/normalRef_allReferences.RDS")
refCombined <- readRDS("PATH/normalRef_allReferences.RDS")

# matching to normal cells of origin
patient_char$tissueComparator <- c("thymus", "thymus", "posterior part of tongue", "lung", "buccal mucosa",
                                   "skin", "skin", "skin", "skin", "lung", 
                                   "lung", "lung", "lung", "lung", "bone",
                                   "bone", "colon", "colon", "colon", "mammary", 
                                   "colon", "colon", "muscle", "cerebellum", "cerebellum")

# biopsy site
patient_char$tissueComparator <- c("lung", "pleura", "posterior part of tongue", "lung", "lung",
                                   "skin", "skin", "skin", "lymph", "lung", 
                                   "pleura", "lung", "skin", "lymph", "spin",
                                   "spin", "spin", "spin", "spin", "skin", 
                                   "skin", "skin", "lung", "cerebellum", "cerebellum")
patient_char$cellTypeComparator <- c("epitheli", "epitheli", "stratified squamous epithelial cell", "epitheli", "epitheli", 
                                     "melanocyte", "melanocyte", "melanocyte", "epitheli", "epithelial", 
                                     "epitheli", "epitheli", "epitheli", "epitheli", "mesenchym|fibro", 
                                     "mesenchym|fibro", "epitheli", "epitheli", "epitheli", "epitheli",
                                     "epitheli", "epitheli", "mesenchym|fibro", "astrocyte", "astrocyte") 
        
subsetNormals_list <- list()
singleMalig_list <- list()

for(i in seq(1, nrow(patient_char))){
  f <- patient_char$sampleName[i]
  tc <- patient_char$tissueComparator[i]
  cc <- patient_char$cellTypeComparator[i]

  tcIndices <- which(grepl(tc, refCombined$tissue))
  ccIndices <- which(grepl(cc, refCombined$cell_type))
  inclIndices <- intersect(tcIndices, ccIndices)
  if(length(inclIndices) < 100){next}
  # print(length(inclIndices))

  subsetNormals_list[[f]] <- subset(refCombined, cells = colnames(refCombined)[inclIndices])
  # print(ncol(subsetNormals_list[[f]]))

  # malignant cells from POG samples
  singleMalig_list[[f]] <- subset(allSamples_list[[f]], identity == "Malignant") %>%
    seurat_rna_pipeline() %>%
    FindNeighbors() %>%
    FindClusters(resolution = 0.15)
  # print(DimPlot(singleMalig_list[[f]], reduction = "rna.umap", group.by = "seurat_clusters") +
  #   ggtitle(f))
}

compare_clusters_vs_normal <- function(seurat_obj,
                                       normal_obj,
                                       genes = relevant_genes,
                                       defaultAssay = "RNA") {
  DefaultAssay(seurat_obj) <- defaultAssay
  DefaultAssay(normal_obj) <- defaultAssay
  
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) stop("None of the supplied genes were found in the data set.")
  
  clusters <- unique(Idents(seurat_obj))
  expr_norm <- GetAssayData(normal_obj, layer = "data")[genes, , drop = FALSE]
  
  all_results <- list()
  
  for (cl in clusters) {
    cluster_cells <- WhichCells(seurat_obj, idents = cl)
    expr_cluster <- GetAssayData(seurat_obj, layer = "data")[genes, cluster_cells, drop = FALSE]
    
    res <- lapply(genes, function(g) {
      v_mal <- expr_cluster[g, ]
      v_norm <- expr_norm[g, ]
      
      # Breadth of expression
      prop_mal <- mean(v_mal > 0)
      prop_norm <- mean(v_norm > 0)
      tab <- matrix(c(sum(v_mal > 0), sum(v_mal == 0),
                      sum(v_norm > 0), sum(v_norm == 0)), nrow = 2)
      p_prop <- fisher.test(tab, alternative = "greater")$p.value
      
      # Magnitude of expression
      p_expr <- wilcox.test(v_mal, v_norm, alternative = "greater")$p.value

      mean_mal <- mean(v_mal)
      mean_norm <- mean(v_norm)
      
      # Combine p-values
      z_prop <- qnorm(1 - p_prop/2) * sign(prop_mal - prop_norm)
      z_expr <- qnorm(1 - p_expr/2) * sign(mean_mal - mean_norm)
      w_prop <- 0.7; w_expr <- 0.3
      z_comb <- (w_prop * z_prop + w_expr * z_expr) / sqrt(w_prop^2 + w_expr^2)
      p_comb <- 2 * (1 - pnorm(abs(z_comb)))
      
      data.frame(cluster = cl,
                 gene = g,
                 prop_mal = prop_mal,
                 prop_norm = prop_norm,
                 mean_mal = mean_mal,
                 mean_norm = mean_norm,
                 p_prop = p_prop,
                 p_expr = p_expr,
                 p_combined = p_comb,
                 stringsAsFactors = FALSE)
    })
    
    all_results[[cl]] <- bind_rows(res)
  }
  
  final_df <- bind_rows(all_results) %>%
    group_by(cluster) %>%
    mutate(fdr_prop = p.adjust(p_prop, method = "BH"),
           fdr_expr = p.adjust(p_expr, method = "BH"),
           fdr_combined = p.adjust(p_combined, method = "BH")) %>%
    arrange(cluster, fdr_combined)
  
  return(final_df)
}

geneExprComp_scRNA <- list()
for(i in seq(1, length(subsetNormals_list))){
  f <- names(subsetNormals_list)[i]
  normal_obj <- subsetNormals_list[[f]]
  seurat_obj <- singleMalig_list[[f]]
  # print(colnames(seurat_obj)[1:5])

  out <- compare_clusters_vs_normal(seurat_obj = seurat_obj, normal_obj = normal_obj) %>%
    mutate(sample = f)
  geneExprComp_scRNA[[f]] <- out
}
geneExprComp_scRNA <- Reduce(rbind, geneExprComp_scRNA)

