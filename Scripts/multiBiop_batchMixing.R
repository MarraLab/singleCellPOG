source("scFunctions_forSubmission.R")

multiBiop_list <- readRDS("PATH/Data/multiBiopsy_separate_list.RDS")
# this object is a list of seurat objects where each object comprises all samples from a patient

for(i in seq(1, length(multiBiop_list))){
    seu <- multiBiop_list[[i]]
    seu <- AddMetaData(seu, allSamples_merged[[]])

    p1 <- DimPlot(seu, reduction = "rna.umap", group.by = c("sample")) +
        ggtitle("Sample") +
        scale_colour_jama() +
        theme(legend.position = "bottom")
    p2 <- DimPlot(seu, reduction = "rna.umap", group.by = c("identity")) +
        ggtitle("Cell Type") +
        scale_colour_nejm()
    
    plot <- p1 + p2
    ggsave(plot = plot, paste0("../../Figures/POGPaperFigures/UpdatedSet/", names(multiBiop_list)[i], "_batchCorr.jpg"), width = 6, height = 4)

    # calculating NMI and ARI on normal cells
    obj <- subset(seu, identity != "Malignant")

    nmi <- NMI(obj@meta.data[["sample"]], obj@meta.data$seurat_clusters, variant = "max")
    ari <- ARI(obj@meta.data[["sample"]], obj@meta.data$seurat_clusters)
    print(paste(nmi, ari))
}


