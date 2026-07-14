source("scFunctions_forSubmission.R")

# genome annotations
annotations <- build_annotations(
    genome = "hg38",
    annotations = c(
        "hg38_basicgenes",
        "hg38_genes_intergenic"
    )
)

regionsList <- list()
for(f in files){
    # get peaks
    seu <- readRDS(paste0(path, f, "/separate_ProcessedSeurat.RDS"))
    peaks <- data.frame(
        peak = rownames(GetAssay(seu, "ATAC"))
    ) %>%
        separate(col = "peak", into = c("chr", "start", "end"), sep = "-")
    peaks$start <- as.numeric(peaks$start)
    peaks$end <- as.numeric(peaks$end)

    # create granges obj
    peaks_granges <- GRanges(
        seqnames = peaks$chr,
        ranges = IRanges(start = peaks$start, end = peaks$end)
    )

    # get overlap with detected SNVs
    snvs <- fread(paste0(path, f, "/self_cbSniffer_AllCounts.tsv")) %>%
        dplyr::filter(UB_ALT > 0)
    snvs_granges <- GRanges(
        seqnames = snvs$chr,
        ranges = IRanges(start = snvs$start, width = 1)
    )

    overlaps <- findOverlaps(peaks_granges, snvs_granges, type = "any") %>%
        as.data.frame()
    
    # number in peaks vs out
    num_inPeaks <- length(unique(overlaps$subjectHits)) / nrow(snvs)

    # annotating where detected SNVs are
    annotated <- annotate_regions(
        regions = snvs_granges,
        annotations = annotations
    ) %>%
        as.data.frame() %>%
        dplyr::filter(!duplicated(cbind(seqnames, start, end))) %>%
        mutate(region = gsub("hg38_genes_", "", annot.type))

    regionsList[[f]] <- dplyr::select(annotated, seqnames, start, end, region) %>%
        mutate(sample = f, prop_inPeaks = num_inPeaks)
}
regionsAgg <- Reduce(rbind, regionsList)
write_csv(regionsAgg, "detectedSNVs_regionAnnot.csv")

regionsAgg <- fread("detectedSNVs_regionAnnot.csv")
