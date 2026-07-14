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
write_csv(regionsAgg, "../Data/detectedSNVs_regionAnnot.csv")

regionsAgg <- fread("../Data/detectedSNVs_regionAnnot.csv")
summRegions <- regionsAgg %>%
    group_by(sample) %>%
    mutate(numVar = n()) %>%
    group_by(region, sample, prop_inPeaks) %>%
    summarize(prop = n()/numVar) %>%
    dplyr::filter(!duplicated(cbind(region, sample, prop)))

# summary statistics for read coverage
snvCov <- data.frame(
    sample = character(),
    median = numeric(),
    quartile1 = numeric(),
    quartile3 = numeric(),
    prop2 = numeric()
)

for(f in files){
    snvs <- fread(paste0(path, f, "/self_cbSniffer_AllCounts.tsv")) %>%
        dplyr::filter(UB_ALT > 0)
    snvs2 <- dplyr::filter(snvs, UB_ALT > 1)
    
    pogSNVs <- fread(paste0(path, f, "/self_POGSNVcalls.tsv"))

    out <- data.frame(
        sample = f,
        median = median(snvs$UB_ALT),
        quartile1 = quantile(snvs$UB_ALT, 0.25),
        quartile3 = quantile(snvs$UB_ALT, 0.75),
        prop2 = nrow(snvs2) / nrow(pogSNVs)
    )
    snvCov <- rbind(snvCov, out)
}
