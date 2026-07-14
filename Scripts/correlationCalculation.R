source("scFunctions_forSubmission.R")

# load all VAFs from cbSniffer
cbSniffer_vafs <- list()
for(f in files){
    cbSniffer_vafs[[f]] <- fread(paste0(path, f, "/self_cbSniffer_AllCounts.tsv")) %>%
        mutate(sample = f)
}
cbSniffer_vafs_df <- Reduce(rbind, cbSniffer_vafs)

# load VAFs from POG
pogVAFs <- list()
for(f in files){
    vafPath <- paste0(path, f, "/VAFS_numCellsContaining.csv")
    if(!file.exists(vafPath)){
        next
    }

    pogVAFs[[f]] <- fread(vafPath) %>%
        mutate(sample = f)
}
pogVAFs_df <- Reduce(rbind, pogVAFs)

# merging
merged <- inner_join(pogVAFs_df, cbSniffer_vafs_df, by = c("sample", "SNV_ID" = "type")) %>%
    dplyr::filter(UB_VAF != 0, UB_DEPTH > 1)

merged %>%
  ggscatter(x = "UB_VAF", y = "AF_t", 
            size = 0.5, alpha = 0.25,
            add = "reg.line", conf.int = F, cor.coef = F,
            xlab = "VAF from cbSniffer",
            ylab = "VAF from bulk WGS") +
    facet_wrap(~sample) +
    geom_smooth(aes(group = 1), method = "lm", color = "black") +
    stat_cor(aes(group = 1), label.y = 1, label.x.npc = "left", size = 3.5, method = "pearson") +
    theme_bw() +
    ylim(0, 1.1) +
    theme(strip.text = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
ggsave("Figures/scBulk_vafCorr.jpg", width = 10, height = 10)