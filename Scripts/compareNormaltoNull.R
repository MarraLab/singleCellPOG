source("scFunctions_forSubmission.R")

nullSNVScores <- list()
for(f in files){
  nullSNVScores[[f]] <- fread(paste0(path, f, "/nullScores.csv")) %>%
    mutate(sample = f,
           barcode = paste(barcode, f, sep = "_"))
}
nullSNVScores <- Reduce(rbind, nullSNVScores)

selfSNVScores <- list()
for(f in files){
  selfSNVScores[[f]] <- fread(paste0(path, f, "/snvScores_perCell.csv")) %>%
    mutate(sample = f,
           barcode = paste(rna_bc, f, sep = "_"))
}
selfSNVScores <- Reduce(rbind, selfSNVScores)

mergedSNVScores <- inner_join(nullSNVScores, selfSNVScores, by = c("barcode", "sample"))

# distribution of SNV scores for normal cells
# cell types
cellTypes <- data.frame(
  barcode = colnames(allSamples_merged),
  cellType = allSamples_merged$identity
)

# distribution 
mergedSNVScores %>%
  dplyr::filter(mutscore.y <= 1, mutscore.x <= 1) %>%
  inner_join(cellTypes, by = "barcode") %>%
  dplyr::filter(cellType != "Malignant") %>% 
  ggplot(aes(x = mutscore.y - mutscore.x)) +
  geom_histogram(binwidth = 0.01)

# effect size
df <- mergedSNVScores %>%
  dplyr::filter(mutscore.y <= 1, mutscore.x <= 1) %>%
  dplyr::select(barcode, sample, mutscore.x, mutscore.y) %>%
  `colnames<-`(c("barcode", "sample", "Null", "Self")) %>%
  pivot_longer(cols = -c(barcode, sample), names_to = "snvCategory", values_to = "score") %>%
  inner_join(cellTypes, by = "barcode") 

fit <- glmmTMB(score ~ snvCategory * cellType + (1 | sample),
               data = df,
               family = gaussian()
              )

eff <- pairs(
  emmeans(fit, ~ snvCategory | cellType)
)

tbl_clean <- eff %>%
  as.data.frame() %>%
  transmute(
    CellType = cellType,
    Contrast = contrast,
    EffectSize = round(estimate, 3),
    Pvalue = signif(p.value, 3)
  )
tbl_clean