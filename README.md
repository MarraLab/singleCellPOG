# A feasibility study for improving the resolution of malignant and immune populations in precision cancer genomic medicine patients through the integration of bulk and single-cell multiomics

## Description of scripts, in the order of the workflow

1. **processing_cellRangerDownstream.R**: Performs QC and data processing on CellRanger outputs to generate Seurat objects containing both snRNA-seq and snATAC-seq information.
2. **processing_snvCNVScoring.R**: Analyses related to CNVs and SNVs
    1. For all samples:
        1. Infer CNVs from snRNA-seq
        2. Calculate SNV scores from snATAC-seq using SNVs detected from matched bulk WGS
        3. Compare SNV scores to a randomly generated background. Run **compareNormaltoNull.R** afterwards to visualize score distributions and calculate effect sizes. 
    2. For longitudinal biopsies:
        1. Calculate SNV scores using SNVs detected from bulk WGS of the later biopsy
        2. Compare changes in the proportion of malignant cells with a given CNV
3. **multiBiop_batchMixing.R**: 
4. **analyses_targetComparisons.R**: Comparisons of expression prevalence for genes associated with therapeutic efficacy between malignant cells and normal cells matched by cell type of origin and biopsy site

## Figure generation

**figures_code.R** uses the following as inputs, which are not uploaded into this repo for patient privacy or file size reasons:

* patient_char_mat_clean: Record of patient data, including cancer type and biopsy site
* allSamples_merged: Seurat object containing all cells from all patients; output of the workflow above
* cnvComp_propCells: Proportion of malignant cells with a given CNV in longitudinal biopsies; output of processing_snvCNVScoring.R
* multiBiopCNV_list: List of Seurat objects containing cells from longitudinal biopsies
* tx_history: History of treatments received by patients

## Analyses not shown in the manuscript

* **correlationCalculation.R**: Correlates VAFs derived from bulk WGS and snATAC-seq across SNVs
* **snvCoverage.R**: Assesses how SNVs detected in snATAC-seq are distributed across peaks and genomic regions

## Session info

R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

locale:
[1] C

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] aricode_1.0.3                           
 [2] annotatr_1.24.0                         
 [3] org.Hs.eg.db_3.16.0                     
 [4] paletteer_1.6.0                         
 [5] circlize_0.4.15                         
 [6] ggsankey_0.0.99999                      
 [7] ggalluvial_0.12.5                       
 [8] qs2_0.1.7                               
 [9] httpgd_2.0.4                            
[10] rstatix_0.7.2                           
[11] emmeans_1.11.2                          
[12] glmmTMB_1.1.11                          
[13] slingshot_2.6.0                         
[14] TrajectoryUtils_1.6.0                   
[15] princurve_2.1.6                         
[16] viridis_0.6.3                           
[17] viridisLite_0.4.2                       
[18] CellChat_2.1.2                          
[19] bigmemory_4.6.4                         
[20] igraph_2.0.3                            
[21] SeuratExtend_1.1.0                      
[22] SeuratExtendData_0.2.1                  
[23] SeuratWrappers_0.3.2                    
[24] ggrepel_0.9.3                           
[25] RColorBrewer_1.1-3                      
[26] ggsci_3.0.0                             
[27] ggpubr_0.6.0                            
[28] harmony_1.2.4                           
[29] Rcpp_1.0.11                             
[30] infercnv_1.3.3                          
[31] scDblFinder_1.12.0                      
[32] BSgenome.Hsapiens.UCSC.hg38_1.4.5       
[33] BSgenome_1.66.3                         
[34] rtracklayer_1.58.0                      
[35] Biostrings_2.66.0                       
[36] XVector_0.38.0                          
[37] TFBSTools_1.36.0                        
[38] JASPAR2020_0.99.10                      
[39] motifmatchr_1.20.0                      
[40] biomaRt_2.54.1                          
[41] TxDb.Hsapiens.UCSC.hg38.knownGene_3.16.0
[42] pheatmap_1.0.12                         
[43] EnsDb.Hsapiens.v86_2.99.0               
[44] ensembldb_2.22.0                        
[45] AnnotationFilter_1.22.0                 
[46] GenomicFeatures_1.50.4                  
[47] AnnotationDbi_1.60.2                    
[48] scran_1.26.2                            
[49] scRNAseq_2.12.0                         
[50] scater_1.26.1                           
[51] scuttle_1.8.4                           
[52] SingleCellExperiment_1.20.1             
[53] SummarizedExperiment_1.28.0             
[54] Biobase_2.58.0                          
[55] GenomicRanges_1.50.2                    
[56] GenomeInfoDb_1.34.9                     
[57] IRanges_2.32.0                          
[58] S4Vectors_0.36.2                        
[59] BiocGenerics_0.44.0                     
[60] MatrixGenerics_1.10.0                   
[61] matrixStats_1.0.0                       
[62] Signac_1.14.0                           
[63] Seurat_5.2.1                            
[64] SeuratObject_5.0.2                      
[65] sp_2.0-0                                
[66] Matrix_1.6-4                            
[67] data.table_1.14.8                       
[68] lubridate_1.9.2                         
[69] forcats_1.0.0                           
[70] stringr_1.5.0                           
[71] dplyr_1.1.2                             
[72] purrr_1.0.1                             
[73] readr_2.1.4                             
[74] tidyr_1.3.0                             
[75] tibble_3.2.1                            
[76] ggplot2_4.0.2                           
[77] tidyverse_2.0.0                         

loaded via a namespace (and not attached):
  [1] TMB_1.9.17                    pbapply_1.7-2                
  [3] lattice_0.20-45               vctrs_0.6.2                  
  [5] mgcv_1.8-41                   blob_1.2.4                   
  [7] survival_3.4-0                nloptr_2.0.3                 
  [9] spatstat.data_3.0-1           later_1.3.1                  
 [11] DBI_1.1.3                     R.utils_2.12.2               
 [13] rappdirs_0.3.3                uwot_0.1.16                  
 [15] dqrng_0.3.0                   zlibbioc_1.44.0              
 [17] htmlwidgets_1.6.2             mvtnorm_1.2-2                
 [19] GlobalOptions_0.1.2           future_1.33.0                
 [21] parallel_4.2.2                irlba_2.3.5.1                
 [23] KernSmooth_2.23-20            promises_1.2.0.1             
 [25] DelayedArray_0.24.0           limma_3.54.2                 
 [27] RcppParallel_5.1.7            RSpectra_0.16-1              
 [29] fastmatch_1.1-3               digest_0.6.33                
 [31] png_0.1-8                     rjags_4-14                   
 [33] bluster_1.8.0                 sctransform_0.4.1            
 [35] cowplot_1.2.0.9000            pkgconfig_2.0.3              
 [37] GO.db_3.16.0                  gridBase_0.4-7               
 [39] spatstat.random_3.1-4         DelayedMatrixStats_1.20.0    
 [41] estimability_1.5.1            ggbeeswarm_0.7.2             
 [43] reformulas_0.4.0              iterators_1.0.14             
 [45] minqa_1.2.5                   statnet.common_4.9.0         
 [47] reticulate_1.44.1             network_1.18.2               
 [49] spam_2.9-1                    beeswarm_0.4.0               
 [51] modeltools_0.2-23             GetoptLong_1.0.5             
 [53] bslib_0.5.0                   zoo_1.8-12                   
 [55] tidyselect_1.2.0              reshape2_1.4.4               
 [57] ica_1.0-3                     rlang_1.1.0                  
 [59] jquerylib_0.1.4               glue_1.6.2                   
 [61] registry_0.5-1                lambda.r_1.2.4               
 [63] CNEr_1.34.0                   ggsignif_0.6.4               
 [65] httpuv_1.6.11                 BiocNeighbors_1.16.0         
 [67] TH.data_1.1-2                 seqLogo_1.64.0               
 [69] annotate_1.76.0               jsonlite_1.8.7               
 [71] bit_4.0.5                     mime_0.12                    
 [73] systemfonts_1.0.4             gridExtra_2.3                
 [75] gplots_3.1.3                  Rsamtools_2.14.0             
 [77] stringi_1.7.12                RcppRoll_0.3.0               
 [79] spatstat.sparse_3.0-1         rbibutils_2.2.16             
 [81] scattermore_1.2               spatstat.explore_3.1-0       
 [83] Rdpack_2.6                    bitops_1.0-7                 
 [85] cli_3.6.1                     RSQLite_2.3.1                
 [87] bigmemory.sri_0.1.8           libcoin_1.0-9                
 [89] timechange_0.2.0              GenomicAlignments_1.34.1     
 [91] nlme_3.1-160                  fastcluster_1.2.3            
 [93] locfit_1.5-9.8                listenv_0.9.0                
 [95] miniUI_0.1.1.1                R.oo_1.25.0                  
 [97] ggnetwork_0.5.13              dbplyr_2.3.3                 
 [99] lifecycle_1.0.3               ExperimentHub_2.6.0          
[101] R.methodsS3_1.8.2             caTools_1.18.2               
[103] codetools_0.2-18              coda_0.19-4                  
[105] vipor_0.4.5                   lmtest_0.9-40                
[107] xtable_1.8-4                  ROCR_1.0-11                  
[109] formatR_1.14                  BiocManager_1.30.27          
[111] abind_1.4-5                   farver_2.1.1                 
[113] FNN_1.1.3.2                   parallelly_1.36.0            
[115] AnnotationHub_3.6.0           RANN_2.6.1                   
[117] poweRlaw_0.70.6               BiocIO_1.8.0                 
[119] RcppAnnoy_0.0.21              goftest_1.2-3                
[121] patchwork_1.3.2.9000          futile.options_1.0.1         
[123] dichromat_2.0-0.1             cluster_2.1.4                
[125] future.apply_1.11.0           ellipsis_0.3.2               
[127] prettyunits_1.1.1             ggridges_0.5.6               
[129] remotes_2.4.2                 unigd_0.1.3                  
[131] argparse_2.2.2                spatstat.utils_3.1-3         
[133] htmltools_0.5.5               BiocFileCache_2.6.1          
[135] yaml_2.3.7                    NMF_0.26                     
[137] utf8_1.2.3                    plotly_4.10.2                
[139] interactiveDisplayBase_1.36.0 XML_3.99-0.13                
[141] withr_2.5.0                   fitdistrplus_1.1-11          
[143] BiocParallel_1.32.6           bit64_4.0.5                  
[145] xgboost_1.7.7.1               rngtools_1.5.2               
[147] multcomp_1.4-25               foreach_1.5.2                
[149] ProtGenerics_1.30.0           progressr_0.18.0             
[151] rsvd_1.0.5                    ScaledMatrix_1.6.0           
[153] memoise_2.0.1                 tzdb_0.4.0                   
[155] curl_5.0.1                    fansi_1.0.4                  
[157] fastDummies_1.7.3             tensor_1.5                   
[159] edgeR_3.40.2                  regioneR_1.30.0              
[161] cachem_1.0.8                  deldir_1.0-9                 
[163] metapod_1.6.0                 rjson_0.2.21                 
[165] clue_0.3-64                   tools_4.2.2                  
[167] sass_0.4.7                    sandwich_3.0-2               
[169] magrittr_2.0.3                RCurl_1.98-1.12              
[171] car_3.1-2                     TFMPvalue_0.0.9              
[173] ape_5.7-1                     xml2_1.3.5                   
[175] httr_1.4.6                    boot_1.3-28                  
[177] globals_0.16.2                R6_2.5.1                     
[179] RcppHNSW_0.4.1                DirichletMultinomial_1.40.0  
[181] progress_1.2.2                KEGGREST_1.38.0              
[183] gtools_3.9.5                  shape_1.4.6                  
[185] statmod_1.5.0                 coin_1.4-2                   
[187] beachmat_2.14.2               rematch2_2.1.2               
[189] BiocVersion_3.16.0            sna_2.7-2                    
[191] BiocSingular_1.14.0           splines_4.2.2                
[193] carData_3.0-5                 colorspace_2.1-0             
[195] generics_0.1.3                pracma_2.4.2                 
[197] pillar_1.9.0                  S7_0.2.0                     
[199] uuid_1.1-0                    GenomeInfoDbData_1.2.9       
[201] plyr_1.8.8                    dotCall64_1.0-2              
[203] gtable_0.3.6                  futile.logger_1.4.3          
[205] stringfish_0.18.0             restfulr_0.0.15              
[207] ComplexHeatmap_2.15.1         fastmap_1.1.1                
[209] doParallel_1.0.17             broom_1.0.5                  
[211] scales_1.4.0                  filelock_1.0.2               
[213] backports_1.4.1               lme4_1.1-34                  
[215] hms_1.1.3                     Rtsne_0.16                   
[217] shiny_1.7.4.1                 numDeriv_2016.8-1.1          
[219] polyclip_1.10-4               grid_4.2.2                   
[221] lazyeval_0.2.2                crayon_1.5.2                 
[223] MASS_7.3-58.1                 sparseMatrixStats_1.10.0     
[225] reshape_0.8.9                 svglite_2.1.1                
[227] compiler_4.2.2                spatstat.geom_3.1-0     