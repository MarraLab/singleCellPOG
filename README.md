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
3. **analyses_targetComparisons.R**: Comparisons of expression prevalence for genes associated with therapeutic efficacy between malignant cells and normal cells matched by cell type of origin and biopsy site

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

