# TFM_Francisco_Montesinos_Exposito

This repository aims to store the code and data used for the co-expression network analysis analysis so as to demonstrate pervasive occurrence of two negatively correlated gene modules in cancer.

This work has been done using a R package called **WGCNA** (Langfelder & Horvath, 2008)(Langfelder & Horvath, 2012).

## NetworkAnalysis_for_10powers.R

In this file every file of the folder with the RNAseq data for the tumors of TCGA will be processed:

  - Removing non expressed genes with *goodSamplesGenes()*.
  - Building the sample clusting tree to detect outliers with *hclust()*.
  - Power detection with *pickSoftThreshold()*.
  - Analysis with the first 10 powers over the threshold using adjacency matrices and TOM (topological overlapping matrices). 

## NetworkAnalysis_with_power.R

Once an analysis has been performed for the 10 powers for each tumor and a power has been chosen to create modules with fewer than 1000 genes, this script will test the selected power for each tumor. The following steps will be executed:

  - The data from the previous list will be loaded.
  - Non-expressed genes will be removed using the function *goodSamplesGenes()*.
  - Analysis will be performed using *blockwiseModules()*.
  - A relationship between each gene and its corresponding module will be established.
  - A heatmap will be generated to identify anticorrelated modules.
  - Anticorrelated modules and their genes will be extracted.

## GO_AnticorrelatedModules_Analysis.R

GO enrichment analysis for anticorrelated modules was performed using enrichGO() of the R package: **clusterProfiler** (Yu et al., 2012). This has not been done in the previous step because of an incompatibility between enrichplot (Yu et al., 2023) and org.Hs.eg.db (Carlson et al., 2019), with WGCNA.
## References

  - Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted correlation network analysis.” BMC Bioinformatics, 559. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559.
  - Langfelder P, Horvath S (2012). “Fast R Functions for Robust Correlations and Hierarchical Clustering.” Journal of Statistical Software, 46(11), 1–17. https://www.jstatsoft.org/v46/i11/.
  - Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. Omics : a journal of integrative biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118.
  - Yu G (2023). enrichplot: Visualization of Functional Enrichment Result. R package version 1.20.0, https://yulab-smu.top/biomedical-knowledge-mining-book/.
  - Carlson M (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2.



