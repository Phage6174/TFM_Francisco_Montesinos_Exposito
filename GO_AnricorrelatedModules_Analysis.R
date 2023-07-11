setwd("C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/")

library(WGCNA)
library(tidyverse)


options(stringsAsFactors = FALSE)
enableWGCNAThreads()

tumor_names <- list.files("input_files_processed/")
tumor_names <- tumor_names[-c(2,3,6,7,8,9,10,14,16,17,20,21,22)]
soft_powers <-
  c(10, 13, 11, 9, 13, 17, 13, 13, 8, 15)  # set soft powers for each tumor



for (i in seq_along(tumor_names)) {
  tumor_name <- tumor_names[i]
  tumor_name <- gsub("\\.txt$", "", tumor_name)
  softPower <- soft_powers[i]
  
  message <- sprintf("Analyzing tumor %s...\n", tumor_name)
  cat(message)
  
  # read data
  message("Step 1: Loading data.")
  
  datExpr0 <- as.data.frame(t(read.csv(
    paste0(
      "C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/input_files_processed/",
      tumor_name,
      ".txt"
    ),
    sep = "\t"
  )))
  
  # filter genes and samples
  message("Step 2: Removing non-expressed genes.")
  
  gsg <- goodSamplesGenes(datExpr0, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      print(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    }
    if (sum(!gsg$goodSamples) > 0) {
      print(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    }
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Create the output folder
  outFolder0 <- file.path("C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/tumor_analysisv3")
  if (!dir.exists(outFolder0)) {
    dir.create(outFolder0, recursive = TRUE)
    message("The output folder was created successfully:", outFolder0)
  }
  
  outFolder <- file.path(outFolder0, tumor_name)
  if (!dir.exists(outFolder)) {
    dir.create(outFolder, recursive = TRUE)
    message("The output folder was created successfully:", outFolder)
  }
  
  outFolderPower <- file.path(outFolder, paste0("power", softPower))
  if (!dir.exists(outFolderPower)) {
    dir.create(outFolderPower, recursive = TRUE)
    message("The output folder was created successfully:", outFolderPower)
  }
  
  setwd(outFolderPower)
  files <- list.files(pattern = "*.RData")
  if (length(files) > 0) {
    for (file in files) {
      load(file)
    }
  }
  
  GOFolder <- file.path(outFolderPower, "GO/")
  if (!dir.exists(GOFolder)) {
    dir.create(GOFolder, recursive = TRUE)
    message("The output folder was created successfully:", GOFolder)
  }
  
  message("Step 6: Eigengene adjecency.")
  
  # Heatmap  shows the eigengene adjacency AIJ = (1 + cor(EI , EJ ))/2
  eigenMat <- matrix(data = NA, nrow = ncol(MEs), ncol = ncol(MEs))
  for(i in 1:nrow(eigenMat)) {
    
    for(j in 1:ncol(eigenMat)) { 
      eigenMat[i,j] <- (1 + cor(MEs[,i], MEs[,j]))/2
    }
  }
  
  ## Where is NDUF family?
  
  geneModules <- bwnet$colors
  names <- as.character(sapply(names(geneModules), function(x) unlist(strsplit(x, split = "\\|", fixed = F))[[1]]))
  names(geneModules) <- names
  genesdf <- cbind(names(geneModules), geneModules)
  moduleList <- split(genesdf,genesdf[,2])
  mitocIndex <- unlist(lapply(moduleList, function(x) {
    sum(grepl("NDUF",x))
  }))
  mitocIndex <- mitocIndex[-1]
  mitocIndex <- names(mitocIndex)[which.max(mitocIndex)]
  mitocIndex <- as.numeric(mitocIndex)
  
  if(mitocIndex == 0) {
    next()
  } else {
    mitocModule <- names(bwnet$colors[bwnet$colors == mitocIndex])
    mitocModule <- as.character(sapply(mitocModule, function(x) unlist(strsplit(x, split = "\\|", fixed = F))[[1]]))
    outMitoc <- file.path(GOFolder, paste0(tumor_name, "_", softPower, "_module", mitocIndex, "_M.csv"))
    write.csv(mitocModule, outMitoc, row.names = FALSE)
    save(mitocModule, file = file.path(GOFolder, paste0(tumor_name, "_", mitocIndex, "_power", softPower,"_M_genes.RData")))
    ## Anticorrelated modules
    
    eigenMat2 <- eigenMat[-1, -1]
    eigenRest <- eigenMat2[mitocIndex,] 
    minIndex <- which(eigenRest < 0.4)
    
    if(length(minIndex) == 0) {
      next()
    } else {
      for(m in 1:length(minIndex)) {
        Tindex <- minIndex[m]
        t_module <- names(bwnet$colors[bwnet$colors == Tindex])
        t_module <- as.character(sapply(t_module, function(x) unlist(strsplit(x, split = "\\|", fixed = F))[[1]]))
        Tmodule <- file.path(GOFolder, paste0(tumor_name, "_", softPower, "_module", Tindex, "_T.csv"))
        write.csv(t_module, Tmodule, row.names = FALSE)
        save(t_module, file = file.path(GOFolder, paste0(tumor_name, "_", Tindex, "_power", softPower,"_T_genes.RData")))
      }
    }
  }
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  # Perform GO enrichment analysis using enrichGO
  ##For module M:
  
  mitocGO <- enrichGO(
    gene = mitocModule,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.2,
    qvalueCutoff = 0.2,
    readable = FALSE,
    keyType = "SYMBOL")
  
  png(filename = file.path(GOFolder, paste0(tumor_name, "_module", mitocIndex, "_barplot.png")),
      width = 800, height = 600)
  Mbarplot <- barplot(mitocGO, showCategory = 15, textsize = 8, main = paste0(
    "Barplot Enrichment Analysis-Module_",mitocIndex,"-", tumor_name))
  print(Mbarplot)
  dev.off()
  
  ## For module T:
  ## Check the number of CSV files with "_T" suffix
  T_file_count <- length(list.files(GOFolder, pattern = "_T.csv", full.names = TRUE))
  
  if (T_file_count == 1) {
    tGO <- enrichGO(
      gene = t_module,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.2,
      qvalueCutoff = 0.2,
      readable = FALSE,
      keyType = "SYMBOL")
    
    png(filename = file.path(GOFolder, paste0(tumor_name, "_module", Tindex, "_barplot.png")),
        width = 800, height = 600)
    Tbarplot <- barplot(tGO, showCategory = 15, textsize = 8, main = paste0(
      "Barplot Enrichment Analysis-Module_",Tindex,"-", tumor_name))
    print(Tbarplot)
    dev.off()
  } else {
    # Multiple CSV files with "_T" suffix, breaking the loop
    print("Multiple CSV files related to T module found. Hence, do it manually.")
  }
} 