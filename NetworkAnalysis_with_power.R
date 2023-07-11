setwd("C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/")

library(WGCNA)
library(tidyverse)


options(stringsAsFactors = FALSE)
enableWGCNAThreads()

tumor_names <- list.files("input_files_processed/")
tumor_names <- tumor_names[-c(4,12,13)]
#para LUSC no hay
soft_powers <-
  c(7, 10, 10, 6, 9, 11, 12, 5, 10, 4, 14, 13, 21, 9, 8, 4, 17, 29, 7)  # set soft powers for each tumor
tumors_power <- as.data.frame(cbind(tumor_names, soft_powers))


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
  
  # Load or create network and module eigengenes (MEs)
  message("Step 3: Creating network and MEs.")
  
  files <- list.files(pattern = "*.RData")
  if (length(files) > 0) {
    for (file in files) {
      load(file)
    }
  } else {
    print("There's no .RData file.")
    # If no existing network and MEs are found, create them
    maxSize <- round(ncol(datExpr0) / 4, digits = 0) + 150
    outBlock <- file.path(outFolderPower, paste0(tumor_name, "_power", softPower))
    bwnet <- blockwiseModules(datExpr0, maxBlockSize = maxSize, power = softPower, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs = TRUE, saveTOMFileBase = outBlock, verbose = 3)
    outName <- file.path(outFolderPower, paste0(tumor_name, "_power", softPower, ".csv"))
    write.csv(table(bwnet$colors), outName)
    outRData <- file.path(outFolderPower, paste0(tumor_name, "_power", softPower, ".RData"))
    
    # Recalculate MEs with color labels
    MEs <- moduleEigengenes(datExpr0, bwnet$colors)$eigengenes
    
    # Save variables in .RData file
    save.image(file = outRData)
  }
  
  message("Step 4: Tree and SoftThreshold plots.")
  
  print("Running the hierachical clustering diagram:")
  pdf(paste0(outFolderPower, tumor_name, "_dendrogram.pdf"))
  mergedColors <-  labels2colors(bwnet$colors)
  plotDendroAndColors(
    bwnet$dendrograms[[1]],
    mergedColors[bwnet$blockGenes[[1]]],
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  message("Step 5: Gene-Module Relationship.")
  
  GOFolder <- file.path(outFolderPower, "GO/")
  if (!dir.exists(GOFolder)) {
    dir.create(GOFolder, recursive = TRUE)
    message("The output folder was created successfully:", GOFolder)
  }
  # create GOFolder
  setwd(GOFolder)
  
  print("Module mapping:")
  module.gene.mapping <- as.data.frame(bwnet$colors)
  colnames(module.gene.mapping) <- c("Modules")
  
  # Convert the rownames to a new column called 'Gene'
  module.gene.mapping <- rownames_to_column(module.gene.mapping, var = "Genes")
  
  # Split the gene names by module and store in a list
  gene_lists <- split(module.gene.mapping$Genes, module.gene.mapping$Modules)
  
  # Create a data frame to store the genes associated with each module label
  gene_df <- data.frame(Module = character(), Gene = character(), stringsAsFactors = FALSE)
  
  # Iterate through the unique module labels and append the corresponding list of genes to the data frame
  for (i in unique(bwnet$colors)) {
    gene_df <- rbind(
      gene_df,
      data.frame(
        Module = paste0("ME", i),
        Gene = paste(gene_lists[[as.character(i)]], collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
  }
  
  modules_genes <- gene_df %>%
    separate_rows(Gene, sep = ";") %>%
    group_by(Module) %>%
    mutate(row = row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = Module, values_from = Gene) %>%
    as.data.frame()
  
  modules_genes <- modules_genes[, -1]
  
  # Save the data frame to a CSV file
  write.csv(
    modules_genes,
    file = paste0("module_genes_power_", softPower, ".csv"),
    row.names = FALSE
  )
  
  
  message("Step 5: Plotting heatmap and dendrogram.")
  
  setwd(outFolderPower)
  outHmap <- file.path(outFolderPower, paste0(tumor_name, "_power", softPower, "_heatmap.pdf"))
  pdf(outHmap)
  
  par(cex = 1.0, xpd = TRUE) # xpd = TRUE permite que el texto se dibuje fuera de los límites del gráfico
  plotEigengeneNetworks(
    MEs,
    "Eigengene adjacency heatmap",
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = FALSE,
    xLabelsAngle = 90, 
  )
  
  text(
    x = 0.11,
    y = 1.045,
    labels = tumor_name,
    adj = 0,
    cex = 1.2
  ) # ajusta la posición y y el tamaño del texto
  
  dev.off()

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
    minIndex <- which(eigenRest < 0.3)
    
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
}