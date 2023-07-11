setwd("/home/usuario/Desktop/WGCNA/")

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

inputFolder <- "/input_files_processed/"
outputFolder <- "/home/usuario/Desktop/WGCNA/output/"

# Get the list of files in the input folder
fileList <- list.files(inputFolder, pattern = ".txt$", full.names = TRUE)

# Loop over each file
for (file in fileList) {
  # Extract the file name without extension
  fileName <- tools::file_path_sans_ext(basename(file))
  
  # Read the file
  inputData <- read.csv(file, sep = "\t")
  
  processWGCNA(inputData, fileName)
}

processWGCNA <- function(inputData, fileName) {
  outFolder <- paste0(outputFolder, fileName, "/")
  dir.create(outFolder, showWarnings = FALSE)
  
  inputFileName <- paste0(inputFolder, fileName, ".csv")
  inputFile <- read.csv(inputFileName, row.names = 1)
  datExpr0 <- as.data.frame(t(inputFile))
  
  gsg <- goodSamplesGenes(datExpr0, verbose = 3)
  
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  ## CLUSTERING SAMPLES TO LOOK FOR OUTLIERS
  sampleTree <- hclust(dist(datExpr0), method = "average")
  
  pdf(paste0(outFolder, fileName, "_sampleTree.pdf"), width = 12, height = 9)
  sizeGrWindow(12, 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
  
  ### POWER DETECTION
  # Choose a set of soft-thresholding powers
  powers <- c(1:30)
  sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
  
  pdf(paste0(outFolder, fileName, "_thresholding_powers.pdf"), width = 9, height = 5)
  sizeGrWindow(9, 5)
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
  dev.off()
  
  minSign <- sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > 0.9][1]
  powers <- c(minSign:(minSign + 10))
  powers <- powers[!(powers > 30)]
  
  ##### LOOP TO PERFORM THE ANALYSIS WITH 10 DIFFERENT POWERS
  for (i in 1:length(powers)) {
    softPower <- powers[i]
    outFolderPower <- paste0(outFolder, "power", softPower, "/")
    dir.create(outFolderPower)
    
    adjacency <- adjacency(datExpr0, power = softPower)
    TOM <- TOMsimilarity(adjacency)
    rm(adjacency)
    dissTOM <- 1 - TOM
    rm(TOM)
    gc()
    
    ## Clustering using TOM
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    
    ## Identificacion de modulos
    minModuleSize <- 30
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
    dynamicColors <- labels2colors(dynamicMods)
    
    ## Calculate eigengenes
    MEList <- moduleEigengenes(datExpr0, colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- 1 - cor(MEs)
    METree <- hclust(as.dist(MEDiss), method = "average")
    
    # Call an automatic merging function
    MEDissThres <- 0.25
    merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    mergedColors <- merge$colors
    mergedMEs <- merge$newMEs
    
    pdf(paste0(outFolderPower, fileName, "_power", softPower, "_dendrogram_and_clusters.pdf"), width = 12, height = 9)
    sizeGrWindow(12, 9)
    plotDendroAndColors(geneTree, mergedColors, "Merged dynamic",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
    
    outName <- paste0(outFolderPower, fileName, "_power", softPower, ".csv")
    write.csv(table(mergedColors), outName)
    
    save.image(file = paste0(outFolderPower, fileName, "_power", softPower, ".RData"))
  }
}
