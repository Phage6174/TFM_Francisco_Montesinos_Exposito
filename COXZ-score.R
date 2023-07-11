setwd("C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/tumor_analysisv3/")
zscore <- file.path("C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/tumor_analysisv3/z-score/")
if (!dir.exists("zscore")) dir.create("zscore")
setwd("zscore")
zscore_t <- read.csv("../../Gene expression-Z-scores.csv", sep = ";")
Amodule <- read.csv("../M_module_tumors2.csv", sep = ";")
Bmodule <- read.csv("../T_module_tumors2.csv", sep = ";")

rownames(zscore_t) <- zscore_t[,1]
zscore_t <- zscore_t[,-1]
library(ggplot2)
library(writexl)


# Function to search for a gene in zscore_t and get the corresponding value
get_gene_value <- function(gene, tumor_name) {
  if (gene %in% rownames(zscore_t)) {
    return(zscore_t[gene, tumor_name])
  } else {
    warning(paste("Gene", gene, "not found in zscore_t. Value set to NA."))
    return(NA)
  }
}

# Process Amodule and Bmodule for each tumor
tumor_names <- "THYMOMA"
results_list_A <- list()
results_list_B <- list()

for (tumor_name in tumor_names) {
  module_genes_A <- Amodule[[tumor_name]]
  module_genes_A <- module_genes_A[nzchar(module_genes_A) & !grepl("\\?", module_genes_A)]
  
  module_genes_B <- Bmodule[[tumor_name]]
  module_genes_B <- module_genes_B[nzchar(module_genes_B) & !grepl("\\?", module_genes_B)]
  
  values_A <- vector()
  values_B <- vector()
  
  cat(paste0("Analyzing the tumor: ", tumor_name))
  
  for (gen in module_genes_A) {
    value <- get_gene_value(gen, tumor_name)
    values_A <- c(values_A, value)
  }
  
  for (gen in module_genes_B) {
    value <- get_gene_value(gen, tumor_name)
    values_B <- c(values_B, value)
  }
  num_rows_A <- length(values_A)
  num_rows_B <- length(values_B)
  
  max_rows <- max(num_rows_A, num_rows_B)
  
  values_A <- c(values_A, rep(NA, max_rows - num_rows_A))
  values_B <- c(values_B, rep(NA, max_rows - num_rows_B))
  
  module_results_A <- cbind(Gene = module_genes_A, Value = values_A, Module = "A")
  module_results_B <- cbind(Gene = module_genes_B, Value = values_B, Module = "B")
  module_results_to_graph <- data.frame(A_module = values_A, B_module = values_B)
  
  module_filename_A <- paste0(tumor_name, "_module_A_results.csv")
  module_filename_B <- paste0(tumor_name, "_module_B_results.csv")
  module_filename_to_graph <- paste0(tumor_name, "_modules_to_graph.xlsx")
  
  module_path <- file.path("C:/Users/Usuario/Desktop/MADOBIS/TFM/TFMv5/tumor_analysisv3/zscore")
  
  if (!dir.exists(module_path)) {
    dir.create(module_path, recursive = TRUE)
    message("The output folder was created successfully:", module_path)
  }
  
  write.csv(module_results_A, file = paste0(module_path, "/", module_filename_A), row.names = FALSE)
  write.csv(module_results_B, file = paste0(module_path, "/", module_filename_B), row.names = FALSE)
  write_xlsx(module_results_to_graph, path = paste0(module_path, "/", module_filename_to_graph))
  
  results_list_A[[tumor_name]] <- module_results_A
  results_list_B[[tumor_name]] <- module_results_B
}


