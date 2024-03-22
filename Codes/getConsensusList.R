#' @description
#' This function takes as input a list containing a list of genes,
#' one for each Classifier used, and a threshold value indicating the
#' minimum number of algorithms in which the gene should be found
#' @param lista_geni list of list containing important genes for each Algorithm
#' used in previous simulations
#' @param n_min threshold number indicating the minimum number of algorithms 
#' in which each gene should be contemporary important
#' @returns consensus.list <- list of genes found by consensus 
getConsensus.list <- function(lista_geni, n_min = 5) {

  library(tidyverse)
  
  # Get all unique genes from all algorithms
  geni_unici <- sort(unique(unlist(lista_geni)))
  # Init a df
  df <- data.frame(RowName = names(lista_geni))
  
  # Add new cols for each unique gene and assign TRUE/FALSE based on its presence
  for (gene in geni_unici) {
    df[[gene]] <- map_lgl(df$RowName, ~ gene %in% lista_geni[[.x]])
  }
  
  # Set rownames as the Algorithm name
  rownames(df) <- df$RowName
  df$RowName <- NULL
  df <- as.data.frame(t(df))
  df <- as.data.frame(ifelse(df == "TRUE", 1, 0))
  
  # Set consensus condition
  ds <- df[rowSums(df) > (n_min - 1), ]
  consensus.list <- list(consensusGenes = rownames(ds))
  
  return(consensus.list)
  
}

# # Example usage function:
# consensus.geneList <- getConsensus.list(lista_geni = 
#                                           readRDS(file = "../Results/list_genes_1.rds"), 
#                                         n_min = 7)
