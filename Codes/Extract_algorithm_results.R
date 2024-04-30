setwd("../Desktop/giorgio/ensembleBP/Codes/")


library(readr)
library(tidyverse)
library(dplyr)




#' @description
#' This function extracts algorithm results from files containing accuracy tables and gene lists.
#' It merges the data, filters by algorithm name, selects specified statistic result , and computes desired statistics.

extract_algorithm_results <- function(path_output ="../Results/",
                                      algorithm_name="AdaBag",
                                      column_name="Balanced Accuracy",
                                      value_type = "max") {
  
  #' @input:
  #' @param Path to the directory containing result files. Default is "../Results/".
  #' @param algorithm_name: Name of the algorithm to filter by. Default is "AdaBag".
  #' @param column_name: Name of the column containing the desired accuracy values. Default is "Balanced Accuracy".
  #' @param value_type: Type of value or statistic to compute ("max", "median", or "mean"). Default is "max".
  
  #' @Output:
  #' @List: 
  #' @return List of top genes associated with the selected algorithm and statistic.
  #' @return df_tot: Merged dataframe containing all genes and accuracy information.
  
  # Read accuracy table filenames and extract algorithm names
  table_accuracy <- list.files(path = path_output,
                               recursive = TRUE,
                               pattern = "^AccuracyTable_",
                               full.names = TRUE)
  df_accuracy <- readr::read_delim(table_accuracy, id = "file_name")
  cv_number <- gsub("\\D", "", basename(df_accuracy$file_name))
  df_accuracy$name <- paste(df_accuracy$...1, cv_number, sep = "_")
  
  # Read gene list filenames and load data
  table_list_genes <- list.files(path = path_output,
                                 recursive = TRUE,
                                 pattern = "^list_genes",
                                 full.names = TRUE)
  dati <- list()
  for (file in table_list_genes) {
    dati[[file]] <- readRDS(file)
  }
  
  # Create dataframes from gene lists
  for (i in seq_along(dati)) {
    current_list <- dati[[i]]
    max_genes <- max(sapply(current_list, length))
    matrix_list <- matrix(NA, nrow = length(current_list), ncol = max_genes)
    rownames(matrix_list) <- names(current_list)
    for (n in 1:length(current_list)) {
      n_genes <- length(current_list[[n]])
      matrix_list[n, 1:n_genes] <- current_list[[n]]
    }
    genes_df <- as.data.frame(matrix_list)
    assign(paste0("listcv_", i), genes_df)
  }
  
  # Combine dataframes
  df_tot <- data.frame()
  for (i in 1:length(dati)) {
    nome_df <- paste0("listcv_", i)
    df_temp <- get(nome_df)
    algorith_names <- paste0(rownames(df_temp), "_", i)
    rownames(df_temp) <- algorith_names
    df_tot <- bind_rows(df_tot, df_temp)
    rownames(df_tot) <- gsub("genes_", "", rownames(df_tot))
    rownames(df_tot) <- gsub("_.*?_", "_", rownames(df_tot))
    df_tot$name <- rownames(df_tot)
    df_merge <- merge(df_tot, df_accuracy, by = "name")
  }
  
  #Filter merged dataframe based on algorithm name
  filtered_df <- df_merge[grep(algorithm_name, df_merge$name), , drop = FALSE]
  
  #Select specified column of statistic result  and convert to numeric
  selected_column <- filtered_df[[column_name]]
  selected_column <- na.omit(selected_column)
  selected_column <- gsub(",", ".", selected_column)
  selected_column <- as.numeric(selected_column)
  
  # Calculate statistic
  if (value_type == "max") {
    selected_value <- (max((selected_column), na.rm = TRUE))
    selected_row <- filtered_df[selected_column == selected_value, ]
  } else if (value_type == "median") {
    mediana <- median(selected_column, na.rm = TRUE)
    selected_value <- median(selected_column, na.rm = TRUE)
    selected_row <- filtered_df[selected_column == selected_value, ]
  } else if (value_type == "mean") {
    mean_value <- mean(selected_column, na.rm = TRUE)
    selected_value <- which.min(abs(selected_column - mean_value))
    selected_row <- filtered_df[selected_value, ]
  }
  
  #Prepare output - Top genes associated with selected algorithm and statistic
  rownames(selected_row) <- selected_row$name
  top_genes <- selected_row[, grep("^V", colnames(selected_row))]
  top_genes <- top_genes[, complete.cases(t(top_genes))]
  if (length(top_genes) == 1) {
    top_genes <- list(top_genes)
  } else if (length(top_genes) > 1) {
    top_genes <- top_genes[1, ]
    top_genes <- list(top_genes)
  }
  
  # Return output as a list
  return(list(top_genes = top_genes, df_tot = df_tot))
  
}


#Example
#result <- extract_algorithm_results(path_output = "../Results/",algorithm_name="svmRadial",
                                   #column_name="Accuracy",
                                   #value_type = "median")