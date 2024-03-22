
# funzione per fare l'arricchimento con String
# gene_list è una lista di geni
# uso algoritmo edge.betweenees per clusterizzare

## TODO: Aggiungiamo il path del grafico come input della funzione
## TODO: cambiare i commenti in Inglese
## TODO: verificare se bisogna installare le librerie prima di chiamarle
## TODO: verificare che il nome dei geni sia nel formato necessario per STRING, 
##       altrimenti cambiare i nomi!!
analyze_gene_network <- function(gene_list, algorithm = "fastgreedy") {
  
  library(STRINGdb)
  library(ggplot2)
  
  # Change type of gene_list from list to df
  gene_data <- as.data.frame(gene_list)
  
  # Set colnames properly
  colnames(gene_data) <- "Gene"
  
  # Initialize STRING database
  ## TODO: questi parametri devono essere tutti richiesti in input dalla funzione
  ##       è comunque da considerare che non tutti lavorano con gli Umani
  string_db <- STRINGdb$new(version = "12", species = 9606, network_type = "full", input_directory = "")
  
  # Map genes to STRING database and remove unmapped ones
  mapped_genes <- string_db$map(gene_data, "Gene", removeUnmappedRows = TRUE)
  
  # Plot and save STRING Network with link and summary
  ## TODO: il path del salvataggio deve essere chiesto da input 
  plot_network_general_path <- file.path(getwd(), "../Results/plot_network_general.pdf")
  pdf(plot_network_general_path)
  string_db$plot_network(mapped_genes, add_link = TRUE, add_summary = TRUE)
  dev.off()
  
  # Compute GO enrichment of general ppi Network
  enrichmentGO <- string_db$get_enrichment(mapped_genes, category = "Process", iea = TRUE)
  # Compute KEGG enrichment of general ppi Network
  enrichmentKEGG <- string_db$get_enrichment(mapped_genes, category = "KEGG", iea = TRUE)
  
  # Compute Clusters
  ## TODO: questo step va fatto solo se richiesto da input
  clustersList <- string_db$get_clusters(mapped_genes$STRING_id, algorithm = algorithm)
  ## TODO: farsi dare il parametro 5 da input
  clustersList <- clustersList[lengths(clustersList) > 5]
  
  # Set path to save PDF plots of clusters
  ## TODO: il path del salvataggio va chiesto da input
  clusters_plot_path <- file.path(getwd(), "../Results/clusters_plot.pdf")
  
  # Initialize pdf plot of ppi clusters
  pdf(clusters_plot_path)
  
  # Build one plot for each cluster
  for (i in 1:length(clustersList)) {
    string_db$plot_network(clustersList[[i]])
  }
  dev.off()
  
  # Build list containing all results
  result <- list(
    plot_network_general = plot_network_general_path,
    enrichmentGO = enrichmentGO,
    enrichmentKEGG = enrichmentKEGG,
    plot_clusters = clusters_plot_path  # Cambiato per includere il percorso del file PDF dei cluster
  )
  
  return(result)
}


# ## Utilizzo della funzione
# list <- readRDS("../Results/list_genes_2.rds")
# gene_list <- list[[4]]
# #
# ## Applico funzione
# result <-  analyze_gene_network(gene_list)









