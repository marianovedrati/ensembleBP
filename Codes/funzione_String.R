
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
  
  # Trasformo la lista di geni in un dataframe
  gene_data <- as.data.frame(gene_list)
  
  # Imposta i nomi delle colonne
  colnames(gene_data) <- "Gene"
  
  # Inizializza il database STRING
  ## TODO: questi parametri devono essere tutti richiesti in input dalla funzione
  string_db <- STRINGdb$new(version = "12", species = 9606, network_type = "full", input_directory = "")
  
  # Mappa i geni al database STRING e rimuovi quelli non mappati
  mapped_genes <- string_db$map(gene_data, "Gene", removeUnmappedRows = TRUE)
  
  # Plot della rete STRING con link e summary
  # Salvataggio del grafico della rete STRING 
  plot_network_general_path <- file.path(getwd(), "../Results/plot_network_general.pdf")
  pdf(plot_network_general_path)
  string_db$plot_network(mapped_genes, add_link = TRUE, add_summary = TRUE)
  dev.off()
  
  # Calcola l'arricchimento nelle annotazioni GO
  enrichmentGO <- string_db$get_enrichment(mapped_genes, category = "Process", iea = TRUE)
  # Calcola l'arricchimento nelle annotazioni KEGG
  enrichmentKEGG <- string_db$get_enrichment(mapped_genes, category = "KEGG", iea = TRUE)
  
  # Ottieni i clusters
  clustersList <- string_db$get_clusters(mapped_genes$STRING_id, algorithm = algorithm)
  ## TODO: farsi dare il parametro 3 da input
  clustersList <- clustersList[lengths(clustersList) > 5]
  
  # Imposta il percorso del file PDF per i plot dei cluster
  clusters_plot_path <- file.path(getwd(), "../Results/clusters_plot.pdf")
  
  # Inizializza il dispositivo grafico PDF per i plot dei cluster
  pdf(clusters_plot_path)
  
  # Crea i plot per tutti i clusters
  for (i in 1:length(clustersList)) {
    string_db$plot_network(clustersList[[i]])
  }
  
  # Chiudi il dispositivo grafico PDF per i plot dei cluster
  dev.off()
  
  # Restituzione dei risultati
  result <- list(
    plot_network_general = plot_network_general_path,
    enrichmentGO = enrichmentGO,
    enrichmentKEGG = enrichmentKEGG,
    plot_clusters = clusters_plot_path  # Cambiato per includere il percorso del file PDF dei cluster
  )
  
  return(result)
}


# # Utilizzo della funzione
list <- readRDS("../Results/list_genes_2.rds")
gene_list <- list[[4]]
#
# # Applico funzione
result <-  analyze_gene_network(gene_list)



result$enrichmentKEGG$description
result$enrichmentGO$description






