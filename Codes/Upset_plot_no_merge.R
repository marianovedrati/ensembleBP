
# funzione per creare Upsetplot senza merge quindi senza ragguppare 
# i 30 algoritmi di predizioni in famiglie

 
#Questa funzione crea un dataframe in cui le righe sono gli algoritmi utilizzati
#e le colonne rappresentano i geni unici predetti dai diversi algoritmi
#Ogni cella del dataframe 
#indica se il gene associato è presente tra le diverse famiglie di algoritmo.
# e mi restituisce un Upset plot

#lista geni: insieme di tutte le liste di geni di tutti gli algoritmi


Upset_plot_no_merge <- function(lista_geni) {
  
  # Carica il pacchetto tidyverse
  library(tidyverse)
  
  #Carico il pacchetto Upset
  library(UpSetR)
  
  # Estrai tutti i geni unici di tutti gli algoritmi
  geni_unici <- sort(unique(unlist(lista_geni)))
  
  # Crea un dataframe vuoto
  df <- data.frame(RowName = names(lista_geni))
  
  # Aggiungi colonne per ogni gene unico e assegna TRUE o FALSE in base alla presenza del gene
  for (gene in geni_unici) {
    df[[gene]] <- map_lgl(df$RowName, ~ gene %in% lista_geni[[.x]])
  }
  
  # Imposta i nomi delle righe uguali al nome dell'algoritmo
  rownames(df) <- df$RowName
  df$RowName <- NULL
  df <- as.data.frame(t(df))
  df <- as.data.frame(ifelse(df == "TRUE", 1, 0))
  print(upset(df))
  
}











