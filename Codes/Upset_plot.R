

### @ Funzione per raggruppare ed unire i diverisi metodi in base all'algoritmo utilizzato

#La funzione combine_genes accetta una lista di liste contenente i nomi dei geni associati
#a ciascun metodo e un metodo di combinazione ("union" o "intersection"). 
#Questa funzione crea un dataframe in cui le righe rappresentano i gruppi di algoitmi 
#e le colonne rappresentano i geni unici presenti tra titti i metodi. 
#Ogni cella del dataframe 
#indica se il gene associato è presente tra le diverse famiglie di algoritmo.

#lista geni: insieme di tutte le liste di gen dai diversi metodi
# method: unione o intersezione dei geni tra i diversi algoritmi con stesso approccio

combine_genes <- function(lista_geni, method) {
  
  # Carica il pacchetto tidyverse
  library(tidyverse)
  
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
  
  df$approccio <- rownames(df)
  df$approccio <- as.factor(gsub(".*_(\\w+)_based", "\\1", df$app))
  
  # Raggruppamento e unione o intersezione dei geni
  grouped <- df %>%
    group_by(approccio) %>%
    summarise_all(if (method == "union") any else if (method == "intersection") all else stop("Metodo non valido. Scegli tra 'union' o 'intersection'."))
  
  group_t <- as.data.frame(t(grouped))
  colnames(group_t) <- group_t[1,]
  group_t <- group_t[-1,]
  return(group_t)
}



# @  Funzione per fare UpsetPlot
# il mio df deve avere i nomi delle dievrse famiglie di algoritmi sulle colonne
# mentre come righe i geni derivati da tutti gli algoritmi
# il df da utilizzare è il df risultatante dallea funzione combine_genes




crea_dataframe_upset <- function(df) {
  
  library(UpSetR)
  df <- as.data.frame(ifelse(df == "TRUE", 1, 0))
  print(upset(df))
  
}


