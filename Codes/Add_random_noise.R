setwd("/Users/giorgiomontesi/Desktop/Universita_di_Siena/A_PhD_Project/Biomarker_Prediction/ensembleBP/Codes")

add_random_noise <- function(count_table, sd_noise) {
  # Ottieni le dimensioni della count table
  num_rows <- nrow(count_table)
  num_cols <- ncol(count_table)
  
  # Genera rumore casuale da una distribuzione normale con deviazione standard specificata
  noise <- matrix(round(rnorm(num_rows * num_cols, mean = 0, sd = sd_noise)), nrow = num_rows, ncol = num_cols)
  
  # Aggiungi il rumore casuale alla count table originale
  noisy_count_table <- count_table + noise
  
  # Assicurati che i valori non diventino negativi
  noisy_count_table[noisy_count_table < 0] <- 0
  
  return(noisy_count_table)
}

# Utilizzo della funzione
# Supponiamo che "original_count_table" sia la tua count table originale
# e "sd_noise" sia la deviazione standard del rumore da aggiungere
noisy_count_table <- add_random_noise(original_count_table, sd_noise = 0.1)




add_proportional_noise <- function(count_table, sd_noise) {
  # Ottieni le dimensioni della count table
  num_rows <- nrow(count_table)
  num_cols <- ncol(count_table)
  
  # Genera rumore proporzionale da una distribuzione normale con deviazione standard specificata
  noise <- matrix(round(rnorm(num_rows * num_cols, mean = 1, sd = sd_noise)), nrow = num_rows, ncol = num_cols)
  
  # Moltiplica ciascun valore per un fattore casuale
  noisy_count_table <- count_table * noise
  
  # Arrotonda i valori e assicurati che non diventino negativi
  noisy_count_table <- round(noisy_count_table)
  noisy_count_table[noisy_count_table < 0] <- 0
  
  return(noisy_count_table)
}

# Utilizzo della funzione
# Supponiamo che "original_count_table" sia la tua count table originale
# e "sd_noise" sia la deviazione standard del rumore da aggiungere
noisy_count_table <- add_proportional_noise(original_count_table, sd_noise = 0.1)



