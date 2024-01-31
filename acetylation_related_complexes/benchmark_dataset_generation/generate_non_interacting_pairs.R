library(tidyverse)
library(readxl)
library(openxlsx)

# Load complex components data
complex_components <- readxl::read_xlsx("./complex_components.xlsx") %>% 
  separate_rows(Components, IDs, sep = ", ")

# Load complex pairs data
non_interacting_complex_pairs <- readxl::read_xlsx("./non_interacting_complexes_pairs.xlsx")


generate_random_pairs <- function(pairs_number,
                                  complex_components,
                                  non_interacting_complex_pairs) {
  
  results <- data.frame(ID1 = character(),
                        ID2 = character(),
                        Complex_1 = character(),
                        Complex_2 = character(),
                        stringsAsFactors = FALSE)
  
  for (i in 1:pairs_number) {
    
    is_in_results <- TRUE
      
    while ( is_in_results == TRUE) {
      
      # Select a random pair of non interacting complexes
      random_complex_pair <- non_interacting_complex_pairs %>% sample_n(1)
      
      # Select a random protein as ID1
      comp_1 <- random_complex_pair$Complex_1
      ID_1 <- complex_components %>% filter(Complex_Name == comp_1) %>% sample_n(1)
      
      
      # Select a random protein as ID2
      comp_2 <- random_complex_pair$Complex_2
      ID_2 <- complex_components %>% filter(Complex_Name == comp_2) %>% sample_n(1)
      
      # Random Pair
      pair <- data.frame(ID1 = c(ID_1$IDs),
                         ID2 = c(ID_2$IDs),
                         Complex_1 = comp_1,
                         Complex_2 = comp_2,
                         stringsAsFactors = FALSE)
      
      # Check if pair is in results
      if (any(apply(results, 1, function(x) identical(x, pair)))) {
        
        print("Pair already in results. Trying again.")
        
      } else {
        
        results <- rbind(results, pair)
        is_in_results <- FALSE
        
      }
        
    }
    
  }
  
  return(results)
  
}

pairs_number <- 200

negative_interactors <-  generate_random_pairs(pairs_number,
                                               complex_components,
                                               non_interacting_complex_pairs)

write.xlsx(negative_interactors, "./negative_interactors.xlsx")
