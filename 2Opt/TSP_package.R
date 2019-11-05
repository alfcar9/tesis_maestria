library(TSP)
library(tidyverse)
set.seed(130912)
num_extract <- function(string){ 
  s <- str_extract(string, "\\d+$") %>% as.double() 
  return(s)	
}	

dist_extract <- function(string){
  s <- str_match(string, "\\S+$")
  return(s)
}

name_instance <- "ulysses22.tsp"
name_instance_path <- paste0("../TSPLIB/TSPLIB_original/", name_instance)
instance_vector <- readLines(name_instance_path)[4:5]
n <- num_extract(instance_vector[1])
edge_type <- dist_extract(instance_vector[2])
if(edge_type == "EUC_2D"){
  cities_pos_df  <- read_table2(name_instance_path, skip = 6, col_names = c("V1", "V2"), n_max = n, col_types = cols("-", "d", "d")) %>% as.data.frame() %>% distinct()
  n <- nrow(cities_pos_df)
} else if (edge_type == "GEO"){
  cities_pos_df  <- read_table2(name_instance_path, skip = 7, col_names = c("V1", "V2"), n_max = n, col_types = cols("-", "d", "d")) %>% as.data.frame() %>% distinct()
  n <- nrow(cities_pos_df)
} else{
  cities_pos_df <- NA
}

cities_pos <- as.matrix(cities_pos_df)
etsp <- ETSP(cities_pos_df)
tour <- solve_TSP(etsp)
plot(etsp, tour, tour_col = "red")
tour
