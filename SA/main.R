############### MAIN ####################################################################

set.seed(130912)
source("libraries.R")
source("required_functions.R")
source("SA_function.R")

name_instance <- "kroA100.tsp"
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
start.time <- Sys.time()
output_tour <- SA_function(position_index = 1, immediate_value = 3, proportion_edges = 0.05, distance_type = edge_type)
name_tour <- paste0("tour", 1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

assign(name_tour, output_tour)
for(j in 2:10){
  output_tour <- SA_function(immediate_value = ceiling(4*runif(1)), proportion_edges = runif(1,1/20,0.7), initial_solution = output_tour[[1]], distance_type = edge_type)
  name_tour <- paste0("tour", j)
  assign(name_tour, output_tour)
}

tour_plot1 <- tour1
tour_plot2 <- tour10

g1 <- TSP_plot(cities_pos_df, tour_plot2[[1]], tour_plot1[[1]]) #+ ggtitle("Iter 1") + labs (x = round(tour_plot2[[2]], 2), y = "")
g1
