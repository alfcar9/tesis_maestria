############### MAIN ####################################################################

set.seed(130912)
source("libraries.R")
source("required_functions.R")
source("SA_function_presentation.R")

name_instance <- "tsp225.tsp"
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

g0 <- ggplot(cities_pos_df, aes(x = V1, y = V2)) + geom_point() + theme_bw() + geom_text(aes(label=c(1:n)), hjust=0, vjust=0) + xlab("") + ylab("")

assign(name_tour, output_tour)
for(j in 2:15){
  output_tour <- SA_function(immediate_value = ceiling(4*runif(1)), proportion_edges = runif(1,1/20,0.7), initial_solution = output_tour[[1]], distance_type = edge_type)
  name_tour <- paste0("tour", j)
  assign(name_tour, output_tour)
}

tour1 <- tour1[[3]]
paths_print <- tour1

for(j in 2:15){
  name_tour <- paste0("tour", j)
  assign(name_tour, get(name_tour)[[3]])
}

for(j in 2:15){
  name_tour <- paste0("tour", j)
  paths_print <- rbind(paths_print, get(name_tour))
}


# tour_plot1 <- tour1
# #tour_plot2 <- tour15

p <- TSP_plot(cities_pos_df, paths_print[1,], paths_print[1,]) + theme(axis.text.y = element_blank(), axis.text.x = element_blank()) + xlab("") + ylab("")
png(paste0("/home/toto/plot_01.png"), width=600, height=500, res=120)
print(p)
dev.off()

for (i in 2:nrow(paths_print)) {
  p <- TSP_plot(cities_pos_df, paths_print[i,], paths_print[i-1,]) + theme(axis.text.y = element_blank(), axis.text.x = element_blank()) + xlab("") + ylab("")
  if(i < 10)
    png(paste0("/home/toto/plot_0", i, ".png"), width=600, height=500, res=120)
  else
    png(paste0("/home/toto/plot_", i, ".png"), width=600, height=500, res=120)
  print(p)
  dev.off()
}


# make.mov <- function(){
#   unlink("plot.mpg")
#   system("convert -delay 0.5 /home/toto/plot*.png /home/toto/plot.mpg")
# }
