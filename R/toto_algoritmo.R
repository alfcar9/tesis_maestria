set.seed(11212)

#######################################################################################
# Parameters
n <- 100 # Number of Nodes in Plot
number_neighbors <- floor(n/2) # Number of the first kth neighbors, the rest edges will get eliminated
#######################################################################################

# Functions
# Define Euclidean distance funcion
dist_euclidean <- function(i1, i2, j1,j2){ 
  return(sqrt((i1-j1)^2+(i2-j2)^2))
}

# Given set of coordiantes calculates euclidean distance and adds noice.
matrix_euclidean <- function(ciudades_pos){
  ciudades_dist <- matrix(integer(n^2), ncol = n)
  for(i in 1:(n-1)){
    i1 <- ciudades_pos[i,1] 
    i2 <- ciudades_pos[i,2]
    for(j in (i+1):n){
      j1 <- ciudades_pos[j,1] 
      j2 <- ciudades_pos[j,2] 
      ciudades_dist[i,j] <- (dist_euclidean(i1, i2, j1, j2) + ifelse(i==j,0, (rnorm(1, sd = 0)^2))) %>% round()
      ciudades_dist[j,i] <- ciudades_dist[i,j]
    }
  }
  return(ciudades_dist)
}

# Updates distance matrix when triangle Inequality is not satissfied
desigualdad_triangulo <- function(ciudades_dist)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      distancia_minima <- ciudades_dist[i,j]
      set <- setdiff(1:n, c(i,j))
      for(k in set){
        distancia_indirecta <- ciudades_dist[i,k] + ciudades_dist[k,j]
        if(distancia_indirecta < distancia_minima){
          distancia_minima <- distancia_minima
        }
      }
      ciudades_dist[i,j] <- distancia_minima
      ciudades_dist[j,i] <- distancia_minima
    }
    return(ciudades_dist)
  }

###################################################################################
# Main

ciudades_pos <- matrix(runif(2*n), ncol = 2) # Generate coordinates of n neighbors randomly
ciudades_pos <- round(100*(ciudades_pos)) # Round to integers

#Generate distance matrix
ciudades_dist <- matrix_euclidean(ciudades_pos)

# For last iteration, save original values
ciudades_dist_sin_mod <- ciudades_dist

ciudades_dist <- desigualdad_triangulo(ciudades_dist)

# Eliminate all but the first kth nearest edges for each node.
vecinos_lista <- list()
for(i in 1:n){
  vecinos_lista[[i]] <- c(0)
  kth_near <- (ciudades_dist[i,] %>% sort())[number_neighbors]
  for(j in 1:n){
    if(ciudades_dist[i,j] > kth_near){
      ciudades_dist[i,j] <- Inf
      ciudades_dist[j,i] <- Inf
    }
    else{
      if(i!=j){
      vecinos_lista[[i]] <- c(vecinos_lista[[i]],j)
      }
    }
  }
  vecinos_lista[[i]] <- tail(vecinos_lista[[i]],-1)
}
 
# Estricturas a usar
caminos <- matrix(integer((n+1)*n), nrow = n)
caminos[1:n,1] <- 1:n
caminos[1:n,(n+1)] <- 1:n
costos <- matrix(integer(n^2), nrow = n)

# Matriz con numero de vecinos
ciudades_dist>0

#Iteraciones
for(num_iter in 1:(n-1)){
  for(i in 1:n){
    camino_min <- Inf
    ciudades_iter <- setdiff(1:n, caminos[i,1:num_iter])
    ciudad_actual <- caminos[i, num_iter]
    for(j in ciudades_iter){
        ciudad_cand <- ciudades_dist[ciudad_actual, j]
        if(ciudad_cand < camino_min){
          camino_min <- ciudad_cand
          ciudad_min <- j
        }
    }
      costos[i, num_iter] <- camino_min
      caminos[i, num_iter+1] <- ciudad_min
   }
}

costos[1:n, n]  <- sapply(1:n, function(i){ciudades_dist_sin_mod[i, caminos[i,n]]})
costos_ciclos <- sapply(1:n, function(i){sum(costos[i,1:n])})
ciclo_index <- which.min(costos_ciclos)
camino_heuritica <- caminos[ciclo_index,]


ciudades_pos_df <- as.data.frame(ciudades_pos)
g1 <- ggplot(ciudades_pos_df, aes(x = V1, y = V2)) + geom_point() + theme_bw()
for(i in 1:n){
  i1 <- ciudades_pos[camino_heuritica[i],1]
  i2 <- ciudades_pos[camino_heuritica[i],2]
  j1 <- ciudades_pos[camino_heuritica[i+1],1]
  j2 <- ciudades_pos[camino_heuritica[i+1],2]
  g1 <- g1 + geom_segment(x = i1 , xend = j1, y = i2, yend = j2)
}
g1
costos_ciclos[ciclo_index]