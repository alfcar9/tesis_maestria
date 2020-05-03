##################################################################
################## FUNCTIONS #####################################
##################################################################

# Function for plotting solutions.
# Takes as input the dataframe with the coordinates of the cities and two cicles for comparisson. 
# The second argument, path_heuristic_original, is the solution before the swapping procedure.
TSP_plot <- function(cities_pos_df, path_heuristic, path_heuristic_original){
  n <- nrow(cities_pos_df)
  g <- ggplot(cities_pos_df, aes(x = V1, y = V2)) +
   geom_point() + theme_bw()  #+ geom_text(aes(label=c(1:n)), hjust=0, vjust=0)
  # Allows plotting if the whole cicle is given or if the cicle is not finished. The second case is mainly for debugging purposes. 
  if(all(!path_heuristic==0)){
    m <- n
  }
  else{
    m <- min(which(path_heuristic==0))-2 
  }
  # Detects which edges appear in the original heuristic so they are plotted in black. Otherwise are plotting in red.
  for(i in 1:m){
    edge <- c(path_heuristic[i], path_heuristic[i+1])
    i_vec <- flatten_dbl(cities_pos_df[edge[1],])
    j_vec <- flatten_dbl(cities_pos_df[edge[2],])
    indexes <- which(path_heuristic_original %in% edge) %>% sort()
    
    if( abs(indexes[2] - indexes[1]) == 1 || (length(indexes) == 3 && indexes[3] - indexes[2] == 1 ) ){
      g <- g + geom_segment(x = i_vec[1], xend = j_vec[1], y = i_vec[2], yend = j_vec[2], color = "black")
    }
    else{
      g <- g + geom_segment(x = i_vec[1] , xend = j_vec[1], y = i_vec[2], yend = j_vec[2], color = "red") 
    }
  }
  # Returns a ggplot
  return(g)
}

# Function for extracting number from string.   
# The names of the instances in TSLIB are given with the number of nodes at the end.  
# The input string in the function of TSP is the directory of the instance.   
# This function takes only the number at the end of the name of the instance directory. 
num_extract <- function(string){ 
  s <- str_extract(string, "\\d+$") %>% as.double() 
  return(s) 
} 

# Extracts the distance function of R
dist_extract <- function(string){
  s <- str_match(string, "\\S+$")
  return(s)
}

# Rounds to nearest integer
nint <- function(value){
  if(value >= 0)
    y <- floor(value)
  else
    y <- -floor(-1*value)
  return(y)
}

lat_long <- function(value){
  PI <- 3.141592
  deg <- nint(value)
  min <- value - deg
  y <- PI*(deg + 5.0 * min / 3.0) / 180.0 
  return(y)
}


# Function for calculating the euclidean distance where the inputs are vectors.
dist_euclidean <- function(i, j){
  result <- (i-j)^2 %>% sum() %>% sqrt() %>% round()
  return(result)
}

# Given a matrix of two columnos and n rows, which set the coordinates of cities. This function
# calculates the matrix of distances.
matrix_cost <- function(matrix_pos, n, distance_type){
  matrix_dist <- matrix(integer(n^2), ncol = n)
  if(distance_type == "EUC_2D"){
    for(i in 1:(n-1)){
      i_vec <- matrix_pos[i,]
      for(j in (i+1):n){
        j_vec <- matrix_pos[j,]
        dij <- dist_euclidean(i_vec, j_vec)
        matrix_dist[i,j] <- dij
        matrix_dist[j,i] <- dij
      }
    }
  } else if (distance_type == "GEO"){
    PI <- 3.141592
    RRR <- 6378.388
    for(i in 1:(n-1)){
      xi <- matrix_pos[i,1]
      latitude_i <- lat_long(xi)
      yi <- matrix_pos[i,2]
      longitude_i <- lat_long(yi)
      for(j in (i+1):n){
        xj <- matrix_pos[j,1]
        latitude_j <- lat_long(xj)
        yj <- matrix_pos[j,2]
        longitude_j <- lat_long(yj)
        q1 <- cos( longitude_i - longitude_j )
        q2 <- cos( latitude_i - latitude_j )
        q3 <- cos( latitude_i + latitude_j )
        dij <-  floor( RRR * acos( 0.5*((1.0+q1)*q2 - (1.0 - q1)*q3)) + 1.0)
        matrix_dist[i,j] <- dij
        matrix_dist[j,i] <- dij
      }
    }
  } else{
    break
  }  
  return(matrix_dist)
}

# Makes symmetric a matrix by checking if matrix M[i][j] = M[j][i]. 
# If the values differ then it replaces the greater by the smaller value.
# This function is useful after we take only the first kth nearest neighbors but because there is no symmetry 
# that is, if a is the nearest neighbor of b, not necessarily b is the closest neighbor of a. Thus, we check 
# if M[i,j] = inifity, but M[j,i] = 100, we make M[i,j] = 100. We add an edge instead of removing one.
matrix_symmetric <- function(matrix, n){
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(matrix[i,j] != matrix[j,i]){
        if(matrix[i,j] < matrix[j,i]){
          matrix[j,i] <- matrix[i,j]
        }
        else{
          matrix[i,j] <- matrix[j,i]
        }
      }
    }
  }
  return(matrix)
}

index_function <- function(paths_vec,j){
  return(which(paths_vec==j)[1])
}


# Detects crossing and replaces them.
crossing_procedure <- function(city_current, city_min, num_iter, i, paths_vec, costs_vec, cities_pos, cities_dist_original){
  epsilon <- 1e-4
  city21 <- city_current
  x21 <- cities_pos[city21,1]
  y21 <- cities_pos[city21,2]
  city22 <- city_min
  x22 <- cities_pos[city22,1]
  y22 <- cities_pos[city22,2]
  if(x22 - x21 != 0){
    m2 <- (y22-y21)/(x22-x21)
  }
  else{
    m2 <- Inf
  }
  b2 <- y21-m2*x21
  k <- num_iter-1
  history_cities <- c()
  while(k >= 1){
    city11 <- paths_vec[k]
    x11 <- cities_pos[city11,1]
    y11 <- cities_pos[city11,2]
    city12 <- paths_vec[k+1]
    x12 <- cities_pos[city12,1]
    y12 <- cities_pos[city12,2]
    if(x12-x11 != 0){
      m1 <- (y12-y11)/(x12-x11)
    }
    else{
      m1 <- Inf
    }
    b1 <- y11 - m1 * x11
    
    if(m1 < Inf & m2 < Inf){
      if(m1 != m2){
        x_intersect <- -(b2 - b1)/(m2 - m1)
        bool1 <- (x11 + epsilon < x_intersect & x_intersect < x12 - epsilon) | (x12 + epsilon < x_intersect & x_intersect < x11 - epsilon)
        bool2 <- (x21 + epsilon < x_intersect & x_intersect < x22 - epsilon) | (x22 + epsilon < x_intersect & x_intersect < x21 - epsilon)
      }
      else{
        if(b1 != b2){
          bool1 <- FALSE
          bool2 <- FALSE
        }
        else{
          x11_aux <- min(x11, x12)
          x12_aux <- max(x11, x12)
          x21_aux <- min(x21, x22)
          x22_aux <- max(x21, x22)
          bool1 <- (x11_aux + epsilon < x21_aux & x21_aux < x12_aux - epsilon) | (x11_aux + epsilon < x22_aux & x22_aux < x12_aux - epsilon) | (x11_aux == x21_aux) | (x22_aux == x12_aux) |
            (x22_aux > x12_aux & x21_aux < x11_aux) | (x12_aux > x22_aux & x11_aux < x21_aux)
          bool2 <- TRUE
        }
      }
    }
    else if( m1 == Inf & m2 < Inf ){
      x_intersect <- x11
      y_intersect <- m2 * x_intersect + b2
      bool1 <- (x21 + epsilon < x_intersect & x_intersect < x22 - epsilon) | (x22 + epsilon < x_intersect & x_intersect < x21 - epsilon)
      bool2 <- (y11 + epsilon < y_intersect & y_intersect < y12 - epsilon) | (y12 + epsilon < y_intersect & y_intersect < y11 - epsilon)
    }
    else if( m1 < Inf & m2 == Inf ){
      x_intersect <- x21
      y_intersect <- m1 * x_intersect + b1
      bool1 <- (x11 + epsilon < x_intersect & x_intersect < x12 - epsilon ) | (x12 + epsilon < x_intersect & x_intersect < x11 - epsilon)
      bool2 <- (y21 + epsilon < y_intersect & y_intersect < y22 - epsilon ) | (y22 + epsilon < y_intersect & y_intersect < y21 - epsilon)
      }
    else{
      # Checks if two vertical segments of line intersect
      bool1 <- TRUE
      if(x11 == x21){
        y11_aux <- min(y11, y12)
        y12_aux <- max(y11, y12)
        y21_aux <- min(y21, y22)
        y22_aux <- max(y21, y22)
        bool2 <- (y11_aux + epsilon < y21_aux & y21_aux < y12_aux - epsilon) | (y11_aux + epsilon < y22_aux & y22_aux < y12_aux - epsilon) | (y11_aux == y21_aux) | (y22_aux == y12_aux) |
          (y22_aux > y12_aux & y21_aux < y11_aux) | (y12_aux > y22_aux & y11_aux < y21_aux)
      }
      else{
        bool2 <- FALSE
      }
    }
    if(bool1 & bool2){ # condition for finding crossing
      if(city12 %in% history_cities)
        break
      if(k != num_iter-1){ # We distinguish the case when the swap is between the last edge, and all the other edges.
        replace_crossing <- rev(paths_vec[(k+1):num_iter])
        paths_vec[(k+1):num_iter] <- replace_crossing
        dist1_replace <- cities_dist_original[city11,city21]
        dist2_replace <- cities_dist_original[city12,city22]
        costs_vec[k:num_iter] <- c(dist1_replace, rev(costs_vec[(k+1):(num_iter-1)]), dist2_replace)
        num_iter_aux <- index_function(paths_vec, city11)
        if(num_iter_aux > 2){
          list_crossing <- crossing_procedure(city11, city21, num_iter_aux, i, paths_vec, costs_vec, cities_pos, cities_dist_original) # Recursive crossing
          paths_vec <- list_crossing[[1]]
          costs_vec <- list_crossing[[2]]
        }
        city_current <- paths_vec[num_iter]
        city21 <- city_current
      }
      else{
        replace_crossing <- rev(paths_vec[k:num_iter])
        paths_vec[k:num_iter] <- replace_crossing
        dist1_replace <- cities_dist_original[paths_vec[(num_iter-2)], city12]
        dist2_replace <- cities_dist_original[city11,city22]
        costs_vec[c(num_iter-2,num_iter)] <- c(dist1_replace, dist2_replace)
        city_current <- paths_vec[num_iter]
        city21 <- city_current
      }
      x21 <- cities_pos[city21,1]
      y21 <- cities_pos[city21,2]
      if(x22 - x21 != 0){
        m2 <- (y22-y21)/(x22-x21)
      }
      else{
        m2 <- Inf
      }
      b2 <- y21-m2*x21
      k <- num_iter-1
      history_cities <- c(history_cities, city21)
    } # end of if crossing exist
    k <- k-1
  } # end of while possible crossing
  list_crossing <- list(paths_vec, costs_vec)
}

