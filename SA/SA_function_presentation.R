SA_function <- function(immediate_value = 2, proportion_edges = 0.2,  position_index = 1, initial_solution, distance_type){
  
  paths_print <- c()
  cities_dist_original <- matrix_cost(cities_pos, n, distance_type) # Generate distance matrix
  cities_dist <- cities_dist_original # For last iteration, save original values
  # Eliminate all but the first kth nearest edges for each node.
  total_neighbors <- floor(n*proportion_edges) # Number of the first kth neighbors, the rest edges will get eliminated
  # Eliminate all but the first kth nearest edges for each node.
  for(i in 1:n){
    kth_near <- (cities_dist_original[i,] %>% sort())[total_neighbors]
    for(j in 1:n){
      if(cities_dist_original[i,j] > kth_near){
        cities_dist[i,j] <- Inf
      }
    }
  }
  # It is possible that the resulting matrix is not symmetric because there is no symmetry in the fact that 
  # the kth neighbor v2 for a v1 node might not be the kth neighbor v1 for v2.
  cities_dist <- matrix_symmetric(cities_dist, n)
  # We define a list that has the neighbors of each node. Because each node can have different amount of neighbors a matrix
  # is not an appropriate path.
  neighbors_list <- list()
  for(i in 1:n){
    neighbors_list[[i]] <- which(as.logical((cities_dist[i,] < Inf)*(cities_dist[i,] > 0) %>% as.integer()))
  }
  # Instanciate structures to use throught the algorithm.
  paths <- matrix(integer((n+1)*n), nrow = n) # The Paths starting at each node. In each row we try to form a cycle.
  paths[,c(1,(n+1))] <- 1:n
  costs <- matrix(integer(n^2), nrow = n) # The cost of traveling from [i,j] to [i,j+1]
  # Vector with number of neighbors
  degree_neighbors <- sapply(1:n, function(i){neighbors_list[[i]] %>% length()}) # Number of neighbors of each node
  degree_neighbors_list <- list() # For each cycle we need to keep control of how many neighbors are left. We define a list.
  for(i in 1:n){
    degree_neighbors_list[[i]] <- degree_neighbors
  }
  
  # First we generate a hamiltonian cycle starting at each node.
  m <- n-1 # Si se quiere ver como se ve el camino en la m-ésima iteracion se puede ajustar, sino dejar m = n-1
  for(num_iter in 1:m){
    for(i in 1:n){
      city_current <- paths[i, num_iter]
      # Se calculan las ciudades vecinas que no se han visitado
      cities_neighbors_to_go <- setdiff(neighbors_list[[city_current]], paths[i, 1:num_iter])
      number_neighbors <- length(cities_neighbors_to_go)
      # Si no hay ciudades vecinas entonces se hace el calculo sobre todas las posibles ciudades
      if(number_neighbors == 0){
        cities_neighbors_to_go <- setdiff(1:n, paths[i, 1:num_iter])
      }
      # Se guardan en un vector solo los grados de los correspondientes nodos por visitar
      degree_neighbors_to_go <- degree_neighbors_list[[i]][cities_neighbors_to_go]
      # Se calcula un valor booleano que previene que nos olvidemos de un nodo. Si existe un nodo con pocos
      # vecinos hay que visitarlo inmediatamente, sino, conviene visitar el más proximo.
      immediate_neighbors <- degree_neighbors_to_go <= immediate_value
      if(all(!immediate_neighbors)){
        # En caso de que todos los nodos tengan multiples vecinos calcular el mas proximo
        city_distances <- cities_dist_original[city_current, cities_neighbors_to_go]
        dist_min <- min(city_distances)
        city_min <- cities_neighbors_to_go[which.min(city_distances)]
      }
      else{
        # En caso de que haya que visitar uno inmediatamente entonces añadirlo
        city_min <- cities_neighbors_to_go[which.min(degree_neighbors_to_go)]
        dist_min <- cities_dist_original[city_current, city_min]
      }
      # Todos los vecinos que tengan asociado a la ciudad que se visita se les resta un grado
      if(number_neighbors!=0){
        degree_neighbors_list[[i]][cities_neighbors_to_go] <- degree_neighbors_list[[i]][cities_neighbors_to_go] - 1
      }
      # Se actualiza el grado de la ciudad actual a 0.
      degree_neighbors_list[[i]][city_current] <- 0
      # Se actualiza la matriz de costo.
      costs[i, num_iter] <- dist_min
      paths[i, num_iter+1] <- city_min
      # Se verifica que no haya crossing
      if(num_iter>2){
        list_crossing <- crossing_procedure(city_current, city_min, num_iter, i, paths[i,], costs[i,], cities_pos, cities_dist_original)
        paths[i,] <- list_crossing[[1]]
        costs[i,] <- list_crossing[[2]]
      } #end of if iter > 2
    } # end of construction of path i
  } # end of iterations
  # Se verifica que la ultima arista no hace crossing
  if(num_iter==(n-1)){
    num_iter <- num_iter+1
    for(i in 1:n){
      city_current <- paths[i, n]
      list_crossing <- crossing_procedure(city_current, i, num_iter, i, paths[i,], costs[i,], cities_pos, cities_dist_original)
      paths[i,] <- list_crossing[[1]]
      costs[i,] <- list_crossing[[2]]
    }
  } #end of if iter > 2
  # Se hace el costo de cerrar el ciclo.
  if(m == n-1){
    costs[1:n, n]  <- sapply(1:n, function(i){cities_dist_original[paths[i, n], i]})
  }
  # The cost of each cycle is calculated
  costs_cycles <- sapply(1:n, function(i){sum(costs[i,])})
  
  # In case we do not want to select the best cycle we can choose second or third best.
  costs_cycles_aux <- costs_cycles
  cycle_indexes <- c() 
  for(i in 1:position_index){
    cycle_index <- which.min(costs_cycles_aux)
    cycle_indexes[i] <- cycle_index
    costs_cycles_aux[cycle_index] <- Inf
  }

  if(missing(initial_solution)){
    # Es el indice del mejor ciclo
    cycle_index <- cycle_indexes[position_index]
    # Se guarda el mejor ciclo y su costo
    path_min_noswap <- paths[cycle_index,]
    cost_min_noswap <- costs[cycle_index,]
    total_cost_min_noswap <- sum(cost_min_noswap)
  }
  else{
    path_min_noswap <- initial_solution
    cost_min_noswap <- c()
    for(i in 1:n){
      cost_min_noswap[i] <- cities_dist_original[path_min_noswap[i], path_min_noswap[i+1]]
    }
    total_cost_min_noswap <- sum(cost_min_noswap)
    if(total_cost_min_noswap > min(costs_cycles)){
      cycle_index <- which.min(costs_cycles)
      path_min_noswap <- paths[cycle_index,]
      cost_min_noswap <- costs[cycle_index,]
      total_cost_min_noswap <- sum(cost_min_noswap)
    }
  }

  # Se inicializan el camino que ha de swapearse.
  path_min_swap <- path_min_noswap
  cost_min_swap <- cost_min_noswap

  # Matriz que a cada fila corresponde un ciclo. Cada celda guarda que numero de visita es cada ciudad. Si en [1,2] = 13, entonces la 13ava ciudad visitada
  # en el ciclo 1 es la ciudad 2.
  index_matrix <- matrix(0L, nrow = n, ncol = n)
  # Se generan los caminos extendidos para poder considerar en una sola matriz todos los caminos posibles leidos como secuencia
  paths_ext <- matrix(0L, nrow = n, ncol = 2*n)
  # Lo análogo con la funcion de costos
  costs_ext <- matrix(0L, nrow = n, ncol = 2*n-1)
  # Se generan estas estructuras
  for(i in 1:n){
    index_matrix[i,] <- sapply(1:n, function(j) index_function(paths[i,], j))
    paths_ext[i,] <-  c(paths[i,], paths[i,2:n])
    costs_ext[i,] <-  c(costs[i,], costs[i, 1:(n-1)])
  }

  path_min_swap_ext <- c(path_min_swap, path_min_swap[2:n])
  cost_min_swap_ext <- c(cost_min_swap, cost_min_swap[1:(n-1)])
  if(missing(initial_solution)){
    index_min_swap <- index_matrix[cycle_index,]
  }
  else{
    index_min_swap <- sapply(1:n, function(j) index_function(path_min_swap, j))
  }

  success <- 1 # Indica que se encontró alguna estructura que swapear, se define 1 para que entre al while al menos una vez
  while(success>0){
    success <- 0  # no se han encontrado estructuras swapeables
    path_list <- list() # se guardan las estructuras swapeables
    path_index_list <- list()  # se guarda donde inicia la estrucutra
    path_cost_vector <- c() # se guarda el costo total de cada estructura
    for(len in 3:(n-1)){  # se fija una longitud
      ini <- 1            # comenzando en el primer nodo
      end <- len+1        # hasta el nodo len + 1 se hace lo siguiente
      for(j in 1:n){   # para cada longitud existen siempre n paths posibles por revisar en el ciclo inicial
        path <- path_min_swap_ext[ini:end] # se guarda la secuencia de ciudades
        path_cost <- sum(cost_min_swap_ext[ini:(end-1)]) # el costo de las transiciones
        path_len <- len + 1 # la longitud del vector que esta dada por len + 1
        extreme1 <- path[1] # se guardan los extremos, la cabeza
        extreme2 <- tail(path, 1) # y la cola
        for(k in 1:n){# vamos a buscar este camino en cada uno de los n-1 caminos generados
          index_begin <- index_matrix[k, extreme1] # obtenemos el indice donde se visita al extremo 1
          index_end <- index_matrix[k, extreme2] # y el indice donde se visita al extremo 2
          if(index_begin > index_end){ # si los indices estan al reves se intercambian
            index_aux <- index_begin
            index_begin <- index_end
            index_end <- index_aux
          }
          path1 <- paths_ext[k, index_begin:index_end] # se busca usando los extremos esa estructura en el ciclo k
          path1_len <- length(path1) # se toma su longitud
          path2 <- paths_ext[k, index_end:(index_begin+n)] # se considera el camino complemento del primero
          path2_len <- length(path2) # se toma su longitud
          if(path_len==path1_len && length(setdiff(path, path1))==0){ # si tienen longitudes iguales y tienen las mismas ciudades
            path1_cost <- sum(costs_ext[k, index_begin:(index_end-1)]) # considera el costo
            if(path1_cost < path_cost){ # entonces si es mas barato es posible hacer un swap
              if(length(path_list)>0){ # evitamos ser redundantes en los caminos reemplazados
                not_in_list <- !any(sapply(1:success, function(i) length(setdiff(path_list[[i]], path1))==0 || length(setdiff(path1, path_list[[i]])) == 0))
              }
              else{
                not_in_list <- TRUE
              }
              if(not_in_list){
                success <- success + 1 # se incrementa success
                path_list[[success]] <- path1 # se agrega a la lista
                path_cost_vector[success] <- path1_cost - path_cost # agrega la utilidad
                path_index_list[[success]] <- c(k, index_begin, index_end) # informacion de donde empieza y donde acaba
              }
            }
          }
          if(path_len==path2_len && length(setdiff(path, path2))==0){ # se hace lo análogo con el camino complemento
            path2_cost <- sum(costs_ext[k,index_end:(index_begin+n-1)])
            if(path2_cost < path_cost){
              if(length(path_list)>0){
                not_in_list <- !any(sapply(1:success, function(i) length(setdiff(path_list[[i]], path2))==0 || length(setdiff(path2, path_list[[i]])) == 0))
              }
              else{
                not_in_list <- TRUE
              }
              if(not_in_list){
                success <- success + 1
                path_list[[success]] <- path2
                path_cost_vector[success] <- path2_cost - path_cost
                path_index_list[[success]] <- c(k, index_end, index_begin+n)
              } # termina if not_in_list
            } # termina if que el costo es menor
          } # termina if de que se contienen las mismas ciudades
        } # termina for de busqueda de ese camino especifico
        # se incrementan el incio y fin y se revisan asi los n posibles caminos
        ini <- ini + 1
        end <- end + 1
      }
    }
    if(success>0){
      m <- 0
      path_cost_vector_aux <- path_cost_vector
      path_index_vector_refined <- c()
      path_list_refined <- list()
      for(k in 1:success){
        index <- which.min(path_cost_vector_aux)
        path_cand <- path_list[[index]]
        if(m > 0){
          k1 <- sum(sapply(1:m, function(i) length(intersect(path_list_refined[[i]], path_cand)))) # paths ajenos
          if(k1==0){
            m <-  m+1
            path_list_refined[[m]] <- path_cand
            path_index_vector_refined[m] <- index
          }
        }
        else{
          m <- m+1
          path_list_refined[[m]] <- path_cand
          path_index_vector_refined[m] <- index
        }
        path_cost_vector_aux <- path_cost_vector_aux[-index]
      }
      for(k in 1:m){
        index_min <- path_index_vector_refined[k] # el numero de success que sera reemplazado
        index_replace <- path_index_list[[index_min]][1] # corresponde al ciclo numero k
        index_begin <- path_index_list[[index_min]][2] # que comienza en la ciudad index_begin
        index_end <- path_index_list[[index_min]][3] # termina en la ciudad index_end
        path_replace <- path_list[[index_min]] # la estructura a remplazarse
        costs_replace <- costs_ext[index_replace, index_begin:(index_end-1)]
        index_head <- index_min_swap[path_replace[1]] # los indices head y tail del ciclo inicial
        index_tail <- index_min_swap[tail(path_replace, 1)] # que va a reemplazarse
        length_half_boolean <- length(path_replace) == n/2+1
        if(length_half_boolean){
          if(index_head < index_tail)
            next_node <- path_min_swap[index_head + 1]
          else
            next_node <- path_min_swap[index_tail + 1]
          if(next_node %in% path_replace)
            no_cut_boolean <- TRUE
          else
            no_cut_boolean <- FALSE
        }
        else{
          lengths_equal_boolean <- length(index_head:index_tail) == length(path_replace)
          if(lengths_equal_boolean)
            no_cut_boolean <- TRUE
          else
            no_cut_boolean <- FALSE
        }
        if( no_cut_boolean ){
          # En este caso se puede sustituir la estructura en el camino sin necesidad de cortarlo
          # Puede ser que haya que ponerlo tal y como esta o en reversa
          if(index_head > index_tail){
            path_replace <- rev(path_replace)
            costs_replace <- rev(costs_replace)
            aux <- index_head
            index_head <- index_tail
            index_tail <- aux
          }
          path_min_swap[index_head:index_tail] <- path_replace
          cost_min_swap[index_head:(index_tail-1)] <- costs_replace
          paths_print <- rbind(paths_print, path_min_swap)
        }
        else{ # En este caso es necesario cortar el camino en partes para el reemplazo
          if(index_head < index_tail){ #nuevamente,pero si head es menor que tail
            path_replace <- rev(path_replace)
            costs_replace <- rev(costs_replace)
            aux <- index_head
            index_head <- index_tail
            index_tail <- aux
          }
          index_cut <- which(path_replace == path_min_swap[1]) #encuentra donde hay que hacer el corte de la estructura
          path_head <- path_replace[1:index_cut] # se parte en dos, head
          path_tail <- path_replace[-1:-(index_cut-1)] # y tail
          costs_head <- costs_replace[1:(index_cut-1)]
          costs_tail <- costs_replace[-1:-(index_cut-1)]
          index_begin <- index_tail + 1 # El primer indice que no es de la estructura
          index_end <- index_head-1 # El ultimo indice que no es de la estructura
          path_min_swap <- c(path_tail, path_min_swap[index_begin:index_end], path_head) # se reemplaza el camino
          cost_min_swap <- c(costs_tail, cost_min_swap[(index_begin-1):index_end], costs_head) # y el costo
          paths_print <- rbind(paths_print, path_min_swap)
        }
        index_min_swap <- sapply(1:n, function(j) {which(path_min_swap==j)[1]})
        path_min_swap_ext <-  c(path_min_swap, path_min_swap[2:n])
        cost_min_swap_ext <-  c(cost_min_swap, cost_min_swap[1:(n-1)])
      }
    }
  }
  total_cost_min_swap <- sum(cost_min_swap)
  list_output <- list(path_min_swap, total_cost_min_swap, paths_print)
  print(total_cost_min_swap)
  return(list_output)
}
