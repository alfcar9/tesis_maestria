set.seed(11212)
n <- 10
A <- matrix(rnorm(n^2,15,4),n)
A <- round(100*t(A)*A)

# Son ceros la diagonal
for(i in 1:n){
  A[i,i] <- 0
}

# Se supone la desigualdad del triangulo
for(i in 1:(n-1)){
  for(j in (i+1):n){
    cand <- A[i,j]
    set <- setdiff(1:n, c(i,j))
    for(k in set){
      test <- A[i,k] + A[k,j]
      if(test < cand){
        cand <- test
      }
    }
    A[i,j] <- cand
    A[j,i] <- cand
  }
}

# Estricturas a usar
camino_lista <- list()
costos_lista <- list()
for(i in 1:(n+1)){
  camino_lista[[i]] <- matrix((1:(i*n))*0, nrow = n)
  camino_lista[[i]][1:n,1] <- 1:n
  costos_lista[[i]] <- (1:n)*0
}

camino_lista[[(n+1)]] <- matrix((1:(n*(n+1)))*0, nrow = n)
camino_lista[[(n+1)]][1:n,1] <- 1:n
camino_lista[[(n+1)]][1:n,(n+1)] <- 1:n

#Iteraciones
for(num_iter in 1:(n-1)){
   for(i in 1:n){
     min_camino <- Inf
     ciudades_iter <- setdiff(1:n, camino_lista[[num_iter]][i,])
     ciudad_actual <- camino_lista[[num_iter]][i,num_iter]
     for(j in ciudades_iter){
       cand <- A[ciudad_actual,j]
       if(cand < min_camino){
         min_camino <- cand
         ciudad <- j
       }
     }
     costos_lista[[num_iter]][i] <- min_camino
     for(k in num_iter:(n-1)){
       camino_lista[[k+1]][i, num_iter] <- ciudad
     }
   }
}
