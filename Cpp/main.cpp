/**********************************************
* Este programa es la heuristica para TSP *
***********************************************/

#include <stdio.h>  //standard input
#include <math.h> // math functions
#include <iostream>
#include <string> // std::string, std::stod
#include <cassert>
#include <cmath> // floor function
#include "tsp_functions.h"
#include <limits> // define infinity
#include <vector> // for vectors
#include <bits/stdc++.h>  // for vectors

using namespace std; //shorthand for standard input

int main(){
	unsigned n; // Declare size of instance, number of city
	double prop_edges, immediate_value;
	unsigned **ppaths = 0; // Matrix that contains the n paths
	double **ppos = 0; // Matrix that contains coordinates
	double **pdist_complete = 0; // Complete distance matrix
	double **pdist_prune = 0; // Pruned with infinites distance matrix
	double **pcosts = 0; // Costs of edges matrix
	unsigned *pdegree_neighbors=0;
	const std::string input = "../Instancias/berlin52.tsp"; // Name of instance file

	n = tsp_functions::finds_number(input); // Gets instance size
	tsp_functions::reads_instance(input, &ppos, n); // Read and saves input file in matrix
	tsp_functions::dist_matrix(ppos, &pdist_complete, n); // Calculates distance matrix
	
	// HYPERPARAMETERS
	prop_edges = 1.0/3.0; // Proportion of neighbors with respect the n-1 neighbors
	immediate_value = 3; // Number of neighbors remaining so it adds to path regardless of the cost

	// Copies and prunes distance matrix with infitiy values
	tsp_functions::prune_matrix(pdist_complete, &pdist_prune, n, prop_edges);

	// It is possible that the resulting matrix is not symmetric because there is no symmetry in the fact that 
	// the kth neighbor v2 for a v1 node might not be the kth neighbor v1 for v2.
	tsp_functions::symmetric_matrix(pdist_prune, n);

	// We define a list that has the neighbors of each node. Because each node can have different amount of neighbors a matrix
	// is not an appropriate path.
	vector<unsigned> **pneighbors_list = 0;
	unsigned **pdegree_neighbors_list = 0;
	tsp_functions::initiate_neighborslist(pdist_prune, &pdegree_neighbors, &pneighbors_list, &pdegree_neighbors_list, n);
	// Initiate matrix of costs and paths
	tsp_functions::initiate_paths(&ppaths, n);
	tsp_functions::initiate_costs(&pcosts, n);
	// for (unsigned i = 0; i < n; ++i)
	// {
	// 		for (auto j = pdegree_neighbors_list[i] -> begin(); j != pdegree_neighbors_list[i] -> end(); ++j){
	// 			cout << *j << " ";
	// 		}
	// 	printf("\n");	
	// }
	// for (int i = 0; i < 1; ++i)
	// {
	// 	for (auto j = pneighbors_list[i] -> begin(); j != pneighbors_list[i] -> end(); ++j) 
	// 	{
	// 		cout << *j << " ";
	// 	}
	// 	printf("\n");
	// }

	// for (auto j = ppaths[0]; j != ppaths[0] + num_iter; ++j) 
	// {
	// 	cout << *j << " ";
	// }
	// printf("\n");
 //    printf("A estas ciudades se puede ir en la iteracion %u \n", num_iter);
 //    for (auto j = city_neighbors_to_go.begin(); j !=  city_neighbors_to_go.end(); ++j)
 //    {
 //    	printf("%u, ", *j);
 //    }
 //    printf("\n");
	// printf("%u\n", city_current);

	unsigned city_current, number_neighbors, city_min, dist_min_index, candidate_city, candidate_city_degree, *candidate_city_ptr;
	double dist_min;
	bool immediate_value_bool;
	vector<unsigned> city_neighbors_to_go, degree_neighbors_to_go;
	vector<double> city_distances;
	std::vector<bool> immediate_neighbors;
	for(unsigned num_iter = 1; num_iter < n; num_iter++){
		for(unsigned i=0; i < 1; ++i){
		    city_current = ppaths[i][num_iter-1];
	 	    number_neighbors = 0;
		    // Se encuentran las ciudades vecinas que no se han visitado así como su respectivo grado
	 	    city_neighbors_to_go.clear();
	 	    degree_neighbors_to_go.clear();
	 	    city_distances.clear();
	 	    immediate_neighbors.clear();

	 	    for(unsigned j = 0; j < pneighbors_list[city_current-1] -> size(); j++){ // Sobre cada posible vecino de current_city
 	    		candidate_city = *(pneighbors_list[city_current-1] -> begin() + j );
 	    		candidate_city_ptr = find(ppaths[i], ppaths[i] + num_iter -1, candidate_city);
 	    		if( candidate_city_ptr == ppaths[i] + num_iter -1) {
 	    			city_neighbors_to_go.push_back(candidate_city);
 	    			candidate_city_degree = *(pdegree_neighbors_list[i] + j );
 	    			degree_neighbors_to_go.push_back(candidate_city_degree);
 	    			immediate_neighbors.push_back(candidate_city_degree <= immediate_value);
 	    			number_neighbors++;
 	    		}
	 	    }
			// Si no hay ciudades vecinas entonces se hace el calculo sobre todas las posibles ciudades
			if(number_neighbors == 0){
				 for (unsigned j = 0; j < n; j++){
				 	for(unsigned k = 0; k < num_iter; k++){
					 	if( j+1 != ppaths[i][k]){
							city_neighbors_to_go.push_back(j+1);
							number_neighbors++;
						}
					}
				 } 
			}

	 		//printf("%u\n", number_neighbors);	
	 	    // Se guardan en un vector solo los grados de los correspondientes nodos por visitar.
	 	    // Se calcula un valor booleano que previene que nos olvidemos de un nodo. Si existe un nodo con pocos
		    // vecinos hay que visitarlo inmediatamente, sino, conviene visitar el más proximo.

	 		immediate_value_bool = std::none_of(immediate_neighbors.begin(), immediate_neighbors.end(), [](bool value) { return value; });
			
		    if(immediate_value_bool){
		  	 	// En caso de que todos los nodos tengan multiples vecinos calcular el mas proximo
	 	     	for (auto j = city_neighbors_to_go.begin(); j != city_neighbors_to_go.end(); ++j){
	 		   		city_distances.push_back(*(pdist_complete[city_current-1] + *j-1));
	 	     	}
	 	     	dist_min_index = std::min_element(city_distances.begin(), city_distances.end()) - city_distances.begin();
	 	     }
	 	 	else{
	 	    	// En caso de que haya que visitar uno inmediatamente entonces añadirlo
	 	    	dist_min_index = std::min_element(degree_neighbors_to_go.begin(), degree_neighbors_to_go.end()) - degree_neighbors_to_go.begin();
	 	    	}
	 	    // for(auto j = city_neighbors_to_go.begin(); j != city_neighbors_to_go.end(); j++){
	 	    // 	printf("%u ", *j);
	 	    // }
	 	    // printf("\n");
	 	    // for(auto j = city_distances.begin(); j != city_distances.end(); j++){
	 	    // 	printf("%f ", *j);
	 	    // }	
	 	    city_min = city_neighbors_to_go[dist_min_index];
	 	    dist_min = city_distances[dist_min_index];
	 	    
	 	    // Todos los vecinos que tengan asociado a la ciudad que se visita se les resta un grado
		    if(number_neighbors != 0){
		    	for(auto j = city_neighbors_to_go.begin(); j != city_neighbors_to_go.end(); j++){
		    		*(pdegree_neighbors_list[i] + *j) = *(pdegree_neighbors_list[i] + *j) -1;
		    	}
		    }		    
		    // Se actualiza el grado de la ciudad actual a 0.
		    *(pdegree_neighbors_list[i] + city_current - 1 ) = 0;
		    // Se actualiza la matriz de costo.
		    pcosts[i][num_iter-1] = dist_min;
		    ppaths[i][num_iter] = city_min;
		    // // Se verifica que no haya crossing
		    // if(num_iter>2){
			   //  city21=city_current
			   //  x21=city_pos[city21,1]
			   //  y21=city_pos[city21,2]
			   //  city22=city_min
			   //  x22=city_pos[city22,1]
			   //  y22=city_pos[city22,2]
			   //  if(x22-x21 != 0){
			   //  	m2 = (y22-y21)/(x22-x21)
			   //  }
		    // 	else{
		    // 		m2 =(y22-y21)/(x22-x21+10e-5)
		    //   	}
		    //   	b2 = y21-m2*x21
		    //   	for(k in (num_iter-2):1){ // ... las restantes aristas
			   //      city11 = paths[i,k]
			   //      x11 = city_pos[city11,1]
			   //      y11 = city_pos[city11,2]
			   //      city12 = paths[i,k+1]
			   //      x12 = city_pos[city12,1]
			   //      y12 = city_pos[city12,2]
			   //      if(x12-x11 != 0){
			   //        m1 = (y12-y11)/(x12-x11)
			   //      }
		    //     	else{
		    //       		m1 = (y12-y11)/(x12-x11+10e-5)
		    //     	}
		    //     	b1 = y11-m1*x11
		    //     	if(m1 != m2){
		    //       		x_intersect = -(b2-b1)/(m2-m1)
		    //       		bool1 = (x11 < x_intersect & x_intersect < x12) | (x12 < x_intersect & x_intersect < x11)
		    //       		bool2 = (x21 < x_intersect & x_intersect < x22) | (x22 < x_intersect & x_intersect < x21)
			   //      	if(bool1 & bool2){ // condition for finding crossing
				  //           replace_crossing = rev(paths[i,(k+1):num_iter])
				  //           paths[i,(k+1):num_iter] = replace_crossing
				  //           dist1_replace = city_dist_original[city11,city21]
				  //           dist2_replace = city_dist_original[city12,city22]
				  //           costs[i,k:num_iter] = c(dist1_replace, rev(costs[i,(k+1):(num_iter-1)]), dist2_replace)
				  //           city_current = paths[i, num_iter]
				  //           city21=city_current
				  //           x21=city_pos[city21,1]
				  //           y21=city_pos[city21,2]
			   //     		} // end of if crossing exist
		    //    		} // end of if slopes are different
		    // 	} // end of for possible crossing
		    // } //end of if iter > 2
  		} // end of construction of path i
	} // end of iterations

	printf("\n");
	tsp_functions::print_matrix(pcosts, 1, n);
	tsp_functions::print_matrix(ppaths, 1, n+1);
	

}
