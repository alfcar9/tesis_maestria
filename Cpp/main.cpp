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

using namespace std; //shorthand for standard input

int main(){
	int n, nneighbors, kth_near; //declare size of instance, number of cities
	double prop_edges, immediate;
	double **pmatrix_pos=0; //declare and initialize pointer
	double **pmatrix_dist_comp=0; //declare and initialize pointer
	double **pmatrix_dist_trun=0; //declare and initialize pointer
	const std::string input = "../Instancias/berlin52.tsp"; // Name of instance file

	n = tsp_functions::finds_number(input); // Gets instance size  
	pmatrix_pos = tsp_functions::reads_instance(input, n); // Read and saves input file in matrix
	pmatrix_dist_comp = tsp_functions::matrix_euclidean(pmatrix_pos, n); // Calculates distance matrix
	tsp_functions::print_matrix(pmatrix_pos, 3, 2); // prints matrix
	tsp_functions::print_matrix(pmatrix_dist_comp, 3, 3); // prints matrix

	// HYPERPARAMETERS
	prop_edges = 1/3; // Proportion of neighbors with respect the n-1 neighbors
	immediate = 3; // Number of neighbors remaining so it adds to path regardless of the cost

	// Copies matrix
	pmatrix_dist_trun = tsp_functions::copy_matrix(pmatrix_dist_comp, n, n);

	// Eliminate all but the first kth nearest edges for each node.
	nneighbors = floor(n*prop_edges); // Number of neighbors we keep
	for(int i; i<n; i++){
		kth_near <- 5; //(cities_dist_original[i,] %>% sort())[total_neighbors]
		for(int j; j<n; j++){
			if(pmatrix_dist_comp[i][j] > kth_near){
				pmatrix_dist_trun[i][j] = numeric_limits<double>::infinity();
			}
		}
	}
}
