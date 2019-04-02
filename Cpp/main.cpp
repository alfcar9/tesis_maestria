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
	// Determines size of instance
	unsigned n; // Declare size of instance, number of city
	const std::string input = "../Instancias/berlin52.tsp"; // Name of instance file
	n = tsp_functions::finds_number(input); // Gets instance size

	// DATA DECLARATIONS
	unsigned immediate_value;
	double prop_edges;

	double **ppos = (double **)malloc(n * sizeof(double*)); // Matrix that contains coordinates
	for(unsigned i = 0; i < n; i++) ppos[i] = (double *)malloc(2 * sizeof(double));

	double **pdist_complete = (double **)malloc(n * sizeof(double*)); // Complete distance matrix
	for(unsigned i = 0; i < n; i++) pdist_complete[i] = (double *)malloc(n * sizeof(double));

	double **pdist_prune = (double **)malloc(n * sizeof(double*)); // Pruned with infinites distance matrix
	for(unsigned i = 0; i < n; i++) pdist_prune[i] = (double *)malloc(n * sizeof(double));

	unsigned **ppaths = (unsigned **)malloc(n * sizeof(unsigned*)); // Matrix that contains the n paths
	for(unsigned i = 0; i < n; i++) ppaths[i] = (unsigned *)malloc( (n+1) * sizeof(unsigned));

	double **pcosts = (double **)malloc(n * sizeof(double*));; // Costs of edges matrix
	for(unsigned i = 0; i < n; i++) pcosts[i] = (double *)malloc( n * sizeof(double));

	unsigned *pneighbors_degree = (unsigned *)malloc(n * sizeof(unsigned));
	
	unsigned **pneighbors_degree_list = (unsigned **)malloc(n * sizeof(unsigned*));
	for(unsigned i = 0; i < n; i++) pneighbors_degree_list[i] = (unsigned *)malloc( n * sizeof(unsigned));

	// We define a list that has the neighbors of each node. Because each node can have different amount of neighbors a matrix
	// is not an appropriate path.
	vector<unsigned> **pneighbors_list = (vector<unsigned> **)malloc(n * sizeof(unsigned*));
	for(unsigned i = 0;	i<n; i++){ pneighbors_list[i] = new vector<unsigned>; }

	// HYPERPARAMETERS
	prop_edges = 1.0/3.0; // Proportion of neighbors with respect the n-1 neighbors
	immediate_value = 3; // Number of neighbors remaining so it adds to path regardless of the cost

	// BEGINNING OF TSP
	tsp_functions::reads_instance(input, &ppos, n); // Read and saves input file in matrix
	tsp_functions::dist_matrix(ppos, &pdist_complete, n); // Calculates distance matrix

	// Copies and prunes distance matrix with infitiy values
	tsp_functions::prune_matrix(pdist_complete, &pdist_prune, n, prop_edges);

	// It is possible that the resulting matrix is not symmetric because there is no symmetry in the fact that 
	// the kth neighbor v2 for a v1 node might not be the kth neighbor v1 for v2.
	tsp_functions::symmetric_matrix(&pdist_prune, n);

	tsp_functions::initiate_neighborslist(pdist_prune, &pneighbors_degree, &pneighbors_list, &pneighbors_degree_list, n);
	// Initiate matrix of costs and paths
	tsp_functions::initiate_paths(&ppaths, n);
	tsp_functions::initiate_costs(&pcosts, n);

	tsp_functions::greedy_paths(immediate_value, n, ppos, &ppaths, &pcosts, pdist_complete, pneighbors_degree, pneighbors_list, &pneighbors_degree_list);
	
	// We calculate cost of last iteration	
	for(unsigned i = 0; i < n; ++i){
		unsigned city_current;
		city_current = ppaths[i][n-1];
		pcosts[i][n-1] = pdist_complete[city_current - 1][i];	
	}

	// We calculate the best greedy cycle
	double costs_cycles[n];
	for(unsigned i = 0; i < n; ++i){
		costs_cycles[i] = 0.0;
		for(unsigned j = 0; j < n; ++j){
		 	costs_cycles[i] += pcosts[i][j];
		}
	}
	for(unsigned j = 0; j < n; j++){
		printf("%.2f ,", pcosts[1][j]);
	}
	printf("\n");
	for(unsigned j = 0; j < n; j++){
		printf("%u ,", ppaths[1][j]);
	}

	//unsigned cycle_index = std::min_element(costs_cycles, costs_cycles + n) - costs_cycles;
	
	printf("El mejor costo es: %f \n \n", costs_cycles[1]);
	//tsp_functions::print_matrix(ppaths, 1, n+1);
	//tsp_functions::print_matrix(pcosts, cycle_index, n);
	//Se guarda el mejor ciclo y su costo
	
	
}