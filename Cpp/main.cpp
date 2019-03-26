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
	int n; // Declare size of instance, number of cities
	double prop_edges; //, immediate;
	int **ppaths = 0; // Matrix that contains the n paths
	double **ppos = 0; // Matrix that contains coordinates
	double **pdist_complete = 0; // Complete distance matrix
	double **pdist_prune = 0; // Pruned with infinites distance matrix
	double **pcosts = 0; // Costs of edges matrix
	int *degree_neighbors=0;
	const std::string input = "../Instancias/berlin52.tsp"; // Name of instance file

	n = tsp_functions::finds_number(input); // Gets instance size
	tsp_functions::reads_instance(input, &ppos, n); // Read and saves input file in matrix
	tsp_functions::dist_matrix(ppos, &pdist_complete, n); // Calculates distance matrix
	
	// HYPERPARAMETERS
	prop_edges = 1.0/3.0; // Proportion of neighbors with respect the n-1 neighbors
	//immediate = 3; // Number of neighbors remaining so it adds to path regardless of the cost

	// Copies and prunes distance matrix with infitiy values
	tsp_functions::prune_matrix(pdist_complete, &pdist_prune, n, prop_edges);

	// It is possible that the resulting matrix is not symmetric because there is no symmetry in the fact that 
	// the kth neighbor v2 for a v1 node might not be the kth neighbor v1 for v2.
	tsp_functions::symmetric_matrix(pdist_prune, n);

	// We define a list that has the neighbors of each node. Because each node can have different amount of neighbors a matrix
	// is not an appropriate path.
	vector<int> **pneighbors_list = 0;
	tsp_functions::initiate_neighborslist(pdist_prune, &degree_neighbors, &pneighbors_list, n);
	// Initiate matrix of costs and paths
	tsp_functions::initiate_paths(&ppaths, n);
	tsp_functions::initiate_costs(&pcosts, n);
	for (int i = 0; i < n; ++i)
	{
		printf("%i, ", degree_neighbors[i]);
	}
	printf("\n");

	for (int i = 0; i < 2; ++i)
	{
		for (auto j = pneighbors_list[i] -> begin(); j != pneighbors_list[i] -> end(); ++j) 
		{
			cout << *j << " ";
		}
		printf("\n");
	}
	
	//tsp_functions::print_matrix(pcosts, 1, n);
	//tsp_functions::print_matrix(ppaths, 1, 1);
	

}
