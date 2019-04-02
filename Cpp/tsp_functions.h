#include <vector> // for vectors
#include <bits/stdc++.h>  // for vectors

using namespace std; //shorthand for standard input
class tsp_functions
{
public:
	//Define Euclidean distance funcion
	static double dist_euclidean(double *pcity1, double *pcity2);
	// Given set of coordiantes calculates euclidean distance matrix.
	static void dist_matrix(double **ppos, double ***pdist_complete, unsigned n);
	// Assures the distance matrix symemetric
	static void symmetric_matrix(double ***pdist, unsigned n);
	// Updates distance matrix when triangle Inequality is not satissfied
	static double **triangle_inequality(double **pdist, unsigned n);
	// Returns the index of thee jth ctiy in the ith path 
	static unsigned index_function(unsigned i, unsigned j, unsigned **paths, unsigned n);
	// Reads instance and defines a matrix
	static void reads_instance(const std::string &input, double ***ppos, unsigned n);
	// Extracts size of instance which is in the name of instance of TSPLIB
	static unsigned finds_number(const std::string &input);
	// Prunsigneds Matrix
	static void print_matrix(double **pmatrix, unsigned m, unsigned n);
	// Prunsigneds Matrix
	static void print_matrix(unsigned **pmatrix, unsigned m, unsigned n);
	// Makes inifity the distance of a proportion of neighbors
	static void prune_matrix(double **pdist_complete, double ***pdist_prune, unsigned n, double prop_edges);
	// Initiates paths matrix
	static void initiate_paths(unsigned ***ppaths, unsigned n);
	// Initiates costs matrix
	static void initiate_costs(double ***pcosts, unsigned n);
	// Initates neighbors list
	static void initiate_neighborslist(double **pdist_prune, unsigned **pdegree_neighbors, 
		vector<unsigned> ***pneighbors_list, unsigned ***pdegree_neighbors_list, unsigned n);
	
	static int greedy_paths(unsigned immediate_value, unsigned n, double **ppos, unsigned ***ppaths, double ***pcosts, double **pdist_complete,
		unsigned *pdegree_neighbors, vector<unsigned> **pneighbors_list, unsigned ***pdegree_neighbors_list);

	static void crossing_procedure(unsigned city_current, unsigned city_min, unsigned num_iter, 
		unsigned i, double **ppos, unsigned ***ppaths, double ***pcosts, double **pdist_complete);
};