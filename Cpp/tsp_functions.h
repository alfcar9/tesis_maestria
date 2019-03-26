#include <vector> // for vectors
#include <bits/stdc++.h>  // for vectors

using namespace std; //shorthand for standard input
class tsp_functions
{
public:
	//Define Euclidean distance funcion
	static double dist_euclidean(double *pcity1, double *pcity2);
	// Given set of coordiantes calculates euclidean distance matrix.
	static void dist_matrix(double **ppos, double ***pdist_complete, int n);
	// Assures the distance matrix symemetric
	static void symmetric_matrix(double **pdist, int n);
	// Updates distance matrix when triangle Inequality is not satissfied
	static double **triangle_inequality(double **pdist, int n);
	// Returns the index of thee jth ctiy in the ith path 
	static int index_function(int i, int j, int **paths, int n);
	// Reads instance and defines a matrix
	static void reads_instance(const std::string &input, double ***ppos, int n);
	// Extracts size of instance which is in the name of instance of TSPLIB
	static int finds_number(const std::string &input);
	// Prints Matrix
	static void print_matrix(double **pmatrix, int m, int n);
	// Prints Matrix
	static void print_matrix(int **pmatrix, int m, int n);
	// Makes inifity the distance of a proportion of neighbors
	static void prune_matrix(double **pdist_complete, double ***pdist_prune, int n, double prop_edges);
	// Initiates paths matrix
	static void initiate_paths(int ***ppaths, int n);
	// Initiates costs matrix
	static void initiate_costs(double ***pcosts, int n);
	// Initates neighbors list
	static void initiate_neighborslist(double **pdist_prune, int **degree_neighbors, vector<int> ***pneighbors_list, int n);
};