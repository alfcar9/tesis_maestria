/**********************************************
* FUNCTIONS
***********************************************/

#include <iostream> // header in standard library
#include <math.h>
#include "tsp_functions.h" //header in local directory
#include <regex> // Regular Expresions
#include <string> // std::string, std::stod
#include <fstream> // read files
#include <sstream>
#include <bits/stdc++.h> // sort method
#include <vector> // neighbors list changes size in runtime
using namespace std; //shorthand for standard input

double tsp_functions::dist_euclidean(double *pcity1, double *pcity2){
  double result;
  result = 0.0;
  for(unsigned i=0; i<2; i++){
      result += pow(pcity1[i]-pcity2[i], 2);
  }
  result = sqrt(result);
  return result;
}

// Given set of coordiantes calculates euclidean distance matrix.
void tsp_functions::dist_matrix(double **ppos, double ***pdist_complete, unsigned n){
	double *pcity1_pos=0;
	double *pcity2_pos=0;
	pcity1_pos = new double[2];
	pcity2_pos = new double[2];

	double **pmatrix=0; 
	pmatrix = new double*[n];
	for(unsigned i=0; i<n; i++){
	    pmatrix[i] = new double[n];
	}
	for(unsigned i=0; i<(n-1); i++){
		pcity1_pos = &ppos[i][0];
		for(unsigned j=(i+1); j<n; j++){
			for(unsigned k=0; k<2;k++){
			pcity2_pos = &ppos[j][0];
		}
			pmatrix[i][j] = dist_euclidean(pcity1_pos, pcity2_pos);
			pmatrix[j][i] = pmatrix[i][j]; 
		}
	}
	*pdist_complete = pmatrix;
}

// Makes a matrix symemetric
void tsp_functions::symmetric_matrix(double **pdist, unsigned n){
  for(unsigned i = 0; i<n-1; i++){
    for(unsigned j = i+1; j<n; j++){
        if(pdist[i][j] != pdist[j][i]){
          if(pdist[i][j] < pdist[j][i]){
            pdist[j][i] = pdist[i][j];
          }
          else{
            pdist[i][j] = pdist[j][i];
          }
        }
     }
  }
}

// Returns the index of thee jth ctiy in the ith path 
unsigned tsp_functions::index_function(unsigned i, unsigned j, unsigned **paths, unsigned n){
	unsigned index_ans, index;
	index = 0;
	index_ans = -1;
	while(index<n){
		if(paths[i][index] == j){
			index_ans = index;	
		}
		index++;
	}
	return index_ans;
}

unsigned tsp_functions::finds_number(const std::string& input){
	unsigned n;
	std::regex rgx("\\-*\\d+\\.*\\d*");
	std::smatch match;
	n = 1;
	if (std::regex_search(input.begin(), input.end(), match, rgx))
	{
	    for (auto m : match)
	    n = stod(m);
	}
	else{
	    n = -1;
	}
	return n;
}

  
// Reads instance of TSPLIB
void tsp_functions::reads_instance(const std::string &input, double ***ppos, unsigned n){
	//declare and initiate matrix of coordunsignedate positions

	double **pmatrix = new double*[n]; 
  	for(unsigned i=0; i<n; i++){
      pmatrix[i] = new double[2];
  	}

  	// reads input file
	ifstream in_file(input);
	string line;

	// exits if file does not exists
	if(!in_file){
	cout << "File does not exist.";
	exit(1);
	}

	// skips the first 6 lines
	for(unsigned i=0; i<6;i++){
	getline(in_file, line);
	}

	unsigned i = 0;
	while (getline(in_file, line) && line != "EOF") {
	stringstream ss;
	ss << line;
	string x, y, z;
	ss >> x >> y >> z;
	double yd = stod(y);
	double zd = stod(z);
	pmatrix[i][0] = yd;
	pmatrix[i][1] = zd;
	i++;
	}
	in_file.close();
	*ppos = pmatrix;
}

void tsp_functions::print_matrix(double **pmatrix, unsigned m, unsigned n){
	for (unsigned i = 0; i < m; i++){
		for( unsigned j = 0; j < n; j++){
		 	printf("%f, ", pmatrix[i][j] );
		}
		printf("\n");
	}
}

void tsp_functions::print_matrix(unsigned **pmatrix, unsigned m, unsigned n){
	for (unsigned i = 0; i < m; i++){
		for(unsigned j = 0; j < n; j++){
		 	printf("%i, ", pmatrix[i][j] );
		}
		printf("\n");
	}
}

// Makes a copy of a matrix and prunes it with infitiny values
void tsp_functions::prune_matrix(double **pdist_complete, double ***pdist_prune, unsigned n, double prop_edges){
	unsigned nneighbors;
	double kth_near, diff;
	// Declare and initiate vector
	double *prow_dist_complete=0;
	prow_dist_complete = new double[n];

	// Declare and initiate row
	double **pmatrix=0;
	pmatrix = new double*[n];
	for(unsigned i=0; i<n; i++){
      pmatrix[i] = new double[n];
  	}

  	// Copies Matrix
	for (unsigned i = 0; i < n; i++){
		for(unsigned j = 0; j < n; j++){
		 	pmatrix[i][j] = pdist_complete[i][j];
		}
	}

	// Eliminate all but the first kth nearest edges for each node.
	nneighbors = floor(n*prop_edges); // Number of neighbors we keep
	for(unsigned i=0; i<n; i++){
		for(unsigned j=0; j<n; j++){
			prow_dist_complete[j] = pdist_complete[i][j];
		} 
		sort(prow_dist_complete, prow_dist_complete+n);
		kth_near = prow_dist_complete[nneighbors];
		
		for(unsigned j=0; j<n; j++){
			diff = pdist_complete[i][j] - kth_near; 
			if(diff > -numeric_limits<double>::epsilon()){
				pmatrix[i][j] = numeric_limits<double>::infinity();
			}
		}
	}
	*pdist_prune = pmatrix;
}

void tsp_functions::initiate_paths(unsigned ***ppaths, unsigned n){
	unsigned **pmatrix=0;
	pmatrix = new unsigned*[n];
	for(unsigned i=0; i<n; i++){
      pmatrix[i] = new unsigned[n+1];
      for(unsigned j=0; j<n; j++){
      	pmatrix[i][j] = 0;
      }
      pmatrix[i][0] = i+1;
      pmatrix[i][n] = i+1;
  	}
  	*ppaths = pmatrix;
}

void tsp_functions::initiate_costs(double ***pcosts, unsigned n){
	double **pmatrix=0;
	pmatrix = new double*[n];
	for(unsigned i=0; i<n; i++){
      pmatrix[i] = new double[n];
      for(unsigned j=0; j<n; j++){
      	pmatrix[i][j] = 0;
      }
  	}
  	*pcosts = pmatrix;
}

void tsp_functions::initiate_neighborslist(double **pdist_prune, unsigned **pdegree_neighbors, vector<unsigned> ***pneighbors_list,  
	unsigned ***pdegree_neighbors_list, unsigned n){
	vector<unsigned> **pneighbors_list_aux = new vector<unsigned>*[n]; // Create neighbors list
	for(unsigned i = 0;	i<n; i++){
		pneighbors_list_aux[i] = new vector<unsigned>;
	}
	
	unsigned *pdegree_neighbors_aux=0;
	pdegree_neighbors_aux = new unsigned[n];

	unsigned **pdegree_neighbors_list_aux = new unsigned*[n]; // Create neighbors degree list

	unsigned nneighbors;
	for(unsigned i=0; i<n; i++){
		nneighbors = 0;
		for (unsigned j = 0; j < n; ++j)
		{
			if (pdist_prune[i][j] < numeric_limits<double>::infinity() && pdist_prune[i][j] > numeric_limits<double>::epsilon())
			{
				pneighbors_list_aux[i] -> push_back(j+1); // en la posición end se inserta j+1
				nneighbors++;
			}
			pdegree_neighbors_aux[i] = nneighbors;
		}
	}

	for(unsigned i = 0;	i<n; i++){
		pdegree_neighbors_list_aux[i] = new unsigned[n];
		for(unsigned j=0; j<n;j++){
			pdegree_neighbors_list_aux[i][j] = pdegree_neighbors_aux[j];
		}
	}

	*pdegree_neighbors = pdegree_neighbors_aux;
	*pdegree_neighbors_list = pdegree_neighbors_list_aux; 
	*pneighbors_list = pneighbors_list_aux;
}