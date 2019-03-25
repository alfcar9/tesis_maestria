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
using namespace std; //shorthand for standard input
#include <vector> // neighbors list changes size in runtime

double tsp_functions::dist_euclidean(double *pcity1, double *pcity2){
  double result;
  result = 0.0;
  for(int i=0; i<2; i++){
      result += pow(pcity1[i]-pcity2[i], 2);
  }
  result = sqrt(result);
  return result;
}

// Given set of coordiantes calculates euclidean distance matrix.
double **tsp_functions::dist_matrix(double **ppos, int n){
	double *pcity1_pos=0;
	double *pcity2_pos=0;
	double **pdist=0;
	pcity1_pos = new double[2];
	pcity2_pos = new double[2]; 
	pdist = new double*[n];
	for(int i=0; i<n; i++){
	    pdist[i] = new double[n];
	}
	for(int i=0; i<(n-1); i++){
		pcity1_pos = &ppos[i][0];
		for(int j=(i+1); j<n; j++){
			for(int k=0; k<2;k++){
			pcity2_pos = &ppos[j][0];
		}
			pdist[i][j] = dist_euclidean(pcity1_pos, pcity2_pos);
			pdist[j][i] = pdist[i][j]; 
		}
	}
	return pdist;
}

// Makes a matrix symemetric
void tsp_functions::symmetric_matrix(double **pdist, int n){
  for(int i = 0; i<n-1; i++){
    for(int j = i+1; j<n; j++){
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

// Updates distance matrix when triangle Inequality is not satissfied
double **tsp_functions::triangle_inequality(double **pdist, int n){
	double dist_min, dist_indir;
	int city_distinct, set[n-2];
	for(int i = 0; i<n-1; i++){
	    for(int j = i+1; j<n; j++){
	     	dist_min = pdist[i][j];
	    	for(int city=0; city<n;city++){
	    		city_distinct = 1;
	    		if(city != i && city != j){
	    			set[city] = city_distinct;
	    			city_distinct ++;
	    		}
	     	}	      
		    for(int city=0; city < n-2; city++){
		    	city_distinct = set[city];
	        	dist_indir = pdist[i][city_distinct] + pdist[city_distinct][j];
		        if(dist_indir < dist_min){
		          dist_min = dist_indir;
		        }
	    	}
	    	pdist[i][j] = dist_min;
	      	pdist[j][i] = dist_min;
	    }
	}
	return pdist;	
}

// Returns the index of thee jth ctiy in the ith path 
int tsp_functions::index_function(int i, int j, int **paths, int n){
	int index_ans, index;
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

int tsp_functions::finds_number(const std::string& input){
	int n;
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
double **tsp_functions::reads_instance(const std::string &input, int n){
	//declare and initiate matrix of coordintate positions
	double **ppos=0;
	ppos = new double*[n];
  	for(int i=0; i<n; i++){
      ppos[i] = new double[2];
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
	for(int i=0; i<6;i++){
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
	ppos[i][0] = yd;
	ppos[i][1] = zd;
	i++;
	}
	in_file.close();
	return ppos;
}

void tsp_functions::print_matrix(double **pmatrix, int m, int n){
	for (int i = 0; i < m; i++){
		for( int j = 0; j < n; j++){
		 	printf("%f, ", pmatrix[i][j] );
		}
		printf("\n");
	}
}

void tsp_functions::print_matrix(int **pmatrix, int m, int n){
	for (int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
		 	printf("%i, ", pmatrix[i][j] );
		}
		printf("\n");
	}
}

// Makes a copy of a matrix and prunes it with infitiny values
double **tsp_functions::prune_matrix(double **pdist_complete, int n, double prop_edges){
	int nneighbors;
	double kth_near;
	// Declare and initiate vector
	double *prow_dist_complete=0;
	prow_dist_complete = new double[n];

	// Declare and initiate row
	double **pdist_prune=0;
	pdist_prune = new double*[n];
	for(int i=0; i<n; i++){
      pdist_prune[i] = new double[n];
  	}

  	// Copies Matrix
	for (int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
		 	pdist_prune[i][j] = pdist_complete[i][j];
		}
	}

	// Eliminate all but the first kth nearest edges for each node.
	nneighbors = floor(n*prop_edges); // Number of neighbors we keep
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			prow_dist_complete[j] = pdist_complete[i][j];
		} 
		sort(prow_dist_complete, prow_dist_complete+n);
		kth_near = prow_dist_complete[nneighbors];
		
		for(int j=0; j<n; j++){
			if(pdist_complete[i][j] > kth_near){
				pdist_prune[i][j] = numeric_limits<double>::infinity();
			}
		}
	}
	return pdist_prune;
}

int **tsp_functions::initiate_paths(int n){
	int **ppaths=0;
	ppaths = new int*[n];
	for(int i=0; i<n; i++){
      ppaths[i] = new int[n+1];
      for(int j=0; j<n; j++){
      	ppaths[i][j] = 0;
      }
      ppaths[i][0] = i+1;
      ppaths[i][n] = i+1;
  	}
  	return ppaths;
}

double **tsp_functions::initiate_costs(int n){
	double **pcosts=0;
	pcosts = new double*[n];
	for(int i=0; i<n; i++){
      pcosts[i] = new double[n];
      for(int j=0; j<n; j++){
      	pcosts[i][j] = 0;
      }
  	}
  	return pcosts;
}

// int **tsp_functions::initiate_neighborslist(double **pdist_prune, int n){
// 	vector<int> pneighbors_list[n];
// 	int number_neighbors;
// 	for(int i=0; i<n; i++){
// 		number_neighbors =0;
// 		for (int j = 0; j < n; ++j)
// 		{
// 			if (pdist_prune[i][j] < numeric_limits<double>::infinity() && pdist_prune[i][j] > 0)
// 			{
// 				number_neighbors++;
// 			}
// 		}
// 	}
// 	return pneighbors_list;
// }

