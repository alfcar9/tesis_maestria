/**********************************************
* FUNCTIONS
***********************************************/

#include <iostream> // header in standard library
#include <math.h>
#include "tsp_functions.h" //header in local directory
#include <regex> // Regular Expresions

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
double **tsp_functions::matrix_euclidean(double **pmatrix_pos, int n){
	double *pcity1_pos=0;
	double *pcity2_pos=0;
	double **pmatrix_dist=0;
	pcity1_pos = new double[2];
	pcity2_pos = new double[2]; 
	pmatrix_dist = new double*[n];
	for(int i=0; i<n; i++){
	    pmatrix_dist[i] = new double[n];
	}
	for(int i=0; i<(n-1); i++){
		pcity1_pos = &pmatrix_pos[i][0];
		for(int j=(i+1); j<n; j++){
			for(int k=0; k<2;k++){
			pcity2_pos = &pmatrix_pos[j][0];
		}
			pmatrix_dist[i][j] = dist_euclidean(pcity1_pos, pcity2_pos);
			pmatrix_dist[j][i] = pmatrix_dist[i][j]; 
		}
	}
	return pmatrix_dist;
}

// Makes a matrix symemetric
double **tsp_functions::matrix_symmetric(double **pmatrix_dist, int n){
  for(int i = 0; i<n-1; i++){
    for(int j = i+1; j<n; j++){
        if(pmatrix_dist[i][j] != pmatrix_dist[j][i]){
          if(pmatrix_dist[i][j] < pmatrix_dist[j][i]){
            pmatrix_dist[j][i] = pmatrix_dist[i][j];
          }
          else{
            pmatrix_dist[i][j] = pmatrix_dist[j][i];
          }
        }
     }
  }
  return pmatrix_dist;
}

// Updates distance matrix when triangle Inequality is not satissfied
double **tsp_functions::triangle_inequality(double **pmatrix_dist, int n){
	double dist_min, dist_indir;
	int city_distinct, set[n-2];
	for(int i = 0; i<n-1; i++){
	    for(int j = i+1; j<n; j++){
	     	dist_min = pmatrix_dist[i][j];
	    	for(int city=0; city<n;city++){
	    		city_distinct = 1;
	    		if(city != i && city != j){
	    			set[city] = city_distinct;
	    			city_distinct ++;
	    		}
	     	}	      
		    for(int city=0; city < n-2; city++){
		    	city_distinct = set[city];
	        	dist_indir = pmatrix_dist[i][city_distinct] + pmatrix_dist[city_distinct][j];
		        if(dist_indir < dist_min){
		          dist_min = dist_indir;
		        }
	    	}
	    	pmatrix_dist[i][j] = dist_min;
	      	pmatrix_dist[j][i] = dist_min;
	    }
	}
	return pmatrix_dist;	
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
// double **reads_instance(string s, int n){

// }