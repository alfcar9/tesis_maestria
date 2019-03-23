/**********************************************
* FUNCTIONS
***********************************************/

#include <iostream> // header in standard library
#include <math.h>
#include "tsp_functions.h" //header in local directory

double tsp_functions::dist_euclidean(double *city1_ptr, double *city2_ptr){
  double result;
  result = 0.0;
  for(int i=0; i<2; i++){
      result += pow(city1_ptr[i]-city2_ptr[i], 2);
  }
  result = sqrt(result);
  return result;
}

// Given set of coordiantes calculates euclidean distance matrix.
double **tsp_functions::matrix_euclidean(double **matrix_pos_ptr, int n){
	double *city1_pos_ptr, *city2_pos_ptr, **matrix_dist_ptr;
	city1_pos_ptr = new double[2];
	city2_pos_ptr = new double[2]; 
	matrix_dist_ptr = new double*[n];
	for(int i=0; i<n; i++){
	    matrix_dist_ptr[i] = new double[n];
	}
	for(int i=0; i<(n-1); i++){
		city1_pos_ptr = &matrix_pos_ptr[i][0];
		for(int j=(i+1); j<n; j++){
			for(int k=0; k<2;k++){
			city2_pos_ptr = &matrix_pos_ptr[j][0];
		}
			matrix_dist_ptr[i][j] = dist_euclidean(city1_pos_ptr, city2_pos_ptr);
			matrix_dist_ptr[j][i] = matrix_dist_ptr[i][j]; 
		}
	}
	return matrix_dist_ptr;
}

// Makes a matrix symemetric
double **tsp_functions::matrix_symmetric(double **matrix_dist_ptr, int n){
  for(int i = 0; i<n-1; i++){
    for(int j = i+1; j<n; j++){
        if(matrix_dist_ptr[i][j] != matrix_dist_ptr[j][i]){
          if(matrix_dist_ptr[i][j] < matrix_dist_ptr[j][i]){
            matrix_dist_ptr[j][i] = matrix_dist_ptr[i][j];
          }
          else{
            matrix_dist_ptr[i][j] = matrix_dist_ptr[j][i];
          }
        }
     }
  }
  return matrix_dist_ptr;
}

// Updates distance matrix when triangle Inequality is not satissfied
double **tsp_functions::triangle_inequality(double **matrix_dist_ptr, int n){
	double dist_min, dist_indir;
	int city, city_distinct, set[n-2];
	for(int i = 0; i<n-1; i++){
	    for(int j = i+1; j<n; j++){
	     	dist_min = matrix_dist_ptr[i][j];
	    	for(int city=0; city<n;city++){
	    		city_distinct = 1;
	    		if(city != i && city != j){
	    			set[city] = city_distinct;
	    			city_distinct += 1;
	    		}
	     	}	      
		    for(int city_distinct = 0; city_distinct < n-2; city_distinct++){
	        	dist_indir = matrix_dist_ptr[i][city_distinct] + matrix_dist_ptr[city_distinct][j];
		        if(dist_indir < dist_min){
		          dist_min = dist_indir;
		        }
	    	}
	    	matrix_dist_ptr[i][j] = dist_min;
	      	matrix_dist_ptr[j][i] = dist_min;
	    }
	}
	return matrix_dist_ptr;	
}
 