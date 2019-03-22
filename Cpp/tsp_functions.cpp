/**********************************************
* FUNCTIONS
***********************************************/

#include <iostream> // header in standard library
#include <math.h>
#include "tsp_functions.h" //header in local directory

double tsp_functions::dist_euclidean(double city1[2], double city2[2]){
  double result;
  result = 0.0;
  for(int i=0; i<=1; i++){
      result += pow(city1[i]-city2[i], 2);
  }
  result = sqrt(result);
  return result;
}

// Given set of coordiantes calculates euclidean distance matrix.
// int *TSP_functions::matrix_euclidean(){
//   int n, array_dist[n];
//   n = 5; 
//   //n = sizeof(matrix_pos)/sizeof(int);
//   for(int i=0; i<n;i++){
//       array_dist[i] = i;
//   }
//   //for(int i=0; i<=(n-2)){
//   //  i_vec = matrix_pos[i,];
//   //   for(j in (i+1):n){
//   //     j_vec = matrix_pos[j,];
//   //     matrix_dist[i,j] = round(dist_euclidean(i_vec, j_vec),2);
//   //     matrix_dist[j,i] = matrix_dist[i,j];
//   //   }
//   // }
//   return array_dist;
// }