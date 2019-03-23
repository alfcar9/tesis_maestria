/**********************************************
* Este programa es la heuristica para TSP *
***********************************************/

#include <stdio.h>  //standard input
#include <math.h> // math functions
#include <iostream>
#include <string> // std::string, std::stod
#include <cassert>
#include <fstream> // read files
#include <sstream>
#include "tsp_functions.h"

using namespace std; //shorthand for standard input

int main(){
  int n;
  double **matrix_pos_ptr, **matrix_dist_ptr;

  matrix_pos_ptr = new double*[n];
  for(int i=0; i<n; i++){
      matrix_pos_ptr[i] = new double[n];
  }

  ifstream in_file("../Instancias/berlin52.tsp");
  //ifstream in_file("example_file.txt");
  string line;
  unsigned i = 0;

  if(!in_file){
    cout << "File does not exist.";
    exit(1);
  }

  n = 52;
  while (getline(in_file, line)) {
  stringstream ss;
  ss << line;
  string x, y, z;
  ss >> x >> y >> z;
  double xd = stod(x);
  double yd = stod(y);
  double zd = stod(z);
  matrix_pos_ptr[i][0] = yd;
  matrix_pos_ptr[i][1] = zd;
  i++;
  }

  in_file.close();
  
  // Read File
  
  matrix_dist_ptr = tsp_functions::matrix_euclidean(matrix_pos_ptr, n);
  //Print Matrix Pos
  for ( int i = 0; i < n; i++){
    for( int j = 0; j < 2; j++){
      printf("%f, ", matrix_pos_ptr[i][j] );
    }
    printf("\n");
  }
  // Print Matrix Dist
  for ( int i = 0; i < n; i++){
    for( int j = 0; j < n; j++){
      printf("%f, ", matrix_dist_ptr[i][j] );
    }
    printf("\n");
  }

  return 0;
}
