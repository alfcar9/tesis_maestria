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
#include <eigen3/Eigen/Dense> //library for linear algebra
#include "tsp_functions.h"

using namespace std; //shorthand for standard input
using namespace Eigen; //shorthand for eigen
int main(){
  int n; //declare size of instance, number of cities
  double **pmatrix_pos=0; //declare and initialize pointer
  double **pmatrix_dist=0; //declare and initialize pointer
  const std::string input = "../Instancias/berlin52.tsp"; // Name of instance file
  n = tsp_functions::finds_number(input);

  pmatrix_pos = new double*[n];
  for(int i=0; i<n; i++){
      pmatrix_pos[i] = new double[n];
  }
  

  ifstream in_file(input);
  //ifstream in_file("example_file.txt");
  string line;
  unsigned i = 0;

  if(!in_file){
    cout << "File does not exist.";
    exit(1);
  }

  for(int i=0; i<6;i++){
    getline(in_file, line);
  }
  while (getline(in_file, line) && line != "EOF") {
  stringstream ss;
  ss << line;
  string x, y, z;
  ss >> x >> y >> z;
  double yd = stod(y);
  double zd = stod(z);
  pmatrix_pos[i][0] = yd;
  pmatrix_pos[i][1] = zd;
  i++;
  }

  in_file.close();
  
  // Read File
  
  pmatrix_dist = tsp_functions::matrix_euclidean(pmatrix_pos, n);
  //Print Matrix Pos
  for ( int i = 0; i < 1; i++){
    for( int j = 0; j < 2; j++){
      printf("%f, ", pmatrix_pos[i][j] );
    }
    printf("\n");
  }
  // Print Matrix Dist
  for ( int i = 0; i < 1; i++){
    for( int j = 0; j < 1; j++){
      printf("%f, ", pmatrix_dist[i][j] );
    }
    printf("\n");
  }

  return 0;
}
