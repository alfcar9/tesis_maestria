/**********************************************
* Este programa es la heuristica para TSP *
***********************************************/

#include <stdio.h>  //standard input
#include <iostream>
#include <string>
#include <math.h> // math functions
#include <cassert>
#include "tsp_functions.h"

using namespace std; //shorthand for standard input
using namespace N;

int main(){
  tsp_functions mc;
  int n;
  //double a[][2] = {{2.0,4.0}, {2.0,4.0}, {1.0,3.0}};
  //double b[][2] = {{3.0,1.0}, {2.0,1.0}, {3.0,5.0}};
  double a[2] = {2.0, 3.0};
  double b[2] = {1.0, 3.0};
  double result;
  //n = sizeof(b)/sizeof(*b);
  result = mc.dist_euclidean(a, b);
  printf("%f \n", result);
  //printf("%i \n", n);
  return 0;
}