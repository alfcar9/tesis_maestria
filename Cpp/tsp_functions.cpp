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
void tsp_functions::dist_matrix(double **ppos, double ***ppdist_complete, unsigned n){
	double **pdist_complete = *ppdist_complete; 
	for(unsigned i=0; i<(n-1); i++){
		for(unsigned j=(i+1); j<n; j++){			
			pdist_complete[i][j] = dist_euclidean(ppos[i], ppos[j]);
			pdist_complete[j][i] = pdist_complete[i][j]; 
		}
	}
}

// Makes a matrix symemetric
void tsp_functions::symmetric_matrix(double ***ppdist, unsigned n){
	double **pdist = *ppdist;
	for(unsigned i = 0; i < n-1; ++i){
		for(unsigned j = i+1; j < n; ++j){
    		if( pdist[i][j] != pdist[j][i] ){
    		    if( pdist[i][j] < pdist[j][i] ){
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
	while( index < n){
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
void tsp_functions::reads_instance(const std::string &input, double ***pppos, unsigned n){
	//declare and initiate matrix of coordunsignedate positions

	double **ppos = *pppos;

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
	ppos[i][0] = yd;
	ppos[i][1] = zd;
	i++;
	}
	in_file.close();
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
void tsp_functions::prune_matrix(double **pdist_complete, double ***ppdist_prune, unsigned n, double prop_edges){
	unsigned nneighbors;
	double kth_near, diff;

	// Declare and initiate array
	double *prow_dist_complete = 0;
	prow_dist_complete = new double[n];

	double **pdist_prune = *ppdist_prune;

	// Copies Distance matrix
	for (unsigned i = 0; i < n; i++){
		for(unsigned j = 0; j < n; j++){
		 	pdist_prune[i][j] = pdist_complete[i][j];
		}
	}

	// Eliminate all but the first kth nearest edges for each node.
	nneighbors = floor(n * prop_edges); // Number of neighbors we keep
	for(unsigned i = 0; i < n; i++){
		for(unsigned j = 0; j < n; j++){
			prow_dist_complete[j] = pdist_complete[i][j];
		} 
		sort(prow_dist_complete, prow_dist_complete + n);
		kth_near = prow_dist_complete[nneighbors];
		for(unsigned j = 0; j < n; j++){
			diff = pdist_complete[i][j] - kth_near; 
			if(diff > -numeric_limits<double>::epsilon()){
				pdist_prune[i][j] = numeric_limits<double>::infinity();
			}
		}
	}
}

void tsp_functions::initiate_paths(unsigned ***pppaths, unsigned n){
	unsigned **ppaths = *pppaths;
	for(unsigned i = 0; i < n; i++){
      for(unsigned j = 0; j < n; j++){
      	ppaths[i][j] = 0;
      }
      ppaths[i][0] = i + 1;
      ppaths[i][n] = i + 1;
  	}
}

void tsp_functions::initiate_costs(double ***ppcosts, unsigned n){
	double **pcosts = *ppcosts;
	for(unsigned i = 0; i < n; i++){
      for(unsigned j = 0; j < n; j++){
      	pcosts[i][j] = 0;
      }
  	}
 }

void tsp_functions::initiate_neighborslist(double **pdist_prune, unsigned **ppneighbors_degree, vector<unsigned> ***ppneighbors_list,  
	unsigned ***ppneighbors_degree_list, unsigned n){
	
	vector<unsigned> **pneighbors_list = *ppneighbors_list;
	for(unsigned i = 0;	i<n; i++){
		pneighbors_list[i] = new vector<unsigned>;
	}
	unsigned *pneighbors_degree = *ppneighbors_degree;	
	unsigned **pneighbors_degree_list = *ppneighbors_degree_list; // Create neighbors degree list

	unsigned nneighbors;
	for(unsigned i=0; i<n; i++){
		nneighbors = 0;
		for (unsigned j = 0; j < n; ++j)
		{
			if (pdist_prune[i][j] < numeric_limits<double>::infinity() && pdist_prune[i][j] > numeric_limits<double>::epsilon())
			{
				pneighbors_list[i] -> push_back(j+1); // en la posición end se inserta j+1
				nneighbors++;
			}
			pneighbors_degree[i] = nneighbors;
		}
	}

	for(unsigned i = 0;	i < n; i++){
		for(unsigned j=0; j < n;j++){
			pneighbors_degree_list[i][j] = pneighbors_degree[j];
		}
	}
}

void tsp_functions::greedy_paths(unsigned immediate_value, unsigned n, double **ppos, unsigned ***pppaths, double ***ppcosts, double **pdist_complete,
 unsigned *pneighbors_degree, vector<unsigned> **pneighbors_list, unsigned ***ppneighbors_degree_list){
	unsigned city_current, number_neighbors, city_min, dist_min_index, candidate_city, candidate_city_degree, *candidate_city_ptr;
	double dist_min;
	bool immediate_value_bool;
	vector<unsigned> city_neighbors_to_go, degree_neighbors_to_go;
	vector<double> city_distances;
	std::vector<bool> immediate_neighbors;
	unsigned **ppaths = *pppaths;
	double **pcosts = *ppcosts;
	unsigned **pneighbors_degree_list = *ppneighbors_degree_list;

	for(unsigned num_iter = 1; num_iter < n; ++num_iter){
		for(unsigned i = 0; i < n; ++i){
		    city_current = ppaths[i][num_iter-1];
	 	    number_neighbors = 0;
		    // Se encuentran las ciudades vecinas que no se han visitado así como su respectivo grado
	 	    city_neighbors_to_go.clear();
	 	    degree_neighbors_to_go.clear();
	 	    city_distances.clear();
	 	    immediate_neighbors.clear();
	 	    for(unsigned j = 0; j < pneighbors_list[city_current-1] -> size(); ++j){ // Sobre cada posible vecino de current_city
 	    		candidate_city = *(pneighbors_list[city_current-1] -> begin() + j );
 	    		candidate_city_ptr = find(ppaths[i], ppaths[i] + num_iter, candidate_city);
 	    		if( candidate_city_ptr == ppaths[i] + num_iter) { // no se ha visitado tal ciudad
 	    			city_neighbors_to_go.push_back(candidate_city);
 	    			candidate_city_degree = *(pneighbors_degree_list[i] + candidate_city- 1 );
 	    			degree_neighbors_to_go.push_back(candidate_city_degree);
 	    			immediate_neighbors.push_back(candidate_city_degree <= immediate_value);
 	    			++number_neighbors;
 	    		}
	 	    }

			// Si no hay ciudades vecinas entonces se hace el calculo sobre todas las posibles ciudades
			if(number_neighbors == 0){
				 for (unsigned j = 0; j < n; ++j){
				 	for(unsigned k = 0; k < num_iter; ++k){
					 	if( j+1 != ppaths[i][k]){
							city_neighbors_to_go.push_back(j+1);
							candidate_city_degree = *(pneighbors_degree_list[i] + j );
							degree_neighbors_to_go.push_back(candidate_city_degree);							
						}
					}
				 } 
			}

	 	    // Se guardan en un vector solo los grados de los correspondientes nodos por visitar.
	 	    // Se calcula un valor booleano que previene que nos olvidemos de un nodo. Si existe un nodo con pocos
		    // vecinos hay que visitarlo inmediatamente, sino, conviene visitar el más proximo.

	 		immediate_value_bool = std::none_of(immediate_neighbors.begin(), immediate_neighbors.end(), [](bool value) { return value; });
		    if(immediate_value_bool){
		  	 	// En caso de que todos los nodos tengan multiples vecinos calcular el mas proximo
	 	     	for (auto j = city_neighbors_to_go.begin(); j != city_neighbors_to_go.end(); ++j){
	 		   		city_distances.push_back(*(pdist_complete[city_current-1] + *j-1));
	 	     	}
	 	     	dist_min_index = std::min_element(city_distances.begin(), city_distances.end()) - city_distances.begin();
	 	     	city_min = city_neighbors_to_go[dist_min_index];
	 	     	dist_min = city_distances[dist_min_index];
	 	     }
	 	 	else{
	 	    	// En caso de que haya que visitar uno inmediatamente entonces añadirlo
	 	    	dist_min_index = std::min_element(degree_neighbors_to_go.begin(), degree_neighbors_to_go.end()) - degree_neighbors_to_go.begin();
	 	    	city_min = city_neighbors_to_go[dist_min_index];
	 	    	dist_min = pdist_complete[city_current-1][city_min-1];
	 	    }	 	
	 	    
	 	    // Todos los vecinos que tengan asociado a la ciudad que se visita se les resta un grado
		    if(number_neighbors != 0){
		    	for(auto j = city_neighbors_to_go.begin(); j != city_neighbors_to_go.end(); ++j){
		    		*(pneighbors_degree_list[i] + *j -1) = *(pneighbors_degree_list[i] + *j -1) -1;
		    	}
		    }		    
		    // Se actualiza el grado de la ciudad actual a 0.
		    *(pneighbors_degree_list[i] + city_current - 1 ) = 0;
		    // Se actualiza la matriz de costo.
		    pcosts[i][num_iter-1] = dist_min;
		    ppaths[i][num_iter] = city_min;
		    // Se verifica que no haya crossing
		    if( num_iter > 2){
			    tsp_functions::crossing_procedure(city_current, city_min, num_iter, i, ppos, pppaths, ppcosts, pdist_complete);
		    } //end of if iter > 2
  		} // end of construction of path i
	} // end of iterations

	// Se verifica que la ultima arista no hace crossing
  	for(unsigned i = 0; i < n; ++i){
		city_current = ppaths[i][n-1];
		tsp_functions::crossing_procedure(city_current, i+1, n, i, ppos, pppaths, ppcosts, pdist_complete);
	}
}

void tsp_functions::crossing_procedure(unsigned city_current, unsigned city_min, unsigned num_iter, unsigned i, 
	double **ppos, unsigned ***pppaths, double ***ppcosts, double **pdist_complete){
	unsigned city11, city12, city21, city22;
	double b1, b2, m1, m2, x11, x12, x21, x22, y11, y12, y21, y22, x_intersect, dist1_replace, dist2_replace;
	bool bool1, bool2;
	vector<unsigned> paths_replace;
	vector<double> costs_replace;
	unsigned **ppaths = *pppaths;
	double **pcosts = *ppcosts;

	city21 = city_current;
    x21 = ppos[city21 - 1][0];
    y21 = ppos[city21 - 1][1];
    city22 = city_min;
    x22 = ppos[city22 - 1][0];
    y22 = ppos[city22 - 1][1];
    if(x22-x21 != 0){
    	m2 = (y22 - y21)/(x22-x21);
    }
	else{
		m2 =(y22-y21)/(x22 - x21 + 10e-5);
  	}
   	b2 = y21 - m2 * x21;
   	for(auto j = num_iter-2; j > 0; --j){ // ... las restantes aristas
        paths_replace.clear();
    	costs_replace.clear();
        city11 = ppaths[i][j-1];
        x11 = ppos[city11 - 1][0];
        y11 = ppos[city11 - 1][1];
        city12 = ppaths[i][j];
        x12 = ppos[city12 -1][0];
        y12 = ppos[city12 -1][1];
        if(x12 - x11 != 0){
          m1 = (y12 - y11)/(x12 - x11);
        }
    	else{
      		m1 = (y12 - y11)/(x12 - x11 + 10e-5);
    	}
     	b1 = y11 - m1 * x11;
     	if(m1 != m2){ //end of if crossing exist
 	  		x_intersect = -(b2 - b1)/(m2 - m1);
       		bool1 = (x11 < x_intersect && x_intersect < x12) || (x12 < x_intersect && x_intersect < x11);
       		bool2 = (x21 < x_intersect && x_intersect < x22) || (x22 < x_intersect && x_intersect < x21);
         	if(bool1 && bool2){ // condition for finding crossing
         		for(unsigned k = j; k < num_iter; ++k){
         			paths_replace.push_back(ppaths[i][k]);  
         		}
         		for(unsigned k = 0; k < num_iter - j; ++k){
         			ppaths[i][num_iter-k-1] = paths_replace[k];  
         		}
				dist1_replace = pdist_complete[city11 - 1][city21 - 1];
	        	dist2_replace = pdist_complete[city12 - 1][city22 - 1];
	        	costs_replace.push_back(dist1_replace);
	        	for(unsigned k = num_iter-2; k >= j; --k){
         			costs_replace.push_back(pcosts[i][k]);  
         		}
         		costs_replace.push_back(dist2_replace);
         		for(unsigned k = 0; k < num_iter - j+1; k++){
         			pcosts[i][j+k-1] = costs_replace[k];  
         		}
	            city_current = ppaths[i][num_iter -1];
	            city21 = city_current;
	            x21 = ppos[city21 -1][0];
	   		    y21 = ppos[city21 -1][1];
       		} // end of if crossing exist
   		} // end of if slopes are different
	} // end of for possible crossing
}