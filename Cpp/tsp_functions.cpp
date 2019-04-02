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
void tsp_functions::symmetric_matrix(double ***pdist, unsigned n){
	for(unsigned i = 0; i < n-1; ++i){
		for(unsigned j = i+1; j < n; ++j){
    		if( (*pdist)[i][j] != (*pdist)[j][i] ){
    		    if( (*pdist)[i][j] < (*pdist)[j][i] ){
            		(*pdist)[j][i] = (*pdist)[i][j];
         		}
         		else{
            		(*pdist)[i][j] = (*pdist)[j][i];
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

int tsp_functions::greedy_paths(unsigned immediate_value, unsigned n, double **ppos, unsigned ***ppaths, double ***pcosts, double **pdist_complete,
 unsigned *pdegree_neighbors, vector<unsigned> **pneighbors_list, unsigned ***ppdegree_neighbors_list){
	unsigned city_current, number_neighbors, city_min, dist_min_index, candidate_city, candidate_city_degree, *candidate_city_ptr;
	double dist_min;
	bool immediate_value_bool;
	vector<unsigned> city_neighbors_to_go, degree_neighbors_to_go;
	vector<double> city_distances;
	std::vector<bool> immediate_neighbors;
	for(unsigned num_iter = 1; num_iter < n; ++num_iter){
		for(unsigned i = 0; i < n; ++i){
		    city_current = (*ppaths)[i][num_iter-1];
	 	    number_neighbors = 0;
		    // Se encuentran las ciudades vecinas que no se han visitado así como su respectivo grado
	 	    city_neighbors_to_go.clear();
	 	    degree_neighbors_to_go.clear();
	 	    city_distances.clear();
	 	    immediate_neighbors.clear();
	 	    for(unsigned j = 0; j < pneighbors_list[city_current-1] -> size(); ++j){ // Sobre cada posible vecino de current_city
 	    		candidate_city = *(pneighbors_list[city_current-1] -> begin() + j );
 	    		candidate_city_ptr = find((*ppaths)[i], (*ppaths)[i] + num_iter, candidate_city);
 	    		if( candidate_city_ptr == (*ppaths)[i] + num_iter) { // no se ha visitado tal ciudad
 	    			city_neighbors_to_go.push_back(candidate_city);
 	    			candidate_city_degree = *((*ppdegree_neighbors_list)[i] + candidate_city- 1 );
 	    			degree_neighbors_to_go.push_back(candidate_city_degree);
 	    			immediate_neighbors.push_back(candidate_city_degree <= immediate_value);
 	    			++number_neighbors;
 	    		}
	 	    }

			// Si no hay ciudades vecinas entonces se hace el calculo sobre todas las posibles ciudades
			if(number_neighbors == 0){
				 for (unsigned j = 0; j < n; ++j){
				 	for(unsigned k = 0; k < num_iter; ++k){
					 	if( j+1 != (*ppaths)[i][k]){
							city_neighbors_to_go.push_back(j+1);
							candidate_city_degree = *((*ppdegree_neighbors_list)[i] + j );
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
		    		*((*ppdegree_neighbors_list)[i] + *j -1) = *((*ppdegree_neighbors_list)[i] + *j -1) -1;
		    	}
		    }		    
		    // Se actualiza el grado de la ciudad actual a 0.
		    *((*ppdegree_neighbors_list)[i] + city_current - 1 ) = 0;
		    // Se actualiza la matriz de costo.
		    (*pcosts)[i][num_iter-1] = dist_min;
		    (*ppaths)[i][num_iter] = city_min;
		    // Se verifica que no haya crossing
		    if( num_iter > 2){
			    tsp_functions::crossing_procedure(city_current, city_min, num_iter, i, ppos, ppaths, pcosts, pdist_complete);
		    } //end of if iter > 2
  		} // end of construction of path i
	} // end of iterations
	return 0;
}

void tsp_functions::crossing_procedure(unsigned city_current, unsigned city_min, unsigned num_iter, unsigned i, 
	double **ppos, unsigned ***ppaths, double ***pcosts, double **pdist_complete){
	unsigned city11, city12, city21, city22;
	double b1, b2, m1, m2, x11, x12, x21, x22, y11, y12, y21, y22, x_intersect, dist1_replace, dist2_replace;
	bool bool1, bool2;
	vector<unsigned> paths_replace;
	vector<double> costs_replace;
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
        city11 = (*ppaths)[i][j-1];
        x11 = ppos[city11 - 1][0];
        y11 = ppos[city11 - 1][1];
        city12 = (*ppaths)[i][j];
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
         			paths_replace.push_back((*ppaths)[i][k]);  
         		}
         		for(unsigned k = 0; k < num_iter - j; ++k){
         			(*ppaths)[i][num_iter-k-1] = paths_replace[k];  
         		}
				dist1_replace = pdist_complete[city11 - 1][city21 - 1];
	        	dist2_replace = pdist_complete[city12 - 1][city22 - 1];
	        	costs_replace.push_back(dist1_replace);
	        	for(unsigned k = num_iter-2; k >= j; --k){
         			costs_replace.push_back((*pcosts)[i][k]);  
         		}
         		costs_replace.push_back(dist2_replace);
         		for(unsigned k = 0; k < num_iter - j+1; k++){
         			(*pcosts)[i][j+k-1] = costs_replace[k];  
         		}
	            city_current = (*ppaths)[i][num_iter -1];
	            city21 = city_current;
	            x21 = ppos[city21 -1][0];
	   		    y21 = ppos[city21 -1][1];
       		} // end of if crossing exist
   		} // end of if slopes are different
	} // end of for possible crossing
}