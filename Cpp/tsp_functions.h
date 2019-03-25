class tsp_functions
{
public:
	//Define Euclidean distance funcion
	static double dist_euclidean(double *pcity1, double *pcity2);
	// Given set of coordiantes calculates euclidean distance matrix.
	static double **dist_matrix(double **ppos, int n);
	// Assures the distance matrix symemetric
	static void symmetric_matrix(double **pdist, int n);
	// Updates distance matrix when triangle Inequality is not satissfied
	static double **triangle_inequality(double **pdist, int n);
	// Returns the index of thee jth ctiy in the ith path 
	static int index_function(int i, int j, int **paths, int n);
	// Reads instance and defines a matrix
	static double **reads_instance(const std::string &input, int n);
	// Extracts size of instance which is in the name of instance of TSPLIB
	static int finds_number(const std::string &input);
	// Prints Matrix
	static void print_matrix(double **pmatrix, int m, int n);
	// Prints Matrix
	static void print_matrix(int **pmatrix, int m, int n);
	// Makes inifity the distance of a proportion of neighbors
	static double **prune_matrix(double **pdist_complete, int n, double prop_edges);
	// Initiates paths matrix
	static int **initiate_paths(int n);
	// Initiates costs matrix
	static double **initiate_costs(int n);
	// Initates neighbors list
	//static vector<int> initiate_neighborslist(double **pdist_prune, int n);
};