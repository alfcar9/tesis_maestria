class tsp_functions
{
public:
	//Define Euclidean distance funcion
	static double dist_euclidean(double *pcity1, double *pcity2);
	// Given set of coordiantes calculates euclidean distance matrix.
	static double **matrix_euclidean(double **pmatrix_pos, int n);
	// Assures the distance matrix symemetric
	static double **matrix_symmetric(double **pmatrix_dist, int n);
	// Updates distance matrix when triangle Inequality is not satissfied
	static double **triangle_inequality(double **pmatrix_dist, int n);
	// Returns the index of thee jth ctiy in the ith path 
	static int index_function(int i, int j, int **paths, int n);
	// Reads instance and defines a matrix
	static double **reads_instance(const std::string &input, int n);
	// Extracts size of instance which is in the name of instance of TSPLIB
	static int finds_number(const std::string &input);
	// Prints Matrix
	static void print_matrix(double **pmatrix, int m, int n);
	// Copies square matrix
	static double **copy_matrix(double **pmatrix, int m, int n);
};