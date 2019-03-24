class tsp_functions
{
public:
	//Define Euclidean distance funcion
	static double dist_euclidean(double *city1_ptr, double *city2_ptr);
	// Given set of coordiantes calculates euclidean distance matrix.
	static double **matrix_euclidean(double **matrix_pos_ptr, int n);
	// Assures the distance matrix symemetric
	static double **matrix_symmetric(double **matrix_dist_ptr, int n);
	// Updates distance matrix when triangle Inequality is not satissfied
	static double **triangle_inequality(double **matrix_dist_ptr, int n);
	// Returns the index of thee jth ctiy in the ith path 
	static int index_function(int i, int j, int **paths, int n);
	// Reads instance and defines a matrix
	//static double **reads_instance(string s, int n);
	static int finds_number(const std::string& input);
};