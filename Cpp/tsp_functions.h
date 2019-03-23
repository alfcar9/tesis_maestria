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
};