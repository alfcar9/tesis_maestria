# Run following code to obtain results for Christofides.

make
instance_name="berlin52"
./tsp ../TSPLIB/TSPLIB_clean/TXTs/$instance_name.txt > Results/$instance_name
rm ../TSPLIB/TSPLIB_clean/TXTs/$instance_name".txt.tour"

# To obtain euclidean or geodesical results comment and uncomment the function 
# TSP::get_distance in tsp.cpp (lines 77-108). 
