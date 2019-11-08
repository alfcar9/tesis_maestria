# Master's Thesis: Swapping Algorithm
This repository is for solving the symmetrical TSP with either the Euclidean distance or Geometrical distance for instances of the TSPLIB.
The contents of the folders are as follows:

## 2Opt.
Stands for two-opt replacement. The script for running the two-opt method is contained in this folders. For this, we downloaded the "TSP" library for R. This library is much more extensive than only two-opt replacement, but yet this was the only one we used for contrasting the SA. Here two files are contained:
###TSP_package.R
Is the script that runs the instances of TSPLIB using the two-opt method from the library TSP. 
### modified_dist.txt
In case we want to input an instance with Geometrical distances, a function of the base package of R needs to be modified. This folders contained how the function dist of R must be rewritten. Simply copy and paste what is in this text file. After running the instances it is recommendable to change back to the original function of dist.

## Christofides for a C++ Christofides implementation.
Is the repository downloaded from https://github.com/sth144/christofides-algorithm-cpp. The file $run_instance.sh$ was added to run from terminal the instances for TSLIB. Again, in case the instances with geometrical distances are tried out, the file tsp.cc needs to be modified. What is required to do is written in the file $run_instance.sh$.

## GA:
Contains files to run a genetic algorithm in MATLAB for the TSP.
### GEO_CSVS
Folder that contains CSVs that must be inputed to the $run_TSPs.m$ file, for those instances that have geometrical distance. These CSVs were generated the cost_matrix function, in the folder SA, in the $required_functions.R$ file.
### mtspf_ga.m
Is the genetic algorithm download for Matlab.
### run_TSPS.m
Matlab file that runs the instances from the folder TSPLIB.

## SA:
Contains the correspodent files for the Swapping Algorithm.
### libraries.R
Contains the libraries we use throught the programs.
### main.R
Contains the file that calls all the files from this folder and executes the SA algorithm for an instance given and plots the results.
### SA_function.R
Is the Swapping Algorithm implementation.
### required_functions.R
Auxiliary functions used during the SA.
## TSPLIB
Contains the instances that we ran for my master's thesis.
### TSPLIB_original:
Contains the files exactly as they were donwloaded from TSPLIB
### TSPLIB_clean:
CSVs contains the files that have only the coordinates of the nodes. TXTs, similary to the CSVs folder but for txt type. 
## Cpp 
Contains an incomplete implementation of the Swapping Algorithm for C++.
