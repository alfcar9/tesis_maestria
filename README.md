# Master's Thesis: Swapping Algorithm

Title: Swapping Algorithm

Author: Alfredo Carrillo del Campo

Date: November 2019

This repository is for solving the symmetrical TSP. It accepts the Euclidean distance or Geometrical distance for instances of the TSPLIB. The contents of the folders are as follows:

## 2Opt/
It stands for two-opt replacement. The script for running the two-opt method is contained in this folder.  It uses the "TSP" library for R. This library is much more extensive than only two-opt replacement, but yet this was the only one we used for contrasting the SA. Here two files are contained:
###TSP_package.R
It is the script that runs the instances of TSPLIB using the two-opt method from the library TSP. 
### modified_dist.txt
In case we want to input an instance with Geometrical distances, a function of the base package of R needs to be modified. This file contains the instructions on how the function dist of R must be rewritten.  Copy and paste what is in this text file when editing the function with the command *trace(dist, edit = TRUE)*. After running the instances, it is recommendable to change back to the original function of dist.

## Christofides/
A C++ Christofides implementation. Is the repository downloaded from [this repository]https://github.com/sth144/christofides-algorithm-cpp. The file *run_instance.sh* was added to run from the terminal the instances for TSLIB. Again, in case the instances with geometrical distances are tried out, the file tsp.cc needs to be modified. The instructions for this modification are written in the file *run_instance.sh*.

## GA/
It contains files to run a genetic algorithm in MATLAB for the TSP.
### GEO_CSVS/
Is the folder that contains CSVs that must be inputted to the *run_TSPs.m* file, for those instances that have geometrical distance. These CSVs were generated the cost_matrix function, in the folder SA, in the *required_functions.R* file.
### mtspf_ga.m
Is the genetic algorithm downloaded for Matlab.
### run_TSPS.m
Matlab file that runs the instances from the folder TSPLIB.

## SA/
It contains the correspondent files for the Swapping Algorithm.
### libraries.R
It contains the libraries we use throughout the programs in this folder.
### main.R
It contains the file that calls all the files from this folder and executes the SA algorithm for an instance given and plots the results.
### SA_function.R
It is the Swapping Algorithm implementation.
### required_functions.R
Auxiliary functions used during the SA.

## TSPLIB/
It Contains the instances that we ran for the master's thesis.
### TSPLIB_original/
It Contains the files exactly as they were downloaded from TSPLIB.
### TSPLIB_clean/
CSVs contains the files that have only the coordinates of the nodes. TXTs, similarly to the CSVs folder but for txt type. 
## Cpp: 
It contains an incomplete implementation of the Swapping Algorithm for C++.
