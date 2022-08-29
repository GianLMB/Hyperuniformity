#!/bin/sh
# In the first part an LxL square lattice is used, and the scaled variance is computed for different values of the radius R.
 
Ls="10 100 200" # different values of L are used.
d=0  # no shuffle of the grid is used in the first part.

OMP_NUM_THREADS=8  # all the virtual processors are used to parallelize the code.
g++ homework1.cpp -fopenmp -o homework1  # the code is compiled in parallel.

for L in $Ls
do
start=`date +%s.%N`
./homework1 <<< $L, $d > scaledvariance_L${L}_d${d}.dat   # the code is executed with L and d from input. The output is saved in the .dat file.
end=`date +%s.%N`

runtime=$( echo "$end - $start" | bc -l )
echo "Run time for L = $L, delta = $d : ${runtime}s"
done
