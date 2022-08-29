#!/bin/sh
L=100
deltas="2"

OMP_NUM_THREADS=8
g++ homework1.cpp -fopenmp -o homework1

for d in $deltas
do
start=`date +%s.%N`
./homework1 > scaledvariance_L${L}_d${d}.dat <<EOF
$L
$d
EOF
end=`date +%s.%N`

runtime=$( echo "$end - $start" | bc -l )
echo "Run time for L = $L, delta = $d : ${runtime}s"
done
