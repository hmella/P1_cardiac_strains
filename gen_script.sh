#!/bin/bash
step=$1
data=$2
for ((i=0; i<=$((data-step)); i=$((i+step))))
  do
     screen -d -m mpirun --allow-run-as-root -n 4 python3 test_base_generation.py -gp 1 -Nd $data -Ns 4000000 -nf 1 -ini $i -fin $((i + step))
 done