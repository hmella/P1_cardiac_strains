#!/bin/bash
step=$1
data=$2
for ((i=0; i<=$((data-step)); i=$((i+step))))
  do
     screen -d -m /usr/local/MATLAB/R2020a/bin/matlab -nodesktop -nosplash -r "ini=$i+1; fin=$((i + step)); test_1; quit"
done
