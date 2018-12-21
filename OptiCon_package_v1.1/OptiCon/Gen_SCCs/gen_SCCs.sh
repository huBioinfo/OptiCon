#!/bin/bash

#$ -o gen_SCCs_log.log

#$ -j y

#$ -l h_vmem=24G

#$ -l mem_free=24G

#$ -cwd

g++ gen_diffSCCs.cpp struct.cpp  

./a.out

/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r gen_RNgeneID  #one should specify the absolute path of the executable Matlab program in your HPC.

mv ./CF_*.mat ../
mv ./RNgeneID.mat ../


