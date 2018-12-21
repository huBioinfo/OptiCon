#!/bin/bash

#$ -o OptiCon_output_log.log

#$ -j y

#$ -l h_vmem=16G

#$ -l mem_free=16G

#$ -cwd

/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r gen_output  #one should specify the absolute path of the executable Matlab program in your HPC.
