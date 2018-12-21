#!/bin/bash

#$ -o GreedySearch_input_log.log

#$ -j y

#$ -l h_vmem=24G

#$ -l mem_free=24G

#$ -cwd


/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r gen_inputs  #one should specify the absolute path of the executable Matlab program in your HPC.

/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r gen_inputs_Files  #one should specify the absolute path of the executable Matlab program in your HPC.

python gen_m_sh_qsubFiles.py

mv ./mFiles/* ./
mv ./shFiles/* ./
chmod +x QSUB.sh
