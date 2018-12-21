#!/bin/bash

/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r cosmic_validation_0  #one should specify the absolute path of the executable Matlab program in your HPC.

mkdir mFiles_randNull
mkdir shFiles_randNull

python gen_m_sh_qsubFiles_randNull.py

mv ./mFiles_randNull/* ./
mv ./shFiles_randNull/* ./
chmod +x QSUB_randNull.sh

./QSUB_randNull.sh
