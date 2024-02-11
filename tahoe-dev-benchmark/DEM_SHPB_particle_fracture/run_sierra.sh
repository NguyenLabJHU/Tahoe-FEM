#!/bin/csh
#MSUB -A ees 
#MSUB -l nodes=1
#MSUB -l partition=sierra
#MSUB -l walltime=16:00:00
#MSUB -q pbatch
#MSUB -m be
#MSUB -V
#MSUB -o ./com_out

##### These are shell commands
date
cd ./
srun -n 1 ./ellip3d  
date
