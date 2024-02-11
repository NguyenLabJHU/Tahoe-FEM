#!/bin/sh
#PBS -m abe
#PBS -M mightzbn@gmail.com
#PBS -N ellip3d-cavity_expansion_retessellate_first
#PBS -l nodes=1:ppn=10
#PBS -l walltime=64:00:00
cd $PBS_O_WORKDIR
./ellip3d 10 > exp_out 