#!/bin/sh
#PBS -m abe
#PBS -M mightzbn@gmail.com
#PBS -N ellip3d-cavity_expansion_retessellate_first
#PBS -l nodes=1:ppn=16
#PBS -l walltime=480:00:00
cd $PBS_O_WORKDIR
./ellip3d 16 > dep_out 
