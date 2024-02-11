#!/bin/sh
#PBS -m abe
#PBS -M yourEmail@gmail.com
#PBS -N ellip3d-cavity_expansion
#PBS -l nodes=1:ppn=8
#PBS -l walltime=60:00:00
cd $PBS_O_WORKDIR
./ellip3d 8 > exp_out 
