#!/bin/sh
#PBS -m abe
#PBS -M mightzbn@gmail.com
#PBS -N ellip3dSPH_dam_break3D
#PBS -l nodes=1:ppn=4
#PBS -l walltime=480:00:00
cd $PBS_O_WORKDIR
./ellip3dSPH 4> bursting_out 
