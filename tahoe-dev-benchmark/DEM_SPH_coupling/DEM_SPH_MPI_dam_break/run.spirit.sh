#!/bin/bash

 

## optional

#PBS -m abe

#PBS -M mightzbn@gmail.com

#PBS -j oe

#PBS -N run

 

## required

#PBS -l select=128:ncpus=16:mpiprocs=16

#PBS -l walltime=60:00:00

#PBS -l place=scatter:excl

#PBS -q standard

#PBS -A ONRDC34502393

 

## cd to work dir

cd $PBS_O_WORKDIR

 

## export additional environment variables

export LD_LIBRARY_PATH=/work4/projects/openmpi/boost-1.57.0_mpt-2.11_intel-compilers-13.0.1/lib:$LD_LIBRARY_PATH

 

module unload mpt mpt

module unload intel-compilers intel-compilers

module load intel-compilers/13.0.1

module load mpt/2.11

 

mpiexec_mpt -np 2048 ./paraEllip3d input.txt

