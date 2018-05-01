#!/bin/sh
#
#SBATCH --job-name=dsmFFT
#SBATCH --ntasks=1
#SBATCH --partition=cpushortprio,cpushort,datanew,data
#SBATCH --exclusive

module purge
module load cmake intel/compiler intel/mpi hdf5/intel/1.10.1

~/git/DSM_package_to_nobu/dsmFFT/dsmFFT.x
