#!/bin/sh
#
#SBATCH --job-name=dsmSH
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cpushortprio,cpushort,datanew,data
#SBATCH --exclusive

module purge
module load slurm cmake intel/compiler intel/mpi hdf5/intel/1.10.1

mpirun ~/git/DSM_package_to_nobu/dsmSH/dsmSHmpi.x
