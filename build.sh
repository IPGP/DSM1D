#!/bin/sh

# intel/compiler/64/15.0/1.133 intel/mpi/64/2015/5.0.2.044 hdf5/intel
# HDF5_ROOT=/cm/shared/apps.BCM6/hdf5/intel/serial/1.8.14

module purge
module load shared slurm
module load cmake intel/compiler intel/mpi hdf5/intel/1.10.1

export CC=icc
export FC=ifort
export I_MPI_F90=ifort
export HDF5_ROOT=/cm/shared/apps/hdf5/intel/1.10.1

cmake ./

make #VERBOSE=1
