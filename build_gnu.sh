#!/bin/sh

module purge
module load cmake gcc openmpi/gcc hdf5/intel

export CC=gcc
export FC=gfortran
#export HDF5_ROOT=/cm/shared/apps/hdf5/intel/1.10.1

#cmake ~/prjMRRS/DSM/myDSM/
cmake ~/git/DSM_package_to_nobu/myDSM/

make #VERBOSE=1
