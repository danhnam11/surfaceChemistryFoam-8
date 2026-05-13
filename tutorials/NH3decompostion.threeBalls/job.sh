#!/bin/bash

# this is comment line

#parallel environment request

#$ -pe mpi_32 32

# our Job name 
#$ -N test.tutorial

##$ -q cpl3.q
#$ -q all-ib.q
#$ -S /bin/bash
##$ -l h=(cpl2-01)

#$ -cwd

source $HOME/OpenFOAM/OpenFOAM-8/etc/bashrc WM_COMPILER_TYPE=ThirdParty WM_COMPILER=Gcc48 WM_LABEL_SIZE=64 WM_MPLIB=OPENMPI FOAMY_HEX_MESH=yes

export MPI_EXEC=~/OpenFOAM/ThirdParty-8/platforms/linux64Gcc48/openmpi-2.1.1/bin/mpirun

solver=catalystFoam

$MPI_EXEC -np $NSLOTS $solver -parallel > run_out.log

