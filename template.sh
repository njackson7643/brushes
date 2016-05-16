#!/bin/bash

#SBATCH --job-name=***
#SBATCH --output=***.out
#SBATCH --error=***.err
#SBATCH --time=***
#SBATCH --qos=***
#SBATCH --partition=***
#SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --ntasks-per-node=***

module load lammps/10Aug2015+intelmpi-5.0+intel-15.0

export OMP_NUM_THREADS = ***

mpirun lmp_intelmpi < ***.in

