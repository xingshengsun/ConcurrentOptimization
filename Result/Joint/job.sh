#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=96:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
# SBATCH --nodes=2   # number of nodes
#SBATCH --mem-per-cpu=6G   # memory per CPU core
#SBATCH -J "Joint"   # job name
#SBATCH --mail-user=xssun10@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module purge
module load openmpi/1.10.7 dakota/6.9.0 ortiz/ls-dyna/7.0.0 gmsh/4.5.4

mpirun -n $SLURM_NPROCS dakota -i dakota.in -write_restart dakota.rst -o dakota.out -e dakota.err
