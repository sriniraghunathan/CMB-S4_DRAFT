#!/bin/bash

##SBATCH -C haswell
##SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem 10000
#SBATCH --time=10:00:00
#SBATCH --output=batch_jobs/myjob.o%j
##SBATCH --mail-user=srinirag@illinois.edu
##SBATCH --mail-type=ALL

##modules
##module load Python/2.7.9-GCC-4.9.2-bare
##module load xcb-proto/1.11-intel-2016.u3-Python-2.7.9
##module use /home/sri/modulefiles/
##module load anaconda
##module load python

# Your script content goes here...
