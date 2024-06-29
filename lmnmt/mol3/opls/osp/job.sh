#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=7MR
#SBATCH --partition=gpu2080
#SBATCH --gres=gpu:1 
#SBATCH --ntasks=1
##SBATCH --tasks-per-node=10
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G 
#SBATCH --time=100:00:00
##SBATCH --output=slurm.out
#SBATCH --array=1-5
#SBATCH --export=ALL

#===============================================
ml load openmm/7.6.0 #gcc/10.2.0 fftw/3.3.8 cuda/11.2
ml load openmpi/3.1.2-gcc-10.2.0

#ml load openmpi/3.1.2-gcc-10.2.0 openmm/7.6.0 gcc/10.2.0 fftw/3.3.8 cuda/11.2
#-----------------------------------------------------

export i=$SLURM_ARRAY_TASK_ID

./run.sh

./analysis.sh

