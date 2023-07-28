#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 1            # number of cores 
#SBATCH -t 2-00:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

module purge
module add r-4.2.2-gcc-11.2.0

make results/params/$1
