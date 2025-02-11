#!/bin/bash
#SBATCH --job-name=5models
#SBATCH --account=mcmanus
#SBATCH --partition=mcmanus
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=0-1:00:00 ## time format is DD-HH:MM:SS

#SBATCH --nodes=10
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    

#SBATCH --error=5models-%A.err ## %A - filled with jobid
#SBATCH --output=5models-%A.out ## %A - filled with jobid

## Useful for remote notification
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=zquinlan@hawaii.edu

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html
# =============== Start of commands to execute ===============

module purge
module load lang/Anaconda3

cd /home/zquinlan/ccaModeling

python simulateManualAdditions.py
