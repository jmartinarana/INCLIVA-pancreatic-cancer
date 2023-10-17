#! /bin/bash 

#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --job-name dnds
#SBATCH --partition long
#SBATCH --cpus-per-task 10 # “-c”
#SBATCH --mem 90G
#SBATCH --output=normal_pancreas_%j.out
#SBATCH --array=0
