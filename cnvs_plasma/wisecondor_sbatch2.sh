#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --job-name wisecondor_2
#SBATCH --partition long
#SBATCH --cpus-per-task 16# “-c”
#SBATCH --mem 40G
#SBATCH --output=wisecondor2_%j.out


export SLURM_TMPDIR="/home/mibarrola/UMP215_MM_PDAC/tmp/slurm_tmpdir"
export TMPDIR="/home/mibarrola/UMP215_MM_PDAC/tmp/slurm_tmpdir"
export SINGULARITYENV_TMPDIR="/home/mibarrola/UMP215_MM_PDAC/tmp/slurm_tmpdir"
module load CONTRA/2.0.8-foss-2016b
module load Singularity/3.6.4-GCC-5.4.0-2.26
module load BEDTools/2.26.0-foss-2016b
module load SAMtools/1.10-foss-2016b
module load R/4.0.2-foss-2016b-X11-20160819

/home/mibarrola/info/software/packages/in_house_scripts_/CNVs/exome/plasmas/run_wisecondor_ONE_parte2.sh \
-f /home/mibarrola/UMP444_exomas_pancreas_1108/conf/cnv_plasma.conf 


