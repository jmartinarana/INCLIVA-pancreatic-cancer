#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --job-name contra4
#SBATCH --partition long
#SBATCH --cpus-per-task 1# “-c”
#SBATCH --mem 4G
#SBATCH --output=contra_4_%j.out
#SBATCH --array=29-30


wbc=( "PDAC-DNAg-016" "PDAC-DNAg-017"  "PDAC-DNAg-018" "PDAC-DNAg-019" "PDAC-DNAg-021" "PDAC-DNAg-022" "PDAC-DNAg-023" "PDAC-DNAg-024" "PDAC-DNAg-025" "PDAC-DNAg-027" "PDAC-DNAg-028" "PDAC-DNAg-029" "PDAC-DNAg-030" "PDAC-DNAg-031" "PDAC-DNAg-032" "PDAC-DNAg-033" "PDAC-DNAg-034" "PDAC-DNAg-036" "PDAC-DNAg-037" "PDAC-DNAg-038" "PDAC-DNAg-039" "PDAC-DNAg-040" "PDAC-DNAg-041" "PDAC-DNAg-043" "PDAC-DNAg-044" "PDAC-DNAg-045" "PDAC-DNAg-047" "PDAC-DNAg-048" "PDAC-DNAg-050" "PDAC-DNAg-051")
plasma=( "1-PDAC-016_L1"  "2-PDAC-017_L1" "3-PDAC-018_L1" "4-PDAC-019_L1" "6-PDAC-021_L1" "7-PDAC-022_L1" "8-PDAC-023_L1" "9-PDAC-024_L1" "10-PDAC-025_L1" "12-PDAC-027_L1" "13-PDAC-028_L1" "14-PDAC-029_L1" "15-PDAC-030_L1" "16-PDAC-031_L1" "17-PDAC-032_L1" "18-PDAC-033_L1" "19-PDAC-034_L1" "21-PDAC-036_L1" "22-PDAC-037_L1" "23-PDAC-038_L1" "24-PDAC-039_L1" "25-PDAC-040_L1" "26-PDAC-041_L1" "28-PDAC-043_L1" "29-PDAC-044_L1"  "30-PDAC-045_L1" "32-PDAC-047_L1" "33-PDAC-048_L1" "35-PDAC-050_L1" "36-PDAC-051_L1")

t=${plasma[$SLURM_ARRAY_TASK_ID]}
n=${wbc[$SLURM_ARRAY_TASK_ID]}

export SLURM_TMPDIR="/home/mibarrola/UMP215_MM_PDAC/tmp/slurm_tmpdir"
export TMPDIR="/home/mibarrola/UMP215_MM_PDAC/tmp/slurm_tmpdir"
export SINGULARITYENV_TMPDIR="/home/mibarrola/UMP215_MM_PDAC/tmp/slurm_tmpdir"
module load CONTRA/2.0.8-foss-2016b
module load Singularity/3.6.4-GCC-5.4.0-2.26
module load BEDTools/2.26.0-foss-2016b
module load SAMtools/1.10-foss-2016b
module load R/4.0.2-foss-2016b-X11-20160819

/home/mibarrola/info/software/packages/in_house_scripts_/CNVs/exome/plasmas/run_contra_plasma_parte4_v2.0.sh \
-f /home/mibarrola/UMP444_exomas_pancreas_1108/conf/cnv_plasma.conf \
-p ${t} \
-n ${n}



