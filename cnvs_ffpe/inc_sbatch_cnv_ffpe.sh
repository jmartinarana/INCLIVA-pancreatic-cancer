#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --job-name CNV_FFPE_PDAC
#SBATCH --partition bigmem
#SBATCH --cpus-per-task 4 # “-c”
#SBATCH --mem 32G
#SBATCH --output=tumor_FFPE_cnv_%j.out
#SBATCH --array=[0-30]


tt=( "TT-PDAC-018"     "TT-PDAC-021" "TT-PDAC-022"   "TT-PDAC-023"   "TT-PDAC-024"   "TT-PDAC-025"   "TT-PDAC-029"   "TT-PDAC-032"   "TT-PDAC-034"   "TT-PDAC-037" "TT-PDAC-039"  "TT-PDAC-044" "TT-PDAC-050" )
nn=( "PDAC-DNAg-018" "PDAC-DNAg-021" "PDAC-DNAg-022" "PDAC-DNAg-023" "PDAC-DNAg-024" "PDAC-DNAg-025" "PDAC-DNAg-029" "PDAC-DNAg-032" "PDAC-DNAg-034" "PDAC-DNAg-037" "PDAC-DNAg-039" "PDAC-DNAg-044" "PDAC-DNAg-050" )

normal=${nn[$SLURM_ARRAY_TASK_ID]}
tumor=${tt[$SLURM_ARRAY_TASK_ID]}



export SLURM_TMPDIR="/home/mibarrola/UMP444_exomas_pancreas_1108/scripts/tmp_slurm"
export TMPDIR="/home/mibarrola/UMP444_exomas_pancreas_1108/scripts/tmp_slurm"
export SINGULARITYENV_TMPDIR="/home/mibarrola/UMP444_exomas_pancreas_1108/scripts/tmp_slurm"

module load Singularity/3.6.4-GCC-5.4.0-2.26
module load BEDTools/2.29.2-foss-2016b

/home/mibarrola/UMP444_exomas_pancreas_1108/scripts/cnvs_ffpe/cnv.sh \
 -f "/home/mibarrola/UMP444_exomas_pancreas_1108/conf/cnv_FFPE.conf" \
 -t ${tumor} \
 -n ${normal}


