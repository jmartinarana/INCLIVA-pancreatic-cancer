#! /bin/bash

sbatch ./contra_sbatch.sh  --wait
sbatch ./contra_sbatch2.sh --wait
sbatch ./contra_sbatch3.sh --wait
sbatch ./contra_sbatch4.sh --wait
sbatch ./wisecondor_sbatch.sh  --wait
sbatch ./wisecondor_sbatch2.sh --wait
sbatch ./wisecondor_sbatch3.sh --wait

