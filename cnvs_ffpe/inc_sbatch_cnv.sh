#!/bin/bash 

function checkDirectories(){
        for dir in ${@}
        do
            if [ ! -d "${dir}" ]
            then
                mkdir -p ${dir}
            fi
        done
}
function usage() {
    echo -e "Usage: $0"
    echo -e "This script runs the pipeline for target seq data"
    echo -e "Mandatory parameters:"
    echo -e "-f configuration file"
    1>&2; exit 1;
}

while getopts "f:t:n:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        *)
                usage ;;
    esac
done

if [ -z "$f" ];
then
        echo -e "ERROR: -f is a mandatory parameter"
        echo -e "ERROR: -t is a mandatory parameter"
        echo -e "ERROR: -n is a mandatory parameter"
        usage
        exit
fi

date
echo -e "Your command: "$@
source ${f}

module load Singularity/3.6.4-GCC-5.4.0-2.26

echo "bam_dir_FFPE $bam_dir_FFPE"
echo "bam_dir_WBC $bam_dir_WBC"

set -x

unset display
$sing $cnvkit batch $bam_dir_FFPE/*.bam --normal $bam_dir_WBC/*.bam --targets $panel/target.bed --fasta $genome --output-reference cnvkit_output --output-dir $output_dir 

sbatch ./inc_sbatch_cnv_ffpe.sh

set +x
