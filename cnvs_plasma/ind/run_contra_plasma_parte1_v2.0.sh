#!/bin/bash

function usage() {
    echo -e "Usage: $0"
    echo -e "This script runs the pipeline for target seq data"
    echo -e "Mandatory parameters:"
    echo -e "-f configuration file"
    echo -e "-p tumor sample"
    echo -e "-n normal sample"
    1>&2; exit 1;
}

while getopts "f:p:n:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        p)
                p=${OPTARG} ;; 
        n)
                n=${OPTARG} ;;
        *)
                usage ;;
    esac
done

if [ -z "$f" ] || [ -z "$p" ] || [ -z "$n" ];
then
        echo -e "ERROR: -f is a mandatory parameter"
        echo -e "ERROR: -p is a mandatory parameter"
        echo -e "ERROR: -n is a mandatory parameter"
        usage
        exit
fi



date

echo -e "Your command: "$@
source ${f}
plasma=${p}
normal=${n}
echo "${output_dir}"

#CONTRA_out_list.txt es una lista de todos los directorios de resultados que se han ido generando, por muestra
echo "plasma  is ${plasma}"
#mkdir -p ${output_dir}/${plasma}
echo "COMMAND--- $contra -t ${target} -s ${plasma_bam}/${plasma}.bam -c ${WBC_bam}/${normal}.bam -f ${genome} --minExon=5000 --maxRegionSize=100000 --targetRegionSize=100000 --sampleName ${plasma} -p -l -o ${output_dir}/${plasma}"
$contra -t ${target} -s ${plasma_bam}/${plasma}.bam -c ${WBC_bam}/${normal}.bam -f ${genome} --minExon=5000 --maxRegionSize=100000 --targetRegionSize=100000 --sampleName ${plasma} -p -l -o ${output_dir}/${plasma}/


