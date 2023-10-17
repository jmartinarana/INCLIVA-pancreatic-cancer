#!/bin/bash

function usage() {
    echo -e "Usage: $0"
    echo -e "This script runs the pipeline for target seq data"
    echo -e "Mandatory parameters:"
    echo -e "-f configuration file"
    1>&2; exit 1;
}

while getopts "f:" opt; do
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
        usage
        exit
fi



date

echo -e "Your command: "$@
source ${f}
echo "${output_dir}"

#CONTRA_out_list.txt es una lista de todos los directorios de resultados que se han ido generando, por muestra

mkdir -p ${output_dir}/NUllD/
echo "COMMAND--- $nde_wrapper ${output_dir}/CONTRA_out_list.txt ${output_dir}/NUllD/ ${target} T F"	
$nde_wrapper ${output_dir}/CONTRA_out_list.txt ${output_dir}/NUllD/ ${target} T F


