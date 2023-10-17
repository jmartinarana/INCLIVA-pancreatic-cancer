#!/bin/bash

function usage() {
    echo -e "Usage: $0"
    echo -e "This script runs the pipeline for target seq data"
    echo -e "Mandatory parameters:"
    echo -e "-f configuration file"
    echo -e "-t tumor sample"
    1>&2; exit 1;
}
while getopts "f:p:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        p)
                p=${OPTARG} ;;
        *)
                usage ;;
    esac
done

if [ -z "$f" ] || [ -z "$p" ] ;
then
        echo -e "ERROR: -f is a mandatory parameter"
        echo -e "ERROR: -p is a mandatory parameter"
        usage
        exit
fi



date
echo -e "Your command: "$@
source ${f}

plasma=${p}


echo "plasma  is ${plasma}"
mkdir -p ${output_100000}/${plasma}
echo "COMMAND--- ${wisecondorXR} predict ${output_wise_plasma}/${plasma}_100000.npz ${output_dir}/${project}_control_PON_100000.npz ${output_100000}/${plasma}/${plasma} --bed --plot --blacklist ${centromere}"
${wisecondorXR} predict ${output_wise_plasma}/${plasma}_100000.npz ${output_dir}/${project}_100000.npz ${output_100000}/${plasma}/${plasma} --bed --plot --blacklist ${centromere}
