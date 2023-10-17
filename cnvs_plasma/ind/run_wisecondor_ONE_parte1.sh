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


echo "plasma  is ${plasma}"
echo "COMMAND --- ${wisecondorX} convert ${plasma_bam}/${plasma}.bam ${output_wise_plasma}/${plasma}_100000.npz -r ${genome} --binsize 100000"
${wisecondorX} convert ${plasma_bam}/${plasma}.bam ${output_wise_plasma}/${plasma}_100000.npz -r ${genome} --binsize 100000

echo "COMMAND --- ${wisecondorX} convert ${WBC_bam}/${normal}.bam ${output_wise_WBC}/${normal}_100000.npz -r ${genome} --binsize 100000"
${wisecondorX} convert ${WBC_bam}/${normal}.bam ${output_wise_WBC}/${normal}_100000.npz -r ${genome} --binsize 100000




