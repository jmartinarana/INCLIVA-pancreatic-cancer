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

if [ -z "$f" ] ;
then
        echo -e "ERROR: -f is a mandatory parameter"
        usage
        exit
fi



date
echo -e "Your command: "$@
source ${f}


echo "COMMAND---${wisecondorX} newref ${output_wise_WBC}/*100000.npz ${output_dir}/${project}_100000.npz --nipt --binsize 100000 --cpus 40"
${wisecondorX} newref ${output_wise_WBC}/*100000.npz ${output_dir}/${project}_100000.npz --nipt --binsize 100000 --cpus 40


#mkdir -p ${output_100000}/${plasma}

#echo "COMMAND--- ${wisecondorX} predict ${output_wise_plasma}/${plasma}_100000.npz ${output_dir}/${project}_control_PON_100000.npz ${output_100000}/${plasma}/${plasma} --bed --plot --blacklist ${centromere}"
#${wisecondorX} predict ${output_wise_plasma}/${plasma}_100000.npz ${output_dir}/${project}_control_PON_100000.npz ${output_100000}/${plasma}/${plasma} --bed --plot --blacklist ${centromere}

