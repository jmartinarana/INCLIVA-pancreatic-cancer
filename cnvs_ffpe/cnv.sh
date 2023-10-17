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
    echo -e "-t tumor sample"
    echo -e "-n normal sample file"
    1>&2; exit 1;
}

while getopts "f:t:n:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        t)
                t=${OPTARG} ;;
        n)
                n=${OPTARG} ;;
        *)
                usage ;;
    esac
done

if [ -z "$f" ] || [ -z "$t" ] || [ -z "$n" ] ;
then
        echo -e "ERROR: -f is a mandatory parameter"
        echo -e "ERROR: -t is a mandatory parameter"
        echo -e "ERROR: -n is a mandatory parameter"
        usage
        exit
fi


date

echo -e "Your command: "$@
conf=${f}
source ${f}
tumor=${t}
normal=${n}


chmod -R 776 $output_dir
echo "COMMAND --- run_facets_v2.0.sh -f $conf -t $tumor -n $normal"
/home/mibarrola/info/software/packages/in_house_scripts_/CNVs/exome/run_facets_v2.0.sh -f  $conf -t $tumor -n $normal
echo "fin run_facets"
chmod -R 776 $output_dir
echo "COMMAND ---varscan2_pipeline_ONE_v2.0.sh -f  $conf -t $tumor -n $normal"
/home/mibarrola/info/software/packages/in_house_scripts_/CNVs/exome/varscan2_pipeline_ONE_v2.0.sh -f $conf -t $tumor -n $normal
echo "fin_varscan2_pipeline_ONE"
chmod -R 776 $output_dir
echo "COMMAND ---run_cvf_csv_facets_v2.1.sh -f $conf -t $tumor"
/home/mibarrola/info/software/packages/in_house_scripts_/CNVs/exome/run_cvf_csv_facets_v2.1.sh -f $conf -t $tumor
echo "fin run_cvf_csv_facets_v2.1"
chmod -R 776 $output_dir
echo "COMMAND --- merge_results_genes_v2.2.sh -f  $conf -t $tumor"
/home/mibarrola/info/software/packages/in_house_scripts_/CNVs/exome/merge/merge_results_genes_v2.2.sh -f  $conf -t $tumor
echo "fin merge_results_genes"
