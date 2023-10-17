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
    echo -e "-t sample name"
    1>&2; exit 1;
}

while getopts "f:t:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        t)
                t=${OPTARG} ;;
        *)
                usage ;;
    esac
done
if [ -z "$f" ] || [ -z "$t" ] ;
then
        echo -e "ERROR: -f is a mandatory parameter"
        echo -e "ERROR: -t is a mandatory parameter"
        usage
        exit
fi



date
echo -e "Your command: "$@
source ${f}
tumor=${t}
echo "Sample  is ${tumor}"
echo "COMMAND --  $from_facets_vcf_to_csv  $output_dir_facets/${tumor}.vcf.gz $output_dir/tmp_${tumor}.bed"
$from_facets_vcf_to_csv  $output_dir_facets/${tumor}.vcf.gz $output_dir/tmp_${tumor}.bed
echo "COMMAND -- $bedtools sort -i $output_dir/tmp_${tumor}.bed > $output_dir/${tumor}.bed"
$bedtools sort -i $output_dir/tmp_${tumor}.bed > $output_dir/${tumor}.bed




#for i in ${!ffpe[@]}; do
#   echo "ffpe  is ${ffpe[$i]}"
#   echo "se va a ejecutar from_facets_vcf_to_csv"
#   echo "COMMAND --  $from_facets_vcf_to_csv  $output_dir_facets/${ffpe[$i]}.vcf.gz $output_dir/tmp_${ffpe[$i]}.bed"
#   $from_facets_vcf_to_csv  $output_dir_facets/${ffpe[$i]}.vcf.gz $output_dir/tmp_${ffpe[$i]}.bed
#   echo "COMMAND -- $bedtools sort -i $output_dir/tmp_${ffpe[$i]}.bed > $output_dir/${ffpe[$i]}.bed"
#   $bedtools sort -i $output_dir/tmp_${ffpe[$i]}.bed > $output_dir/${ffpe[$i]}.bed   
  
#done
