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
    echo -e "-t tumor bam file"
    echo -e "-n normal bam file"
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
source ${f}
tumor=${t}
normal=${n}

echo "COMMAND--  $cnv_facets -t $bam_dir_FFPE/${tumor}.bam -n $bam_dir_WBC/${normal}.bam -N 40 -T $panel/target.bed -vcf $vcf_dir_FFPE/${tumor}.vcf.gz -g hg38 -o  $output_dir_facets/${tumor}"
$cnv_facets -t $bam_dir_FFPE/${tumor}.bam -n $bam_dir_WBC/${normal}.bam -N 40 -T $panel/target.bed -vcf $vcf_dir_FFPE/${tumor}.vcf.gz -g hg38 -o $output_dir_facets/${tumor}


#for i in ${!ffpe[@]}; do
#  echo "ffpe  is ${ffpe[$i]}"
#  echo "wbc  is ${wbc[$i]}"

#  echo "COMMAND-- $sing $cnv_facets -t $bam_dir_FFPE/${ffpe[$i]}.bam -n $bam_dir_WBC/${wbc[$i]}.bam -N 40 -T $panel/target.bed -vcf $vcf_dir_FFPE/${ffpe[$i]}.vcf.gz -g hg38 -o  $output_dir_facets/${ffpe[$i]}"
#  $cnv_facets -t $bam_dir_FFPE/${ffpe[$i]}.bam -n $bam_dir_WBC/${wbc[$i]}.bam -N 40 -T $panel/target.bed -vcf $vcf_dir_FFPE/${ffpe[$i]}.vcf.gz -g hg38 -o $output_dir_facets/${ffpe[$i]}

#done


