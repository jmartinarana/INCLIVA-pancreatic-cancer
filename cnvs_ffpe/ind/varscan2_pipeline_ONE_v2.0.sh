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
    echo -e "-n normal sample"
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

if [ -z "$f" ] || [ -z "$t" ] || [ -z "$n" ];
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
echo "inicio varscan para $tumor y $normal"
# calculo del Data-Ratio 
echo "COMMAND--  $samtools stats $bam_dir_WBC/${normal}.bam > ${output_dir}/${normal}_stats.txt"
$samtools stats $bam_dir_WBC/${normal}.bam > ${output_dir}/${normal}_stats.txt
echo "COMMAND-- $samtools stats $bam_dir_FFPE/${tumor}.bam > ${output_dir}/${tumor}_stats.txt"
$samtools stats $bam_dir_FFPE/${tumor}.bam > ${output_dir}/${tumor}_stats.txt
   
a=`cat  ${output_dir}/${normal}_stats.txt | grep "bases mapped" | head -n 2 | tail -n +2 | awk '{print $5}'`
echo "a=$a"
b=`cat  ${output_dir}/${tumor}_stats.txt | grep "bases mapped" | head -n 2 | tail -n +2 | awk '{print $5}'`
echo "b=$b"
dataratio=`echo "scale=4; ${a}/${b}" | bc`
echo "dataratio= $dataratio"
# calculo del copynumber
echo "COMMAND--  $samtools mpileup -q 1 -f ${genome} $bam_dir_WBC/${normal}.bam $bam_dir_FFPE/${tumor}.bam | $varscan copynumber varScan ${tumor} --output $output_dir_varscan/ --mpileup 1 --min-coverage 100 --min-segment-size 100 --data-ratio ${dataratio}" 
$samtools mpileup -q 1 -f ${genome} $bam_dir_WBC/${normal}.bam $bam_dir_FFPE/${tumor}.bam | $varscan copynumber varScan ${tumor} --output $output_dir_varscan/ --mpileup 1 --min-coverage 100 --min-segment-size 100 --data-ratio ${dataratio} > $output_dir_varscan/${tumor}.copynumber
mv ./${tumor}.copynumber  $output_dir_varscan
#Copycaller y CBS binary segmentation
echo "COMMAND--  $varscan copyCaller $output_dir_varscan/${tumor}.copynumber --output-file $output_dir_varscan/${tumor}.copynumber.called --output-homdel-file $output_dir_varscan/${tumor}.copynumber.called.homdel "
$varscan copyCaller $output_dir_varscan/${tumor}.copynumber --output-file $output_dir_varscan/${tumor}.copynumber.called --output-homdel-file $output_dir_varscan/${tumor}.copynumber.called.homdel
echo "COMMAND--  $sing2 $dnacopy_cnvs $output_dir_varscan ${tumor}"
$sing2 $dnacopy_cnvs $output_dir_varscan/ ${tumor}
echo "fin varscan"


