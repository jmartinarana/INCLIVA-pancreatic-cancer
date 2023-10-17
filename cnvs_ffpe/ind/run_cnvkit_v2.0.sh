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

echo "bam_dir_FFPE $bam_dir_FFPE"
echo "bam_dir_WBC $bam_dir_WBC"
echo "--command $sing $cnvkit batch $bam_dir_FFPE/${tumor}.bam --normal $bam_dir_WBC/${normal}.bam --targets $panel/target.bed --fasta $genome --output-reference my_reference_FIS_Exomas_t.cnn --output-dir $output_dir --diagram --scatter "
$sing $cnvkit batch $bam_dir_FFPE/${tumor}.bam --normal $bam_dir_WBC/${normal}.bam --targets $panel/target.bed --fasta $genome --output-reference my_reference_FIS_Exomas_t.cnn --output-dir $output_dir --diagram --scatter 

#python3 /nfs/home/software/packages/cnvkit/cnvkit.py batch /media/scratch2/FIS_exomas/CNVs/bams/plasmas/*.bam --normal /media/scratch2/FIS_exomas/CNVs/bams/normal/*.bam --targets /nfs/home/panel_designs/HyperExome/target.bed --fasta /nfs/home/references/genomes/human/GRCh38/Homo_sapiens.GRCh38.fa --output-reference my_reference_FIS_Exomas_p.cnn --output-dir ./results/plasmas/ --diagram --scatter -p 40

#python2 /nfs/home/software/packages/cnvkit/cnvkit.py batch /media/scratch2/FIS_exomas/CNVs/bams/organoides/*.bam --normal /media/scratch2/FIS_exomas/CNVs/bams/normal/*.bam --targets /nfs/home/panel_designs/HyperExome/target.bed --fasta /nfs/home/references/genomes/human/GRCh38/Homo_sapiens.GRCh38.fa --output-reference my_reference_FIS_Exomas_o.cnn --output-dir ./results/organoides/ --diagram --scatter -p 40

