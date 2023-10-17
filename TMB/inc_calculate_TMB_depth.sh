#!/bin/bash
source $1

mkdir -p ${tmp_dir}/QC/mapping
mkdir -p ${tmp_dir}/VCFs
mkdir -p ${tmp_dir}/variant_calling

for sample in ${tumor_plasma[@]}
do

	echo "----------------------- STARTING ANALYSIS FOR ${sample} -----------------------"
	echo -e "COMMAND -- samtools depth -@ ${threads} -d 9898989898 -a ${analysis_dir_tumor_plasma}/mapping/${sample}.bam  -b ${panel}/target.bed | awk '{if ($3 >= ${cov_tmb}) print $1"\t"$2-1"\t"$2}' | bedtools merge > ${tmp_dir}/QC/mapping/covered_intervals_${sample}_${cov_tmb}.bed"
	#cat ${tmp_tumor_plasma}/mapping/${sample}_depth |  awk '{if ($3 >= '$cov_tmb') print $1"\t"$2-1"\t"$2}' | bedtools merge > ${tmp_dir}/QC/mapping/covered_intervals_${sample}_${cov_tmb}.bed

	echo "COMMAND -- vcf-subset -e -a -c ${sample} ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_cancer_somatic_annotated_vep_oncokb.vcf > ${tmp_dir}/VCFs/${sample}.vcf"
	vcf-subset -e -a -c ${sample} ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_cancer_somatic_annotated_vep_oncokb.vcf > ${tmp_dir}/VCFs/${sample}.vcf

	echo "COMMAND -- bedtools intersect -a ${tmp_dir}/VCFs/${sample}.vcf -b ${tmp_dir}/QC/mapping/covered_intervals_${sample}_${cov_tmb}.bed -header > ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf"
	bedtools intersect -a ${tmp_dir}/VCFs/${sample}.vcf -b ${tmp_dir}/QC/mapping/covered_intervals_${sample}_${cov_tmb}.bed -header > ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf

	echo "COMMAND -- bgzip -f ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf"
	bgzip -f ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf

	echo "COMMAND -- tabix -p vcf ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf.gz"
	tabix -p vcf ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf.gz

	list_of_vcfs+=" ${tmp_dir}/QC/mapping/${sample}_${cov_tmb}.vcf.gz"
done

echo "COMMAND -- vcf-merge -c none -t ${list_of_vcfs} | bgzip > ${tmp_dir}/variant_calling/all_samples_TMB.vcf.gz"
vcf-merge -c none -t ${list_of_vcfs} | bgzip > ${tmp_dir}/variant_calling/all_samples_TMB.vcf.gz 

echo "COMMAND -- python2.7 /home/mibarrola/UMP444_exomas_pancreas_1108/scripts/inc_vcf_to_csv.py -f ${tmp_dir}/variant_calling/all_samples_TMB.vcf.gz -o ${analysis_dir_tumor_plasma}/${name}_${cancer_sample_is}_cancer_somatic_annotated_vep_oncokb_${cov_tmb}_TMB.csv -m ${af_t} -n ${max_freq_for_heterozygous}"
python2.7 /media/scratch2/FIS_exomas/TANDA3/TMB/Sheila/new/inc_vcf_to_csv.py -f ${tmp_dir}/variant_calling/all_samples_TMB.vcf.gz -o ${analysis_dir_tumor_plasma}/${name}_${cancer_sample_is}_cancer_somatic_annotated_vep_oncokb_${cov_tmb}_TMB.csv -m ${af_t} -n ${max_freq_for_heterozygous}
