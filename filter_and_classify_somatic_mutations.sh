#!/bin/bash
source $1
mkdir -p ${tmp_dir}/variant_calling
mkdir -p ${analysis_dir}

function annotate_variants(){
	echo -e "COMMAND -- ${vep} -i ${1} --use_transcript_ref --force_overwrite --offline --cache --merged --dir ${vep_db} --dir_plugins ${vep_db}/Plugins --fasta ${4} --sift b --polyphen b --gene_phenotype --numbers --vcf_info_field ANN --terms SO --hgvs --shift_hgvs 1 --canonical --biotype --xref_refseq --max_af --af_gnomad --pubmed --minimal --force_overwrite --vcf --af --symbol --fork ${7} --buffer_size 1000 -o ${2}/variant_calling/${6}_annotated_vep.vcf --pick --pick_order rank,tsl,appris,canonical"
	${vep} -i ${1} --use_transcript_ref --force_overwrite --offline --cache --merged --dir ${vep_db} --dir_plugins ${vep_db}/Plugins --fasta ${4} --sift b --polyphen b --gene_phenotype --numbers --vcf_info_field ANN --terms SO --hgvs --shift_hgvs 1 --canonical --biotype --xref_refseq --max_af --af_gnomad --pubmed --minimal --force_overwrite --vcf --af --symbol --fork ${7} --buffer_size 1000 -o ${2}/variant_calling/${6}_annotated_vep.vcf --pick --pick_order rank,tsl,appris,canonical 
	date
	
	echo -e "COMMAND -- ${bcftools} annotate -a ${5}/cancerhospots.bed.gz -c CHROM,FROM,TO,HS -h <(echo -e 'INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">') ${2}/variant_calling/${6}_annotated_vep.vcf > ${2}/variant_calling/${6}_annotated_uniq_h.vcf"
	${bcftools} annotate -a ${5}/cancerhospots.bed.gz -c CHROM,FROM,TO,HS -h <(echo -e '##INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">') ${2}/variant_calling/${6}_annotated_vep.vcf > ${2}/variant_calling/${6}_annotated_uniq_h.vcf
	
	echo -e "COMMAND -- ${snpSift} annotate ${COSMIC} ${2}/variant_calling/${6}_annotated_uniq_h.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf"
	${snpSift} annotate ${COSMIC} ${2}/variant_calling/${6}_annotated_uniq_h.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf
	

	echo "awk \'{if ($5 != \"*\") print $0}\' ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c_tmp0.vcf"
	awk '{if ($5 != "*") print $0}' ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c_tmp0.vcf

	echo -e "COMMAND -- ${inc_mv_cosmic_id_to_annotation} -i ${2}/variant_calling/${6}_annotated_uniq_h_c_tmp0.vcf -o ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf"
	${sing} ${inc_mv_cosmic_id_to_annotation} -i ${2}/variant_calling/${6}_annotated_uniq_h_c_tmp0.vcf -o ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf
	
	echo -e "COMMAND -- ${snpSift} annotate ${5}/mixed_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${2}/variant_calling/${6}_annotated_mixed.vcf"
	${snpSift} annotate ${5}/mixed_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${2}/variant_calling/${6}_annotated_mixed.vcf
	
	echo -e "COMMAND -- ${snpSift} annotate ${5}/inc_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${2}/variant_calling/${6}_annotated_db.vcf"
	${snpSift} annotate ${5}/inc_db.vcf.gz ${2}/variant_calling/${6}_annotated_mixed.vcf > ${2}/variant_calling/${6}_annotated_db.vcf
###############	echo -e "COMMAND -- ${sing} ${inc_add_conservation} ${2}/variant_calling/${6}_annotated_db.vcf ${2}/variant_calling/${6}_annotated_cons.vcf ${5}/conservation.bed"
#	${sing} ${inc_add_conservation} ${2}/variant_calling/${6}_annotated_db.vcf ${2}/variant_calling/${6}_annotated_cons.vcf ${5}/conservation.bed
#	echo -e "COMMAND -- ${vcfuniq} ${2}/variant_calling/${6}_annotated_cons.vcf > ${2}/variant_calling/${6}_annotated_cons_uq.vcf"
#	${vcfuniq} ${2}/variant_calling/${6}_annotated_cons.vcf > ${2}/variant_calling/${6}_annotated_cons_uq.vcf
#	echo -e "COMMAND -- ${vcfuniq} ${2}/variant_calling/${6}_annotated_cons_uq.vcf > ${2}/variant_calling/${6}_annotated_vep.vcf"
#	${vcfuniq} ${2}/variant_calling/${6}_annotated_cons_uq.vcf > ${2}/variant_calling/${6}_annotated_vep.vcf
#	echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${6}_annotated_vep.vcf"
#	${sing} bgzip -f ${2}/variant_calling/${6}_annotated_vep.vcf
#	echo "${snpSift} annotate -info EVIDENCE_LEVEL,ONCOKB_CLASIFICATION,REVIEWED_VARIANT,TREATMENT ${ONCOKB} ${2}/variant_calling/${6}_annotated_vep.vcf.gz > ${2}/variant_calling/${6}_annotated_vep_oncokb.vcf.gz"
###############	${snpSift} annotate -info hSl,Oncogenic,VARIANT_IN_ONCOKB,Treatment ${ONCOKB} ${2}/variant_calling/${6}_annotated_vep.vcf.gz > ${2}/variant_calling/${6}_annotated_vep_oncokb.vcf.gz
	echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${6}_annotated_db.vcf"
	${sing} bgzip -f ${2}/variant_calling/${6}_annotated_db.vcf

	echo -e "COMMAND -- ${snpSift} annotate -info GENEINFO,hSl,Oncogenic,VARIANT_IN_ONCOKB,Treatment ${ONCOKB} ${2}/variant_calling/${6}_annotated_db.vcf.gz | sed 's/Other Biomarkers/Other_Biomarkers/' > ${2}/variant_calling/${6}_annotated_vep_oncokb.vcf"
	${snpSift} annotate -info GENEINFO,hSl,Oncogenic,VARIANT_IN_ONCOKB,Treatment ${ONCOKB} ${2}/variant_calling/${6}_annotated_db.vcf.gz  | sed 's/Other Biomarkers/Other_Biomarkers/' > ${2}/variant_calling/${6}_annotated_vep_oncokb.vcf
	
	echo -e "COMMAND -- ${inc_vcf_to_csv} -f ${2}/variant_calling/${6}_annotated_vep_oncokb.vcf -o ${3}/${6}_annotated_vep_oncokb.csv -m ${af_t} -n 1 "
	${inc_vcf_to_csv} -f ${2}/variant_calling/${6}_annotated_vep_oncokb.vcf -o ${3}/${6}_annotated_vep_oncokb.csv -m ${af_t} -n 1

	echo -e "COMMAND -- ${inc_muts_classification}  ${3}/${6}_annotated_vep_oncokb.csv ${3}/${6}_annotated_and_classified.csv ${oncogenes}"
	${inc_muts_classification}  ${3}/${6}_annotated_vep_oncokb.csv ${3}/${6}_annotated_and_classified.csv ${oncogenes}
}

number_of_samples=`echo ${#names[@]} -1 | bc`
for index in $(seq 0 ${number_of_samples}) 
do
	echo ""
	echo ""
	echo "--------------------------------------------------------------------------------------------------------------------------------------------"
	date
	echo "STARTING ANALYSIS FOR SAMPLE ${names[${index}]}"	
	paired_samples=(${normal[${index}]} ${tumor_plasma[${index}]})
	sample_type=("normal" "tumor_or_plasma")
	analysis_directories=(${analysis_normal} ${analysis_tumor_plasma})
	tmp_directories=(${tmp_normal} ${tmp_tumor_plasma})
	posible_frequencies=(${af_n} ${af_t})
	echo "---- Calling variants with lofreq at a very low frequency ${a_and_b_lofreq_relaxed} in known pathogenic sites ----"
	# This code is repeated twice, once for normal samples and once for either plasma or FFPE samples
	for either_normal_tumor in $(seq 0 1)
	do
		echo "---- Calling lofreq variants in pathogenic regions in ${sample_type[${either_normal_tumor}]} samples ----"
		echo -e "COMMAND -- ${bedtools} intersect -a ${ONCOKB_PATHOGENIC} -b ${analysis_directories[${either_normal_tumor}]}/QC/mapping/${paired_samples[${either_normal_tumor}]}_covered_regions_${cov}x.bed > ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_oncokb_pathogenic_regions_${cov}.bed"
		#${bedtools} intersect -a ${ONCOKB_PATHOGENIC} -b ${analysis_directories[${either_normal_tumor}]}/QC/mapping/${paired_samples[${either_normal_tumor}]}_covered_regions_${cov}x.bed > ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_oncokb_pathogenic_regions_${cov}.bed

		echo -e "COMMAND -- ${lofreq} -f ${genome} -l ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_oncokb_pathogenic_regions_${cov}.bed --force-overwrite -o ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf --call-indels ${a_and_b_lofreq_relaxed} ${analysis_directories[${either_normal_tumor}]}/mapping/${paired_samples[${either_normal_tumor}]}.bam"
		#${lofreq} -f ${genome} -l ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_oncokb_pathogenic_regions_${cov}.bed --force-overwrite -o ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf --call-indels ${a_and_b_lofreq_relaxed} ${analysis_directories[${either_normal_tumor}]}/mapping/${paired_samples[${either_normal_tumor}]}.bam

		echo -e "COMMAND -- ${bgzip} -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf"
		#${bgzip} -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf

		echo -e "COMMAND -- ${tabix} -p vcf -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz"
		${tabix} -p vcf -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz

		#Filtering variants and formatting lofreq output
		echo -e "COMMAND -- ${sing} ${inc_select_variants_lofreq} -i ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz -o ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf -v ${posible_frequencies[${either_normal_tumor}]} -d ${min_depth_to_filter_variant}"
		${sing} ${inc_select_variants_lofreq} -i ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp0_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz -o ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf -v ${posible_frequencies[${either_normal_tumor}]} -d ${min_depth_to_filter_variant}
		
		echo -e "COMMAND -- ${vcfsort} ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf | ${vcfuniq} > ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp2_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf"
		${vcfsort} ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf | ${vcfuniq}  > ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp2_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf

		echo -e "COMMAND -- ${sing} bgzip -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp2_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf"
		${sing} bgzip -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp2_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf
		
		echo -e "COMMAND --${sing} ${inc_anno_format_lofreq} ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp2_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp3_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf"
		${sing} ${inc_anno_format_lofreq} ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp2_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp3_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf

		echo -e "COMMAND -- ${sing} bgzip -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp3_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf"
		${sing} bgzip -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp3_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf
		
		echo -e "COMMAND --${sing} ${inc_anno_bam_info} ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp3_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf Clipped"
		${sing} ${inc_anno_bam_info} ${tmp_directories[${either_normal_tumor}]}/variant_calling/tmp3_${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf Clipped

		echo -e "COMMAND --${sing} bgzip -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf"
                ${sing} bgzip -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf

                echo -e "COMMAND -- ${tabix} -p vcf -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf.gz"
                ${tabix} -p vcf -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf.gz

                echo "${bedtools} intersect -a ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf.gz -b ${ONCOKBVCF} -header | vcf-sort | ${vcfuniq} | bgzip -f > ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz"
                ${bedtools} intersect -a ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site_onetmp.vcf.gz -b ${ONCOKBVCF} -header | vcf-sort |  ${vcfuniq} | bgzip -f > ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz

                echo "${tabix} -p vcf -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz"
                ${tabix} -p vcf -f  ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz


		# Joining lofreq calls with relaxed parameters with previous lofreq calls
		echo -e "COMMAND -- ${vcfisec}  -f -n +1 ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_anno.vcf.gz ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz > ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_final.vcf"
		${vcfisec}  -f -n +1 ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_anno.vcf.gz ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_oncokb_site.vcf.gz > ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_final.vcf

		echo -e "COMMAND -- ${bgzip} -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_final.vcf"
		${bgzip} -f ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_final.vcf

		echo -e "COMMAND -- ${tabix} -f -p vcf ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_final.vcf.gz"
		${tabix} -f -p vcf ${tmp_directories[${either_normal_tumor}]}/variant_calling/${paired_samples[${either_normal_tumor}]}_lofreq_final.vcf.gz	

		echo ""
	done

	echo "---- Lofreq variants. Filtering by depth and strand bias. Classifying variants ----"
	# Selecting germline variants
	echo -e "COMMAND -- ${vcfisec}  -f -n +2 ${tmp_normal}/variant_calling/${normal[${index}]}_lofreq_final.vcf.gz ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_lofreq_final.vcf.gz | ${snpSift} filter \" (DP4[2] != 0) & (DP4[3] != 0) \"  | ${snpSift} filter \" ((DP4[2] + DP4[3]) >= ${sum_of_alt_allele_depths}) \" | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_germline.vcf.gz"
	${vcfisec}  -f -n +2 ${tmp_normal}/variant_calling/${normal[${index}]}_lofreq_final.vcf.gz ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_lofreq_final.vcf.gz | ${snpSift} filter " (DP4[2] != 0) & (DP4[3] != 0) "  | ${snpSift} filter " ((DP4[2] + DP4[3]) >= ${sum_of_alt_allele_depths}) " | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_germline.vcf.gz

	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_germline.vcf.gz"
        tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_germline.vcf.gz

	# Selecting somatic variants
        echo -e "COMMAND -- ${vcfisec} -f -c ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_lofreq_final.vcf.gz ${tmp_normal}variant_calling/${normal[${index}]}_lofreq_final.vcf.gz | ${snpSift} filter \" (DP4[2] != 0) & (DP4[3] != 0) \" | ${snpSift} filter \" ((DP4[2] + DP4[3]) >= ${sum_of_alt_allele_depths}) \" | ${snpSift} filter \" (SBSB =< ${max_SB_lofreq} ) \" | ${snpSift} filter \"( GEN[0].AF > ${af_t} )\" > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_sm.vcf"
	${vcfisec} -f -c ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_lofreq_final.vcf.gz ${tmp_normal}variant_calling/${normal[${index}]}_lofreq_final.vcf.gz | ${snpSift} filter " (DP4[2] != 0) & (DP4[3] != 0) " | ${snpSift} filter " ((DP4[2] + DP4[3]) >= '${sum_of_alt_allele_depths}') " | ${snpSift} filter " (SB <= '${max_SB_lofreq}') " | ${snpSift} filter "( GEN[0].AF > '${af_t}' )"   > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_sm.vcf

	# Select only those variants in plasma/tumor that are within covered areas in normal samples
        echo -e "COMMAND -- ${bedtools} intersect -v -header -a ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_sm.vcf -b ${analysis_normal}/QC/mapping/${normal[${index}]}_non_covered_regions_10x.bed | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_somatic.vcf.gz"
	${bedtools} intersect -v -header -a ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_sm.vcf -b ${analysis_normal}/QC/mapping/${normal[${index}]}_non_covered_regions_10x.bed | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_somatic.vcf.gz

	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_somatic.vcf.gz"
	tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_somatic.vcf.gz

	echo "---- GATK variants. Filtering by depth and strand bias. Classifying variants ----"
	echo -e "COMMAND -- ${vcfisec}  -f -n +2 ${tmp_normal}/variant_calling/${normal[${index}]}_normal.vcf.gz ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_normal.vcf.gz | ${snpSift} filter \" (GEN[0].SB[2] > 0) & (GEN[0].SB[3] > 0)\" | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_germline.vcf.gz"
	${vcfisec}  -f -n +2 ${tmp_normal}/variant_calling/${normal[${index}]}_normal.vcf.gz ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_normal.vcf.gz | ${snpSift} filter "(GEN[0].SB[2] > 0) & (GEN[0].SB[3] > 0)" |  bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_germline.vcf.gz

        echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_germline.vcf.gz"
        tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_germline.vcf.gz

	echo -e "COMMAND -- ${vcfisec} -f -c ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_normal.vcf.gz ${tmp_normal}/variant_calling/${normal[${index}]}_normal.vcf.gz | ${snpSift} filter \"( GEN[0].AF > ${af_t} ) & (GEN[0].SB[2] > 0) & (GEN[0].SB[3] > 0)\"  > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_gatk_sm.vcf"
        ${vcfisec} -f -c ${tmp_tumor_plasma}/variant_calling/${tumor_plasma[${index}]}_normal.vcf.gz ${tmp_normal}/variant_calling/${normal[${index}]}_normal.vcf.gz | ${snpSift} filter "( GEN[0].AF > '${af_t}' ) & (GEN[0].SB[2] > 0) & (GEN[0].SB[3] > 0)"  > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_gatk_sm.vcf

	# Select only those variants in plasma/tumor that are within covered areas in normal samples
        echo -e "COMMAND -- ${bedtools} intersect -v -header -a ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_gatk_sm.vcf -b ${analysis_normal}/QC/mapping/${normal[${index}]}_non_covered_regions_10x.bed | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_somatic.vcf.gz"
        ${bedtools} intersect -v -header -a ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_gatk_sm.vcf -b ${analysis_normal}/QC/mapping/${normal[${index}]}_non_covered_regions_10x.bed | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_somatic.vcf.gz
	
	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_somatic.vcf.gz"
	tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_somatic.vcf.gz

	echo "---- Combining GATK + LOFREQ for germline variants ----"
	echo -e "COMMAND -- ${vcfisec} -f -n +1 ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_germline.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_germline.vcf.gz > ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf"
        ${vcfisec} -f -n +1 ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_germline.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_germline.vcf.gz > ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf

	if [[ "${cancer_sample_is}" == "plasma" ]]
	then
		echo "---- The cancer sample is PLASMA ----"
		echo -e "COMMAND -- ${sing} ${inc_reformat_genotype} -i ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf -o ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_other_somatic.vcf -v ${af_t} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}"
		${sing} ${inc_reformat_genotype} -i ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf -o ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_other_somatic.vcf -v ${af_t} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}

		echo -e "COMMAND -- ${snpSift} filter \"( GEN[0].AF < ${freq_other_somatic})\" ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_other_somatic.vcf  > ${tmp_dir}/${names[${index}]}_other_somatic.vcf"
		${snpSift} filter "( GEN[0].AF < ${freq_other_somatic})" ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_other_somatic.vcf > ${tmp_dir}/${names[${index}]}_other_somatic.vcf

		echo -e "COMMAND -- bgzip -f ${tmp_dir}/${names[${index}]}_other_somatic.vcf"
		bgzip -f ${tmp_dir}/${names[${index}]}_other_somatic.vcf

		echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_other_somatic.vcf.gz"
		tabix -p vcf -f ${tmp_dir}/${names[${index}]}_other_somatic.vcf.gz
		
		list_of_other_somatic+=" ${tmp_dir}/${names[${index}]}_other_somatic.vcf.gz"
		echo ""
	fi

	echo -e "COMMAND -- ${sing} ${inc_reformat_genotype} -i ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf -o ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_germline.vcf -v ${af_n} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}"
	${sing} ${inc_reformat_genotype} -i ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf -o ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_germline.vcf -v ${af_n} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}
	
	echo -e "COMMAND -- ${snpSift} filter \"( GEN[0].AF > ${af_n})\" ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_germline.vcf > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf"
	${snpSift} filter "( GEN[0].AF > ${af_n})" ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_lofreq_gatk_germline.vcf > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf
	
	echo -e "COMMAND -- bgzip -f  ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf"
	bgzip -f  ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf
	
	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf.gz"
	tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf.gz	
	

	echo "---- Combining GATK + LOFREQ for somatic variants ----"
        echo "${vcfisec} -f -n +1 ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_somatic.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_somatic.vcf.gz > ${tmp_dir}/tmp_${tumor_plasma[${index}]}_lofreq_gatk_somatic.vcf"
        ${vcfisec} -f -n +1 ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_somatic.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_gatk_somatic.vcf.gz > ${tmp_dir}/tmp_${tumor_plasma[${index}]}_lofreq_gatk_somatic.vcf

        echo -e "COMMAND -- ${sing} ${inc_reformat_genotype} -i ${tmp_dir}/tmp_${tumor_plasma[${index}]}_lofreq_gatk_somatic.vcf -o ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf -v ${af_t} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}"
        ${sing} ${inc_reformat_genotype} -i ${tmp_dir}/tmp_${tumor_plasma[${index}]}_lofreq_gatk_somatic.vcf -o ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf -v ${af_t} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}

        echo -e "COMMAND -- bgzip -f  ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf"
        bgzip -f  ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf

        echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf.gz"
        tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf.gz


	echo "---- Combining GATK + Lofreq with Haplotype for germline variants ----"
	echo "${bcftools} norm -m-any ${sarek_dir}/${normal[${index}]}/HaplotypeCaller/HaplotypeCaller_${normal[${index}]}.vcf.gz | ${bcftools} norm -Ov --check-ref -w -f ${genome} | bgzip -f > ${tmp_dir}/HaplotypeCaller_${normal[${index}]}_tmp.vcf.gz"
	${bcftools} norm -m-any ${sarek_dir}/${normal[${index}]}/HaplotypeCaller/HaplotypeCaller_${normal[${index}]}.vcf.gz | ${bcftools} norm -Ov --check-ref -w -f ${genome} | bgzip -f > ${tmp_dir}/HaplotypeCaller_${normal[${index}]}_tmp.vcf.gz

	echo "${tabix} -p vcf -f ${tmp_dir}/HaplotypeCaller_${normal[${index}]}_tmp.vcf.gz"
	${tabix} -p vcf -f ${tmp_dir}/HaplotypeCaller_${normal[${index}]}_tmp.vcf.gz

	# Removing weird values from haplotype caller: those with AD=* or with both reference and alternative 0,0
	echo "${snpSift} filter \" ( GEN[0].AD[1] > 0 | GEN[0].AD[1] > 0) & (GEN[0].AD[0] != * & GEN[0].AD[1] != *) \" ${tmp_dir}/HaplotypeCaller_${normal[${index}]}_tmp.vcf.gz | bgzip -c > ${tmp_dir}/HaplotypeCaller_clean_${normal[${index}]}_tmp.vcf.gz"
	${snpSift} filter " ( GEN[0].AD[1] > 0 | GEN[0].AD[1] > 0) & (GEN[0].AD[0] != * & GEN[0].AD[1] != *) " ${tmp_dir}/HaplotypeCaller_${normal[${index}]}_tmp.vcf.gz | bgzip -c > ${tmp_dir}/HaplotypeCaller_clean_${normal[${index}]}_tmp.vcf.gz
	echo "${tabix} -p vcf -f ${tmp_dir}/HaplotypeCaller_clean_${normal[${index}]}_tmp.vcf.gz"
	${tabix} -p vcf -f ${tmp_dir}/HaplotypeCaller_clean_${normal[${index}]}_tmp.vcf.gz

	echo -e "COMMAND -- ${vcfisec} -f -n +1 ${tmp_dir}/HaplotypeCaller_clean_${normal[${index}]}_tmp.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf.gz  | bgzip -f > ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz"
	${vcfisec} -f -n +1 ${tmp_dir}/HaplotypeCaller_clean_${normal[${index}]}_tmp.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_germline.vcf.gz  | bgzip -f > ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz

	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz"
	tabix -p vcf -f ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz
	
	echo "---- Final selection of germline variants ----"
	echo -e "COMMAND -- ${snpSift} filter \"( GEN[0].AF > ${af_n} ) | ( AF > ${af_n} )\" ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz |  bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_germline.vcf.gz"
	${snpSift} filter "( GEN[0].AF > ${af_n}) | ( AF > ${af_n} )" ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_germline.vcf.gz
	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_germline.vcf.gz"
	tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_germline.vcf.gz

	echo "---- Final selection of somatic variants ----"
	echo -e "COMMAND -- ${vcfisec} -f -c ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf.gz ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz | ${snpSift} filter \"( GEN[0].AF > ${af_t} )\" | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz"
        ${vcfisec} -f -c ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_lofreq_gatk_somatic.vcf.gz ${tmp_dir}/tmp_${names[${index}]}_${cancer_sample_is}_all_germline_plus_HC.vcf.gz | ${snpSift} filter "( GEN[0].AF > '${af_t}' )" | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz

	echo -e "COMMAND -- ${tabix} -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz"
        ${tabix} -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz	

	if [[ "${cancer_sample_is}" == "plasma" ]]
        then	
		echo -e "COMMAND -- ${vcfisec} -f -c ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz ${tmp_dir}/${names[${index}]}_other_somatic.vcf.gz | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp1_sm.vcf.gz"
		${vcfisec} -f -c ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz ${tmp_dir}/${names[${index}]}_other_somatic.vcf.gz | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp1_sm.vcf.gz

		echo "COMMAND -- mv ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp1_sm.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz"
		mv ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp1_sm.vcf.gz ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz

		echo -e "COMMAND -- ${tabix} -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz"
	        ${tabix} -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz
	fi

	echo -e "COMMAND -- ${bedtools} intersect -v -header -a ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz -b ${analysis_normal}/QC/mapping/${normal[${index}]}_non_covered_regions_10x.bed | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_somatic.vcf.gz"
	${bedtools} intersect -v -header -a ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_tmp_sm.vcf.gz -b ${analysis_normal}/QC/mapping/${normal[${index}]}_non_covered_regions_10x.bed | bgzip -f > ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_somatic.vcf.gz

	echo -e "COMMAND -- tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_somatic.vcf.gz"
	tabix -p vcf -f ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_somatic.vcf.gz

	list_of_somatic+=" ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_somatic.vcf.gz"
	list_of_germline+=" ${tmp_dir}/${names[${index}]}_${cancer_sample_is}_germline.vcf.gz"
done

echo "---- Merging individual files for all samples ----"
echo -e "COMMAND -- ${vcfmerge} -c none ${list_of_somatic} | ${sing} bgzip -f > ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_somatic_all.vcf.gz"
${vcfmerge} -c none ${list_of_somatic} | ${sing} bgzip -f  > ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_somatic_all.vcf.gz

echo -e "COMMAND -- ${vcfmerge} -c none ${list_of_germline} | ${sing} bgzip -f > ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz"
${vcfmerge} -c none ${list_of_germline} | ${sing} bgzip -f > ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz

echo -e "COMMAND -- ${sing} tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_somatic_all.vcf.gz"
tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_somatic_all.vcf.gz

echo -e "COMMAND -- ${sing} tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz"
tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz

echo "---- Annotating cancer-somatic variants ----"
echo -e "COMMAND -- annotate_variants ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_somatic_all.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} ${name}_${cancer_sample_is}_cancer_somatic ${threads}"
annotate_variants ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_somatic_all.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} "${name}_${cancer_sample_is}_cancer_somatic" ${threads}

if [ -z ${existing_germline_vcf} ]
then
	echo "---- Annotating germline variants ----"
	echo -e "COMMAND -- annotate_variants ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} ${name}_${cancer_sample_is}_germline ${threads}"
	annotate_variants ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} "${name}_${cancer_sample_is}_germline" ${threads}
else
	echo "---- Concatenating germline files from previous analysis in either FFPE or plasma samples ----"
	echo -e "COMMAND -- ${vcfisec}  -f -n +1 ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz ${existing_germline_vcf} > ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf"
	${vcfisec}  -f -n +1 ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_all.vcf.gz ${existing_germline_vcf} > ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf
	
	echo -e "COMMAND -- ${sing} bgzip -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf"	
	${sing} bgzip -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf

	echo -e "COMMAND -- ${sing} tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf.gz"
	${sing} tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf.gz
	
	echo "---- Annotating germline variants ----"
	echo -e "COMMAND -- annotate_variants ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} ${name}_${cancer_sample_is}_germline_FFPE_and_plasma ${threads}"
        annotate_variants ${tmp_dir}/variant_calling/${name}_${cancer_sample_is}_germline_FFPE_and_plasma.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} "${name}_${cancer_sample_is}_germline_FFPE_and_plasma" ${threads}
fi



if [[ "${cancer_sample_is}" == "plasma" ]]
then
        echo -e "COMMAND -- ${vcfmerge} -c none ${list_of_other_somatic} | ${sing} bgzip -f > ${tmp_dir}/variant_calling/${name}_other_somatic_all.vcf.gz"
        ${vcfmerge} -c none ${list_of_other_somatic} | ${sing} bgzip -f > ${tmp_dir}/variant_calling/${name}_other_somatic_all.vcf.gz

	echo -e "COMMAND -- ${sing} tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_other_somatic_all.vcf.gz"
	${sing} tabix -p vcf -f ${tmp_dir}/variant_calling/${name}_other_somatic_all.vcf.gz

	echo "---- Annotating other somatic variants ----"
	echo -e "annotate_variants ${tmp_dir}/variant_calling/${name}_other_somatic_all.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} "${name}_other_somatic" ${threads}"
	annotate_variants ${tmp_dir}/variant_calling/${name}_other_somatic_all.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} "${name}_other_somatic" ${threads}
fi
echo "END OF PIPELINE"
date
