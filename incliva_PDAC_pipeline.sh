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


function fastq_preprocessing(){
    echo "Running fastq preprocessing"
    date
    echo "COMMAND -- ${fastqc} -o ${5}/QC/rawdata -t ${6} ${1} ${2}"
    ${fastqc} -o ${5}/QC/rawdata -t ${6} --nogroup ${1} ${2}
    date
    if [ "${7}" == "yes" ]
    then
        echo "Processing UMIs"
        date
        echo "COMMAND -- $umitools extract --bc-pattern=${ns_in_umi} --stdin=${8} --read2-in=${1} --stdout=${4}/QC/${sample}_umi_R1.fastq.gz --read2-stdout"	
	$umitools extract --bc-pattern=${ns_in_umi} --stdin=${8} --read2-in=${1} --stdout=${4}/QC/${sample}_umi_R1.fastq.gz --read2-stdout
        echo "COMMAND -- $umitools extract --bc-pattern=${ns_in_umi} --stdin=${8} --read2-in=${2} --stdout=${4}/QC/${sample}_umi_R2.fastq.gz --read2-stdout"
        $umitools extract --bc-pattern=${ns_in_umi} --stdin=${8} --read2-in=${2} --stdout=${4}/QC/${sample}_umi_R2.fastq.gz --read2-stdout



        R1=`echo ${4}/QC/${sample}_umi_R1.fastq.gz`
        R2=`echo ${4}/QC/${sample}_umi_R2.fastq.gz`
        echo "COMMAND -- $cutadapt -j 1 -a $illumina_adapter -A $illumina_adapter -m 1 -o ${4}/QC/${3}_cut_R1.fastq -p ${4}/QC/${3}_cut_R2.fastq ${R1} ${R2}"
        $cutadapt -j 1 -a $illumina_adapter -A $illumina_adapter -m 1  -o ${4}/QC/${3}_cut_R1.fastq -p ${4}/QC/${3}_cut_R2.fastq ${R1} ${R2}
    else
        echo "Sample does not have UMIs"

        echo "COMMAND -- $cutadapt -j 1 -a $illumina_adapter -A $illumina_adapter -m 1 -o ${4}/QC/${3}_cut_R1.fastq -p ${4}/QC/${3}_cut_R2.fastq ${R1} ${R2}"
        $cutadapt -j 1 -a $illumina_adapter -A $illumina_adapter -m 1 -o ${4}/QC/${3}_cut_R1.fastq -p ${4}/QC/${3}_cut_R2.fastq ${R1} ${R2}
    fi
    
    echo "COMMAND -- ${prinseq} -fastq ${4}/QC/${3}_cut_R1.fastq -fastq2 ${4}/QC/${3}_cut_R2.fastq -min_qual_mean 30 -out_good ${4}/QC/filtered_${3} -out_bad ${4}/QC/bad_${3}"
    ${prinseq} -fastq ${4}/QC/${3}_cut_R1.fastq -fastq2 ${4}/QC/${3}_cut_R2.fastq -min_qual_mean 30 -out_good ${4}/QC/filtered_${3} -out_bad ${4}/QC/bad_${3}
    echo "COMMAND -- gzip ${4}/QC/filtered_${3}_1.fastq"
    gzip ${4}/QC/filtered_${3}_1.fastq
    echo "COMMAND -- gzip ${4}/QC/filtered_${3}_2.fastq"
    gzip ${4}/QC/filtered_${3}_2.fastq
    
    echo "COMMAND -- ${fastqc} -o ${5}/QC/rawdata -t ${6} --nogroup ${4}/QC/filtered_${3}_1.fastq.gz ${4}/QC/filtered_${3}_2.fastq.gz"
    ${fastqc} -o ${5}/QC/rawdata -t ${6} --nogroup ${4}/QC/filtered_${3}_1.fastq.gz ${4}/QC/filtered_${3}_2.fastq.gz
    date
}

function mapping_and_bampostprocessing(){
    echo -e "COMMAND -- ${bwa} mem -T 0 -M -R "@RG\tID:AgilentXTHS\tSM:${6}\tPL:Illumina_${seq}" -t ${5} ${7} ${1} ${2} > ${3}/mapping/${6}.sam"
    ${bwa} mem -T 0 -M -R "@RG\tID:AgilentXTHS\tSM:${6}\tPL:Illumina_${seq}" -t ${5} ${7} ${1} ${2} > ${3}/mapping/${6}.sam
    date
    echo -e "COMMAND -- ${samtools} view -@ ${threads} -S -q ${min_mapping_quality} -b -h ${3}/mapping/${6}.sam > ${3}/mapping/${6}_raw.bam"
    ${samtools} view -@ ${threads} -S -q ${min_mapping_quality} -b -h ${3}/mapping/${6}.sam > ${3}/mapping/${6}_raw.bam
    
    echo -e "${samtools} flagstat -@ ${threads} ${3}/mapping/${6}_raw.bam | head -5 | tail -1 | awk '{print $1}'"
    bam_all_reads=`${samtools} flagstat -@ ${threads} ${3}/mapping/${6}_raw.bam | head -5 | tail -1 | awk '{print $1}'`
    echo ${bam_all_reads} > ${3}/mapping/stats_${6}_bam_all_reads

    echo -e "COMMAND -- ${samtools} sort --threads ${5} ${3}/mapping/${6}_raw.bam > ${3}/mapping/${6}_sorted.bam"
    ${samtools} sort --threads ${5} ${3}/mapping/${6}_raw.bam > ${3}/mapping/${6}_sorted.bam
    
    echo -e "COMMAND -- Removing sam file:"
    l1=`cat ${3}/mapping/${6}.sam | wc -l`
    l2=`cat ${3}/mapping/${6}_raw.bam | wc -l`
    state1=`${samtools} quickcheck ${3}/mapping/${6}.sam && echo -e 'Ok'`
    state2=`${samtools} quickcheck ${3}/mapping/${6}_raw.bam && echo -e 'Ok'`

    
    if [ "${state1}" == "Ok" ] && [ "${state2}" == "Ok" ] && [ "${l1}" -gt "${l2}" ]
    then
    	echo -e "\n\nCOMMAND -- Integrity of file: ${3}/mapping/${6}.sam and ${3}/mapping/${6}_raw.bam seem to be correct."	
    	echo -e "COMMAND -- In order to save disk space, this file will be deleted:"
    	echo -e "COMMAND -- rm ${3}/mapping/${6}.sam\n\n"
        rm ${3}/mapping/${6}.sam   
        rm ${3}/mapping/${6}_raw.bam 
    else
    	echo -e "\n\nERROR -- Something went wrong, please check integrity of ${3}/mapping/${6}.sam  and/or ${3}/mapping/${6}_raw.bam. Run will be stopped\n\n"
    	exit
    fi
    
    echo -e "COMMAND -- ${samtools} index -@ ${threads} ${3}/mapping/${6}_sorted.bam"
    ${samtools} index -@ ${threads} ${3}/mapping/${6}_sorted.bam

    echo -e "COMMAND -- ${samtools} view  -@ ${threads} -h -L ${8}/capture.bed -b ${3}/mapping/${6}_sorted.bam > ${3}/mapping/ontarget_${6}_sorted.bam"
    ${samtools} view -@ ${threads} -h -L ${8}/capture.bed -b ${3}/mapping/${6}_sorted.bam > ${3}/mapping/ontarget_${6}_sorted.bam
    
    echo -e "COMMAND -- ${samtools} index -@ ${threads} ${3}/mapping/ontarget_${6}_sorted.bam"
    ${samtools} index -@ ${threads} ${3}/mapping/ontarget_${6}_sorted.bam
       
    if [ "${9}" == "yes" ]
    then
        date
        echo -e "COMMAND -- ${umitools} dedup --output-stats ${3}/QC/${sample}_umidedup_stats -I ${3}/mapping/${6}_sorted.bam -S ${3}/mapping/${sample}_processed.bam -L ${3}/QC/${sample}_umi.log --paired"
        ${umitools} dedup --output-stats ${3}/QC/${sample}_umidedup_stats -I ${3}/mapping/${6}_sorted.bam -S ${3}/mapping/${sample}_processed.bam -L ${3}/QC/${sample}_umi.log --paired
        date
    else
        if [ "${10}" == "yes" ]
        then
            echo -e "No duplicates will be removed because the design is with no UMI and amplicon-based"
	    echo "COMMAND -- ${bamclipper} -b ${3}/mapping/${6}_sorted.bam -p ${8}/amplicons.pe.bed -n 20"
            ${bamclipper} -b ${3}/mapping/${6}_sorted.bam -p ${8}/amplicons.pe.bed -n 20
            echo "COMMAND -- mv ${sample}_sorted.primerclipped.bam ${3}/mapping/${sample}_processed.bam"
            mv ${sample}_sorted.primerclipped.bam ${3}/mapping/${sample}_processed.bam
            echo "COMMAND -- mv ${sample}_sorted.primerclipped.bam.bai ${3}/mapping/${sample}_processed.bam.bai"
            mv ${sample}_sorted.primerclipped.bam.bai ${3}/mapping/${sample}_processed.bam.bai
        else
            echo -e "COMMAND -- ${picard} MarkDuplicates I=${3}/mapping/${6}_sorted.bam O=${3}/mapping/${sample}_processed.bam REMOVE_DUPLICATES=True M=${sample}_duplicates_metrics"
            ${picard} MarkDuplicates I=${3}/mapping/${6}_sorted.bam O=${3}/mapping/${sample}_processed.bam REMOVE_DUPLICATES=True M=${3}/mapping/${sample}_duplicates_metrics
	    

        fi
    fi
    
    
    l3=`zcat ${3}/mapping/${6}_sorted.bam| wc -l`
    l4=`zcat ${3}/mapping/ontarget_${6}_sorted.bam| wc -l`

    state3=`${samtools} quickcheck ${3}/mapping/${6}_sorted.bam && echo -e 'Ok'`
    state4=`${samtools} quickcheck ${3}/mapping/ontarget_${6}_sorted.bam && echo -e 'Ok'`
    
    
    echo -e "COMMAND -- Removing bams_1"    
    echo -e "L3: $l3"
    echo -e "L4: $l4"
    echo -e "$state3"
    echo -e "$state4"
    
    if [ "${state3}" == "Ok" ] && [ "${state4}" == "Ok" ] && [ "${l3}" -gt "${l4}" ]
    then
 	echo -e "\n\nCOMMAND -- Integrity of files: ontarget_${6}_sorted.bam  and ${6}_sorted.bam seem to be correct. Running some statistics\n\n"
 	echo -e "COMMAND -- ${picard} CollectAlignmentSummaryMetrics I=${3}/mapping/${6}_sorted.bam O=${3}/mapping/${6}_alignment_metrics R=${7}"
 	${picard} CollectAlignmentSummaryMetrics I=${3}/mapping/${6}_sorted.bam O=${3}/mapping/${6}_alignment_metrics R=${7}
 	mapping_reads=`head ${3}/mapping/${sample}_alignment_metrics | tail -1 | awk '{print $6}'`
	echo ${mapping_reads} > ${3}/mapping/stats_${sample}_mapping_reads
 	echo -e "mapping_reads, $mapping_reads"

 	echo -e "COMMAND -- ${picard} CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${6}_sorted.bam O=${3}/mapping/ontarget_${6}_alignment_metrics R=${7}"
 	${picard} CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${6}_sorted.bam O=${3}/mapping/ontarget_${6}_alignment_metrics R=${7}
 	ontarget_mapping_reads=`head ${3}/mapping/ontarget_${sample}_alignment_metrics | tail -1 | awk '{print $6}'`
        echo ${ontarget_mapping_reads} > ${3}/mapping/stats_${sample}_ontarget_mapping_reads
 	echo -e "ontarget_mapping_reads, $ontarget_mapping_reads"
    
	
 	echo -e "\n\nCOMMAND -- In order to save disk space, these files will be deleted:"
 	echo -e "COMMAND -- rm ${3}/mapping/${6}_sorted.bam"
        rm ${3}/mapping/${6}_sorted.bam   
        echo -e "COMMAND -- rm ${3}/mapping/ontarget_${6}_sorted.bam\n\n"
        rm ${3}/mapping/${6}_sorted.bam
    else
 	echo -e "\n\nCOMMAND -- Something went wrong, please check integrity of files ${3}/mapping/${6}_sorted.bam and/or ${3}/mapping/ontarget_${6}_sorted.bam. Run will be stopped\n\n"
 	exit
    fi

    echo "COMMAND -- mv ${3}/mapping/${sample}_processed.bam ${3}/mapping/${sample}_processed_tmp.bam"
    mv ${3}/mapping/${sample}_processed.bam ${3}/mapping/${sample}_processed_tmp.bam
 
    echo -e "COMMAND -- ${picard} CleanSam I=${3}/mapping/${sample}_processed_tmp.bam O=${3}/mapping/${sample}_processed_tmp2.bam"
    ${picard} CleanSam I=${3}/mapping/${sample}_processed_tmp.bam O=${3}/mapping/${sample}_processed.bam 

    echo "COMMAND -- ${samtools} index -@ ${threads}  ${3}/mapping/${6}_processed.bam"
    ${samtools} index ${3}/mapping/${6}_processed.bam
    echo "COMMAND -- ${samtools} view -@ ${threads}  -b -h -L ${8}/capture.bed ${3}/mapping/${6}_processed.bam > ${3}/mapping/ontarget_${6}_processed.bam"
    ${samtools} view -@ ${threads}  -b -h -L ${8}/capture.bed ${3}/mapping/${6}_processed.bam > ${3}/mapping/ontarget_${6}_processed.bam
   
    echo -e "COMMAND -- ${gatk} --java-options "-Xmx30g" BaseRecalibrator -I ${3}/mapping/${6}_processed.bam -R ${7} --intervals ${8}/capture.bed --known-sites ${8}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites ${8}/All_20200219.vcf.gz --output ${3}/variant_calling/tmp_${6}.table"
    ${gatk} --java-options "-Xmx30g" BaseRecalibrator -I ${3}/mapping/${6}_processed.bam -R ${7} --intervals ${8}/capture.bed --known-sites ${8}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites ${8}/All_20200219.vcf.gz --output ${3}/variant_calling/tmp_${6}.table
    echo -e "COMMAND -- ${gatk} --java-options "-Xmx30g" ApplyBQSR -R ${7} -I ${3}/mapping/${6}_processed.bam --bqsr-recal-file ${3}/variant_calling/tmp_${6}.table -O ${3}/mapping/${6}_gatk.bam"
    ${gatk} --java-options "-Xmx30g" ApplyBQSR -R ${7} -I ${3}/mapping/${6}_processed.bam --bqsr-recal-file ${3}/variant_calling/tmp_${6}.table -O ${3}/mapping/${6}_gatk.bam


    l5=`zcat ${3}/mapping/${6}_processed.bam |wc -l`
    l6=`zcat ${3}/mapping/ontarget_${6}_processed.bam | wc -l`
    state5=`${samtools} quickcheck ${3}/mapping/${6}_processed.bam && echo -e 'Ok'`
    state6=`${samtools} quickcheck ${3}/mapping/ontarget_${6}_processed.bam && echo -e 'Ok'`
    
    echo -e "$l5, $l6, $state5, $state6"

        
    if [ "${state5}" == "Ok" ] && [ "${state6}" == "Ok" ] && [ "${l5}" -gt "${l6}" ]
    then
 	echo -e "\n\nCOMMAND -- Integrity of files ${3}/mapping/${6}_processed.bam  and ${3}/mapping/ontarget_${6}_processed.bam seem to be correct. Running some statistics\n\n"
 	echo -e "COMMAND -- ${picard} CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${6}_processed.bam O=${3}/mapping/ontarget_${6}_processed_alignment_metrics R=${7}"
 	${picard} CollectAlignmentSummaryMetrics I=${3}/mapping/ontarget_${6}_processed.bam O=${3}/mapping/ontarget_${6}_processed_alignment_metrics R=${7}
 	dedup_reads=`head ${3}/mapping/ontarget_${6}_processed_alignment_metrics | tail -1 | awk '{print $6}'`
	echo ${dedup_reads} > ${3}/mapping/stats_${sample}_dedup_reads
 	echo -e "dedup_reads $dedup_reads"
	
	echo "COMMAND -- ${picard} CollectInsertSizeMetrics I=${3}/mapping/ontarget_${6}_processed.bam O=${3}/QC/${6}_insert H=${3}/QC/${6}_hist"
        ${picard} CollectInsertSizeMetrics I=${3}/mapping/ontarget_${6}_processed.bam O=${3}/QC/${6}_insert H=${3}/QC/${6}_hist

 	raw_data=`zcat ${12} | wc -l | awk '{print $1/2}'`
	if [ "${9}" == "yes" ]
 	then
#            raw_data=`zcat ${12} | wc -l | awk '{print $1/2}'`
            tmp_duplicates=`echo -e "${ontarget_mapping_reads} - ${dedup_reads}" | bc`
            duplicates=`echo -e "scale=4; ${tmp_duplicates}/${raw_data}" | bc | awk '{print $1*100}'`
 	else
 #           raw_data=`zcat ${12} | wc -l | awk '{print $1/2}'`
            if [ "${10}" == "yes" ]
            then
 		duplicates="0"
 		echo -e "duplicates $duplicates"
            else
 		duplicates=`head -8 ${3}/mapping/${sample}_duplicates_metrics | tail -1 | awk '{print $10*100}'`
 		echo -e "duplicates $duplicates"
            fi
 	fi
	
	echo ${raw_data} > ${3}/mapping/stats_${sample}_raw_data
	echo ${duplicates} > ${3}/mapping/stats_${sample}_duplicates
 	echo -e "\n\nCOMMAND -- In order to save disk space, these files will be deleted:"
 	echo -e "COMMAND -- rm ${3}/mapping/${6}_processed.bam"
        rm ${3}/mapping/${6}_processed.bam
        rm ${3}/mapping/ontarget_${6}_processed.bam
        
    else
 	echo -e "\n\nERROR -- Something went wrong, please check integrity of files ${3}/mapping/${6}_processed.bam. Run will be stopped\n\n"
 	exit
    fi
    
    echo -e "COMMAND -- lofreq indelqual -u 20,20 -f ${7} -o ${4}/mapping/${6}.bam ${3}/mapping/${6}_gatk.bam"
    ${lofreq} indelqual -u 20,20 -f ${7} -o ${4}/mapping/${6}.bam ${3}/mapping/${6}_gatk.bam
    echo "COMMAND -- ${samtools} index -@ ${threads} ${4}/mapping/${6}.bam"
    ${samtools} index -@ ${threads}  ${4}/mapping/${6}.bam

    echo "COMMAND -- state7=`${samtools} quickcheck ${3}/mapping/${6}_gatk.bam && echo -e 'Ok'`"
    state7=`${samtools} quickcheck ${3}/mapping/${6}_gatk.bam && echo -e 'Ok'`
    state8=`${samtools} quickcheck ${4}/mapping/${6}.bam && echo -e 'Ok'`


    if [ "${state7}" == "Ok" ] && [ "${state8}" == "Ok" ]
    then
        echo -e "\n\nCOMMAND -- Integrity of file: ${3}/mapping/${6}_gatk.bam seems to be correct."
        echo -e "COMMAND -- NO LO BORRO POR AHORA In order to save disk space, this file will be deleted:"
#        echo -e "COMMAND -- rm ${3}/mapping/${6}_gatk.bam\n\n"
 #       rm ${3}/mapping/${6}_gatk.bam

    else
        echo -e "\n\nCOMMAND -- Something went wrong, please check integrity of ${3}/mapping/${6}_gatk.bam or ${3}/mapping/${6}.bam. Run will be stopped\n\n"
        exit
    fi

    
    date
    echo "End: mapping and bam postprocessing"
    

}

function variant_calling(){
    echo -e "COMMAND -- ${bedtools} sort -i ${7}/capture.bed | ${bedtools} merge > ${2}/mapping/target.bed"
    ${bedtools} sort -i ${7}/capture.bed | ${bedtools} merge > ${2}/mapping/target.bed
    echo -e "COMMAND -- ${samtools} depth -d 9898989898 -a -b ${2}/mapping/target.bed ${1} > ${2}/mapping/${5}_depth"
    ${samtools} depth -d 9898989898 -a -b ${2}/mapping/target.bed ${1} > ${2}/mapping/${5}_depth
    echo -e "COMMAND -- awk '{if ($3 >= '$8') print $1"\t"$2-1"\t"$2}' ${2}/mapping/${5}_depth | ${bedtools} merge > ${3}/QC/mapping/${5}_covered_regions_${8}x.bed"
    awk '{if ($3 >= '$8') print $1"\t"$2-1"\t"$2}' ${2}/mapping/${5}_depth | ${bedtools} merge > ${3}/QC/mapping/${5}_covered_regions_${8}x.bed
    echo -e "COMMAND -- ${gatk} --java-options "-Xmx30g" Mutect2 --base-quality-score-threshold 19 --min-base-quality-score 20 --native-pair-hmm-threads ${4} -R ${6} -I ${1} -L ${3}/QC/mapping/${5}_covered_regions_${8}x.bed -tumor ${5} -O ${2}/variant_calling/tmp1_gatk_${5}.vcf --max-reads-per-alignment-start 0"
    ${gatk} --java-options "-Xmx30g" Mutect2 --native-pair-hmm-threads ${4} -R ${6} -I ${1} -L ${3}/QC/mapping/${5}_covered_regions_${8}x.bed -tumor ${5} -O ${2}/variant_calling/tmp1_gatk_${5}.vcf --max-reads-per-alignment-start 0
    date
    echo -e "COMMAND -- ${bedtools} intersect -header -a ${2}/variant_calling/tmp1_gatk_${5}.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_gatk_${5}.vcf"
    ${bedtools} intersect -header -a ${2}/variant_calling/tmp1_gatk_${5}.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_gatk_${5}.vcf
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/tmp2_gatk_${5}.vcf "
    ${sing} bgzip -f ${2}/variant_calling/tmp2_gatk_${5}.vcf 
    echo -e "COMMAND -- ${sing} tabix -p vcf -f ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz"
    ${sing} tabix -p vcf -f ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz
    echo -e "COMMAND -- ${bcftools} norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | ${bcftools} norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${5}.vcf"
    ${bcftools} norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | ${bcftools} norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${5}.vcf #normaliza indels
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/gatk_${5}.vcf"
    ${sing} bgzip -f ${2}/variant_calling/gatk_${5}.vcf
    echo -e "COMMAND -- ${sing} tabix -p vcf -f ${2}/variant_calling/gatk_${5}.vcf.gz"
    ${sing} tabix -p vcf -f ${2}/variant_calling/gatk_${5}.vcf.gz
    date
    echo -e "COMMAND -- ${inc_select_variants} -i ${2}/variant_calling/gatk_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp.vcf -v ${af} -d ${8}"
    ${sing} ${inc_select_variants} -i ${2}/variant_calling/gatk_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp.vcf -v ${af} -d ${8}
    echo -e "COMMAND -- ${vcfsort} ${2}/variant_calling/${5}_tmp.vcf | ${vcfuniq} > ${2}/variant_calling/${5}.vcf"
    ${vcfsort} ${2}/variant_calling/${5}_tmp.vcf | ${vcfuniq} > ${2}/variant_calling/${5}_normal_tmp.vcf
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}_normal_tmp.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_normal_tmp.vcf.gz"
    ${sing} bgzip -f ${2}/variant_calling/${5}_normal_tmp.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_normal_tmp.vcf.gz
    echo -e "COMMAND -- ${sing} ${inc_anno_bam_info} ${2}/variant_calling/${5}_normal_tmp.vcf.gz ${2}/variant_calling/${5}_normal.vcf Clipped"
    ${sing} ${inc_anno_bam_info} ${2}/variant_calling/${5}_normal_tmp.vcf.gz ${2}/variant_calling/${5}_normal.vcf Clipped
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}_normal.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_normal.vcf.gz"
    ${sing} bgzip -f ${2}/variant_calling/${5}_normal.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_normal.vcf.gz
    date
}


function variant_calling3(){
    rm ${2}/variant_calling/tmp1_lofreq_${5}.vcf
    #echo -e "COMMAND -- lofreq call-parallel -f ${6} -l ${3}/QC/mapping/${5}_covered_regions_${8}x.bed --pp-threads ${threads} --force-overwrite -o ${2}/variant_calling/tmp1_lofreq_${5}.vcf  --call-indels ${1}  ${lofreq_a} ${lofreq_b}"
    #lofreq call-parallel -f ${6} -l ${3}/QC/mapping/${5}_covered_regions_${8}x.bed --pp-threads ${threads} --force-overwrite -o ${2}/variant_calling/tmp1_lofreq_${5}.vcf --call-indels ${1}  ${lofreq_a} ${lofreq_b}
    # usamos lofreq call por error de ejecucion en el cluster
    echo -e "COMMAND -- ${lofreq} call -f ${6} -l ${3}/QC/mapping/${5}_covered_regions_${8}x.bed  --force-overwrite -o ${2}/variant_calling/tmp1_lofreq_${5}.vcf  --call-indels ${1}  ${lofreq_a}  ${lofreq_b}"
    ${lofreq} call -f ${6} -l ${3}/QC/mapping/${5}_covered_regions_${8}x.bed --force-overwrite -o ${2}/variant_calling/tmp1_lofreq_${5}.vcf --call-indels ${1}  ${lofreq_a} ${lofreq_b}

    date
    echo -e "COMMAND -- ${bedtools} intersect -header -a ${2}/variant_calling/tmp1_lofreq_${5}.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_lofreq_${5}.vcf"
    ${bedtools} intersect -header -a ${2}/variant_calling/tmp1_lofreq_${5}.vcf -b ${3}/QC/mapping/${5}_covered_regions_${8}x.bed > ${2}/variant_calling/tmp2_lofreq_${5}.vcf
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/tmp2_lofreq_${5}.vcf"
    ${sing} bgzip -f ${2}/variant_calling/tmp2_lofreq_${5}.vcf
    echo -e "COMMAND -- ${sing} tabix -p vcf -f ${2}/variant_calling/tmp2_lofreq_${5}.vcf.gz"
    ${sing} tabix -p vcf -f ${2}/variant_calling/tmp2_lofreq_${5}.vcf.gz
    echo -e "COMMAND -- ${bcftools} norm -m-any ${2}/variant_calling/tmp2_gatk_${5}.vcf.gz | ${bcftools} norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/gatk_${5}.vcf"
    ${bcftools} norm -m-any ${2}/variant_calling/tmp2_lofreq_${5}.vcf.gz | ${bcftools} norm -Ov --check-ref -w -f ${6} > ${2}/variant_calling/lofreq_${5}.vcf
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/lofreq_${5}.vcf"
    ${sing} bgzip -f ${2}/variant_calling/lofreq_${5}.vcf
    echo -e "COMMAND -- ${sing} tabix -p vcf -f ${2}/variant_calling/lofreq_${5}.vcf.gz"
    ${sing} tabix -p vcf -f ${2}/variant_calling/lofreq_${5}.vcf.gz
    date
    echo -e "COMMAND -- ${sing} ${inc_select_variants_lofreq} -i ${2}/variant_calling/lofreq_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp_lofreq.vcf -v ${af} -d ${8}"
    ${sing} ${inc_select_variants_lofreq} -i ${2}/variant_calling/lofreq_${5}.vcf.gz -o ${2}/variant_calling/${5}_tmp_lofreq.vcf -v ${af} -d ${8}
    echo -e "COMMAND -- ${vcfsort} ${2}/variant_calling/${5}_tmp_lofreq.vcf | ${vcfuniq} > ${2}/variant_calling/${5}.vcf"   
    ${vcfsort} ${2}/variant_calling/${5}_tmp_lofreq.vcf | ${vcfuniq} > ${2}/variant_calling/${5}.vcf 
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}.vcf.gz"
    ${sing} bgzip -f ${2}/variant_calling/${5}.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}.vcf.gz
    echo -e "COMMAND --${sing} ${inc_anno_format_lofreq} ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf"
    ${sing} ${inc_anno_format_lofreq} ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf.gz"
    ${sing} bgzip -f ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf.gz
    echo -e "COMMAND --${sing} ${inc_anno_bam_info} ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf.gz ${2}/variant_calling/${5}_lofreq_anno.vcf Clipped"
    ${sing} ${inc_anno_bam_info} ${2}/variant_calling/${5}_lofreq_anno_tmp.vcf.gz ${2}/variant_calling/${5}_lofreq_anno.vcf Clipped
    echo -e "COMMAND --${sing} bgzip -f ${2}/variant_calling/${5}_lofreq_anno.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_lofreq_anno.vcf.gz"
    ${sing} bgzip -f ${2}/variant_calling/${5}_lofreq_anno.vcf && ${sing} tabix -f -p vcf ${2}/variant_calling/${5}_lofreq_anno.vcf.gz
    echo -e "COMMAND --mv ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/${5}_lofreq.vcf.gz"
    mv ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/${5}_lofreq.vcf.gz
    echo -e "COMMAND --rm ${2}/variant_calling/${5}.vcf.gz.tbi"
    rm ${2}/variant_calling/${5}.vcf.gz.tbi
    if [ "${9}" == "yes" ]
    then
        echo -e "COMMAND -- ${inc_gatk_realiable} ${2}/variant_calling/${5}_normal.vcf.gz ${2}/variant_calling/${5}_normal_conf.vcf Low_confidence"
        ${sing} ${inc_gatk_realiable} ${2}/variant_calling/${5}_normal.vcf.gz ${2}/variant_calling/${5}_normal_conf.vcf Low_confidence
    	echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}_normal_conf.vcf"
        ${sing} bgzip -f ${2}/variant_calling/${5}_normal_conf.vcf
    	echo -e "COMMAND -- ${sing} tabix -p vcf -f ${2}/variant_calling/${5}_normal_conf.vcf.gz"
        ${sing} tabix -p vcf -f ${2}/variant_calling/${5}_normal_conf.vcf.gz
        echo -e "COMMAND -- ${inc_gatk_realiable} ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${5}_lofreq_conf.vcf High_confidence"
        ${sing} ${inc_gatk_realiable} ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${5}_lofreq_conf.vcf High_confidence
    	echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}_lofreq_conf.vcf"
        ${sing} bgzip -f ${2}/variant_calling/${5}_lofreq_conf.vcf
    	echo -e "COMMAND -- ${sing} tabix -p vcf -f ${2}/variant_calling/${5}_lofreq_conf.vcf.gz"
        ${sing} tabix -p vcf -f ${2}/variant_calling/${5}_lofreq_conf.vcf.gz
        echo -e "COMMAND -- vcf-isec -f -n +1 ${2}/variant_calling/${5}_lofreq_conf.vcf.gz ${2}/variant_calling/${5}_normal_conf.vcf.gz > ${2}/variant_calling/${5}_gt.vcf"
        ${vcfisec} -f -n +1 ${2}/variant_calling/${5}_lofreq_conf.vcf.gz ${2}/variant_calling/${5}_normal_conf.vcf.gz > ${2}/variant_calling/${5}_gt.vcf
    else
    	echo -e "COMMAND -- ${vcfisec} -f -n +1 ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${5}_normal.vcf.gz > ${2}/variant_calling/${5}_gt.vcf"
        ${vcfisec} -f -n +1 ${2}/variant_calling/${5}_lofreq_anno.vcf.gz ${2}/variant_calling/${5}_normal.vcf.gz > ${2}/variant_calling/${5}_gt.vcf
    fi
    echo -e "COMMAND -- ${inc_reformat_genotype} -i ${2}/variant_calling/${5}_gt.vcf -o ${2}/variant_calling/${5}.vcf -v ${af} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}"
    ${sing} ${inc_reformat_genotype} -i ${2}/variant_calling/${5}_gt.vcf -o ${2}/variant_calling/${5}.vcf -v ${af} -m ${max_freq_for_heterozygous} -d ${min_depth_to_filter_variant}
    echo -e "COMMAND -- ${sing} bgzip -f ${2}/variant_calling/${5}.vcf"
    ${sing} bgzip -f ${2}/variant_calling/${5}.vcf
    echo -e "COMMAND -- ${sing} tabix ${2}/variant_calling/${5}.vcf.gz"
    ${sing} tabix ${2}/variant_calling/${5}.vcf.gz
    echo -e "COMMAND -- ${sing} ${inc_check_empty_vcf} ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/gatk_${5}.vcf.gz"
    ${sing} ${inc_check_empty_vcf} ${2}/variant_calling/${5}.vcf.gz ${2}/variant_calling/gatk_${5}.vcf.gz
    date
}

function statistics(){
    
    echo "COMMAND -- ${samtools} flagstat -@ ${threads} ${1} | ${samtools} flagstat ${1} | awk '{print $6}' | head -9 | tail -1 | sed 's/(//' | sed 's/%//'"
    per_read=`${samtools} flagstat -@ ${threads}  ${1} | awk '{print $6}' | head -9 | tail -1 | sed 's/(//' | sed 's/%//'`
    echo "per_read $per_read"
    echo ${per_read} > ${3}/mapping/stats_${sample}_per_read
    mapping_reads=`cat ${3}/mapping/stats_${sample}_mapping_reads | awk '{print $1}'`
    echo ${mapping_reads} > ${3}/mapping/stats_${sample}_mapping_reads
    raw_data=`cat ${3}/mapping/stats_${sample}_raw_data`
    echo "raw_data $raw_data"
    echo ${raw_data} > ${3}/mapping/stats_${sample}_raw_data
    all_mapped=`echo "scale=4; ${mapping_reads}/${raw_data}" | bc | awk '{print $1*100}'`
    echo "all_mapped $all_mapped"
    echo ${all_mapped} > ${3}/mapping/stats_${sample}_all_mapped
    
    prinseq_reads=`zcat ${3}/QC/filtered_${sample}_1.fastq.gz | wc -l | awk '{print $1/2}'`
    echo "prinseq_reads $prinseq_reads"
    echo ${prinseq_reads} > ${3}/mapping/stats_${sample}_prinseq_reads
 
    prinseq_per=`echo "scale=4; ${prinseq_reads}/${raw_data}" | bc | awk '{print 100-($1*100)}'`
    echo "prinseq_per $prinseq_per"
    echo ${prinseq_per} > ${3}/mapping/stats_${sample}_prinseq_per

    bam_target_reads=`head ${3}/mapping/ontarget_${sample}_alignment_metrics | tail -1 | awk '{print $6}'`
    echo "bam_target_reads $bam_target_reads"
    echo ${bam_target_reads} >  ${3}/mapping/stats_${sample}_bam_target_reads

    dedup_reads=`cat ${3}/mapping/stats_${sample}_dedup_reads | awk '{print $1}'`
    usable_reads=`echo "scale=4; ${dedup_reads}/${raw_data}" | bc | awk '{print $1*100}'`
    echo "usable_reads $usable_reads"
    echo ${usable_reads} >  ${3}/mapping/stats_${sample}_usable_reads   



    specificity=`echo "scale=4; ${bam_target_reads}/${mapping_reads}" | bc | awk '{print $1*100}'`
    echo "specificity $specificity"
    echo ${specificity}  >  ${3}/mapping/stats_${sample}_specificity

    bam_all_reads=`cat ${3}/mapping/stats_${sample}_bam_all_reads`
    echo -e "bam_all_reads, $bam_all_reads"
    #echo ${bam_all_reads}    >  ${3}/mapping/stats_${sample}_bam_all_reads

    all_vars=`zcat ${3}/variant_calling/${5}.vcf.gz | grep -v "#" -c`
    echo -e "all_vars $all_vars"
    echo ${all_vars}  >  ${3}/mapping/stats_${sample}_all_vars

    mean_insert_tmp=`head -8 ${3}/QC/${5}_insert | sed 's/,/./g' | tail -1 | awk '{print $6}'`
    echo -e "mean_insert_tmp $mean_insert_tmp"
    echo ${mean_insert_tmp}  >  ${3}/mapping/stats_${sample}_mean_insert_tmp
 
    mean_insert=`echo -e "scale=2; ${mean_insert_tmp}/1" | bc`
    echo -e "mean_insert, $mean_insert"
    echo ${mean_insert}  >  ${3}/mapping/stats_${sample}_mean_insert

    deviation_tmp=`head -8 ${3}/QC/${5}_insert | sed 's/,/./g' | tail -1 | awk '{print $7}'`
    echo -e "deviation_tmp, $deviation_tmp"
    echo ${deviation_tmp} >  ${3}/mapping/stats_${sample}_deviation_tmp

    deviation=`echo -e "scale=2; ${deviation_tmp}/1" | bc`
    echo -e "deviation, $deviation"
    echo ${deviation} >  ${3}/mapping/stats_${sample}_deviation
    
    design=`awk '{print $3-$2+1}' ${3}/mapping/target.bed | ${sing} suma`
    echo -e "design $design"
    echo ${design}  >  ${3}/mapping/stats_${sample}_design

    echo -e "COMMAND -- bedtools subtract -a ${3}/mapping/target.bed -b ${4}/QC/mapping/${5}_covered_regions_${10}x.bed > ${3}/QC/${5}_non_covered_${10}x.bed"
    ${bedtools} subtract -a ${7}/capture.bed -b ${4}/QC/mapping/${5}_covered_regions_${10}x.bed > ${3}/QC/${5}_non_covered_${10}x.bed

    echo -e "COMMAND -- bedtools intersect -a ${3}/QC/${5}_non_covered_${10}x.bed -b ${7}/inc_db.vcf -wao | awk '{if ($5 ~ /chr/) print $12"\t"$1"\t"$2"\t"$3"\t"$4; else print "--\t"$1"\t"$2"\t"$3"\t"$4}' | sed 's/.*CDS=//' | sed 's/;CNT.*\tchr/\tchr/' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > ${4}/QC/mapping/${5}_non_covered_regions_${10}x.bed"
    ${bedtools} intersect -a ${3}/QC/${5}_non_covered_${10}x.bed -b ${7}/inc_db.vcf -wao | awk '{if ($5 ~ /chr/) print $12"\t"$1"\t"$2"\t"$3"\t"$4; else print "--\t"$1"\t"$2"\t"$3"\t"$4}' | sed 's/.*CDS=//' | sed 's/;CNT.*\tchr/\tchr/' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > ${4}/QC/mapping/${5}_non_covered_regions_${10}x.bed
    #if [ "${12}" != "no" ]
    #    then
    #    samplemed=`awk '{if ($3>=25) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    #else
    #    samplemed=`awk '{if ($3>=30) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    #fi

    sample10=`awk -v depth="$d_10" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    echo -e "sample10 $sample10"
    echo ${sample10}  >  ${3}/mapping/stats_${sample}_sample10

    depth10=`echo -e "scale=4; ${sample10}/${design}" | bc | awk '{print $1*100}'`
    echo -e "depth10 $depth10"
    echo ${depth10}  >  ${3}/mapping/stats_${sample}_depth10

    sample20=`awk -v depth="$d_20" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    echo -e "sample20 $sample20"
    echo ${sample20}  >  ${3}/mapping/stats_${sample}_sample20

    depth20=`echo -e "scale=4; ${sample20}/${design}" | bc | awk '{print $1*100}'`
    echo -e "depth20 $depth20"
    echo ${depth20}  >  ${3}/mapping/stats_${sample}_depth20

    sample50=`awk -v depth="$d_50" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    echo -e "sample50 $sample50"
    echo ${sample50}  >  ${3}/mapping/stats_${sample}_sample50

    depth50=`echo -e "scale=4; ${sample50}/${design}" | bc | awk '{print $1*100}'`
    echo -e "depth50 $depth50"
    echo ${depth50}  >  ${3}/mapping/stats_${sample}_depth50

    sample100=`awk -v depth="$d_100" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    echo -e "sample100 $sample100"
    echo ${sample100}  >  ${3}/mapping/stats_${sample}_sample100

    depth100=`echo -e "scale=4; ${sample100}/${design}" | bc | awk '{print $1*100}'`
    echo -e "depth100 $depth100"
    echo ${depth100}  >  ${3}/mapping/stats_${sample}_depth100

    sample250=`awk -v depth="$d_250" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    echo -e "sample250 $sample250"
    echo ${sample250}  >  ${3}/mapping/stats_${sample}_sample250

    depth250=`echo -e "scale=4; ${sample250}/${design}" | bc | awk '{print $1*100}'`
    echo -e "depth250 $depth250"
    echo ${depth250}  >  ${3}/mapping/stats_${sample}_depth250

    sample500=`awk -v depth="$d_500" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
    echo -e "sample500 $sample500"
    echo ${sample500}  >  ${3}/mapping/stats_${sample}_sample500

    depth500=`echo -e "scale=4; ${sample500}/${design}" | bc | awk '{print $1*100}'`
    echo -e "depth500 $depth500"
    echo ${depth500}  >  ${3}/mapping/stats_${sample}_depth500




#    samplemin=`awk -v depth="$d_min" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
#    echo -e "samplemin $samplemin"
#    echo ${samplemin}  >  ${3}/mapping/stats_${sample}_samplemin

#    depthmin=`echo -e "scale=4; ${samplemin}/${design}" | bc | awk '{print $1*100}'`
#    echo -e "depthmin $depthmin"
#    echo ${depthmin}  >  ${3}/mapping/stats_${sample}_depthmin

#    samplemax=`awk -v depth="$d_max" '{if ($3>=depth) print $1"\t"$2-1"\t"$2}' ${3}/mapping/${5}_depth | ${bedtools} merge | awk '{print $3-$2+1}' | ${sing} suma`
#    echo -e "samplemax $samplemax"
#    echo ${samplemax}  >  ${3}/mapping/stats_${sample}_samplemax

#    depthmax=`echo -e "scale=4; ${samplemax}/${design}" | bc | awk '{print $1*100}'`
#    echo -e "depthmax $depthmax"
#    echo ${depthmax}  >  ${3}/mapping/stats_${sample}_depthmax

    echo -e "COMMAND -- ${sing} ${inc_calculate_IQR_and_coverage_plot} ${3}/mapping/${5}_depth ${3}/QC/${5}_R_coverage_summary ${4}/QC/mapping/${5}_coverage_plot.pdf ${5}"
    ${sing} ${inc_calculate_IQR_and_coverage_plot} ${3}/mapping/${5}_depth ${3}/QC/${5}_R_coverage_summary ${4}/QC/mapping/${5}_coverage_plot.pdf ${5}
    echo -e "COMMAND -- ${sing} ${inc_af} ${3}/variant_calling ${4}/QC/mapping ${5}"
    ${sing} ${inc_af} ${3}/variant_calling ${4}/QC/mapping  ${5}

    first_quantile=`tail -1 ${3}/QC/${5}_R_coverage_summary | awk '{print $2}'`
    echo -e "first_quantile $first_quantile"
    echo ${first_quantile}  >  ${3}/mapping/stats_${sample}_first_quantile

    third_quantile=`tail -1 ${3}/QC/${5}_R_coverage_summary | awk '{print $5}'`
    echo -e "third_quantile $third_quantile"
    echo ${third_quantile}  >  ${3}/mapping/stats_${sample}_third_quantile

    median_coverage=`tail -1 ${3}/QC/${5}_R_coverage_summary | awk '{print $3}'`
    echo -e "median_coverage $median_coverage"
    echo ${median_coverage}  >  ${3}/mapping/stats_${sample}_median_coverage

    rpu=`${sing} ${inc_calculate_RPU} ${3}/QC/${5}_umidedup_stats_per_umi_per_position.tsv`
    echo -e "rpu $rpu"
    echo ${rpu}  >  ${3}/mapping/stats_${sample}_rpu

    echo  "${5},${raw_data},${prinseq_per},${all_mapped},${per_read},${duplicates},${usable_reads},${specificity},${first_quantile},${third_quantile},${median_coverage},${rpu},${depthmin},${depthmed},${depthmax},${mean_insert},${deviation},${all_vars}" >> ${4}/QC/mapping/${5}_global_mapping_stats.csv

    echo -e "COMMAND -- global mapping stats in ${4}/QC/mapping/${5}_global_mapping_stats.csv"
    date
}

function check_list_of_transcripts(){
    while IFS='' read -r line || [[ -n "${line}" ]]; do
        transcript=`echo -e ${line} | awk '{print $3}'`
        gene=`echo -e ${line} | awk '{print $1}'`
        has_final_info=`grep $transcript ${3}`
        has_intermediate_info=`grep $transcript ${2}`
        if [  -z "${has_final_info}" ]
        then
            if [ ! -z "${has_intermediate_info}" ]
            then
                echo -e "-------------ERROR----SOMETHING MAY HAPPENS WITH THE LIST OF TRANSCRIPTS-----------------"
                echo -e "----------------CHECK THAT NO RESULTS ARE OBTAINED FOR ${line}-----------------------"
            fi
        fi
    done < "${1}"
}

function annotate_variants(){
     echo -e "COMMAND -- vep -i ${1} --use_transcript_ref --force_overwrite --offline --cache --merged --dir ${5}/v102 --dir_plugins ${5}/v102/Plugins --fasta ${4} --sift b --polyphen b --gene_phenotype --numbers --vcf_info_field ANN --terms SO --hgvs --shift_hgvs 1 --canonical --biotype --xref_refseq --max_af --af_esp --af_gnomad --pubmed --minimal --force_overwrite --vcf --af --symbol --domains --fork ${7} --buffer_size 200 -o ${2}/variant_calling/${6}_annotated_vep.vcf"
     ${vep} -i ${1} --use_transcript_ref --force_overwrite --offline --cache --merged --dir ${5}/v102 --dir_plugins ${5}/v102/Plugins --fasta ${4} --sift b --polyphen b --gene_phenotype --numbers --vcf_info_field ANN --terms SO --hgvs --shift_hgvs 1 --canonical --biotype --xref_refseq --max_af --af_esp --af_gnomad --pubmed --minimal --force_overwrite --vcf --af --symbol --domains --fork ${7} --buffer_size 200 -o ${2}/variant_calling/${6}_annotated_vep.vcf
     date
     echo -e "COMMAND -- ${bcftools} annotate -a ${5}/cancerhospots.bed.gz -c CHROM,FROM,TO,HS -h <(echo -e 'INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">') ${2}/variant_calling/${6}_annotated_vep.vcf > ${2}/variant_calling/${6}_annotated_uniq_h.vcf"
     ${bcftools} annotate -a ${5}/cancerhospots.bed.gz -c CHROM,FROM,TO,HS -h <(echo -e '##INFO=<ID=HS,Number=1,Type=String,Description="Tag variants in cancer hotspots cancerhotspot.org">') ${2}/variant_calling/${6}_annotated_vep.vcf > ${2}/variant_calling/${6}_annotated_uniq_h.vcf
     echo -e "COMMAND -- ${snpSift} annotate ${5}/AllMutations_COSMIC_v91.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf"
     ${snpSift} annotate ${5}/AllMutations_COSMIC_v91.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h.vcf > ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf
     echo -e "COMMAND -- ${inc_mv_cosmic_id_to_annotation} -i ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf -o ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf"
     ${sing} ${inc_mv_cosmic_id_to_annotation} -i ${2}/variant_calling/${6}_annotated_uniq_h_c.vcf -o ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf
     echo -e "COMMAND -- ${snpSift} annotate ${5}/mixed_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${2}/variant_calling/${6}_annotated_mixed.vcf"
     ${snpSift} annotate ${5}/mixed_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${2}/variant_calling/${6}_annotated_mixed.vcf
     echo -e "COMMAND -- ${snpSift} annotate ${5}/inc_db.vcf.gz ${2}/variant_calling/${6}_annotated_uniq_h_c_cok.vcf > ${2}/variant_calling/${6}_annotated_db.vcf"
     ${snpSift} annotate ${5}/inc_db.vcf.gz ${2}/variant_calling/${6}_annotated_mixed.vcf > ${2}/variant_calling/${6}_annotated_db.vcf
     echo -e "COMMAND -- ${sing} ${inc_add_conservation} ${2}/variant_calling/${6}_annotated_db.vcf ${2}/variant_calling/${6}_annotated_cons.vcf ${5}/conservation.bed"
     ${sing} ${inc_add_conservation} ${2}/variant_calling/${6}_annotated_db.vcf ${2}/variant_calling/${6}_annotated_cons.vcf ${5}/conservation.bed 
     echo -e "COMMAND -- ${vcfuniq} ${2}/variant_calling/${6}_annotated_cons.vcf > ${3}/variant_calling/${6}_annotated.vcf"
     ${vcfuniq} ${2}/variant_calling/${6}_annotated_cons.vcf > ${3}/variant_calling/${6}_annotated.vcf
     echo -e "COMMAND -- ${sing} bgzip -f ${3}/variant_calling/${6}_annotated.vcf"
     ${sing} bgzip -f ${3}/variant_calling/${6}_annotated.vcf
    
    date
     if [ "${9}" == "yes" ]
         then
         echo -e "COMMAND -- ${inc_vcf_to_csv_lofreq_conf} -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m ${af} -n 1 -r ${5}/list_of_transcripts"
         	#${sing} ${inc_vcf_to_csv_lofreq_conf} -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m ${af} -n 1 -r ${5}/list_of_transcripts
		${sing} ${inc_vcf_to_csv_lofreq_conf} -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m ${af} -n 1 
         else
         echo -e "COMMAND -- ${inc_vcf_to_csv_lofreq_mixed} -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m ${af} -n 1 -r ${5}/list_of_transcripts"
         ${sing} ${inc_vcf_to_csv_lofreq_mixed} -f ${3}/variant_calling/${6}_annotated.vcf.gz -o ${3}/variant_calling/${6}_annotated.csv -m ${af} -n 1 -r ${5}/list_of_transcripts
    fi
    echo -e "COMMAND -- ${sing} ${inc_aa_3to1} ${3}/variant_calling/${6}_annotated.csv"
    ${sing} ${inc_aa_3to1} ${3}/variant_calling/${6}_annotated.csv
    echo -e "COMMAND -- ${inc_af_merge} ${2}/variant_calling ${3}/QC/mapping"
    ${sing} ${inc_af_merge} ${2}/variant_calling ${3}/QC/mapping
    echo -e "COMMAND --  ${inc_samples_contamination} ${3} 0.7"
    ${sing} ${inc_samples_contamination} ${3}
    echo -e "COMMAND -- ${sing} ${inc_non_covered_regions} ${3}/QC/mapping/${6}_global_mapping_stats.csv ${3} ${2} ${5}/target.bed ${8}"
    ${sing} ${inc_non_covered_regions} ${3}/QC/mapping/${6}_global_mapping_stats.csv ${3} ${2} ${5}/target.bed ${8}
    #echo -e "COMMAND -- ${sing} ${inc_gene_cover} ${3} ${2} ${5} ${8}"
    #${sing} ${inc_gene_cover} ${3} ${2} ${5} ${8}
    if [ "${9}" == "yes" ]
        then
        
        echo -e ${10}
        if [ "${10}" != "no" ]
            then
            if [ "${BRCA}" == "yes" ]
            then
                echo -e "COMMAND -- ${inc_amplicons_BRCA_universal} ${3} ${10} ${8} ${5}"
                ${sing} ${inc_amplicons_BRCA_universal} ${3} ${10} ${8} ${5}
            else
                echo -e "COMMAND -- ${inc_amplicons_TP53_universal} ${3} ${10} ${8}"
                ${sing} ${inc_amplicons_TP53_universal} ${3} ${10} ${8}
            fi
            else
            echo -e "COMMAND -- ${inc_amplicons_cover} ${3} ${5}/capture_amplicons.bed ${8}"
            ${sing} ${inc_amplicons_cover} ${3} ${5}/capture_amplicons.bed ${8}
            fi        
    fi
    
    if [ "${9}" == "yes" ]
        then
        echo -e "COMMAND -- Rscript ${inc_csv_format_conf} ${3}/variant_calling/${6}_annotated.csv ${3}/variant_calling"
        ${sing} ${inc_csv_format_conf} ${3}/variant_calling/${6}_annotated.csv ${3}/variant_calling
        else
        echo -e "COMMAND -- ${inc_csv_format} ${3}/variant_calling/${6}_annotated.csv ${3}/variant_calling"
        ${sing} ${inc_csv_format} ${3}/variant_calling/${6}_annotated.csv ${3}/variant_calling
    fi
    date
 #   echo -e "check_list_of_transcripts ${5}/list_of_transcripts ${2}/variant_calling/${6}_annotated_vep.vcf ${3}/variant_calling/${6}_annotated.csv ${2}/variant_calling"
#    check_list_of_transcripts ${5}/list_of_transcripts ${2}/variant_calling/${6}_annotated_vep.vcf ${3}/variant_calling/${6}_annotated.csv ${2}/variant_calling
    
}

function usage() {
    echo -e "Usage: $0"
    echo -e "This script runs the pipeline for target seq data"
    echo -e "Mandatory parameters:"
    echo -e "-f configuration file"
    echo -e "-s sample"
    1>&2; exit 1;
}



while getopts "f:s:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        s)
                s=${OPTARG} ;;
        *)
                usage ;;
    esac
done

if [ -z "$f" ] || [ -z "$s" ];
then
        echo -e "ERROR: -f is a mandatory parameter"
        echo -e "ERROR: -s is a mandatory parameter"
        usage
        exit
fi


#Se ejecuta al inicio, crea los enlaces simbólicos y los directorios, si no existen

date
echo -e "Your command: "$@
source ${f}
samples=${s}

checkDirectories ${tmp_dir}/mapping
checkDirectories ${tmp_dir}/variant_calling
checkDirectories ${tmp_dir}/QC
checkDirectories ${analysis_dir}/mapping
checkDirectories ${analysis_dir}/variant_calling
checkDirectories ${analysis_dir}/QC/mapping
checkDirectories ${analysis_dir}/QC/rawdata
ln -s -f ${INC_somatic} ${panel}/
ln -s -f ${INC_somatic}.tbi ${panel}/
ln -s -f ${vep_db} ${panel}/
ln -s -f ${dbSNP} ${panel}/
ln -s -f ${dbSNP}.tbi ${panel}/
ln -s -f ${GATK_indels} ${panel}/
ln -s -f ${GATK_indels}.tbi ${panel}/
ln -s -f ${COSMIC} ${panel}/
ln -s -f ${COSMIC}.tbi ${panel}/
ln -s -f ${HOTSPOTS} ${panel}/
ln -s -f ${HOTSPOTS}.tbi ${panel}/
ln -s -f ${INC_somatic_unzip} ${panel}/
ln -s -f ${MIXED} ${panel}/
ln -s -f ${MIXED}.tbi ${panel}/
${bedtools} intersect -a ${GATK_small} -b ${panel}/capture.bed -header > ${panel}/exac_roche.vcf
${sing} bgzip -f ${panel}/exac_roche.vcf && ${sing} tabix -p vcf -f ${panel}/exac_roche.vcf.gz

#Lo siguiente crea crea un csv con el header por defecto de las mapping stats. Atención, si se lanza el pipe incompleto "machaca" las stats anteriores. Habría que cambiarlo.
echo -e "Sample,Rawdata,%LowQReads,%MappedReads,%Paired_reads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,RPU,Nt_10x,Nt_25x,Nt_50x,Mean_insert,Insert_SD,NumVars" > ${analysis_dir}/QC/mapping/${name}_global_mapping_stats.csv
if [ ${panel} != "/nfs/home/panel_designs/TP53" ]    #Este diseño de panel no está en el cluster del cipf
    then
        echo -e "------NO AMPLICONS------"
        echo -e "NO ES TP53"
        p="no"
    else
        echo -e "Sample,Rawdata,%LowQReads,%MappedReads,%Paired_reads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_10x,Nt_30x,Nt_50x,Mean_insert,Insert_SD,NumVars" > ${analysis_dir}/QC/mapping/${name}_global_mapping_stats.csv   #crea un csv por defectode las mapping stats si es amplicon
        echo -e "ES TP53"
        p=`echo -e ${samples} | awk '{split($0,a,","); print a[1]}'`
    fi
if [ ${BRCA} == "yes" ]
    then
        p=`echo -e ${samples} | awk '{split($0,a,","); print a[1]}'`
    fi
sample_list=`echo -e ${samples} | sed 's/,/ /gi'`
cwd=`pwd`

for sample in ${sample_list[@]}
do
    echo -e "Running analysis pipeline on sample ${sample}"
    R1=`ls ${rawdata}/${sample}_1*fast*`
    R2=`ls ${rawdata}/${sample}_2*fast*`
    if [ ${has_umi} == "no" ]
    then
	echo -e "------$sample DOES NOT HAVE UMIs------"
    else
	echo -e "------$sample HAS UMIs------"
	R1=`ls ${rawdata}/${sample}_1*fast*`    # vigilar esto tema nomenclatura ----------
    	R2=`ls ${rawdata}/${sample}_2*fast*`
	umi=`ls ${rawdata}/${sample}_UMI*fast*`
    fi
    if [ ${amplicon} == "yes" ]
    then
        echo -e "------DESIGN FOR SAMPLE $sample IS AMPLICON-BASED------"
    else
        echo -e "------DESIGN FOR SAMPLE $sample IS NOT AMPLICON-BASED------"
    fi
    echo -e "COMMAND -- checkDirectories ${tmp_dir}/QC/${sample}"
    checkDirectories ${tmp_dir}/QC/${sample}  
    echo "COMMAND -- fastq_preprocessing ${R1} ${R2} ${sample} ${tmp_dir} ${analysis_dir} ${threads} ${has_umi} ${umi}"
    #fastq_preprocessing ${R1} ${R2} ${sample} ${tmp_dir} ${analysis_dir} ${threads} ${has_umi} ${umi}
    echo -e "COMMAND -- mapping_and_bampostprocessing ${tmp_dir}/QC/filtered_${sample}_1.fastq.gz ${tmp_dir}/QC/filtered_${sample}_2.fastq.gz ${tmp_dir} ${analysis_dir} ${threads} ${sample} ${genome} ${panel} ${has_umi} ${amplicon} ${analysis_dir}/mapping/${sample}.bam ${R1} ${name} ${cov} ${p}"
    #mapping_and_bampostprocessing ${tmp_dir}/QC/filtered_${sample}_1.fastq.gz ${tmp_dir}/QC/filtered_${sample}_2.fastq.gz ${tmp_dir} ${analysis_dir} ${threads} ${sample} ${genome} ${panel} ${has_umi} ${amplicon} ${analysis_dir}/mapping/${sample}.bam ${R1} ${name} ${cov} ${p}
    echo -e "COMMAND -- variant_calling ${analysis_dir}/mapping/${sample}.bam ${tmp_dir} ${analysis_dir} ${threads} ${sample} ${genome} ${panel} ${cov} ${amplicon}"
    #variant_calling ${analysis_dir}/mapping/${sample}.bam ${tmp_dir} ${analysis_dir} ${threads} ${sample} ${genome} ${panel} ${cov} ${amplicon}
    echo -e "COMMAND -- variant_calling3 ${analysis_dir}/mapping/${sample}.bam ${tmp_dir} ${analysis_dir} ${threads} ${sample} ${genome} ${panel} ${cov} ${amplicon}"
    variant_calling3 ${analysis_dir}/mapping/${sample}.bam ${tmp_dir} ${analysis_dir} ${threads} ${sample} ${genome} ${panel} ${cov} ${amplicon}   
    statistics ${analysis_dir}/mapping/${sample}.bam ${R1} ${tmp_dir} ${analysis_dir} ${sample} ${name} ${panel} ${genome} ${has_umi} ${cov} ${amplicon} ${p}
    #list_files+=" ${tmp_dir}/variant_calling/${sample}.vcf.gz "
done

date
#echo "Start Combining files"
#echo -e "COMMAND -- $vcfmerge -c none $list_files | ${sing} bgzip > ${analysis_dir}/variant_calling/${name}.vcf.gz"
#$vcfmerge -c none $list_files | ${sing} bgzip > ${analysis_dir}/variant_calling/${name}.vcf.gz     

#echo -e "COMMAND -- ${sing} tabix -p vcf -f ${analysis_dir}/variant_calling/${name}.vcf.gz"
#${sing} tabix -p vcf -f ${analysis_dir}/variant_calling/${name}.vcf.gz

#echo -e "COMMAND -- annotate_variants ${analysis_dir}/variant_calling/${name}.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} ${name} ${threads} ${cov} ${amplicon} ${p}"
#annotate_variants ${analysis_dir}/variant_calling/${name}.vcf.gz ${tmp_dir} ${analysis_dir} ${genome} ${panel} ${name} ${threads} ${cov} ${amplicon} ${p}

#date
#echo -e "Start: Checking for errors"
#echo -e "COMMAND -- ${sing} ${inc_nohup_error} ${analysis_dir}"
#${sing} ${inc_nohup_error} ${analysis_dir}


echo -e "END OF PIPELINE"
date


