#! /bin/bash -x

sample_names=$1

for sample in ${sample_names}; do
	zcat ${sample}.vcf.gz | grep -v "#" | awk -v var="${sample}" '{print var"\t"$1"\t"$2"\t"$4"\t"$5}' >> $2
