#!/usr/bin/bash

# creates the batch file for GATK database array creation

VCF_DIR='/expanse/projects/sebat1/g2mh_data/gvcfs'
OUT_PREFIX='/expanse/projects/sebat1/a1sriniv/g2mh/joint_calling/batches/vcf_batch_'
BATCH_SIZE=200

files=($(ls "$VCF_DIR"/*.gvcf.gz))
total_num_files=${#files[@]}

batch_index=1
line=''

for ((i=0; i<total_num_files; i++)); do
	file="${files[i]}"
	line+="-V ${file} "

	# if batch size is reached
	if (( (i + 1) % BATCH_SIZE == 0 || i + 1 == total_num_files )); then
		echo "$line" > "${OUT_PREFIX}_${batch_index}.txt"
		line=""
		((batch_index++))
	fi
done
