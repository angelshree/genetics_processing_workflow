#!/bin/bash
#SBATCH --job-name=g2mhVEP1
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o logs_bcftools/siteonly_out
#SBATCH -e logs_bcftools/siteonly_err

date

module load cpu/0.17.3b gcc/10.2.0/npcyll4 bcftools/1.12

VCF_DIR='/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs/annotated_vcfs'
#BED_FILE='nimh_fu_genes_vep_query.bed'
OUT_DIR='vcf_filtered'

for vcf in "$VCF_DIR"/*VQSR.vcf.gz; do
	# extract filename
	fname=$(basename "$vcf")
	base="${fname%.vcf.gz}"

	echo "Processing $fname...."
	echo "base: $base"

	# create output filename
	outvcf="${OUT_DIR}/${base}_siteonly.vcf.gz"

	echo "filtering ${fname}...."

	echo "output name: $outvcf"

	# use bcftools to remove genotypes
	bcftools view \
		-G \
		-O z \
		-o "$outvcf" \
		"$vcf" 
	
	#echo "vcf is filtered"

	# index the output vcf
	bcftools index "$outvcf"


	#break

done

echo "filtering completed\n"
date
