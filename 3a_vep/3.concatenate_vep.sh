#!/bin/bash
#SBATCH --job-name=g2mhVEP
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o logs_bcftools/out_%a
#SBATCH -e logs_bcftools/err_%a

date
module load cpu/0.17.3b gcc/10.2.0/npcyll4 bcftools/1.12

# combining chunked vcfs for final processing
for chr in chr{1..22} chrX chrY; do
	echo "merging $chr chunks"

	# first, combine original genome together
	files=$(ls vcf_bgzip/${chr}_chunk_*.vcf.gz 2>/dev/null)
	echo "files: $files"
	if [[ -n "$files" ]]; then
		bcftools concat -a -O z -o vcf_bgzip/${chr}_merged.vcf.gz $files
		bcftools index vcf_bgzip/${chr}_merged.vcf.gz
	else
		# no chunks
		cp vcf_bgzip/${chr}_jointcall_VQSR.filtered.vcf.gz vcf_bgzip/${chr}_merged.vcf.gz
		bcftools index vcf_bgzip/${chr}_merged.vcf.gz
	fi

	files=$(ls vep_annot_out/${chr}_chunk_*.vcf.gz 2>/dev/null)
	echo "files: $files"
        if [[ -n "$files" ]]; then
                bcftools concat -a -O z -o vep_annot_out/${chr}_merged.vcf.gz $files
                bcftools index vep_annot_out/${chr}_merged.vcf.gz
	else
		# no chunks
		cp vep_annot_out/${chr}_jointcall_VQSR.vep.vcf.gz vep_annot_out/${chr}_merged.vcf.gz
		bcftools index vep_annot_out/${chr}_merged.vcf.gz
        fi
done

# now, combine the chromosomes together
echo 'combining chromosmes together'

# genotype vcf
bcftools concat -a -O z -o vcf_bgzip/allchr.vcf.gz vcf_bgzip/chr*_merged.vcf.gz
bcftools index vcf_bgzip/allchr.vcf.gz

# vep vcf
bcftools concat -a -O z -o vep_annot_out/allchr.vep.vcf.gz vep_annot_out/chr*_merged.vcf.gz
bcftools index vep_annot_out/allchr.vep.vcf.gz

echo 'done! :)'
date
