#!/bin/bash
#SBATCH --job-name=rsid_annotate2
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/annotate_rsid_chrY.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/annotate_rsid_chrY.err

## create the array jobs
CHROM='Y'

echo
date
echo
echo "start annotating rsids for chromosome ${CHROM}"

module load cpu/0.17.3b gcc/10.2.0/npcyll4
module load bcftools/1.12

TOPMED_VCF='/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/bravo-dbsnp-all.rsid.vcf.gz'
DIR="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs/annotated_vcfs"

for i in {1..5}; do
	bcftools annotate -a ${TOPMED_VCF} -c ID \
	-Oz -o ${DIR}/rsid/chr${CHROM}_chunk_${i}_jointcall_VQSR.vcf.gz ${DIR}/chr${CHROM}_chunk_${i}_jointcall_VQSR.vcf.gz
	bcftools index ${DIR}/rsid/chr${CHROM}_chunk_${i}_jointcall_VQSR.vcf.gz
	echo "chunk ${i} done"
done

echo "annotating rsids for chromosome ${CHROM} is done"
echo
date
echo
