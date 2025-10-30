#!/bin/bash
#SBATCH --job-name=VQSR5a
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -t 1-10:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/annotate_jointcall_VQSR_chr%a.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/annotate_jointcall_VQSR_chr%a.err
#SBATCH --array=1-6

date

CHR=$SLURM_ARRAY_TASK_ID
module load cpu/0.17.3b gcc/10.2.0/npcyll4
module load bcftools/1.12

INPUT_DIR="/expanse/projects/sebat1/g2mh_data/jointcalled_vcfs"
OUT_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs"

for i in {1..20}; do # this is the number of chunks
	bcftools annotate -a ${OUT_DIR}/recalibrated_sites_snv_only.vcf.gz -c FILTER -Oz -o ${OUT_DIR}/annotated_vcfs/chr${CHR}_chunk_${i}_jointcall_VQSR_SNPs.vcf.gz ${INPUT_DIR}/chr${CHR}_chunk_${i}_jointcall.vcf.gz
	bcftools index ${OUT_DIR}/annotated_vcfs/chr${CHR}_chunk_${i}_jointcall_VQSR_SNPs.vcf.gz
done

for i in {1..20}; do # this is the number of chunks
	bcftools annotate -a ${OUT_DIR}/recalibrated_sites_indel_only.vcf.gz -c FILTER -Oz -o ${OUT_DIR}/annotated_vcfs/chr${CHR}_chunk_${i}_jointcall_VQSR.vcf.gz ${OUT_DIR}/annotated_vcfs/chr${CHR}_chunk_${i}_jointcall_VQSR_SNPs.vcf.gz
        bcftools index ${OUT_DIR}/annotated_vcfs/chr${CHR}_chunk_${i}_jointcall_VQSR.vcf.gz
done

echo "done"
