#!/bin/bash
#SBATCH --job-name=make_bed
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=32GB
#SBATCH -t 1-10:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/make_bed.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/make_bed.err

module load cpu/0.17.3b gcc/10.2.0/npcyll4
module load bcftools/1.12
DIR='/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs/annotated_vcfs/rsid'

date

echo "chr1"
for i in {1..20}; do
	bcftools view -f PASS ${DIR}/chr1_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr1 done"

echo "chr2"
for i in {1..20}; do
        bcftools view -f PASS ${DIR}/chr2_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr2 done"

echo "chr3"
for i in {1..20}; do
        bcftools view -f PASS ${DIR}/chr3_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr3 done"

echo "chr4"
for i in {1..20}; do
        bcftools view -f PASS ${DIR}/chr4_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr4 done"

echo "chr5"
for i in {1..20}; do
        bcftools view -f PASS ${DIR}/chr5_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr5 done"

echo "chr6"
for i in {1..20}; do
        bcftools view -f PASS ${DIR}/chr6_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr6 done"

echo "chr7"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr7_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr7 done"

echo "chr8"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr8_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr8 done"

echo "chr9"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr9_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr9 done"

echo "chr10"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr10_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr10 done"

echo "chr11"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr11_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr11 done"

echo "chr12"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr12_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr12 done"

echo "chr13"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr13_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr13 done"

echo "chr14"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr14_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr14 done"

echo "chr15"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chr15_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr15 done"

echo "chr16"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr16_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr16 done"

echo "chr17"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr17_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr17 done"

echo "chr18"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr18_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr18 done"

echo "chr19"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr19_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr19 done"

echo "chr20"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr20_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr20 done"

echo "chr21"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr21_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr21 done"

echo "chr22"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chr22_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chr22 done"

echo "chrX"
for i in {1..10}; do
        bcftools view -f PASS ${DIR}/chrX_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chrX done"

echo "chrY"
for i in {1..5}; do
        bcftools view -f PASS ${DIR}/chrY_chunk_${i}_jointcall_VQSR.vcf.gz | awk '!/^#/ {gsub(/^chr/, "", $1); print $1 "\t" $3 "\t" 0 "\t" $2 "\t" $4 "\t" $5}' >> ${DIR}/g2mh_freeze1.bim
done
echo "chrY done"

date
