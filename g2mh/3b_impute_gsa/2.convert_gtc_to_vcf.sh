#!/bin/bash
#SBATCH --job-name=GSA_VCF
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/make_VCF.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/make_VCF.err


date

module load gcc/10.2.0/npcyll4
module load samtools/1.13/2cnwok7

GTC_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa"
BCFTOOLS_DIR="/expanse/projects/sebat1/a1sriniv/tools/bcftools-1.22/bin"
PLINK_DIR="/expanse/projects/sebat1/a1sriniv/tools/plink"
OUT_DIR="/expanse/projects/sebat1/g2mh_data/gsa/vcf_format"

echo "start converting gtc files to VCF"
echo ""

${BCFTOOLS_DIR}/bcftools plugin /expanse/projects/sebat1/a1sriniv/tools/gtc2vcf/gtc2vcf.so \
--bpm ${GTC_DIR}/refs/GSA-24v3-0_A2.bpm  \
--csv ${GTC_DIR}/refs/GSA-24v3-0_A2.csv \
--egt ${GTC_DIR}/refs/GSA-24v3-0_A1_ClusterFile.egt \
--gtc ${GTC_DIR}/gtc_out/ \
--fasta-ref /expanse/projects/sebat1/a1sriniv/refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
-Oz -o ${OUT_DIR}/All_GSA_freeze1.vcf.gz

${BCFTOOLS_DIR}/bcftools sort ${OUT_DIR}/All_GSA_freeze1.vcf.gz -Oz -o ${OUT_DIR}/Sorted_All_GSA_freeze1.vcf.gz
${BCFTOOLS_DIR}/bcftools index ${OUT_DIR}/Sorted_All_GSA_freeze1.vcf.gz
echo ""
echo "finished making the vcf file."
echo ""
echo "subsetting each chromosome now:"
echo ""
for i in {1..22}; do
	${BCFTOOLS_DIR}/bcftools view -r chr${i} ${OUT_DIR}/Sorted_All_GSA_freeze1.vcf.gz -Oz -o ${OUT_DIR}/chr${i}_GSA_freeze1.vcf.gz
	${BCFTOOLS_DIR}/bcftools index ${OUT_DIR}/chr${i}_GSA_freeze1.vcf.gz
	${PLINK_DIR}/plink --vcf ${OUT_DIR}/chr${i}_GSA_freeze1.vcf.gz --make-bed --out ${OUT_DIR}/chr${i}_GSA_freeze1_plink
done
echo
echo "subsetting each chromosome done"

date

