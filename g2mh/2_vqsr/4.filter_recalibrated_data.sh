#!bin/bash
module load cpu/0.17.3b  gcc/10.2.0/npcyll4
module load bcftools/1.12

INPUT_PATH="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs"
echo "start filtering snv file"
bcftools view -e 'FILTER="."' ${INPUT_PATH}/recalibrated_sites_snps.vcf.gz -Oz -o ${INPUT_PATH}/recalibrated_sites_snv_only.vcf.gz
echo "filtering snv file done"
echo ""
echo "start indexing snv file"
bcftools index ${INPUT_PATH}/recalibrated_sites_snv_only.vcf.gz
echo "indexing snv file done"
echo ""

echo "start filtering indel file"
bcftools view -e 'FILTER="."' ${INPUT_PATH}/recalibrated_sites_indel.vcf.gz -Oz -o ${INPUT_PATH}/recalibrated_sites_indel_only.vcf.gz
echo "filtering indel file done"
echo ""
echo "start indexing indel file"
bcftools index ${INPUT_PATH}/recalibrated_sites_indel_only.vcf.gz
echo "indexing snv indel done"
echo ""
