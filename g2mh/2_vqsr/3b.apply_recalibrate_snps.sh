#!/bin/bash
#SBATCH --job-name=VQSR-g2mh
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/apply_VQSR_snps.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/apply_VQSR_snps.err


module load singularitypro
GATK_CONTAINER="/expanse/projects/sebat1/j3guevar/CONTAINERS/gatk4:4.6.2.0--py310hdfd78af_0"
VCF_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/vcfs"
REF_GENOME_PATH="/expanse/projects/sebat1/a1sriniv/g2mh/1a_cnv_genotyping/wgs/resources/"
RESOURCE_PATH="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/resources"
VQSR_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/vcfs"
OUT_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs"

# recalibrate snps
singularity exec --bind /expanse/projects/sebat1/ $GATK_CONTAINER gatk ApplyVQSR \
	-V ${VCF_DIR}/All_sites_only_g2mh_sorted.vcf.gz \
	--recal-file ${VQSR_DIR}/All_sites_g2mh_snps.recal \
	--tranches-file ${VQSR_DIR}/All_sites_g2mh_snps.tranches \
	--truth-sensitivity-filter-level 99.8 \
	--create-output-variant-index true \
	-mode SNP \
	-R ${REF_GENOME_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	-O ${OUT_DIR}/recalibrated_sites_snps.vcf.gz

echo "done"


