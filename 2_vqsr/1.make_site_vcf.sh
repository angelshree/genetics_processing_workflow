#!/bin/bash
#SBATCH --job-name=VQSR0
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/logs/make_siteVCF.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/logs/make_siteVCF.err


date
echo ""
echo "start script"
echo ""

# load gatk
module load singularitypro
GATK_CONTAINER="/expanse/projects/sebat1/j3guevar/CONTAINERS/gatk4:4.6.2.0--py310hdfd78af_0"
BCFTOOLS_DIR="/expanse/projects/sebat1/a1sriniv/tools/bcftools-1.22/bin"
VCF_DIR="/expanse/projects/sebat1/g2mh_data/jointcalled_vcfs/"
OUT_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/vqsr/vcfs"

# make a site only vcf to work with
for VCF in ${VCF_DIR}/*vcf.gz; do
	filename=$(basename "${VCF}")
	NAME="${filename%.vcf.gz}"
	singularity exec --bind /expanse/projects/sebat1/ $GATK_CONTAINER gatk MakeSitesOnlyVcf \
	-I ${VCF} \
	-O ${OUT_DIR}/sites_only_${NAME}.vcf.gz
done

find ${OUT_DIR}  -type f -name 'sites_only*' ! -name '*tbi'  > ${OUT_DIR}/input_list
${BCFTOOLS_DIR}/bcftools merge -l ${OUT_DIR}/input_list -Oz -o ${OUT_DIR}/All_sites_only_g2mh.vcf.gz
${BCFTOOLS_DIR}/bcftools sort ${OUT_DIR}/All_sites_only_g2mh.vcf.gz -o ${OUT_DIR}/All_sites_only_g2mh_sorted.vcf.gz

singularity exec --bind /expanse/projects/sebat1/ $GATK_CONTAINER gatk IndexFeatureFile -I ${OUT_DIR}/All_sites_only_g2mh_sorted.vcf.gz

echo ""
echo "script finished"
echo ""

