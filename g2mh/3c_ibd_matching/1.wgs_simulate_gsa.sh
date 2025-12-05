#!/usr/bin/bash
#SBATCH --job-name=wgs_sim
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3c_ibd_matching/logs/wgs_simulate.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3c_ibd_matching/logs/wgs_simulate.err

GSA_VCF_DIR='/expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/preimpute/check_bim_out/vcfs_final/vcfs_upload'
WGS_VCF_DIR='/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs/annotated_vcfs/combined_chrs/with_utrecht'
BCFTOOLS_DIR="/expanse/projects/sebat1/a1sriniv/tools/bcftools-1.22/bin"

for i in {1..22}; do
	${BCFTOOLS_DIR}/bcftools query -f '%CHROM\t%POS\t%POS\n' ${GSA_VCF_DIR}/chr${i}_GSA_freeze1_plink-updated-chr${i}_mod.vcf.gz > gsa_snps.txt	
	${BCFTOOLS_DIR}/bcftools view -R gsa_snps.txt ${WGS_VCF_DIR}/chr${i}_jointcall_VQSR_combined.vcf.gz -Oz -o wgs_simulated/WGS_chr${i}_subset.vcf.gz
	${BCFTOOLS_DIR}/bcftools index wgs_simulated/WGS_chr${i}_subset.vcf.gz
	break
done
