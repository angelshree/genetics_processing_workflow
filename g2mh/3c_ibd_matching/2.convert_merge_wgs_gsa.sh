#!/usr/bin/bash
#SBATCH --job-name=combine_convert_g2mh
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3c_ibd_matching/logs/combine_convert.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3c_ibd_matching/logs/combine_convert.err

PLINK_DIR="/expanse/projects/sebat1/a1sriniv/tools/plink"
GSA_VCF_DIR='/expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/preimpute/check_bim_out/vcfs_final/vcfs_upload'
BCFTOOLS_DIR="/expanse/projects/sebat1/a1sriniv/tools/bcftools-1.22/bin"

# index GSA (if not done already)
#for i in {1..22}; do
#	${BCFTOOLS_DIR}/bcftools index ${GSA_VCF_DIR}/chr${i}_GSA_freeze1_plink-updated-chr${i}_mod.vcf.gz
#done

# first, concatenate the VCFs so we have one VCF for GSA and one VCF for WGS
echo "concatenating WGS"
#${BCFTOOLS_DIR}/bcftools concat -Oz -o merged/WGS_merged_subset.vcf.gz wgs_simulated/WGS_chr{1..22}_subset.vcf.gz
echo " "
echo "concatenating GSA"
#ls ${GSA_VCF_DIR}/chr*_GSA_freeze1_plink-updated-chr*_mod.vcf.gz | sort -V > merge_list.txt
#${BCFTOOLS_DIR}/bcftools concat -f merge_list.txt -Oz -o merged/GSA_merged.vcf.gz
echo " "
echo "indexing"
# index
#${BCFTOOLS_DIR}/bcftools index merged/WGS_merged_subset.vcf.gz
#${BCFTOOLS_DIR}/bcftools index merged/GSA_merged.vcf.gz

# remove multiallelic sites
echo " "
echo "removing multiallelic sites"
#${BCFTOOLS_DIR}/bcftools view -m2 -M2 -v snps merged/WGS_merged_subset.vcf.gz -Oz -o filtered_vcfs/WGS_merged_biallelic_snps.vcf.gz
#${BCFTOOLS_DIR}/bcftools view -m2 -M2 -v snps merged/GSA_merged.vcf.gz -Oz -o filtered_vcfs/GSA_merged_biallelic_snps.vcf.gz
#${BCFTOOLS_DIR}/bcftools index filtered_vcfs/WGS_merged_biallelic_snps.vcf.gz
#${BCFTOOLS_DIR}/bcftools index filtered_vcfs/GSA_merged_biallelic_snps.vcf.gz

# convert to plink format
echo " "
echo "converting to plink format"
#${PLINK_DIR}/plink --vcf filtered_vcfs/WGS_merged_biallelic_snps.vcf.gz --make-bed --out plink_format/WGS_merged_subset
#${PLINK_DIR}/plink --vcf filtered_vcfs/GSA_merged_biallelic_snps.vcf.gz --make-bed --out plink_format/GSA_merged

# rename missing variant IDs
echo " "
echo "handling missing variant IDs"
#${PLINK_DIR}/plink --bfile plink_format/WGS_merged_subset --set-missing-var-ids @:# --make-bed --out plink_format/WGS_fixed
#${PLINK_DIR}/plink --bfile plink_format/GSA_merged --set-missing-var-ids @:# --make-bed --out plink_format/GSA_fixed

# remove locations that would become multiallelic after merging WGS and GSA
echo " "
echo "remove loci that will become multiallelic after merging"
#${PLINK_DIR}/plink --bfile plink_format/WGS_fixed --bmerge plink_format/GSA_fixed --merge-mode 6 --out merge_check  --allow-extra-chr --allow-no-sex
#awk '{print $2}' merge_check.missnp > to_exclude.txt
#${PLINK_DIR}/plink --bfile plink_format/WGS_fixed --exclude merge_check.missnp --make-bed --out plink_format/WGS_clean
#${PLINK_DIR}/plink --bfile plink_format/GSA_fixed --exclude merge_check.missnp --make-bed --out plink_format/GSA_clean

# handle potential allele mismatches
#${PLINK_DIR}/plink --bfile plink_format/WGS_biallelic --bmerge plink_format/GSA_biallelic --flip-scan
#${PLINK_DIR}/plink --bfile plink_format/WGS_biallelic --flip plink.flipscan --make-bed --out plink_format/WGS_flipped

# final merge
echo " "
echo "final combining"
#${PLINK_DIR}/plink --bfile plink_format/WGS_clean --bmerge plink_format/GSA_clean --make-bed --out plink_format/WGS_GSA_combined

echo " "
echo "LD pruning of combined data"
${PLINK_DIR}/plink --bfile plink_format/WGS_GSA_combined --indep-pairwise 50 5 0.2 --out king_prune
${PLINK_DIR}/plink --bfile plink_format/WGS_GSA_combined --extract king_prune.prune.in --make-bed --out plink_format/WGS_GSA_combined_pruned
