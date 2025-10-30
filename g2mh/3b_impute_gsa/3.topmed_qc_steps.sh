#!/usr/bin/bash
#SBATCH --job-name=GSA_VCF
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/vcf_qc.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/vcf_qc.err

PLINK_DIR="/expanse/projects/sebat1/a1sriniv/tools/plink/plink"
BIM_DIR="/expanse/projects/sebat1/g2mh_data/gsa/vcf_format"
FREQ_DIR="/expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/preimpute/freq"

for i in {1..22}; do
	${PLINK_DIR} --freq --bfile ${BIM_DIR}/chr${i}_GSA_freeze1_plink --out ${FREQ_DIR}/chr${i}_GSA_freeze1_freq_plink
	perl /expanse/projects/sebat1/a1sriniv/tools/topmed_impute/HRC-1000G-check-bim.pl -b ${BIM_DIR}/chr${i}_GSA_freeze1_plink.bim -f ${FREQ_DIR}/chr${i}_GSA_freeze1_freq_plink.frq -r /expanse/projects/sebat1/a1sriniv/tools/topmed_impute/PASS.Variantsbravo-dbsnp-all.tab -h -o preimpute/check_bim_out/chr${i}
	sed -i "s|^[[:space:]]*plink|${PLINK_DIR}|" preimpute/check_bim_out/chr${i}/Run-plink.sh
        sh preimpute/check_bim_out/chr${i}/Run-plink.sh	
	rm /expanse/projects/sebat1/g2mh_data/gsa/vcf_format/TEMP*
done


