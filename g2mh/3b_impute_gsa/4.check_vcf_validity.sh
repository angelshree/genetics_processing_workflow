#!/bin/bash
#SBATCH --job-name=GSA_gtc
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 1-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/check_vcf_valid.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/check_vcf_valid.err

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/bin:$PATH"
eval "$($HOME/bin/micromamba shell hook -s bash)"
micromamba activate py2

for i in {1..22}; do
	python check_vcf.py -r /expanse/projects/sebat1/a1sriniv/g2mh/1a_cnv_genotyping/wgs/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa -o preimpute/check_vcf_out/chr${i} preimpute/check_bim_out/chr${i}/chr${i}_GSA_freeze1_plink-updated-chr${i}.vcf  
done
