#!/usr/bin/bash
##SBATCH --job-name=remove_samples
##SBATCH --account=ddp195
##SBATCH --partition=ind-shared
##SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=40
##SBATCH --mem=32GB
##SBATCH -t 1-10:00 #Runtime in D-HH:MM
##SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/sample_remove.log
##SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/logs/sample_remove.err

module load cpu/0.17.3b gcc/10.2.0/npcyll4
module load bcftools/1.12

DIR='/expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs/annotated_vcfs/rsid'

# iterate through each vcf to remove samples that need to be removed
for VCF in ${DIR}/*.vcf.gz; do
	echo $VCF
done




