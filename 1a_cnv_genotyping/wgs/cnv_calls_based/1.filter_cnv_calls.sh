#!/usr/bin/bash
##SBATCH --job-name=g2mh_cnv_filtering
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH -t 1-00:00 #Runtime in D-HH:MM
#SBATCH -o bed_conversion.out
#SBATCH -e bed_conversion.err

date

# load bcftools module from slurm
module load cpu/0.17.3b  gcc/10.2.0/npcyll4 bcftools/1.12

cnv_vcf_dir='/expanse/projects/sebat1/g2mh_data/dragen/cnv'
recurrent_cnv_bed='/expanse/projects/sebat1/a1sriniv/g2mh/cnv_genotyping/wgs/other_recurrent_loci/recurrent_cnv_loci.bed'
outdir='cnv_vcfs_filt'
outdir_bed='cnv_vcf_beds'

for vcf in ${cnv_vcf_dir}/*.cnv.vcf.gz; do
	base=$(basename "$vcf" .cnv.vcf.gz)
	echo ${base}
	bcftools view -R ${recurrent_cnv_bed} "$vcf" -Oz -o "${outdir}/${base}.filtered.vcf.gz"
	bcftools index "${outdir}/${base}.filtered.vcf.gz"
	bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ALT\t%INFO[\t%GT:%SM:%CN:%BC:%PE]\n' "${outdir}/${base}.filtered.vcf.gz" > "${outdir_bed}/${base}.filtered.bed"
done

date
