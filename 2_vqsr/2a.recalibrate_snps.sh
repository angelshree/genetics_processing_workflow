#!/bin/bash
#SBATCH --job-name=VQSR-g2mh
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/logs/VQSR_snv.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/logs/VQSR_snv.err


# start script
echo "start VQSR for snvs"

module load singularitypro
GATK_CONTAINER="/expanse/projects/sebat1/j3guevar/CONTAINERS/gatk4:4.6.2.0--py310hdfd78af_0"

singularity exec --bind /expanse/projects/sebat1 $GATK_CONTAINER gatk VariantRecalibrator \
	-V /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/vcfs/All_sites_only_g2mh_sorted.vcf.gz \
	-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
	--trust-all-polymorphic \
	-an ReadPosRankSum -an MQRankSum -an FS -an QD -an SOR -an DP \
	-mode SNP \
	-R /expanse/projects/sebat1/a1sriniv/g2mh/cnv_genotyping/wgs/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	--resource:hapmap,known=false,training=true,truth=true,prior=15 /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/resources/hapmap_3.3.hg38.vcf.gz \
	--resource:omni,known=false,training=true,truth=true,prior=12 /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/resources/1000G_omni2.5.hg38.vcf.gz \
	--resource:1000G,known=false,training=true,truth=false,prior=10 /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2 /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf \
	-O /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/vcfs/All_sites_g2mh_snps.recal \
	--tranches-file /expanse/projects/sebat1/a1sriniv/g2mh/vqsr/vcfs/All_sites_g2mh_snps.tranches \
	--rscript-file All_sites_g2mh_snps.plots.R \
	--dont-run-rscript
echo "done"
