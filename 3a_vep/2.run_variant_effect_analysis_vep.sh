#!/bin/bash
#SBATCH --job-name=g2mhVEP
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o logs_vep/out_%a
#SBATCH -e logs_vep/err_%a
#SBATCH --array=1-20%20 #--array=1-260%30

date

module load singularitypro

OUT_PATH="vep_annot_out"
VEP_CONTAINER="/expanse/projects/sebat1/mahangari/containers/ensembl-vep:110.1--pl5321h2a3209d_0"
GRCH38_REF="/expanse/projects/sebat1/a1sriniv/g2mh/1a_cnv_genotyping/wgs/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CACHE="/expanse/projects/sebat1/vep_cache/version_110"
PLUGIN="/expanse/projects/sebat1/vep_plugins/V2/VEP_plugins"

VCF=$(realpath /expanse/projects/sebat1/a1sriniv/g2mh/2_vqsr/filtered_vcfs/annotated_vcfs/*VQSR.vcf.gz | grep -v "All" | sort -V | head -n $SLURM_ARRAY_TASK_ID - | tail -n 1)
VCF_NAME=$(echo ${VCF} | rev | cut -d / -f 1 | rev | cut -d . -f 1)
CHROM=$(echo ${VCF_NAME} | cut -d _ -f 1)

echo $VCF
echo $VCF_NAME
echo $CHROM

singularity exec --bind /expanse/projects/sebat1/ ${VEP_CONTAINER} \
	vep \
	--input_file $VCF \
	--format vcf \
	--output_file ${OUT_PATH}/${VCF_NAME}.vep.vcf.gz \
	--vcf \
	--compress_output bgzip \
	--minimal \
    --canonical \
    --symbol \
	--assembly GRCh38 \
	--cache \
	--dir_cache ${CACHE} \
	--offline \
	--fasta ${GRCH38_REF} \
	--allele_number \
	--pick_allele \
	--force_overwrite \
	--fork 4 \
	--stats_text \
	--dir_plugins ${PLUGIN} \
	--plugin pLI \
	--plugin dbNSFP,/expanse/projects/sebat1/vep_resources/dbNSFP4.4a/dbNSFP4.4a_grch38.gz,transcript_match=1,MPC_score \
	--plugin LOEUF,file=/expanse/projects/sebat1/vep_resources/LOEUF/loeuf_dataset_grch38.tsv.gz,match_by=gene \
	--plugin AlphaMissense,file=/expanse/projects/sebat1/vep_resources/AlphaMissense/AlphaMissense_hg38.tsv.gz \
	--custom file=/expanse/projects/sebat1/mahangari/gnomad/gnomad.exomes.v4.0.sites.${CHROM}.vcf.bgz,short_name=gnomADe_V4,format=vcf,type=exact,fields=AF \
	--custom file=/expanse/projects/sebat1/mahangari/gnomad/gnomad.exomes.v4.0.sites.${CHROM}.vcf.bgz,short_name=gnomADe_V4_nfe,format=vcf,type=exact,fields=AF_nfe

date
