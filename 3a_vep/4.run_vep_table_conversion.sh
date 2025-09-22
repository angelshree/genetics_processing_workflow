#!/usr/bin/bash
#SBATCH --job-name=g2mh_table_conv
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o logs_table/out_%a
#SBATCH -e logs_table/err_%a

date

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/bin:$PATH"
eval "$($HOME/bin/micromamba shell hook -s bash)"
micromamba activate python

echo "mamba shell activated"
which python

VCF_DIR='vcf_bgzip'
VEP_DIR='vep_annot_out'

echo "vcf directory: ${VCF_DIR}"
echo "vep directory: ${VEP_DIR}"

for chr in chr{1..22} chrX chrY; do
    # Handle chunked files
    for vcf in ${VCF_DIR}/${chr}_chunk_*_jointcall_VQSR.bgz.vcf.gz; do
        # Skip if no matching files
	#echo "$vcf"
        [[ -e "$vcf" ]] || continue

	base=$(basename "$vcf" .bgz.vcf.gz)

        # Derive the matching VEP file
	vep="${VEP_DIR}/${base}.vep.vcf.gz"

        if [[ -f "$vep" ]]; then
            echo "Running script on: $vcf and $vep, both chunked"
            python convert_vep_results_to_table.py \
		   --vcf_genot ${vcf} \
		   --vcf_annot ${vep} \
		   --psam /expanse/projects/sebat1/j3guevar/autism_exposure_proj/vep/tmp.txt \
		   --regions_file nimh_fu_genes_vep_query.bed \
		   --out table_out/${base}_table.txt 
        else
            echo "⚠️ No matching VEP file for $vcf"
        fi

    done

    # Handle unchunked file
    vcf="${VCF_DIR}/${chr}_jointcall_VQSR.bgz.vcf.gz"
    vep="${VEP_DIR}/${chr}_jointcall_VQSR.vep.vcf.gz"
    base=$(basename "$vcf" .bgz.vcf.gz)
    if [[ -f "$vcf" && -f "$vep" ]]; then
        echo "Running script on: $vcf and $vep, unchunked"
        python convert_vep_results_to_table.py \
		--vcf_genot ${vcf} \
                --vcf_annot ${vep} \
                --psam /expanse/projects/sebat1/j3guevar/autism_exposure_proj/vep/tmp.txt \
                --regions_file nimh_fu_genes_vep_query.bed \
                --out table_out/${base}_table.txt
    elif [[ -f "$vcf" || -f "$vep" ]]; then
        echo "⚠️ One file found but not both for $chr"
    fi

done
