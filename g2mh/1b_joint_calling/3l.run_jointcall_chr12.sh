#!/bin/bash
#SBATCH --job-name=g2mh-jointcall
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH --cpus-per-task=32
#SBATCH --mem=80GB
#SBATCH -o jointcall_logs/jointcall_all_chr12_chunk_%a.log
#SBATCH -e jointcall_logs/jointcall_all_chr12_chunk_%a.err
#SBATCH --array=1-10%5

CHUNK=$SLURM_ARRAY_TASK_ID
CHROM='chr12'
REGION=$(sed -n "${CHUNK}p" /expanse/projects/sebat1/a1sriniv/g2mh/joint_calling/chunks/"${CHROM}"_chunks.txt)

#echo ${CHUNK}
#echo ${CHROM}
echo ${REGION}
#exit 0

echo "start script to jointly call genotypes for ${REGION} of chromosome ${CHROM}"
echo
date
echo

module load singularitypro
GATK_CONTAINER="/expanse/projects/sebat1/j3guevar/CONTAINERS/gatk4:4.6.2.0--py310hdfd78af_0"

singularity exec --bind /expanse/projects/sebat1/ $GATK_CONTAINER gatk GenotypeGVCFs \
	-R /expanse/projects/sebat1/a1sriniv/g2mh/cnv_genotyping/wgs/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	-V gendb:///expanse/projects/sebat1/a1sriniv/g2mh/joint_calling/db/${CHROM} \
	--tmp-dir /expanse/projects/sebat1/a1sriniv/g2mh/joint_calling/tmpdir \
	-L ${REGION} \
	-O /expanse/projects/sebat1/g2mh_data/jointcalled_vcfs/${CHROM}_chunk_${CHUNK}_jointcall.vcf.gz
echo
echo "script to jointly genotype ${CHROM} done."
date
echo
