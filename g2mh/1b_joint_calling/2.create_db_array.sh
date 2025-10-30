#!/bin/bash
#SBATCH --job-name=g2mh-db-chr%a
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o db_logs/batch6/db_chr%a.log
#SBATCH -e db_logs/batch6/db_chr%a.err
#SBATCH --array=9-25%20

# create the array jobs
case $SLURM_ARRAY_TASK_ID in
23)
        CHROM="chrX"
        ;;
24)
        CHROM="chrY"
        ;;
25)
        CHROM="chrM"
        ;;
*)
        CHROM=chr$SLURM_ARRAY_TASK_ID
        ;;
esac

#CHROM='chr22'

# based on batches, decide what to do with the workspace
# if first batch, create. Otherise, update
if [[ "$1" == "batches/vcf_batch__1.txt" ]]; then 
	argument='--genomicsdb-workspace-path '
else
	argument='--genomicsdb-update-workspace-path '
fi

date
echo "start creating db"

export BATCH=$(cat $1)
#echo $BATCH

#echo $argument
#exit 0

module load singularitypro
GATK_CONTAINER="/expanse/projects/sebat1/j3guevar/CONTAINERS/gatk4:4.6.2.0--py310hdfd78af_0"

singularity exec --bind /expanse/projects/sebat1/ $GATK_CONTAINER gatk GenomicsDBImport \
	-R /expanse/projects/sebat1/a1sriniv/g2mh/cnv_genotyping/wgs/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	${BATCH} \
	$argument /expanse/projects/sebat1/a1sriniv/g2mh/joint_calling/db/${CHROM} \
	--reader-threads 50 \
       	--genomicsdb-shared-posixfs-optimizations true \
	--tmp-dir /expanse/projects/sebat1/a1sriniv/g2mh/joint_calling/tmpdir \
	--intervals ${CHROM} 

echo "done processing ${CHROM} db file"
echo
date
echo
