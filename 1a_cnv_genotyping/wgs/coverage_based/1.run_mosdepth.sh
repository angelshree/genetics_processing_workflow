#!/usr/bin/bash
##SBATCH --job-name=g2mh_wgs_1
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH -t 0-02:00 #Runtime in D-HH:MM
#SBATCH -o logs_mosdepth/out_%x_%a.out
#SBATCH -e logs_mosdepth/err_%x_%a.err
#SBATCH --array=0-800%30

# note there are more than 1000 alignment files in G2MH, so need to split this up
# ls /expanse/projects/sebat1/g2mh_data/alignments/*.cram > all_alignments.txt
# split -l 800 all_alignments.txt split_alignments_
# then, run sbatch --job-name <job_batch_identifier> 1.run_mosdepth.sh <batch_file_>

date

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$HOME/bin:$PATH"
eval "$($HOME/bin/micromamba shell hook -s bash)"
micromamba activate mosdepth 

#echo "mamba shell activated"

# note: mosdepth is a conda install

BATCH_FILE=$1

# Skip if the index is too large for this batch
NUM_LINES=$(wc -l < "$BATCH_FILE")
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NUM_LINES" ]; then
    echo "Index $SLURM_ARRAY_TASK_ID exceeds lines in $BATCH_FILE, skipping"
    exit 0
fi

FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$BATCH_FILE")

# filepaths
out='/expanse/projects/sebat1/a1sriniv/g2mh/cnv_genotyping/wgs/mosdepth_out'

window_size=10000
echo "window size: ${window_size}"

# iterating through all files to run mosdepth

echo "sample: ${FILE}"
		
outfile_name=$(echo "$FILE" | cut -d'/' -f7)
echo "$outfile_name"
outfile_name=$(echo "${outfile_name}" | sed 's/.....$//')	
echo "output file: ${out}/${outfile_name}_${window_size}"
#exit 0

mosdepth \
	--threads $(nproc) \
	--by ${window_size} \
	--no-per-base \
	--fast-mode \
	--fasta /expanse/projects/sebat1/j3guevar/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	${out}/${outfile_name}_${window_size} \
	${FILE}
		
echo

date
