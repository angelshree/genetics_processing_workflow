#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -J plink_py
#SBATCH -o logs_py/py_%a.out
#SBATCH -e logs_py/py_%a.err
#SBATCH --array=1-22

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?Array index required}
VCF="/expanse/projects/sebat1/g2mh_data/gsa/vcf_format/imputed_vcfs/chr${CHR}.dose.vcf.gz"

python PLINK_PIPELINE_01.py \
  --chr "${CHR}" \
  --vcf "${VCF}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory 16000 \
  --maf 0.01
