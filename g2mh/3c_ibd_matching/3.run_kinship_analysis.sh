#!/usr/bin/bash
#SBATCH --job-name=run_ibd
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH -t 2-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3c_ibd_matching/logs/run_ibd.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3c_ibd_matching/logs/run_ibd.err

KING_DIR="/expanse/projects/sebat1/a1sriniv/tools/king"

${KING_DIR} -b plink_format/WGS_GSA_combined_pruned.bed --kinship --prefix king_out/king_ibd 
