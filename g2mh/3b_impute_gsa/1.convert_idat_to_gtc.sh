#!/bin/bash
#SBATCH --job-name=GSA_gtc
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH -t 1-00:00 #Runtime in D-HH:MM
#SBATCH -o /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/make_gtc.log
#SBATCH -e /expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/logs/make_gtc.err


date

echo "start"

/expanse/projects/sebat1/a1sriniv/tools/iaap-cli/iaap-cli gencall \
	/expanse/projects/sebat1/mahangari/g2mh/freeze1/gsa/refs/GSA-24v3-0_A2.bpm \
	/expanse/projects/sebat1/mahangari/g2mh/freeze1/gsa/refs/GSA-24v3-0_A1_ClusterFile.egt \
	--idat-folder /expanse/projects/sebat1/g2mh_data/idat/ \
	--output-gtc \
	/expanse/projects/sebat1/a1sriniv/g2mh/3b_impute_gsa/gtc_out


echo "end"

date

