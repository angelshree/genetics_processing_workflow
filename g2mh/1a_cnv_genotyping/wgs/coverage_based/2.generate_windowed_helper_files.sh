#!/usr/bin/bash

# this script takes information about gc content, segdups, simple repeats, gaps, and masked regions and formats them into 10kb bed files that are consistent with the output of mosdepth for downstream analysis of CNV genotypes

module load cpu/0.15.4 gcc/10.2.0 bedtools2/2.27.1

# first, make 10kb windows of full genome
zcat mosdepth_out/884-NDARTB560PVR_10000.regions.bed.gz | awk '{print $1 "\t" $2 "\t" $3}' > resources/window_10kb_hg38.bed

# gc content 
bedtools nuc -fi resources/GRCh38_full_analysis_set_plus_decoy_hla.fa -bed resources/window_10kb_hg38.bed > window_10kb_hg38.gc.bed

# segdup coverage 
bedtools coverage -a resources/window_10kb_hg38.bed -b resources/genomicSuperDups.sorted.bed > resources/segdups_coverage_10kb.hg38.bed

# simple repeat coverage
bedtools coverage -a resources/window_10kb_hg38.bed -b resources/simpleRepeat.bed > resources/simple_repeat_coverage_10kb.hg38.bed

# gap coverage
bedtools coverage -a resources/window_10kb_hg38.bed -b resources/gap.bed > resources/gap_coverage_10kb.hg38.bed

# rmsk (masked) coverage
bedtools coverage -a resources/window_10kb_hg38.bed -b resources/rmsk.bed > resources/rmsk_10kb.hg38.bed
