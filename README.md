## genetics_processing_workflow

# SNV/INDEL Pipeline Overview

## 1. Input Data
- Variant calling begins with GVCFs generated using **DRAGEN**.

## 2. GenomicsDB Construction
- Samples are divided into batches, and a **GenomicsDB** is created for each batch.  
- Ultimately, we maintain one GenomicsDB per chromosome.

## 3. Joint Genotyping
- **GATK Joint Genotyping** is performed using the chromosome-specific GenomicsDB as input.  
- For large chromosomes (and large sample sizes), we process the data in chunks.  
- Chunk sizes vary depending on chromosome length and are determined using `get_chunks.py`.  
- Approximately 24 sbatch scripts are used, each launching arrays. Each chunk typically requires ~1000 minutes of runtime.

## 4. Variant Quality Score Recalibration (VQSR)
- **SNPs** and **indels** are recalibrated separately, producing both a recalibration file and a tranches file.  
- Recalibration is applied to SNPs and indels, with site-only VCFs generated using `make_site_vcf.sh` (which strips genotypes while retaining variant sites).  
- Filtering is applied to recalibrated variants, retaining only those with a `.` in the FILTER column (though an alternative approach would be to retain only `PASS` variants).

## 5. Post-processing
- VQSR annotations are incorporated into the VCFs using **bcftools annotate**.  
- rsIDs are added at this stage as well.
