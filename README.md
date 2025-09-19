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

# Structural Variant (SV) Pipeline Overview

## 1. Recurrent CNVs
- Start with **DRAGEN CNV** call files.  
- **Filter CNV calls** and extract recurrent loci from VCFs (using partial overlap).  
- Format recurrent loci into a **BED file**.  
- Repeat the process with QC filtering (e.g., based on the `FILTER` flag, retaining only `PASS` calls).

### Molly's Script
- Core functionality: Python + shell scripts, essentially a **bedtools overlap** (reciprocal and one-way overlap).  

### QC Checks
- Run **mosdepth** for coverage.  
- Run QC scripts (James' coverage QC scripts).  
- Use the **coverage caller script** for visualization.  
- Assess **CNV prevalence** across recurrent CNVs.

## 2. Generic SVs
- The pipeline is still under development.  
- Regardless of approach, **DRAGEN SV calls** will serve as input to a genotyper.  
- Candidate genotypers:  
  - **sv2** (+ **SURVIVOR**)  
  - **Graphtyper** ([DecodeGenetics/graphtyper](https://github.com/DecodeGenetics/graphtyper))  
  - **GATK-SV** ([broadinstitute/gatk-sv](https://github.com/broadinstitute/gatk-sv)) — uncertain if joint-genotyping is supported.
 
# GSA Genotype Microarray Pipeline Overview

## 1. Input Data
- Begin with **report** and **intensity files** (Illumina platform).  
- Confirm whether all report files are stored in the same folder.

## 2. CNV Workflow
- Current workflow (Omar):  
  - Use **PennCNV** and **iPattern**, followed by **harmonization** and **QC**.  
- Obtain CNV calls and process with **Molly's script** (to identify recurrent CNVs).  
- Retain both **recurrent CNVs** and **other CNVs** for downstream analyses.

## 3. QC Steps
- Assess **prevalence** of CNVs.  
- Use the **coverage caller** for QC evaluation.  

## 4. Imputation
- Discussion with Mamad and Molly regarding imputation strategy.  
- Plan to use **Ricopolli** (to be detailed later).  

## 5. Integration with WGS
- Combine **WGS + GSA** data for downstream analyses:  
  - **Polygenic Scores (PGS)**  
  - **IBD matching**  
  - **Ancestry SNPs and PCA** (scripts available from Molly)  
