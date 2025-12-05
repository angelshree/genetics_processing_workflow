# Convert first fileset
plink2 \
  --pfile /expanse/projects/sebat1/a1sriniv/g2mh/4a_pgs/1_wgs/05_summary/cohort \
  --make-bed \
  --out wgs \
  --memory 16000 \
  --threads $(nproc)

# Convert second fileset
plink2 \
  --pfile /expanse/projects/sebat1/a1sriniv/g2mh/4a_pgs/2_gsa/05_summary/cohort \
  --make-bed \
  --out imputed \
  --memory 16000 \
  --threads $(nproc)

plink --memory 16000 --threads $(nproc) --bfile imputed --bmerge wgs --make-bed --out merged --geno 0.01
plink2 --bfile merged --make-pgen --out mergedp --threads $(nproc) --memory 16000
