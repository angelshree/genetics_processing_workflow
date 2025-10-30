###################################################################################################
# sex checks on GSA
###################################################################################################

###################################################################################################
# imports
library(conflicted)
library(data.table) # For high-performance data manipulation
library(dplyr)
###################################################################################################

###################################################################################################
manifest <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/demographics_distributions/full_sample_manifest.csv')

# filter for wgs samples and remove unwanted cols
manifest <- manifest %>%
  filter(WGS.GSA == 'GSA') %>%
  select(GUID, sex_final)
###################################################################################################

###################################################################################################
gsa_reports <- c(list.files("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/adjusted_LRR", full.names = TRUE))
#gsa_reports <- head(gsa_reports, 5)
# processing portion
results <- lapply(gsa_reports, function(f) {
  dt <- fread(f, na.strings=c("NA","--"))
  
  if(length(names(dt)) == 11) {
    print(f)
  }
  
  # standardize colnames
  setnames(dt, gsub("[ .-]+","_",names(dt)))
  
  # make sure LRR is numeric
  dt[, Log_R_Ratio := as.numeric(Log_R_Ratio)]
  dt$Chr <- as.character(dt$Chr)
  
  # compute averages for X and Y
  mean_X <- dt[Chr %in% c("X","chrX"), mean(Log_R_Ratio, na.rm=TRUE)]
  mean_Y <- dt[Chr %in% c("Y","chrY"), mean(Log_R_Ratio, na.rm=TRUE)]
  mean_autosome <- dt[Chr %in% as.character(1:22), mean(Log_R_Ratio, na.rm=TRUE)]
  
  data.table(
    Sample_ID = gsub("\\.report$","",basename(f)),
    X_LRR = mean_X,
    Y_LRR = mean_Y,
    autosome_LRR = mean_autosome
  )
})
avgXY_wide <- rbindlist(results, fill=TRUE)
###################################################################################################

scatter.smooth(avgXY_wide$X_LRR, avgXY_wide$Y_LRR)

###################################################################################################
avgXY_wide <- avgXY_wide %>% filter(!Sample_ID == '887-NDARFM133AW6.report.adjusted.bak')
avgXY_wide$GUID = substr(avgXY_wide$Sample_ID, 5, nchar(avgXY_wide$Sample_ID)-16)
manifest_gsa_combined <- merge(manifest, avgXY_wide, by='GUID')
manifest_gsa_combined$X_LRR_Ratio = (manifest_gsa_combined$X_LRR - manifest_gsa_combined$autosome_LRR)
manifest_gsa_combined$Y_LRR_Ratio = (manifest_gsa_combined$Y_LRR - manifest_gsa_combined$autosome_LRR)
###################################################################################################

ggplot(manifest_gsa_combined, aes(x = X_LRR, y = Y_LRR, color = sex_final)) +
  geom_point() +
  theme_classic() +
  labs(title = 'avg x-LRR vs. avg y-LRR for G2MH GSA samples')

ggplot(manifest_gsa_combined, aes(x = X_LRR_Ratio, y = Y_LRR_Ratio, color = sex_final)) +
  geom_point() +
  theme_classic() +
  labs(title = 'avg x-LRR vs. avg y-LRR ratios for G2MH GSA samples')








