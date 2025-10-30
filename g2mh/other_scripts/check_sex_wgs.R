###################################################################################################
# sex checks on WGS
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
  filter(WGS.GSA == 'WGS') %>%
  select(GUID, sex_final)
###################################################################################################

###################################################################################################
# function to check sex based on coverage
check_sex <- function(filepath) {
  cov <- read.table(filepath,
                    header=FALSE,
                    sep = "\t",
                    col.names = c("chr", "start", "end", "coverage"))
  cov$len <- cov$end - cov$start
  chr_cov <- aggregate(cbind(cov$coverage * cov$len, cov$len),
                       by = list(cov$chr), FUN = sum)
  names(chr_cov) <- c("chr", "cov_sum", "len_sum")
  chr_cov$mean_cov <- chr_cov$cov_sum / chr_cov$len_sum
  auto_cov <- mean(chr_cov$mean_cov[chr_cov$chr %in% paste0("chr", 1:22)])
  x_cov    <- chr_cov$mean_cov[chr_cov$chr == "chrX"]
  y_cov    <- chr_cov$mean_cov[chr_cov$chr == "chrY"]
  x_ratio <- x_cov / auto_cov
  y_ratio <- y_cov / auto_cov
  sex <- if (x_ratio > 0.75 & y_ratio < 0.1) {
    "XX"
  } else if (x_ratio < 0.75 & y_ratio > 0.1) {
    "XY"
  } else {
    "Ambiguous"
  }
  return(data.frame(filepath, x_ratio, y_ratio, sex))
}
###################################################################################################

###################################################################################################
wgs_bed_files <- c(list.files("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/mosdepth_out/", pattern = "^[0-9]{3}-.*\\.regions.bed.gz$", full.names = TRUE))

# this takes a bit
sex_all_wgs <- wgs_bed_files %>%
  lapply(check_sex) %>%
  bind_rows()

sex_all_wgs <- sex_all_wgs %>%
  mutate(GUID = str_extract(filepath, "(?<=-)[^_]+"))

manifest_wgs_combined <- merge(manifest, sex_all_wgs, by='GUID')

manifest_wgs_combined <- manifest_wgs_combined %>%
  mutate(wgs_sex = ifelse(sex == 'XX', 'Female', ifelse(sex == 'XY', 'Male', 'Ambiguous'))) %>%
  mutate(sex_discrepancy = (sex_final != wgs_sex))

# discrepant samples: NDARDZ613CM4 (Male XXY), NDARZV565FTG (Female XO)

###################################################################################################

###################################################################################################
# plot
ggplot(manifest_wgs_combined, aes(x = x_ratio, y = y_ratio, color = sex_final)) +
  geom_point() +
  theme_classic() +
  labs(title = 'x-ratio vs. y-ratio for G2MH WGS samples')
