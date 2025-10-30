# checking GSA calls - what is the null distribution and where do our 22q samples lie?

library(dplyr)
library(scales)
library(tidyr)

setwd("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/cnv_calls/")

# import raw calls
batch1 <- read.csv("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/gsa_raw_cnv_calls/G2MH_batch1.formatted_Merge_IP_UP.plt", sep='\t')
batch2 <- read.csv("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/gsa_raw_cnv_calls/G2MH_batch2.formatted_Merge_IP_UP.plt", sep='\t')
batch3 <- read.csv("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/gsa_raw_cnv_calls/G2MH_batch3.formatted_Merge_IP_UP.plt", sep='\t')
batch4 <- read.csv("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/gsa_raw_cnv_calls/G2MH_batch4.formatted_Merge_IP_UP.plt", sep='\t')
gsa_batches_combined <- rbind(batch1, batch2, batch3, batch4)
gsa_batches_combined <- gsa_batches_combined %>%
  filter(algorithm == 'IP,UP') %>%
  mutate(cnv_length = as.integer(consEnd) - as.integer(consStart)) %>%
  filter(consCallType != 'DELDUP') %>%
  mutate(guid = sampleId %>%
           str_sub(5) %>%                      # remove first 4 characters
           str_remove("\\.report$") %>%        # remove '.report' at end
           str_replace_all("_", ""))           # remove '_'

gsa_del_sum <- gsa_batches_combined %>% 
  filter(consCallType == 'DEL') %>%
  group_by(guid) %>%
  summarise(del_count = sum(cnv_length, na.rm = TRUE))
gsa_dup_sum <- gsa_batches_combined %>% 
  filter(consCallType == 'DUP') %>%
  group_by(guid) %>%
  summarise(dup_count = sum(cnv_length, na.rm = TRUE))

# import sample manifest (to get people in GSA that we KNOW have a CNV for sure)
manifest <- read.csv('~/sebatlab Dropbox/G2MH/NetWorkPaper/AnalysisByAnjali/analysis/post_data_freeze/demographics_distributions/full_sample_manifest.csv')
manifest <- manifest %>%
  filter(WGS.GSA == 'GSA') %>%
  select(Genomic_ID, guid, redcap_data_access_group, age_final, sex_final, group_class_breakpoints) %>%
  mutate(guid = guid %>% str_replace_all("_", ""))

manifest_with_cnv_length <- manifest %>%
  left_join(gsa_del_sum, by='guid') %>%
  left_join(gsa_dup_sum, by='guid')

samples_with_putative_cnv <- manifest_with_cnv_length %>% filter(group_class_breakpoints != 'no_cnv')
samples_no_cnv <- manifest_with_cnv_length %>% filter(group_class_breakpoints == 'no_cnv')
sample_sd_stats <- samples_no_cnv %>%
  summarise(
    mean_del = mean(del_count, na.rm = TRUE),
    sd_del = sd(del_count, na.rm = TRUE),
    mean_dup = mean(dup_count, na.rm = TRUE),
    sd_dup = sd(dup_count, na.rm = TRUE))
    
mean_val <- sample_sd_stats$mean_del
sd_val <- sample_sd_stats$sd_del
p1 <- ggplot(samples_no_cnv, aes(x = del_count)) +
      geom_histogram(binwidth = 50000, fill = "black", alpha = 0.5) +
      geom_vline(xintercept = mean_val, color = "blue", linetype = "solid", size = 1) +
      geom_vline(xintercept = mean_val + sd_val, color = "red", linetype = "dashed") +
      geom_vline(xintercept = mean_val + 2*sd_val, color = "red", linetype = "dashed") +
      geom_vline(xintercept = mean_val + 3*sd_val, color = "red", linetype = "dashed") +
      geom_point(data = samples_with_putative_cnv,
                 aes(x = del_count, y = 90), 
                 color = "purple", size = 1) +
      geom_text(data = samples_with_putative_cnv,
                aes(x = del_count, y = 100, label = guid,  angle = 90,  ),
                vjust = -0.5, color = "purple", ) +
      scale_x_continuous(labels = comma) +
      labs(
        title = "Total Deleted Bases per Sample with SD Lines",
        x = "Total deletion length (bp)",
        y = "Count"
      ) +
      theme_minimal()
samples_above_del_sd <- samples_no_cnv %>% filter(del_count >= (mean_val + 3*sd_val))
ggsave("del_lengths_histogram.png", plot = p1, width = 8, height = 6, dpi = 300)

mean_val <- sample_sd_stats$mean_dup
sd_val <- sample_sd_stats$sd_dup
p2 <- ggplot(samples_no_cnv, aes(x = dup_count)) +
  geom_histogram(binwidth = 50000, fill = "black", alpha = 0.5) +
  geom_vline(xintercept = mean_val, color = "blue", linetype = "solid", size = 1) +
  geom_vline(xintercept = mean_val + sd_val, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean_val + 2*sd_val, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean_val + 3*sd_val, color = "red", linetype = "dashed") +
  geom_point(data = samples_with_putative_cnv,
             aes(x = dup_count, y = 65), 
             color = "purple", size = 1) +
  geom_text(data = samples_with_putative_cnv,
            aes(x = dup_count, y = 75, label = guid,  angle = 90,  ),
            vjust = -0.5, color = "purple") +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Total Duplicated Bases per Sample with SD Lines",
    x = "Total duplication length (bp)",
    y = "Count"
  ) +
  theme_minimal()
samples_above_dup_sd <- samples_no_cnv %>% filter(dup_count >= (mean_val + 3*sd_val))
ggsave("dup_lengths_histogram.png", plot = p1, width = 8, height = 6, dpi = 300)



# more constrained del and dup lengths to see the additonal samples easier
mean_val <- sample_sd_stats$mean_del
sd_val <- sample_sd_stats$sd_del
p1 <- ggplot(samples_no_cnv %>% filter(del_count < 7500000), aes(x = del_count)) +
  geom_histogram(binwidth = 50000, fill = "black", alpha = 0.5) +
  geom_vline(xintercept = mean_val, color = "blue", linetype = "solid", size = 1) +
  geom_vline(xintercept = mean_val + sd_val, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean_val + 2*sd_val, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean_val + 3*sd_val, color = "red", linetype = "dashed") +
  geom_point(data = samples_with_putative_cnv,
             aes(x = del_count, y = 90), 
             color = "purple", size = 1) +
  geom_text(data = samples_with_putative_cnv,
            aes(x = del_count, y = 100, label = guid,  angle = 90,  ),
            vjust = -0.5, color = "purple", ) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Total Deleted Bases per Sample with SD Lines and Extra Samples",
    x = "Total deletion length (bp)",
    y = "Count"
  ) +
  theme_minimal()
samples_above_del_sd <- samples_no_cnv %>% filter(del_count >= (mean_val + 3*sd_val))
ggsave("del_lengths_histogram_filt.png", plot = p1, width = 16, height = 6, dpi = 300)

mean_val <- sample_sd_stats$mean_dup
sd_val <- sample_sd_stats$sd_dup
p2 <- ggplot(samples_no_cnv %>% filter(dup_count < 4000000), aes(x = dup_count)) +
  geom_histogram(binwidth = 50000, fill = "black", alpha = 0.5) +
  geom_vline(xintercept = mean_val, color = "blue", linetype = "solid", size = 1) +
  geom_vline(xintercept = mean_val + sd_val, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean_val + 2*sd_val, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean_val + 3*sd_val, color = "red", linetype = "dashed") +
  geom_point(data = samples_with_putative_cnv,
             aes(x = dup_count, y = 75), 
             color = "purple", size = 1) +
  geom_text(data = samples_with_putative_cnv,
            aes(x = dup_count, y = 85, label = guid,  angle = 90,  ),
            vjust = -0.5, color = "purple") +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Total Duplicated Bases per Sample with SD Lines and Extra Samples",
    x = "Total duplication length (bp)",
    y = "Count"
  ) +
  theme_minimal()
samples_above_dup_sd <- samples_no_cnv %>% filter(dup_count >= (mean_val + 3*sd_val))
ggsave("dup_lengths_histogram_filt.png", plot = p2, width = 16, height = 6, dpi = 300)

###################################################################################################
# based on the proband's group class, the theoretical distribution will be shifted
# i.e. if the proband has 22q Dup, the 3SD threshold for filtering based on dup length
# will be shifted up by 2,000,000 bp; if the proband has 16pdel, the 3SD threshold for filtering 
# based on del length will be shifted down by 600,000

datafreeze <- read.csv("~/sebatlab Dropbox/G2MH/DataFreeze2025/Download_08_27_2025/Project2_DataRelease/August2025/g2mh_project2_data_release_082025.csv") 
datafreeze <- datafreeze %>%
  select(rarecnv_id, guid, rarecnv_family_id, proband, proband_rarecnv_id) %>%
  left_join(
    datafreeze %>%
      select(proband_rarecnv_id = rarecnv_id, proband_guid = guid),
    by = "proband_rarecnv_id"
  ) %>%
  mutate(guid = guid %>% str_replace_all("_", "")) %>%
  mutate(proband_guid = proband_guid %>% str_replace_all("_", "")) %>%
  filter(is.na(proband))

manifest <- read.csv('~/sebatlab Dropbox/G2MH/NetWorkPaper/AnalysisByAnjali/analysis/post_data_freeze/demographics_distributions/full_sample_manifest.csv')
manifest <- manifest %>%
  select(GUID, group_class_breakpoints) 

nda_uploads <- read.csv("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/demographics_distributions/G2MH_NDA_Upload_20250715.tsv",  sep = '\t')

datafreeze_manifest_combined <- merge(datafreeze, manifest, by.x = 'proband_guid', by.y = 'GUID', all.x = TRUE)
datafreeze_manifest_combined <- datafreeze_manifest_combined %>%
  separate(group_class_breakpoints, into = c("proband_locus", "proband_bp"), sep = " ") %>%
  filter(!is.na(proband_bp)) %>%
  filter(!(proband_bp %in% c('&', 'Unknown'))) %>%
  select(-proband)

# use this to get the lengths to transform each group
group_class_breakpoint_options <- c('BP4-5' = 549933,
                                    'BP1-3' = 563570,
                                    'BP2-3' = 223587, 
                                    'BP1-5' = 1717560, 
                                    'BP2-5' = 1377577, 
                                    'A-D' = 2062618,
                                    'A-B' = 1200181, 
                                    'A-C' = 1615181,
                                    'B-D' = 862437, 
                                    'C-D' = 362437,
                                    'D-E' = 1044162, 
                                    'D-F' = 1745621, 
                                    'D-G' = 2384162, 
                                    'D-H' = 2634162, 
                                    'E-F' = 486459, 
                                    'E-G' = 1125000,
                                    'E-H' = 1375000,
                                    'F-G' = 400000,
                                    'F-H' = 650000,
                                    .default = NA)

cnv_lengths <- manifest_with_cnv_length %>% select(guid, group_class_breakpoints, del_count, dup_count)

# add starting CNV lengths to matrix and manually adjust the 'other' breakpoint individuals
length_transformation <- merge(cnv_lengths, datafreeze_manifest_combined, by='guid', all.y= 'TRUE')
length_transformation <- length_transformation %>%
  mutate(proband_starting_length = recode(as.character(proband_bp), !!!group_class_breakpoint_options)) %>%
  filter(guid %in% nda_uploads$GUID) %>%
  mutate(proband_starting_length = ifelse(guid == 'NDARAU294LVK', 3117563, proband_starting_length)) %>%
  mutate(proband_starting_length = ifelse(guid == 'NDARFM574WAM', 1050000, proband_starting_length)) %>%
  mutate(proband_starting_length = ifelse(guid == 'NDARHR431PAD', 3808239, proband_starting_length)) %>%
  mutate(proband_starting_length = ifelse(guid == 'NDARJY224HN5', 3808239, proband_starting_length)) %>%
  mutate(proband_starting_length = ifelse(guid == 'NDARKZ481KXH', 1050000, proband_starting_length)) %>%
  mutate(proband_starting_length = ifelse(guid == 'NDARMP323GD1', 3117563, proband_starting_length)) %>%
  mutate(del_count = ifelse(is.na(del_count), 0, del_count)) %>%
  mutate(dup_count = ifelse(is.na(dup_count), 0, dup_count))
  
# now, for each person: take mean + starting_val + 3*sd and determine if they fall outside the range
# do this for both dup and del; if the proband has a del add their val; same if dup
# del_filter_out is TRUE if the sample's del count falls outside the 3SD range of mean+starting_val
gsa_sd_calcs_del <- length_transformation %>%
  mutate(deletion = grepl('Deletion', proband_locus)) %>%
  mutate(del_add = ifelse(deletion, proband_starting_length, 0)) %>%
  mutate(del_filter_out = (sample_sd_stats$mean_del + del_add + 3*sample_sd_stats$sd_del) < del_count) 
  
# now, for each person: take mean + starting_val + 3*sd and determine if they fall outside the range
# do this for both dup and del; if the proband has a del add their val; same if dup
# del_filter_out is TRUE if the sample's del count falls outside the 3SD range of mean+starting_val
gsa_sd_calcs_dup <- length_transformation %>%
  mutate(duplication = grepl('Duplication', proband_locus)) %>%
  mutate(dup_add = ifelse(duplication, proband_starting_length, 0)) %>%
  mutate(dup_filter_out = (sample_sd_stats$mean_dup + dup_add + 3*sample_sd_stats$sd_dup) < dup_count) 

write.csv(gsa_sd_calcs_del, '~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/cnv_calls/del_gsa_sd_calcs.csv')
write.csv(gsa_sd_calcs_dup, '~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/cnv_calls/dup_gsa_sd_calcs.csv')

# use the info in gsa_sd_calcs_del to remove people where filter=TRUE