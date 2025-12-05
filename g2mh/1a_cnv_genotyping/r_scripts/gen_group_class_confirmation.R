###################################################################################################
# Load necessary libraries
###################################################################################################
library(conflicted)
library(data.table) # For high-performance data manipulation
library(GenomicRanges)
library(ggnewscale)
library(patchwork)
library(scico)
library(tidyverse)  # This loads dplyr, ggplot2, and other core tidyverse packages

###################################################################################################
# Set up
###################################################################################################
# Specify preferences for conflicting functions
# I use dplyr functions like these ones. When you load other libraries, sometimes these function names are overridden by some other library's function with the same names.
# So I use this to ensure that the dplyr functions are the ones that are called (otherwise you can just write e.g. dplyr::select).
conflict_prefer("filter", "dplyr")
conflict_prefer("first", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
setwd("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/coverage_analysis")

###################################################################################################
# Loading coverage data
###################################################################################################

summary_files <- c(list.files("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/mosdepth_out/", pattern = "^[0-9]{3}-.*\\.mosdepth.summary.txt$", full.names = TRUE))
bed_files <- c(list.files("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/mosdepth_out/", pattern = "^[0-9]{3}-.*\\.regions.bed.gz$", full.names = TRUE))

# Create a named list to store the data
dfs_summary_files = lapply(summary_files, fread)
dfs_bed_files= lapply(bed_files, fread)

merge_summary_files = function(dfs, filepaths, chromosome) {
  dfs2 = lapply(seq_along(dfs), function(i) {
    df = dfs[[i]]
    filename = basename(filepaths[i])
    sample = sub("_.*", "", filename)
    #print(sample)
    df_n = df %>% 
      filter(chrom == chromosome) %>%
      mutate(sample = sample) %>%
      select(-min, -max, -length, -bases)
    #print(table(df_n$chrom))
    return(df_n)
  })
  df_merged_summary_files = Reduce(function(x, y) merge(x, y, by = c("chrom", "mean", "sample"), all = TRUE), dfs2)
  return(df_merged_summary_files)
}

df_merged_summary_chr16 = merge_summary_files(dfs_summary_files, summary_files, "chr16")
df_merged_summary_chr22 = merge_summary_files(dfs_summary_files, summary_files, "chr22")

merge_bed_files = function(dfs, filepaths, chrom) {
  dfs2 = lapply(seq_along(dfs), function(i) {
    df = dfs[[i]]
    filename =  basename(filepaths[i])
    sample = sub("_.*", "", filename)
    colnames(df)[colnames(df) == "V4"] = paste0(sample) # Rename column "V4" to the filepath
    df <- df %>% filter(V1 == chrom)
    #print(table(df$V1))
    return(df)
  })
  df_merged_bed_files = Reduce(function(x, y) merge(x, y, by = c("V1", "V2", "V3"), all = TRUE), dfs2)
  return(df_merged_bed_files)
}

df_merged_chr16_10000 = merge_bed_files(dfs_bed_files, bed_files, "chr16") %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3")
df_merged_chr22_10000 = merge_bed_files(dfs_bed_files, bed_files, "chr22") %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3")

###################################################################################################
# Loading recurrent regions file and reformatting the regions (for 16p.11.2 and 22q11.2)
###################################################################################################
recurrent_cnvs_hg38_sorted.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/recurrent_cnvs_hg38_sorted.bed") %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "cnv" = "V4")
recurrent_cnvs_hg38_sorted.bed_16p11.2 = recurrent_cnvs_hg38_sorted.bed %>%
  filter(grepl("16p11.2", cnv)) %>%
  filter(grepl("BP1_BP2|BP2_BP3|distal|BP4_BP5|proximal", cnv))
recurrent_cnvs_hg38_sorted.bed_22q11.2 = recurrent_cnvs_hg38_sorted.bed %>%
  filter(grepl("22q11.2", cnv)) %>%
  filter(grepl("A_B|B_C|C_D|D_E|E_F|F_G|F_H", cnv))
# Reformatting the recurrent regions a little bit.
recurrent_regions_16p11.2 = recurrent_cnvs_hg38_sorted.bed_16p11.2 %>%
  add_row(`#chrom` = "chr16", start = 27500000, end = 28000000, cnv = "16p11.2_leftFlank") %>%
  add_row(`#chrom` = "chr16", start = 31000000, end = 31500000, cnv = "16p11.2_rightFlank") %>%
  arrange(start) %>%
  filter(cnv != "16p11.2_distal") %>%
  filter(cnv != "16p11.2_proximal")
recurrent_regions_22q11.2 = recurrent_cnvs_hg38_sorted.bed_22q11.2 %>%
  add_row(`#chrom` = "chr22", start = 17500000, end = 18000000, cnv = "22q11.2_leftFlank") %>%
  add_row(`#chrom` = "chr22", start = 24600000, end = 25000000, cnv = "22q11.2_rightFlank") %>%
  add_row(`#chrom` = "chr22", start = 23950000, end = 24600000, cnv = "22q11.2_BP7_BP8_G_H") %>%
  arrange(start) %>%
  filter(cnv != "22q11.2_BP6_BP8_F_H") %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP1_BP2_A_B", "22q11.2_A_B", cnv)) %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP2_BP3_B_C", "22q11.2_B_C", cnv)) %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP3_BP4_C_D", "22q11.2_C_D", cnv)) %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP4_BP5_D_E", "22q11.2_D_E", cnv)) %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP5_BP6_E_F", "22q11.2_E_F", cnv)) %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP6_BP7_F_G", "22q11.2_F_G", cnv)) %>%
  mutate(cnv = if_else(cnv == "22q11.2_BP7_BP8_G_H", "22q11.2_G_H", cnv))

###################################################################################################
# Loading coverage data for problematic regions (gap, rmsk, simpleRepeat, segdups)
###################################################################################################
gap_coverage.chr16.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/gap_coverage.chr16.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "gap_coverage" = "V7")
gap_coverage.chr22.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/gap_coverage.chr22.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "gap_coverage" = "V7")
rmsk_coverage.chr16.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/rmsk_coverage.chr16.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "rmsk_coverage" = "V7")
rmsk_coverage.chr22.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/rmsk_coverage.chr22.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "rmsk_coverage" = "V7")
simpleRepeat_coverage.chr16.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/simpleRepeat_coverage.chr16.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "simpleRepeat_coverage" = "V7")
simpleRepeat_coverage.chr22.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/simpleRepeat_coverage.chr22.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "simpleRepeat_coverage" = "V7")
segdups_coverage.chr16.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/segdups_coverage.chr16.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "segdups_coverage" = "V7")
segdups_coverage.chr22.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/segdups_coverage.chr22.bed") %>%
  select(1, 2, 3, 7) %>%
  rename("#chrom" = "V1", "start" = "V2", "end" = "V3", "segdups_coverage" = "V7")

###################################################################################################
# Loading GC content BED files
###################################################################################################
chr16_window_size_10kb.gc.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/window_beds/chr16_window_size_10kb.gc.bed") %>%
  rename("#chrom" = "#1_usercol", "start" = "2_usercol", "end" = "3_usercol", "AT" = "4_pct_at", "GC" = "5_pct_gc") %>%
  select(`#chrom`, start, end, GC)
chr22_window_size_10kb.gc.bed = fread("~/sebatlab Dropbox/G2MH/cnv_scripts/resources/window_beds/chr22_window_size_10kb.gc.bed") %>%
  rename("#chrom" = "#1_usercol", "start" = "2_usercol", "end" = "3_usercol", "AT" = "4_pct_at", "GC" = "5_pct_gc") %>%
  select(`#chrom`, start, end, GC)

###################################################################################################
# Creating dataframes containing coverage data, GC content, and problematic region coverage 
# Then filtering the dataframes based on problematic region coverage
###################################################################################################
df_chr16 = df_merged_chr16_10000 %>%
  left_join(chr16_window_size_10kb.gc.bed) %>%
  left_join(gap_coverage.chr16.bed) %>%
  left_join(segdups_coverage.chr16.bed) %>%
  left_join(simpleRepeat_coverage.chr16.bed) %>%
  left_join(rmsk_coverage.chr16.bed) %>%
  filter(gap_coverage < 0.20) %>%
  filter(segdups_coverage < 0.20) %>%
  filter(simpleRepeat_coverage < 0.20)
df_chr22 = df_merged_chr22_10000 %>%
  left_join(chr22_window_size_10kb.gc.bed) %>%
  left_join(gap_coverage.chr22.bed) %>%
  left_join(segdups_coverage.chr22.bed) %>%
  left_join(simpleRepeat_coverage.chr22.bed) %>%
  left_join(rmsk_coverage.chr22.bed) %>%
  filter(gap_coverage < 0.20) %>%
  filter(segdups_coverage < 0.20) %>%
  filter(simpleRepeat_coverage < 0.20)

###################################################################################################
# GC correction
# (LOESS takes a while, try not to re-run)
###################################################################################################

# GC correction for df_chr16
sample_cols = names(df_chr16)[4:1089]
# Create a dataframe with identifier columns
df_pred = df_chr16 %>% select(`#chrom`, start, end)
for (sample in sample_cols) {
  # Fit LOESS for this sample
  loess_fit = loess(formula = as.formula(paste0("`", sample, "` ~ GC")), data = df_chr16)
  # Predict expected coverage based on GC
  df_pred[[sample]] <- predict(loess_fit, newdata = df_chr16)
}
# Create a new dataframe to hold adjusted coverage
df_chr16_gc_adjusted <- df_chr16 %>% select(`#chrom`, start, end)
for (sample in sample_cols) {
  baseline = mean(df_chr16[[sample]])
  df_chr16_gc_adjusted[[sample]] = df_chr16[[sample]] - df_pred[[sample]] + baseline
}

# GC correction for df_chr22
sample_cols = names(df_chr22)[4:1089]
# Create a dataframe with identifier columns
df_pred = df_chr22 %>% select(`#chrom`, start, end)
for (sample in sample_cols) {
  # Fit LOESS for this sample
  loess_fit = loess(formula = as.formula(paste0("`", sample, "` ~ GC")), data = df_chr22)
  # Predict expected coverage based on GC
  df_pred[[sample]] <- predict(loess_fit, newdata = df_chr22)
}
# Create a new dataframe to hold adjusted coverage
df_chr22_gc_adjusted <- df_chr22 %>% select(`#chrom`, start, end)
for (sample in sample_cols) {
  baseline = mean(df_chr22[[sample]])
  df_chr22_gc_adjusted[[sample]] = df_chr22[[sample]] - df_pred[[sample]] + baseline
}
# save the environment here so it can be reloaded and we can avoid rerunning long steps...
# save(....) to <>.RData
#save.image(file = "~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/my_environment.RData")

###################################################################################################
# Get coverage of recurrent CNV regions
# (The map_dfr steps also take a while, so try not to re-run them.)
###################################################################################################

# reload here if skipping earlier steps
load("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/my_environment.RData")

get_recurrent_coverage = function(sample, df, recurrent_df, chrom) {
  df_coverage = df %>% select(`#chrom`, start, end, matches(sample))
  df_coverage$coverage = df_coverage[, 4] # The sample's adjusted coverage
  mean_coverage = mean(df_coverage$coverage)
  gr_cnv = GRanges(
    seqnames = recurrent_df$`#chrom`,
    ranges = IRanges(start = recurrent_df$start, end = recurrent_df$end),
    cnv = recurrent_df$cnv
  )
  gr_cov <- GRanges(
    seqnames = df_coverage$`#chrom`,
    ranges = IRanges(start = df_coverage$start, end = df_coverage$end),
    coverage = df_coverage$coverage
  )
  hits <- findOverlaps(gr_cnv, gr_cov)
  # For each CNV region, calculate mean coverage over overlapping coverage bins
  mean_cov <- tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                     queryHits(hits), mean)
  recurrent_df$recurrent_mean_cov = NA
  recurrent_df$recurrent_mean_cov[as.integer(names(mean_cov))] = mean_cov
  # Get p-value
  region_stats <- tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                         queryHits(hits),
                         function(x) {
                           if(length(x) > 1) {
                             test <- t.test(x, mu = mean_coverage)
                             c(mean = mean(x), p_value = test$p.value)
                           } else {
                             c(mean = x, p_value = NA)
                           }
                         })
  region_stats = do.call(rbind, region_stats)
  #recurrent_df$recurrent_mean_cov <- region_stats[,"mean"]
  recurrent_df$p_value    <- region_stats[,"p_value"]
  output = recurrent_df %>%
    mutate(dummy = 1) %>% 
    #pivot_wider(id_cols = dummy, names_from = cnv, values_from = recurrent_mean_cov) %>%
    pivot_wider(id_cols = dummy, names_from = cnv, values_from = c(recurrent_mean_cov, p_value) ) %>%
    select(-dummy) %>%
    mutate(sample = sample) %>%
    mutate(!!sym(chrom) := mean_coverage) %>%
    select(sample, !!sym(chrom), everything())
  return(output)
}

sample_cols = names(df_chr16)[4:1089]
df_recurrent_chr16 = map_dfr(sample_cols, ~ get_recurrent_coverage(.x, df_chr16_gc_adjusted, recurrent_regions_16p11.2, chrom = "chr16"))
df_recurrent_chr22 = map_dfr(sample_cols, ~ get_recurrent_coverage(.x, df_chr22_gc_adjusted, recurrent_regions_22q11.2, chrom = "chr22"))
df_recurrent_wgs = full_join(df_recurrent_chr16, df_recurrent_chr22, by = "sample") %>%
  rename_with(~ str_remove(., "recurrent_mean_cov_"), contains("recurrent_mean_cov_")) %>%
  mutate(across(starts_with("16p11.2"), ~ . / chr16)) %>%
  mutate(across(starts_with("22q11.2"), ~ . / chr22)) %>%
  mutate(across(starts_with("16p11.2"), ~ round(. * 2) / 2 )) %>%
  mutate(across(starts_with("22q11.2"), ~ round(. * 2) / 2 )) %>%
  rowwise() %>%
  mutate(NON_SUBSEQUENT_CNV_16p11.2 = {
    vals = c_across(c(`16p11.2_BP1_BP2`, `16p11.2_BP2_BP3`, `16p11.2_BP4_BP5`))  # adjust your column range
    non_one_indices = which(vals != 1)
    # If 0 or 1 non-1 value, they’re trivially consecutive.
    if(length(non_one_indices) <= 1) {
      FALSE
    } else {
      any(diff(non_one_indices) != 1)
    }
  }) %>%
  mutate(NON_SUBSEQUENT_CNV_22q11.2 = {
    vals = c_across(c(`22q11.2_A_B`,  `22q11.2_B_C`,  `22q11.2_C_D`,  `22q11.2_D_E`,  `22q11.2_E_F`,  `22q11.2_F_G`, `22q11.2_G_H`))
    non_one_indices = which(vals != 1)
    # If 0 or 1 non-1 value, they’re trivially consecutive.
    if(length(non_one_indices) <= 1) {
      FALSE
    } else {
      any(diff(non_one_indices) != 1)
    }
  }) %>%
  mutate(MORE_THAN_TWO_UNIQUE_VALS_16p11.2 = n_distinct(c_across(c(`16p11.2_BP1_BP2`, `16p11.2_BP2_BP3`, `16p11.2_BP4_BP5`))) > 2) %>%
  mutate(MORE_THAN_TWO_UNIQUE_VALS_22q11.2 = n_distinct(c_across(c(`22q11.2_A_B`,  `22q11.2_B_C`,  `22q11.2_C_D`,  `22q11.2_D_E`,  `22q11.2_E_F`,  `22q11.2_F_G`, `22q11.2_G_H`))) > 2) %>%
  ungroup() %>%
  mutate(across(contains("p_value"), ~ -log10(.)))


###################################################################################################
# Stitching CNV calls
###################################################################################################

# Helper function to stitch a single row's CNV events.
stitch_event = function(vals, names) {
  valid_idx = which(vals != 1)
  if (length(valid_idx) == 0) return(NA_character_)
  groups = split(valid_idx, cumsum(c(1, diff(valid_idx) != 1)))
  event_type_map = function(x) {
    if (x == 0.5) "Deletion" else if (x == 1.5) "Duplication" else if (x == 2) "Triplication" else NA_character_
  }
  stitched = sapply(groups, function(idxs) {
    event_value = vals[idxs[1]]
    evt_type = event_type_map(event_value)
    # Assume column name format: prefix_BPstart_BPend
    first_parts = str_split(names[idxs[1]], "_")[[1]]
    last_parts  = str_split(names[idxs[length(idxs)]], "_")[[1]]
    prefix   = first_parts[1]
    start_bp = first_parts[2]
    end_bp   = last_parts[3]
    paste0(prefix, "_", evt_type, ' ', start_bp, "-", end_bp)
  })
  paste(stitched, collapse = ";")
}
stitch_cnv = function(df, cols_to_stitch, skip_column = "NON_SUBSEQUENT_CNV_16p11.2", output_col) {
  df %>%
    rowwise() %>%
    mutate(!!sym(output_col) := if (!.data[[skip_column]]) {
      stitch_event(c_across(all_of(cols_to_stitch)), cols_to_stitch)
    } else {
      NA_character_  # or another placeholder if desired
    }) %>%
    ungroup()
}


# Now you can apply it twice in a dplyr chain:
result_df = df_recurrent_wgs %>%
  stitch_cnv(
    cols_to_stitch = c("16p11.2_BP1_BP2", "16p11.2_BP2_BP3", "16p11.2_BP4_BP5"),
    skip_column = "NON_SUBSEQUENT_CNV_16p11.2",
    output_col = "16p11.2_stitched"
  ) %>%
  stitch_cnv(
    cols_to_stitch = c("22q11.2_A_B", "22q11.2_B_C", "22q11.2_C_D", "22q11.2_D_E", "22q11.2_E_F", "22q11.2_F_G", "22q11.2_G_H"),
    skip_column = "NON_SUBSEQUENT_CNV_22q11.2",
    output_col = "22q11.2_stitched"
  ) %>%
  mutate(GUID = str_replace(sample, ".*-", "")) %>%
  select(GUID, everything())

###################################################################################################
###################################################################################################
###################################################################################################
# importing sample manifest to cross-check
sample_manifest <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/demographics_distributions/full_sample_manifest.csv') %>%
  select(guid, group_class_breakpoints, comment, nda_dataset, redcap_data_access_group) %>%
  mutate(guid_nosep = gsub('_', '', guid))

# removing the duplicate samples that were reprocessed
sample_manifest <- sample_manifest[-16,]
sample_manifest <- sample_manifest[-1927,]
###################################################################################################
###################################################################################################
###################################################################################################

###################################################################################################
# check discrepancies
result_df_1 <- merge(result_df, sample_manifest, by.x='GUID', by.y='guid_nosep', all.x=TRUE)
result_df_1$`16p11.2_stitched` <- gsub("BP([0-9]+)-BP([0-9]+)", "BP\\1-\\2", result_df_1$`16p11.2_stitched`)
result_df_2 <- result_df_1 %>%
  mutate(stitched_breakpoints = ifelse(is.na(`16p11.2_stitched`), `22q11.2_stitched`, `16p11.2_stitched`)) %>%
  mutate(check_both = !is.na(`16p11.2_stitched`) & !is.na(`22q11.2_stitched`))

both_loci <- result_df_2 %>% filter(check_both)
result_df_filt <- result_df_2 %>% filter(!check_both)

# checking what already has gen_group_class
g2mh_project2_data_release_march2025 <- read_csv("~/sebatlab Dropbox/G2MH/DataFreeze2025/SoftRelease 2/g2mh_project2_data_release_march2025.csv") %>%
  select(guid, rarecnv_id, gen_group_class, rarecnv_family_id) %>%
  mutate(checked = ifelse(!is.na(gen_group_class), TRUE, FALSE)) %>%
  select(-gen_group_class)

sample_status <- merge(result_df_2, g2mh_project2_data_release_march2025, by.x = 'GUID', by.y = 'guid', all.x = TRUE)
not_checked <- sample_status %>% filter(nda_dataset == 'Jul_2025_WGS')

write.csv(sample_status, '~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/wgs_results_full_uncorrected.csv')

#discrepancies <- result_df_filt %>% filter(group_class_breakpoints != stitched_breakpoints)
###################################################################################################

###################################################################################################
# Plotting data (specifically discrepant data)
###################################################################################################

plot_data_with_segdups <- function(sample, df, recurrent_df, df_annot, x1 = 0, x2 = 50000000, y_lower_limit = 0, y_upper_limit = 100, site = '', redcap_ident = '') {
  ordered_levels <- recurrent_df[order(start)]$cnv %>% unique()
  df_coverage <- df %>% select(`#chrom`, start, end, contains(sample))
  df_coverage$coverage <- df_coverage[, 4]  # The sample's adjusted coverage
  mean_coverage <- mean(df_coverage$coverage, na.rm = TRUE)
  # Calculate recurrent_df mean coverage over CNV regions if needed
  gr_cnv <- GRanges(
    seqnames = recurrent_df$`#chrom`,
    ranges = IRanges(start = recurrent_df$start, end = recurrent_df$end),
    cnv = recurrent_df$cnv
  )
  gr_cov <- GRanges(
    seqnames = df_coverage$`#chrom`,
    ranges = IRanges(start = df_coverage$start, end = df_coverage$end),
    coverage = df_coverage$coverage
  )
  hits <- findOverlaps(gr_cnv, gr_cov)
  mean_cov <- tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                     queryHits(hits), mean)
  recurrent_df$recurrent_mean_cov <- NA
  recurrent_df$recurrent_mean_cov[as.integer(names(mean_cov))] <- mean_cov
  # Create the plot with two fill scales
  p <- ggplot() +
    # Background gradient using all intervals in df_segdups
    geom_rect(data = df_annot,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = 5, fill = segdups_coverage),
              inherit.aes = FALSE) +
    scale_fill_gradient(low = "transparent", high = "black") +
    ggnewscale::new_scale_fill() +
    # Recurrent regions from recurrent_df
    geom_rect(data = recurrent_df,
              aes(xmin = as.numeric(start), xmax = as.numeric(end),
                  ymin = -Inf, ymax = Inf, fill = cnv),
              inherit.aes = FALSE, alpha = 0.5) +
    scale_fill_brewer(limits = ordered_levels, palette = "Set3") +
    # Coverage points
    geom_point(data = df_coverage, aes(x = start, y = coverage),
               size = 1, alpha = 0.7, na.rm = TRUE) +
    # Mean coverage segments for recurrent regions
    geom_segment(data = recurrent_df,
                 aes(x = as.numeric(start), xend = as.numeric(end),
                     y = recurrent_mean_cov, yend = recurrent_mean_cov),
                 color = "blue", size = 2) +
    scale_x_continuous(
      limits = c(x1, x2),
      breaks = seq(0, max(df_coverage$start, na.rm = TRUE), by = 1000000),
      labels = scales::label_number()
    ) +
    scale_y_continuous(
      limits = c(y_lower_limit, y_upper_limit)
    ) +
    labs(y = "coverage", title = paste0(sample, " ", "(", site, "); redcap: ", redcap_ident)) +
    geom_hline(yintercept = mean_coverage, color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

plot_data <- function(sample, df, recurrent_df, df_annot, x1 = 0, x2 = 50000000, y_lower_limit = 0, y_upper_limit = 100, other_title = "") {
  ordered_levels <- recurrent_df[order(start)]$cnv %>% unique()
  df_coverage <- df %>% select(`#chrom`, start, end, contains(sample))
  df_coverage$coverage <- df_coverage[, 4]  # The sample's adjusted coverage
  mean_coverage <- mean(df_coverage$coverage, na.rm = TRUE)
  # Calculate recurrent_df mean coverage over CNV regions if needed
  gr_cnv <- GRanges(
    seqnames = recurrent_df$`#chrom`,
    ranges = IRanges(start = recurrent_df$start, end = recurrent_df$end),
    cnv = recurrent_df$cnv
  )
  gr_cov <- GRanges(
    seqnames = df_coverage$`#chrom`,
    ranges = IRanges(start = df_coverage$start, end = df_coverage$end),
    coverage = df_coverage$coverage
  )
  hits <- findOverlaps(gr_cnv, gr_cov)
  mean_cov <- tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                     queryHits(hits), mean)
  recurrent_df$recurrent_mean_cov <- NA
  recurrent_df$recurrent_mean_cov[as.integer(names(mean_cov))] <- mean_cov
  # Create the plot with two fill scales
  p <- ggplot() +
    scale_fill_gradient(low = "transparent", high = "black") +
    ggnewscale::new_scale_fill() +
    # Recurrent regions from recurrent_df
    geom_rect(data = recurrent_df,
              aes(xmin = as.numeric(start), xmax = as.numeric(end),
                  ymin = -Inf, ymax = Inf, fill = cnv),
              inherit.aes = FALSE, alpha = 0.5) +
    scale_fill_brewer(limits = ordered_levels, palette = "Set3") +
    # Coverage points
    geom_point(data = df_coverage, aes(x = start, y = coverage),
               size = 1, alpha = 0.7, na.rm = TRUE) +
    # Mean coverage segments for recurrent regions
    geom_segment(data = recurrent_df,
                 aes(x = as.numeric(start), xend = as.numeric(end),
                     y = recurrent_mean_cov, yend = recurrent_mean_cov),
                 color = "blue", size = 2) +
    scale_x_continuous(
      limits = c(x1, x2),
      breaks = seq(0, max(df_coverage$start, na.rm = TRUE), by = 1000000),
      labels = scales::label_number()
    ) +
    scale_y_continuous(
      limits = c(y_lower_limit, y_upper_limit)
    ) +
    labs(y = "coverage", title = paste0(sample, " ", "(", other_title, ")")) +
    geom_hline(yintercept = mean_coverage, color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

# Zoomed-in coordinates
x1_chr16 = 25000000
x2_chr16 = 33000000
x1_chr22 = 16000000
x2_chr22 = 27000000

for (i in 1:nrow(not_checked)) {
  id <- not_checked$guid[i]
  if (id == 'NDARVC394BJM') {next}
  site_id <- not_checked$redcap_data_access_group[i]
  redcap_info <- not_checked$group_class_breakpoints[i]
  #print(paste0(id, '; ', site_id, '; ', redcap_info))
  plot1 <- plot_data_with_segdups(id, df_chr16_gc_adjusted, recurrent_regions_16p11.2, segdups_coverage.chr16.bed, x1 = x1_chr16, x2 = x2_chr16, y_lower_limit = 00, y_upper_limit = 100, site=site_id, redcap_ident=redcap_info)
  plot2 <- plot_data_with_segdups(id, df_chr22_gc_adjusted, recurrent_regions_22q11.2, segdups_coverage.chr22.bed, x1 = x1_chr22, x2 = x2_chr22, y_lower_limit = 00, y_upper_limit = 100, site=site_id, redcap_ident=redcap_info)
  #print(plot1+plot2)
  plot_title <- paste0('discrepancies_plots/', id, '_', site_id, '.png')
  #print(plot_title)
  ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)
}

site_id = 'maastricht'
redcap_info = '22q11.2_Deletion E-F'
plot1 <- plot_data_with_segdups('NDARFV978LTD', df_chr16_gc_adjusted, recurrent_regions_16p11.2, segdups_coverage.chr16.bed, x1 = x1_chr16, x2 = x2_chr16, y_lower_limit = 00, y_upper_limit = 100, site=site_id, redcap_ident=redcap_info)
plot2 <- plot_data_with_segdups('NDARFV978LTD', df_chr22_gc_adjusted, recurrent_regions_22q11.2, segdups_coverage.chr22.bed, x1 = x1_chr22, x2 = x2_chr22, y_lower_limit = 00, y_upper_limit = 100, site=site_id, redcap_ident=redcap_info)
ggsave(paste0('discrepancies_plots/', 'NDARFV978LTD', '_', 'maastricht', '.png'), plot = plot1+plot2, width = 20, height = 8, dpi = 300)

###################################################################################################
# Plotting p-values from coverage data
###################################################################################################

p_value_cols_16p = c(
  "p_value_16p11.2_leftFlank",
  "p_value_16p11.2_BP1_BP2",
  "p_value_16p11.2_BP2_BP3",
  "p_value_16p11.2_BP4_BP5",
  "p_value_16p11.2_rightFlank"
)
p_value_cols_22q = c(
  "p_value_22q11.2_leftFlank",
  "p_value_22q11.2_A_B",
  "p_value_22q11.2_B_C",
  "p_value_22q11.2_C_D",
  "p_value_22q11.2_D_E",
  "p_value_22q11.2_E_F",
  "p_value_22q11.2_F_G",
  "p_value_22q11.2_G_H",
  "p_value_22q11.2_rightFlank"
)

plot_pvalue = function(df, pvalue_col, group_col) {
  p <- ggplot(df, aes(x = 0, y = {{ pvalue_col }}, color = factor({{ group_col }}))) +
    geom_jitter(height = 0) +
    labs(x = NULL, y = "-log10(p_value)", title = "log(p_value) Distributions") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    scale_color_manual(values = c("0.5" = "blue", "1" = "black", "1.5" = "red", "2" = "green")) +
    ggtitle(deparse(substitute(group_col)))
  print(p)
}
# Plotting p-values for 16p regions
plot_pvalue(df_recurrent_wgs, p_value_16p11.2_leftFlank, `16p11.2_leftFlank`)
plot_pvalue(df_recurrent_wgs, p_value_16p11.2_BP1_BP2, `16p11.2_BP1_BP2`)
plot_pvalue(df_recurrent_wgs, p_value_16p11.2_BP2_BP3, `16p11.2_BP2_BP3`)
plot_pvalue(df_recurrent_wgs, p_value_16p11.2_BP4_BP5, `16p11.2_BP4_BP5`)
plot_pvalue(df_recurrent_wgs, p_value_16p11.2_rightFlank, `16p11.2_rightFlank`)
# Plotting p-values for 22q regions
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_leftFlank, `22q11.2_leftFlank`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_A_B, `22q11.2_A_B`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_B_C, `22q11.2_B_C`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_C_D, `22q11.2_C_D`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_D_E, `22q11.2_D_E`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_E_F, `22q11.2_E_F`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_F_G, `22q11.2_F_G`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_G_H, `22q11.2_G_H`)
plot_pvalue(df_recurrent_wgs, p_value_22q11.2_rightFlank, `22q11.2_rightFlank`)

# remove what we won't need anymore
rm(dfs_bed_files)
rm(dfs_summary_files)
rm(loess_fit)
save.image(file = "~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/wgs.RData")

###################################################################################################
###################################################################################################
###################################################################################################
# GSA data
###################################################################################################
###################################################################################################
###################################################################################################

###################################################################################################
# Read and process GSA report files
###################################################################################################
load("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/wgs.RData")

read_sample_report <- function(filepath) {
  sample_name <- sub("\\.report.adjusted$", "", basename(filepath))
  #print(sample_name)
  df <- fread(filepath, header = TRUE, sep = "\t")
  #print(names(df))
  fread(filepath, header = TRUE, sep = "\t") %>%
    filter(`Log R Ratio` != "NaN") %>% 
    filter(Chr %in% c(16, 22)) %>%
    select(Chr, Position, `Log R Ratio`) %>%
    distinct(Chr, Position, .keep_all = TRUE) %>%
    rename(!!sample_name := `Log R Ratio`)
}
read_sample_reports <- function(folder) {
  files <- list.files(folder, pattern = "\\.report.adjusted$", full.names = TRUE)
  df_list <- lapply(files, read_sample_report)
  names(df_list) <- sub("\\.report.adjusted$", "", basename(files))
  df_list
}

report_dfs1 = read_sample_reports("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/adjusted_LRR")

#load("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/gsa1.RData")
df_reports = purrr::reduce(report_dfs1, full_join, by = c("Chr", "Position")) %>%
  mutate(Chr = paste0("chr", Chr)) %>%
  rename("#chrom" = "Chr") %>%
  arrange(`#chrom`, Position)

save.image(file = "~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/gsa1.RData")
# bad_file <- "~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/adjusted_LRR/887-NDARFM133AW6.report.adjusted"
# names(fread(bad_file, nrows = 0))

###################################################################################################
# Find overlaps with problematic regions 
###################################################################################################

df_gap = rbind(gap_coverage.chr16.bed, gap_coverage.chr22.bed)
df_rmsk = rbind(rmsk_coverage.chr16.bed, rmsk_coverage.chr22.bed)
df_segdups = rbind(segdups_coverage.chr16.bed, segdups_coverage.chr22.bed)
df_simpleRepeat = rbind(simpleRepeat_coverage.chr16.bed, simpleRepeat_coverage.chr22.bed)

# Convert the positions dataframe (df1) to a GRanges object, treating each position as a one-base range.
gr_positions <- GRanges(
  seqnames = df_reports$`#chrom`,
  ranges = IRanges(start = df_reports$Position, end = df_reports$Position)
)
# Convert the segdup dataframe (df2) to a GRanges object.
gr_gap <- GRanges(
  seqnames = df_gap$`#chrom`,
  ranges = IRanges(start = df_gap$start, end = df_gap$end),
  gap_coverage = df_gap$gap_coverage
)
# Convert the segdup dataframe (df2) to a GRanges object.
gr_rmsk <- GRanges(
  seqnames = df_rmsk$`#chrom`,
  ranges = IRanges(start = df_rmsk$start, end = df_rmsk$end),
  rmsk_coverage = df_rmsk$rmsk_coverage
)
# Convert the segdup dataframe (df2) to a GRanges object.
gr_segdups <- GRanges(
  seqnames = df_segdups$`#chrom`,
  ranges = IRanges(start = df_segdups$start, end = df_segdups$end),
  segdups_coverage = df_segdups$segdups_coverage
)
# Convert the segdup dataframe (df2) to a GRanges object.
gr_simpleRepeat <- GRanges(
  seqnames = df_simpleRepeat$`#chrom`,
  ranges = IRanges(start = df_simpleRepeat$start, end = df_simpleRepeat$end),
  simpleRepeat_coverage = df_simpleRepeat$simpleRepeat_coverage
)
# Find overlaps: which row in df1 (gr_positions) falls into which interval in df2 (gr_segdup)
gap_hits = findOverlaps(gr_positions, gr_gap)
rmsk_hits = findOverlaps(gr_positions, gr_rmsk)
segdups_hits = findOverlaps(gr_positions, gr_segdups)
simpleRepeat_hits = findOverlaps(gr_positions, gr_simpleRepeat)

# Create a new column in df1 with the corresponding segdup_coverage.
df_reports$gap_coverage = NA
df_reports$gap_coverage[queryHits(gap_hits)] = mcols(gr_gap)$gap_coverage[subjectHits(gap_hits)]
df_reports$rmsk_coverage = NA
df_reports$rmsk_coverage[queryHits(rmsk_hits)] <- mcols(gr_rmsk)$rmsk_coverage[subjectHits(rmsk_hits)]
df_reports$segdups_coverage <- NA
df_reports$segdups_coverage[queryHits(segdups_hits)] <- mcols(gr_segdups)$segdups_coverage[subjectHits(segdups_hits)]
df_reports$simpleRepeat_coverage <- NA
df_reports$simpleRepeat_coverage[queryHits(simpleRepeat_hits)] <- mcols(gr_simpleRepeat)$simpleRepeat_coverage[subjectHits(simpleRepeat_hits)]

###################################################################################################
# Filter problematic regions 
###################################################################################################

df_reports.filtered = df_reports %>%
  filter(gap_coverage < 0.20) %>%
  filter(segdups_coverage < 0.20) %>%
  filter(simpleRepeat_coverage < 0.20) %>%
  filter(! (if_any(everything(), is.na)) )

df_reports.filtered.chr16 = df_reports.filtered %>% filter(`#chrom` == "chr16")
df_reports.filtered.chr22 = df_reports.filtered %>% filter(`#chrom` == "chr22")

###################################################################################################
# Get recurrent CNV calls from GSA data
###################################################################################################

get_recurrent_coverage_report = function(sample, df, recurrent_df, chrom) {
  df = drop_na(df)
  df$start = df$Position
  df$end = df$Position
  df_coverage = df %>% select(`#chrom`, start, end, matches(sample))
  df_coverage$coverage = df_coverage[, 4] # The sample's adjusted coverage
  mean_coverage = mean(df_coverage$coverage)
  gr_cnv = GRanges(
    seqnames = recurrent_df$`#chrom`,
    ranges = IRanges(start = recurrent_df$start, end = recurrent_df$end),
    cnv = recurrent_df$cnv
  )
  gr_cov <- GRanges(
    seqnames = df_coverage$`#chrom`,
    ranges = IRanges(start = df_coverage$start, end = df_coverage$end),
    coverage = df_coverage$coverage
  )
  hits <- findOverlaps(gr_cnv, gr_cov)
  # For each CNV region, calculate mean coverage over overlapping coverage bins
  mean_cov <- tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                     queryHits(hits), mean)
  recurrent_df$recurrent_mean_cov = NA
  recurrent_df$recurrent_mean_cov[as.integer(names(mean_cov))] = mean_cov
  # Get p-value
  region_stats <- tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                         queryHits(hits),
                         function(x) {
                           if(length(x) > 1) {
                             test <- t.test(x, mu = mean_coverage)
                             c(mean = mean(x), p_value = test$p.value)
                           } else {
                             c(mean = x, p_value = NA)
                           }
                         })
  region_stats = do.call(rbind, region_stats)
  #recurrent_df$recurrent_mean_cov <- region_stats[,"mean"]
  recurrent_df$p_value    <- region_stats[,"p_value"]
  output = recurrent_df %>%
    mutate(dummy = 1) %>% 
    #pivot_wider(id_cols = dummy, names_from = cnv, values_from = recurrent_mean_cov) %>%
    pivot_wider(id_cols = dummy, names_from = cnv, values_from = c(recurrent_mean_cov, p_value) ) %>%
    select(-dummy) %>%
    mutate(sample = sample) %>%
    mutate(!!sym(chrom) := mean_coverage) %>%
    select(sample, !!sym(chrom), everything())
  return(output)
}

reports_sample_cols = names(df_reports.filtered)[3:587]
df_reports_chr16 = map_dfr(reports_sample_cols, ~ get_recurrent_coverage_report(.x, df_reports.filtered.chr16, recurrent_regions_16p11.2, chrom = "chr16"))
df_reports_chr22 = map_dfr(reports_sample_cols, ~ get_recurrent_coverage_report(.x, df_reports.filtered.chr22, recurrent_regions_22q11.2, chrom = "chr22"))
df_recurrent_gsa = full_join(df_reports_chr16, df_reports_chr22, by = "sample") %>%
  rename_with(~ str_remove(., "recurrent_mean_cov_"), contains("recurrent_mean_cov_")) %>%
  mutate(across(contains("p_value"), ~ -log10(.))) %>%
  mutate(PUTATIVE_CNV_16p11.2 = if_any(c("p_value_16p11.2_BP1_BP2", "p_value_16p11.2_BP2_BP3", "p_value_16p11.2_BP4_BP5"), ~ . > 10)) %>%
  mutate(PUTATIVE_CNV_22q11.2 = if_any(c("p_value_22q11.2_A_B", "p_value_22q11.2_B_C", "p_value_22q11.2_C_D", "p_value_22q11.2_D_E", "p_value_22q11.2_E_F", "p_value_22q11.2_F_G", "p_value_22q11.2_G_H"), ~ . > 10))

df_recurrent_gsa.calls = df_recurrent_gsa %>%
  mutate(
    `16p11.2_leftFlank` = case_when(
      p_value_16p11.2_leftFlank > 7 & `16p11.2_leftFlank` > 0  ~ 1.5,
      p_value_16p11.2_leftFlank > 7 & `16p11.2_leftFlank` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `16p11.2_BP1_BP2` = case_when(
      p_value_16p11.2_BP1_BP2 > 7 & `16p11.2_BP1_BP2` > 0  ~ 1.5,
      p_value_16p11.2_BP1_BP2 > 7 & `16p11.2_BP1_BP2` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `16p11.2_BP2_BP3` = case_when(
      p_value_16p11.2_BP2_BP3 > 7 & `16p11.2_BP2_BP3` > 0  ~ 1.5,
      p_value_16p11.2_BP2_BP3 > 7 & `16p11.2_BP2_BP3` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `16p11.2_BP4_BP5` = case_when(
      p_value_16p11.2_BP4_BP5 > 7 & `16p11.2_BP4_BP5` > 0  ~ 1.5,
      p_value_16p11.2_BP4_BP5 > 7 & `16p11.2_BP4_BP5` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `16p11.2_buffer_after` = case_when(
      p_value_16p11.2_rightFlank > 7 & `16p11.2_rightFlank` > 0  ~ 1.5,
      p_value_16p11.2_rightFlank > 7 & `16p11.2_rightFlank` < 0  ~ 0.5,
      TRUE ~ 1
    )
  ) %>%
  mutate(
    `22q11.2_buffer_before` = case_when(
      p_value_22q11.2_leftFlank > 7 & `22q11.2_leftFlank` > 0  ~ 1.5,
      p_value_22q11.2_leftFlank > 7 & `22q11.2_leftFlank` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_A_B` = case_when(
      p_value_22q11.2_A_B > 7 & `22q11.2_A_B` > 0  ~ 1.5,
      p_value_22q11.2_A_B > 7 & `22q11.2_A_B` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_B_C` = case_when(
      p_value_22q11.2_B_C > 7 & `22q11.2_B_C` > 0  ~ 1.5,
      p_value_22q11.2_B_C > 7 & `22q11.2_B_C` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_C_D` = case_when(
      p_value_22q11.2_C_D > 7 & `22q11.2_C_D` > 0  ~ 1.5,
      p_value_22q11.2_C_D > 7 & `22q11.2_C_D` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_D_E` = case_when(
      p_value_22q11.2_D_E > 7 & `22q11.2_D_E` > 0  ~ 1.5,
      p_value_22q11.2_D_E > 7 & `22q11.2_D_E` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_E_F` = case_when(
      p_value_22q11.2_E_F > 7 & `22q11.2_E_F` > 0  ~ 1.5,
      p_value_22q11.2_E_F > 7 & `22q11.2_E_F` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_F_G` = case_when(
      p_value_22q11.2_F_G > 7 & `22q11.2_F_G` > 0  ~ 1.5,
      p_value_22q11.2_F_G > 7 & `22q11.2_F_G` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_G_H` = case_when(
      p_value_22q11.2_G_H > 7 & `22q11.2_G_H` > 0  ~ 1.5,
      p_value_22q11.2_G_H > 7 & `22q11.2_G_H` < 0  ~ 0.5,
      TRUE ~ 1
    ),
    `22q11.2_buffer_after` = case_when(
      p_value_22q11.2_rightFlank > 7 & `22q11.2_rightFlank` > 0  ~ 1.5,
      p_value_22q11.2_rightFlank > 7 & `22q11.2_rightFlank` < 0  ~ 0.5,
      TRUE ~ 1
    ),
  ) %>%
  rowwise() %>%
  mutate(NON_SUBSEQUENT_CNV_16p11.2 = {
    vals = c_across(c(`16p11.2_BP1_BP2`, `16p11.2_BP2_BP3`, `16p11.2_BP4_BP5`))  # adjust your column range
    non_one_indices = which(vals != 1)
    # If 0 or 1 non-1 value, they’re trivially consecutive.
    if(length(non_one_indices) <= 1) {
      FALSE
    } else {
      any(diff(non_one_indices) != 1)
    }
  }) %>%
  mutate(NON_SUBSEQUENT_CNV_22q11.2 = {
    vals = c_across(c(`22q11.2_A_B`,  `22q11.2_B_C`,  `22q11.2_C_D`,  `22q11.2_D_E`,  `22q11.2_E_F`,  `22q11.2_F_G`, `22q11.2_G_H`))
    non_one_indices = which(vals != 1)
    # If 0 or 1 non-1 value, they’re trivially consecutive.
    if(length(non_one_indices) <= 1) {
      FALSE
    } else {
      any(diff(non_one_indices) != 1)
    }
  }) %>%
  mutate(MORE_THAN_TWO_UNIQUE_VALS_16p11.2 = n_distinct(c_across(c(`16p11.2_BP1_BP2`, `16p11.2_BP2_BP3`, `16p11.2_BP4_BP5`))) > 2) %>%
  mutate(MORE_THAN_TWO_UNIQUE_VALS_22q11.2 = n_distinct(c_across(c(`22q11.2_A_B`,  `22q11.2_B_C`,  `22q11.2_C_D`,  `22q11.2_D_E`,  `22q11.2_E_F`,  `22q11.2_F_G`, `22q11.2_G_H`))) > 2) %>%
  ungroup()

###################################################################################################
# Stitch recurrent CNV calls from GSA data
###################################################################################################

# Now you can apply it twice in a dplyr chain:
result_df_gsa = df_recurrent_gsa.calls %>%
  stitch_cnv(
    cols_to_stitch = c("16p11.2_BP1_BP2", "16p11.2_BP2_BP3", "16p11.2_BP4_BP5"),
    skip_column = "NON_SUBSEQUENT_CNV_16p11.2",
    output_col = "16p11.2_stitched"
  ) %>%
  stitch_cnv(
    cols_to_stitch = c("22q11.2_A_B", "22q11.2_B_C", "22q11.2_C_D", "22q11.2_D_E", "22q11.2_E_F", "22q11.2_F_G", "22q11.2_G_H"),
    skip_column = "NON_SUBSEQUENT_CNV_22q11.2",
    output_col = "22q11.2_stitched"
  ) %>%
  mutate(GUID = str_replace(sample, ".*-", "")) %>%
  select(GUID, everything())

###################################################################################################
# Plot p-values from GSA data (from the log R values)
###################################################################################################

# Plotting p-values for 16p11.2
plot_pvalue(df_recurrent_gsa, p_value_16p11.2_leftFlank, `16p11.2_leftFlank`)
plot_pvalue(df_recurrent_gsa, p_value_16p11.2_BP1_BP2, `16p11.2_BP1_BP2`)
plot_pvalue(df_recurrent_gsa, p_value_16p11.2_BP2_BP3, `16p11.2_BP2_BP3`)
plot_pvalue(df_recurrent_gsa, p_value_16p11.2_BP4_BP5, `16p11.2_BP4_BP5`)
plot_pvalue(df_recurrent_gsa, p_value_16p11.2_rightFlank, `16p11.2_rightFlank`)
# Plotting p-values for 22q11.2
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_leftFlank, `22q11.2_leftFlank`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_A_B, `22q11.2_A_B`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_B_C, `22q11.2_B_C`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_C_D, `22q11.2_C_D`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_D_E, `22q11.2_D_E`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_E_F, `22q11.2_E_F`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_F_G, `22q11.2_F_G`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_G_H, `22q11.2_G_H`)
plot_pvalue(df_recurrent_gsa, p_value_22q11.2_rightFlank, `22q11.2_rightFlank`)

###################################################################################################
# Plot report data
###################################################################################################

###################################################################################################
# check discrepancies
result_df_11 <- merge(result_df_gsa, sample_manifest, by.x='GUID', by.y='guid_nosep')
result_df_11$`16p11.2_stitched` <- gsub("BP([0-9]+)-BP([0-9]+)", "BP\\1-\\2", result_df_11$`16p11.2_stitched`)
result_df_21 <- result_df_11 %>%
  mutate(stitched_breakpoints = ifelse(is.na(`16p11.2_stitched`), `22q11.2_stitched`, `16p11.2_stitched`)) %>%
  mutate(check_both = !is.na(`16p11.2_stitched`) & !is.na(`22q11.2_stitched`))

both_loci <- result_df_21 %>% filter(check_both)
result_df_filt1 <- result_df_21 %>% filter(!check_both)

sample_status_gsa <- merge(result_df_21, g2mh_project2_data_release_march2025, by.x = 'GUID', by.y = 'guid')
not_checked_gsa <- sample_status_gsa %>% filter(!checked) %>% filter(nda_dataset == 'Jul_2025_GSA')

#save.image(file = "~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/processed_new_data_env.RData")
write.csv(sample_status_gsa, '~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/gsa_results_full_uncorrected.csv')

#discrepancies <- result_df_filt %>% filter(group_class_breakpoints != stitched_breakpoints)
###################################################################################################

###################################################################################################
# Plotting data (specifically discrepant data)
###################################################################################################

plot_report_data = function(sample, df, recurrent_df, x1 = 0, x2 = 50000000, chrom, site, redcap_ident) {
  ordered_levels = recurrent_df[order(start)]$cnv %>% unique()
  df_coverage = df %>% 
    filter(`#chrom` == chrom) %>%
    rename("start" = "Position") %>%
    mutate(end = start) %>%
    select(`#chrom`, start, end, contains(sample))
  df_coverage$coverage = df_coverage[, 4] # The sample's adjusted coverage
  mean_coverage = mean(df_coverage$coverage)
  gr_cov <- GRanges(
    seqnames = df_coverage$`#chrom`,
    ranges = IRanges(start = df_coverage$start, end = df_coverage$end),
    coverage = df_coverage$coverage
  )
  gr_cnv = GRanges(
    seqnames = recurrent_df$`#chrom`,
    ranges = IRanges(start = recurrent_df$start, end = recurrent_df$end),
    cnv = recurrent_df$cnv
  )
  hits = findOverlaps(gr_cnv, gr_cov)
  # For each CNV region, calculate mean coverage over overlapping coverage bins
  mean_cov = tapply(mcols(gr_cov)$coverage[subjectHits(hits)],
                    queryHits(hits), mean)
  recurrent_df$recurrent_mean_cov = NA
  recurrent_df$recurrent_mean_cov[as.integer(names(mean_cov))] = mean_cov
  ggplot(df_coverage, aes(x = start, y = coverage)) + 
    geom_point(size = 1, alpha = 0.9, na.rm = TRUE) +
    geom_rect(data = recurrent_df, aes(xmin = as.numeric(start), xmax = as.numeric(end), ymin = -Inf, ymax = Inf, fill = cnv), inherit.aes = FALSE, alpha = 0.5) +
    geom_segment(data = recurrent_df, aes(x = as.numeric(start), xend = as.numeric(end), y = recurrent_mean_cov, yend = recurrent_mean_cov), color = "blue", size = 1) +
    scale_fill_brewer(limits = ordered_levels, palette = "Set3") +
    scale_x_continuous(
      limits = c(x1, x2),
      breaks = seq(0, max(df_coverage$start, na.rm = TRUE), by = 1000000),
      labels = scales::label_number()) +  # Format x-axis labels
    scale_y_continuous(
      limits = c(-2, 2)
    ) +
    labs(y = "LRR", title = paste0(sample, " ", "(", site, "); redcap: ", redcap_ident)) +
    geom_hline(yintercept = mean_coverage, color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (i in 1:nrow(not_checked_gsa)) {
  id <- not_checked_gsa$guid[i]
  site_id <- not_checked_gsa$redcap_data_access_group[i]
  redcap_info <- not_checked_gsa$group_class_breakpoints[i]
  #print(paste0(id, '; ', site_id, '; ', redcap_info))
  plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
  plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
  #print(plot1+plot2)
  plot_title <- paste0('discrepancies_plots_gsa_adjusted/', id, '_', site_id, '.png')
  #print(plot_title)
  ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)
  #break
}

id = 'NDAREF566FJP'
site_id = 'ucla'
redcap_info = 'no_cnv'
plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
plot_title <- paste0('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/weird_samples/', id, '_', site_id, '.png')
ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)

id = 'NDARCJ441XEP'
site_id = 'ucla'
redcap_info = 'no_cnv'
plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
plot_title <- paste0('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/weird_samples/', id, '_', site_id, '.png')
ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)

id = 'NDARZK179ZPX'
site_id = 'montreal'
redcap_info = 'no_cnv'
plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
plot_title <- paste0('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/weird_samples/', id, '_', site_id, '.png')
ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)

id = 'NDARHP361ZMZ'
site_id = 'pennchop'
redcap_info = '22q11.2_Deletion B-D'
plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
plot_title <- paste0('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/weird_samples/', id, '_', site_id, '.png')
ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)

id = 'NDARYW554GDN'
site_id = 'cardiff'
redcap_info = '22q11.2_Duplication A-D'
plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
plot_title <- paste0('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/weird_samples/', id, '_', site_id, '.png')
ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)

###################################################################################################
# plot ids that are above 3SD del length in GSA or 3SD dup length in GSA
id = 'NDARWF808APK'
site_id = 'cardiff'
redcap_info = 'no_cnv'
plot1 <- plot_report_data(id, df_reports.filtered, recurrent_regions_16p11.2, x1 = x1_chr16, x2 = x2_chr16, chrom = "chr16", site=site_id, redcap_ident=redcap_info)
plot2 <- plot_report_data(id, df_reports.filtered, recurrent_regions_22q11.2, x1 = x1_chr22, x2 = x2_chr22, chrom = "chr22", site=site_id, redcap_ident=redcap_info)
plot_title <- paste0('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/gsa_put_filter/', id, '_', site_id, '.png')
ggsave(plot_title, plot = plot1+plot2, width = 20, height = 8, dpi = 300)
###################################################################################################
# combine wgs and gsa
combined_all_samples <- rbind(sample_status %>% select(GUID, nda_dataset, redcap_data_access_group, group_class_breakpoints, rarecnv_family_id, stitched_breakpoints),
                              sample_status_gsa %>% select(GUID, nda_dataset, redcap_data_access_group, group_class_breakpoints, rarecnv_family_id, stitched_breakpoints))

combined_new_samples <- rbind(not_checked %>% select(GUID, nda_dataset, redcap_data_access_group, group_class_breakpoints, rarecnv_family_id, stitched_breakpoints),
                              not_checked_gsa %>% select(GUID, nda_dataset, redcap_data_access_group, group_class_breakpoints, rarecnv_family_id, stitched_breakpoints))
#combined_new_samples <- not_checked_gsa %>% select(GUID, nda_dataset, redcap_data_access_group, group_class_breakpoints, rarecnv_family_id, stitched_breakpoints)
combined_new_samples <- combined_new_samples %>%
  mutate(same_locus = (grepl("16p11.2", combined_new_samples$group_class_breakpoints) &
                         grepl("16p11.2", combined_new_samples$stitched_breakpoints)) |
           (grepl("22q11.2", combined_new_samples$group_class_breakpoints) &
              grepl("22q11.2", combined_new_samples$stitched_breakpoints)) |
           (combined_new_samples$group_class_breakpoints == "no_cnv" &
              is.na(combined_new_samples$stitched_breakpoints))
  ) %>%
  mutate(prev_unknown = grepl("Unknown", combined_new_samples$group_class_breakpoints)) %>%
  mutate(fully_same = (combined_new_samples$group_class_breakpoints == combined_new_samples$stitched_breakpoints)) %>%
  mutate(fully_same = ifelse(is.na(fully_same), TRUE, fully_same)) %>%
  mutate(same_locus_diff_bp = same_locus & !prev_unknown & !fully_same) %>%
  mutate(diff_locus = !same_locus) 
###################################################################################################

###################################################################################################
wgs_calls <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/wgs_results_full_uncorrected.csv') %>%
  select(rarecnv_id, nda_dataset, GUID, stitched_breakpoints) %>%
  #mutate(guid = gsub('_', '', GUID)) %>%
  filter(nda_dataset == 'Jul_2025_WGS') 
gsa_calls <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/gsa_results_full_uncorrected.csv') %>%
  select(rarecnv_id, nda_dataset, GUID, stitched_breakpoints) %>%
  #mutate(guid = gsub('_', '', GUID)) %>%
  filter(nda_dataset == 'Jul_2025_GSA')

###################################################################################################

###################################################################################################
df_wgs_gsa <- wgs_calls %>%
  rbind(gsa_calls) %>%
  mutate(has_16p11.2_deletion = grepl('16p11.2_Deletion', stitched_breakpoints)) %>%
  mutate(has_16p11.2_duplication = grepl('16p11.2_Duplication', stitched_breakpoints)) %>%
  mutate(has_22q11.2_deletion = grepl('22q11.2_Deletion', stitched_breakpoints)) %>%
  mutate(has_22q11.2_duplication = grepl('22q11.2_Duplication', stitched_breakpoints)) %>%
  mutate(breakpoints_16p = ifelse(has_16p11.2_deletion | has_16p11.2_duplication, sub(".*\\s+", "", stitched_breakpoints), NA)) %>%
  mutate(breakpoints_22q = ifelse(has_22q11.2_deletion | has_22q11.2_duplication, sub(".*\\s+", "", stitched_breakpoints), NA)) %>%
  filter(nda_dataset %in% c('Jul_2025_WGS', 'Jul_2025_GSA')) %>%
  # changing incorrect calls and flagging discrepant samples
  mutate(discrepancy_with_RedCap = ifelse(GUID %in% c('NDARBE330WPB', 'NDARMB155JTZ', 'NDARNH098VWR', 
                                                      'NDARHL050KY2', 'NDARZK179ZPX', 'NDARVU809TJM',
                                                      'NDARUL801BZ4', 'NDARME032MUJ', 'NDARAC013KU9', 
                                                      'NDARXE398XMK', 'NDARNT837CL1', 'NDARVC394BJM'), TRUE, FALSE)) %>%
  mutate(gen_group_class = NA) %>%
  mutate(gen_group_nonstandard_type = NA) %>%
  mutate(gen_group_nonstandard_type_breakpoints = NA) %>%
  mutate(gen_group_class = ifelse(has_22q11.2_deletion, 1, gen_group_class)) %>%
  mutate(gen_group_class = ifelse(has_22q11.2_duplication, 2, gen_group_class)) %>%
  mutate(gen_group_class = ifelse(has_16p11.2_deletion, 3, gen_group_class)) %>%
  mutate(gen_group_class = ifelse(has_16p11.2_duplication, 4, gen_group_class)) %>%
  # NDARNH098VWR
  mutate(has_22q11.2_deletion = ifelse(GUID == 'NDARNH098VWR', TRUE, has_22q11.2_deletion)) %>%
  mutate(gen_group_class = ifelse(GUID == 'NDARNH098VWR', 1, gen_group_class)) %>%
  mutate(breakpoints_22q = ifelse(GUID == 'NDARNH098VWR', 'oth', breakpoints_22q)) %>%
  mutate(gen_group_nonstandard_type = ifelse(GUID == 'NDARNH098VWR', 'chr22:21960000-22220000 Deletion', gen_group_nonstandard_type)) %>%
  mutate(gen_group_nonstandard_type_breakpoints = ifelse(GUID == 'NDARNH098VWR', 'chr22:21960000-22220000', gen_group_nonstandard_type_breakpoints)) %>%
  # NDARVU809TJM
  mutate(has_22q11.2_deletion = ifelse(GUID == 'NDARVU809TJM', TRUE, has_22q11.2_deletion)) %>%
  mutate(gen_group_class = ifelse(GUID == 'NDARVU809TJM', 1, gen_group_class)) %>%
  mutate(breakpoints_22q = ifelse(GUID == 'NDARVU809TJM', 'oth', breakpoints_22q)) %>%
  mutate(gen_group_nonstandard_type = ifelse(GUID == 'NDARVU809TJM', 'chr22:21960000-22220000 Deletion', gen_group_nonstandard_type)) %>%
  mutate(gen_group_nonstandard_type_breakpoints = ifelse(GUID == 'NDARVU809TJM', 'chr22:21960000-22220000', gen_group_nonstandard_type_breakpoints)) %>%
  # NDARUL801BZ4
  mutate(breakpoints_22q = ifelse(GUID == 'NDARUL801BZ4', 'A-D', breakpoints_22q)) %>%
  #NDARXD988NAU
  mutate(has_22q11.2_deletion = ifelse(GUID == 'NDARXD988NAU', TRUE, has_22q11.2_deletion)) %>%
  mutate(breakpoints_22q = ifelse(GUID == 'NDARXD988NAU', 'oth', breakpoints_22q)) %>%
  mutate(gen_group_nonstandard_type = ifelse(GUID == 'NDARXD988NAU', 'chr22:21960000-22220000 Deletion', gen_group_nonstandard_type)) %>%
  mutate(gen_group_nonstandard_type_breakpoints = ifelse(GUID == 'NDARXD988NAU', 'chr22:21960000-22220000', gen_group_nonstandard_type_breakpoints)) %>%
  mutate(gen_group_class = ifelse(GUID == 'NDARXD988NAU', 12, gen_group_class)) %>%
  # NDARCK212TF3
  mutate(has_22q11.2_duplication = ifelse(GUID == 'NDARCK212TF3', TRUE, has_22q11.2_deletion)) %>%
  mutate(breakpoints_22q = ifelse(GUID == 'NDARCK212TF3', 'A-D', breakpoints_22q)) %>%
  mutate(gen_affected = ifelse(is.na(gen_group_class), 0, 1)) %>%
  # NDARHJ855WUV
  mutate(breakpoints_22q = ifelse(GUID == 'NDARHJ855WUV', 'oth', breakpoints_22q)) %>%
  mutate(gen_group_nonstandard_type = ifelse(GUID == 'NDARHJ855WUV', 'chr22:B-F Deletion', gen_group_nonstandard_type)) %>%
  mutate(gen_group_nonstandard_type_breakpoints = ifelse(GUID == 'NDARHJ855WUV', 'B-F', gen_group_nonstandard_type_breakpoints)) 
  

df_redcap <- df_wgs_gsa %>%
  select(GUID, rarecnv_id, gen_affected, discrepancy_with_RedCap, gen_group_class, breakpoints_16p, breakpoints_22q, gen_group_nonstandard_type, gen_group_nonstandard_type_breakpoints) %>%
  mutate(gen_group_16p_type = case_when(
    str_detect(breakpoints_16p, "BP4-5") ~ "1",
    str_detect(breakpoints_16p, "BP1-3") ~ "2",
    str_detect(breakpoints_16p, "BP2-3") ~ "3",
    str_detect(breakpoints_16p, "BP1-5") ~ "4",
    str_detect(breakpoints_16p, "BP2-5") ~ "5",
    str_detect(breakpoints_16p, "oth") ~ "oth",
    TRUE ~ NA_character_
  )) %>%
  mutate(gen_group_22q_type = case_when(
    str_detect(breakpoints_22q, "A-D") ~ "1",
    str_detect(breakpoints_22q, "A-B") ~ "2",
    str_detect(breakpoints_22q, "A-C") ~ "3",
    str_detect(breakpoints_22q, "B-D") ~ "4",
    str_detect(breakpoints_22q, "C-D") ~ "5",
    str_detect(breakpoints_22q, "D-E") ~ "6",
    str_detect(breakpoints_22q, "D-F") ~ "7",
    str_detect(breakpoints_22q, "D-G") ~ "8",
    str_detect(breakpoints_22q, "D-H") ~ "9",
    str_detect(breakpoints_22q, "E-F") ~ "10",
    str_detect(breakpoints_22q, "E-G") ~ "11",
    str_detect(breakpoints_22q, "E-H") ~ "12",
    str_detect(breakpoints_22q, "F-G") ~ "13",
    str_detect(breakpoints_22q, "F-H") ~ "14",
    str_detect(breakpoints_22q, "Other") ~ "oth",
    TRUE ~ NA_character_
  )) %>%
  select(GUID, rarecnv_id, gen_affected, gen_group_class, gen_group_16p_type, gen_group_22q_type, gen_group_nonstandard_type, gen_group_nonstandard_type_breakpoints, discrepancy_with_RedCap) %>%
  mutate(across(everything(), ~replace(., is.na(.), ""))) %>%
  mutate(discrepancy_with_RedCap = as.logical(discrepancy_with_RedCap)) 

# integrate with datafreeze genetics now
phenotypes <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/demographics_distributions/RareCNVPhenotypicCol_DATA_2025-07-22_1927.csv') %>%
  filter(!(guid == '')) %>%
  filter(!(rarecnv_id %in% c(1, 2, 3))) %>%
  mutate(GUID_nosep = gsub('_', '', guid)) %>%
  mutate(GUID = guid) %>%
  select(GUID, GUID_nosep, rarecnv_id, gen_affected, gen_group_class, gen_group_16p_type, gen_group_22q_type, gen_group_nonstandard_type, gen_group_nonstandard_type_breakpoints) %>%
  mutate(discrepancy_with_RedCap = FALSE)
  

phenotypes$GUID <- as.character(phenotypes$GUID)
phenotypes$GUID_nosep <- as.character(phenotypes$GUID_nosep)
df_redcap$GUID  <- as.character(df_redcap$GUID)
phenotypes$gen_affected <- as.character(phenotypes$gen_affected)
df_redcap$gen_affected  <- as.character(df_redcap$gen_affected)

id_nums <- phenotypes %>% select(GUID, GUID_nosep, rarecnv_id)
phenotypes <- phenotypes %>% 
  select(-rarecnv_id, -GUID) %>%
  mutate(GUID = GUID_nosep) %>%
  select(-GUID_nosep)
df_redcap <- df_redcap %>% select(-rarecnv_id)

datafreeze_genetics_new_comb <- phenotypes %>%
  rows_update(df_redcap, by = "GUID")  %>%
  merge(id_nums, by.x = 'GUID', by.y='GUID_nosep') %>%
  mutate(GUID = GUID.y) %>%
  select(GUID, rarecnv_id, gen_affected, gen_group_class, gen_group_16p_type, gen_group_22q_type, gen_group_nonstandard_type, gen_group_nonstandard_type_breakpoints, discrepancy_with_RedCap) %>%
  filter(!GUID %in% c('nadlknalnf')) %>%
  mutate(affected = ifelse(GUID == 'NDARMB155JTZ', 0, gen_affected)) %>%
  mutate(affected = ifelse(GUID == 'NDARBE330WPB', 0, gen_affected))

write.csv(datafreeze_genetics_new_comb, '~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/gen_group_classes_with_new_samples.csv')

# test
wgs_calls <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/wgs_results_full_uncorrected.csv')
gsa_calls <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/cnv_genotyping/gsa_results_full_uncorrected.csv')
test <- datafreeze_genetics_new_comb %>%
  mutate(GUID = gsub('_', '', GUID)) %>%
  filter(gen_affected == 1)
test2 <- datafreeze_genetics_new_comb %>% filter(GUID %in% setdiff(wgs_calls$GUID, test$GUID))
###################################################################################################
# SEX CHECKS FOR WGS

wgs_bed_files <- c(list.files("~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/mosdepth_out/", pattern = "^[0-9]{3}-.*\\.regions.bed.gz$", full.names = TRUE))

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

sex_all_wgs <- wgs_bed_files %>%
  lapply(check_sex) %>%
  bind_rows()

check_sex('/Users/athabasca/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/mosdepth_out/896-NDARBE330WPB_10000.regions.bed.gz')
check_sex('/Users/athabasca/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/mosdepth_out/896-NDARMB155JTZ_10000.regions.bed.gz')
