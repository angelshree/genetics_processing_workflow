library(dplyr)
library(ggplot2)

################################################################################
# get phenotypes (IQ) and major group class for all subjects
################################################################################
# set up phenotypes
phenotypes <- read.csv('~/sebatlab Dropbox/G2MH/DataFreeze2025/Download_08_27_2025/Project2_DataRelease/August2025/g2mh_project2_data_release_082025.csv')
group_class_mapping <- c("1" = '22q11.2_Deletion', "2" = '22q11.2_Duplication',
                         "3" = '16p11.2_Deletion', "4" = '16p11.2_Duplication',
                         "5" = '22q11.2_Duplication & 16p11.2_Deletion', 
                         "6" = '22q11.2_Duplication & 16p11.2_Duplication',
                         "7" = '22q11.2_Deletion & 16p11.2_Deletion',
                         "8" = '22q11.2_Deletion & 16p11.2_Duplication',
                         "11" = '22q11.2_Triplication',
                         "20" = '22q11.2_Other')
group_16p_mapping <- c("0" = 'Unknown', "1" = 'BP4-5', "2" = 'BP1-3', "3" = 'BP2-3',
                       "4" = 'BP1-5', "5" = 'BP2-5')
group_22q_mapping <- c('0' = 'Unknown', '1' = 'A-D', '2' = 'A-B', '3' = 'A-C',
                       '4' = 'B-D', '5' = 'C-D', '6' = 'D-E', '7' = 'D-F',
                       '8' = 'D-G', '9' = 'D-H', '10' = 'E-F', '11' = 'E-G',
                       '12' = 'E-H', '13' = 'F-G', '14' = 'F-H', '15' = 'TBX1_variant',
                       'oth' = 'Other')
phenotypes <- phenotypes %>%
  mutate(across(names(phenotypes), ~ ifelse(. == '', NA, .))) %>%
  mutate(across(c('group_16p_type', 'group_22q_type'), ~ ifelse(is.na(group_class), NA, .))) %>%
  mutate(final_group_class = ifelse(!is.na(gen_group_class), gen_group_class, group_class)) %>%
  mutate(final_group_22q_type = ifelse(!is.na(gen_group_22q_type), gen_group_22q_type, group_22q_type)) %>%
  mutate(final_group_16p_type = ifelse(!is.na(gen_group_16p_type), gen_group_16p_type, group_16p_type)) %>%
  mutate(final_group_class = recode(as.character(final_group_class), !!!group_class_mapping)) %>%
  mutate(final_group_16p_type = recode(as.character(final_group_16p_type), !!!group_16p_mapping)) %>%
  mutate(final_group_22q_type = recode(as.character(final_group_22q_type), !!!group_22q_mapping)) %>%
  mutate(final_group_class = ifelse(is.na(final_group_class), 'no_cnv', final_group_class)) %>%
  filter(final_group_class %in% c('22q11.2_Deletion', '22q11.2_Duplication', '16p11.2_Deletion', '16p11.2_Duplication', 'no_cnv')) %>%
  mutate(final_group_16p_type_1 = ifelse((final_group_class %in% c('16p11.2_Deletion', '16p11.2_Duplication') & is.na(final_group_16p_type)), 'Unknown', final_group_16p_type)) %>%
  mutate(final_group_22q_type_1 = ifelse((final_group_class %in% c('22q11.2_Deletion', '22q11.2_Duplication') & is.na(final_group_22q_type)), 'Unknown', final_group_22q_type)) %>%
  filter(!(final_group_22q_type_1 %in% c('Unknown', 'Other'))) %>%
  filter(!(final_group_16p_type_1 %in% c('Unknown', 'Other', 'oth'))) %>%
  mutate(group_class_breakpoints = ifelse(final_group_class %in% c('22q11.2_Deletion', '22q11.2_Duplication'), paste0(final_group_class, ' ', final_group_22q_type_1),
                                          ifelse(final_group_class %in% c('16p11.2_Deletion', '16p11.2_Duplication'), paste0(final_group_class, ' ', final_group_16p_type_1), final_group_class)))
phenotypes = phenotypes %>%
  select(guid, dem_sex, dem_age, smry_iq_age, group_class_breakpoints, smry_iq_fsiq, smry_iq_viq, smry_iq_piq, redcap_data_access_group) %>%
  mutate(dem_age2 = dem_age*dem_age) %>%
  mutate(smry_iq_age2 = smry_iq_age*smry_iq_age) %>%
  filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5')) %>%
  mutate(age_final = ifelse(!is.na(smry_iq_age), smry_iq_age, dem_age)) %>%
  mutate(age2 = age_final**2)



################################################################################
# import all files and filter for only the rare variants
################################################################################

# read in files
vep_gene_variants_files <- list.files('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/vep_gene_annotations', pattern = "\\.txt$", full.names=TRUE)
vep_gene_variants <- lapply(vep_gene_variants_files, fread)
vep_gene_variants <- do.call(rbind, vep_gene_variants)

# filter out common gnomAD variants
vep_gene_variants <- vep_gene_variants %>%
  filter(gnomADe_V4_AF <= 0.05)

# filter out variants that are too common in this cohort specifically
# take out any specific variants in more than 5% of the population
variant_counts <- vep_gene_variants %>%
  group_by(chrom, start, end, consequence, ref, alt) %>%
  summarise(count = n(), .groups = "drop")
n_samples <- 1070
variant_counts <- variant_counts %>%
  mutate(cohort_freq = count / n_samples) %>%
  filter(cohort_freq <= 0.05)
variants_filtered <- vep_gene_variants %>%
  inner_join(variant_counts, by = c("chrom", "start", "end", "consequence", "ref", "alt")) %>%
  mutate(guid = substr(offspring, 5, nchar(offspring))) %>%
  filter(is.na(MPC_score) | MPC_score > 2)

################################################################################
# combine variants with group class and IQ scores
################################################################################

# gene lists of interest
fu_genes_matrix <- read_xlsx('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/individual_gene_variants/Fu_ASD_DD_NDD_Genes.xlsx')
nimh_genes_matrix <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/individual_gene_variants/NIMH_gene_list.csv')

variants_class_combined <- merge(phenotypes, variants_filtered, by='guid', all.x=TRUE)

# examine LOF and missense variants in NIMH genes
variants_nimh <- variants_class_combined %>%
  filter(symbol %in% nimh_genes_matrix$gene)
table(variants_nimh$consequence)
table(variants_nimh$symbol)

# separate missense variants from the rest (since they function differently)
variants_nimh_lof <- variants_nimh %>% filter(!(consequence == 'missense_variant'))
variants_nimh_missense <- variants_nimh %>%
  filter(consequence == 'missense_variant') %>%
  filter(!is.na(MPC_score))

# plot gene counts
gene_counts <- rbind(variants_nimh_missense, variants_nimh_lof) %>%
  count(symbol, name = "freq") %>%
  arrange(desc(freq))
top_genes <- gene_counts %>% slice_max(freq, n = 40)
ggplot(top_genes, aes(x = reorder(symbol, freq), y = freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "counts for top 40 genes mutated in G2MH cohort",
    x = "Gene",
    y = "Variant Count"
  ) +
  theme_minimal(base_size = 14)

# add 'has lof' and 'has missense' marker to genomatrix
phenotypes <- phenotypes %>%
  mutate(has_lof_mutation = ifelse(guid %in% variants_nimh_lof$guid, TRUE, FALSE)) %>%
  mutate(has_missense_mutation = ifelse(guid %in% variants_nimh_missense$guid, TRUE, FALSE))

# train models on IQ ~ lof mutation
results_list <- list()
genomatrix1 <- phenotypes %>%
  drop_na(smry_iq_fsiq, has_missense_mutation, age_final, age2, dem_sex) %>%
  filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5'))
for (cnv_class in c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5') ) {
  print(cnv_class)
  sample_subset <- genomatrix1 %>% filter(group_class_breakpoints == cnv_class)
  sample_subset <- sample_subset[complete.cases(sample_subset[, c("smry_iq_fsiq", "has_missense_mutation", "age_final", "age2", "dem_sex", "redcap_data_access_group")]), ]
  model <- lm(smry_iq_fsiq ~ has_missense_mutation + age_final + age2 + dem_sex + redcap_data_access_group, data = sample_subset)
  tidy_fit <- broom::tidy(model)
  eff <- tidy_fit[tidy_fit$term == "has_missense_mutationTRUE", ]
  results_list[[cnv_class]] <- data.frame(
    group_class_breakpoints = cnv_class,
    estimate = eff$estimate,
    p.value = eff$p.value,
    direction = ifelse(eff$estimate > 0, "Higher IQ", "Lower IQ")
  )
}
results <- do.call(rbind, results_list)
df_plot <- genomatrix1 %>%
  left_join(results %>% select(group_class_breakpoints, p.value, direction), by = 'group_class_breakpoints')
ggplot(df_plot, aes(x = has_missense_mutation, y = smry_iq_fsiq, fill = has_missense_mutation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6,
               position = position_dodge(width = 0.8)) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~ group_class_breakpoints, nrow = 1) +  # one panel per CNV group
  theme_minimal() +
  labs(
    x = "has secondary CNV",
    y = "FSIQ scores",
    title = "FSIQ scores across G2MH participants with and without missense mutation",
    fill = "Secondary Mutation"
  ) +
  geom_text(
    data = results,
    aes(x = 1.5, y = max(genomatrix1$smry_iq_fsiq)*0.95,
        label = ifelse(is.na(p.value), "NA", paste0("p = ", signif(p.value, 2)))),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# use group class as a covariate
genomatrix1 <- phenotypes %>% filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5'))
sample_subset <- genomatrix1[complete.cases(genomatrix1[, c("smry_iq_fsiq", "group_class_breakpoints", "has_missense_mutation", "age_final", "age2", "dem_sex", "redcap_data_access_group")]), ]
genomatrix1$group_class_breakpoints <- as.factor(genomatrix1$group_class_breakpoints)
genomatrix1$group_class_breakpoints <- relevel(genomatrix1$group_class_breakpoints, ref='22q11.2_Deletion A-D')
genomatrix1$group_class_breakpoints <- droplevels(genomatrix1$group_class_breakpoints)
model <- lm(smry_iq_fsiq ~ has_missense_mutation + group_class_breakpoints + age_final + age2 + dem_sex + redcap_data_access_group, data = sample_subset)
pval <- tidy(model) %>%
  filter(term == "has_missense_mutationTRUE") %>%
  pull(p.value)
ggplot(genomatrix1, aes(x = has_missense_mutation, y = smry_iq_fsiq, fill = has_missense_mutation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  theme_minimal() +
  labs(
    x = "Has Secondary Mutation",
    y = "FSIQ",
    title = "FSIQ scores across G2MH participants with and without missense mutation"
  ) +
  geom_text(
    aes(x = 1.5, y = 100,
        label = paste0("p = ", signif(pval, 2))),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme(legend.position = "none")


