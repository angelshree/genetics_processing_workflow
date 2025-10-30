# main effects of CNV on each core psych variable
library(car)
library(ggplot2)
library(MASS)
library(ggpubr)
library(broom)
conflicts_prefer(dplyr::recode)

################################################################################
# setup
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
  mutate(smry_iq_age2 = smry_iq_age*smry_iq_age)

# set up CNV calls
wgs_cnv_calls <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/cnv_calls/wgs_cnv_calls.csv')
rownames(wgs_cnv_calls) = as.character(wgs_cnv_calls$GUID)
loci <- grep("^X[0-9]+", names(wgs_cnv_calls), value = TRUE)
loci <- loci[!grepl("oneway$", loci)]   
loci <- setdiff(loci, c("X22q11.21_AD", "X22q11.21_DH", "X16p11.2", "X15q35", "X15q13"))
loci <- gsub("^X", "locus.", loci)
colnames(wgs_cnv_calls) <- gsub("^X", "locus.", colnames(wgs_cnv_calls))

# combine with phenotypes
wgs_cnv_calls <- wgs_cnv_calls[as.character(phenotypes$guid),c('identifier',loci)]
wgs_cnv_calls <- wgs_cnv_calls %>% filter(!is.na(identifier))
phenotypes <- phenotypes %>% filter(guid %in% rownames(wgs_cnv_calls))
genomatrix <- cbind(phenotypes, wgs_cnv_calls)

# make sure sex is a factor with male as ref
genomatrix$dem_sex <- relevel(as.factor(genomatrix$dem_sex), ref='1')

# getting rid of weird vals
genomatrix[genomatrix == "NI"] <- NA
genomatrix[genomatrix == "M"] <- NA
genomatrix[genomatrix == "UNK"] <- NA

# calculate who has a secondary CNV
genomatrix <- genomatrix %>%
  select(-all_of(grep('locus.22q11.2', names(genomatrix))), -all_of(grep('locus.16p11.2', names(genomatrix)))) %>%
  mutate(num_secondary_cnvs = rowSums(across(starts_with("locus."), ~ .x != 2, .names = NULL), na.rm = TRUE)) %>%
  mutate(has_secondary_cnv = (num_secondary_cnvs > 0)) %>%
  filter(!(group_class_breakpoints == 'no_cnv')) %>%
  group_by(group_class_breakpoints) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  mutate(age_final = ifelse(!is.na(smry_iq_age), smry_iq_age, dem_age)) %>%
  mutate(age2 = age_final**2)
genomatrix$has_secondary_cnv <- as.factor(genomatrix$has_secondary_cnv)

# raw IQ dists across group class
plot <- genomatrix %>% filter((group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5')))
ggplot(plot, aes(x = smry_iq_fsiq, color = group_class_breakpoints, fill = group_class_breakpoints)) +
  geom_density(alpha = 0.3) +   # alpha makes the fill transparent
  theme_minimal() +
  labs(title = "Density Plot of FSIQ across groups",
       x = "Value", y = "Density")
ggplot(plot, aes(x = age_final, color = group_class_breakpoints, fill = group_class_breakpoints)) +
  geom_density(alpha = 0.3) +   # alpha makes the fill transparent
  theme_minimal() +
  labs(title = "Density Plot of age across groups",
       x = "Value", y = "Density")
plot(density((genomatrix %>% filter(group_class_breakpoints == '22q11.2_Deletion A-D'))$smry_iq_fsiq, na.rm = TRUE), main = "IQ Density across 22q11.2_Deletion A-D")
plot(density((genomatrix %>% filter(group_class_breakpoints == '22q11.2_Duplication A-D'))$smry_iq_fsiq, na.rm = TRUE), main = "IQ Density  across 22q11.2_Duplication A-D")
plot(density((genomatrix %>% filter(group_class_breakpoints == '16p11.2_Deletion BP4-5'))$smry_iq_fsiq, na.rm = TRUE), main = "IQ Density across 16p11.2_Deletion BP4-5")
plot(density((genomatrix %>% filter(group_class_breakpoints == '16p11.2_Duplication BP4-5'))$smry_iq_fsiq, na.rm = TRUE), main = "IQ Density across 16p11.2_Duplication BP4-5")

# raw IQ dists across group class
plot(density((genomatrix %>% filter(group_class_breakpoints == '22q11.2_Deletion A-D'))$age_final, na.rm = TRUE), main = "age Density across 22q11.2_Deletion A-D")
plot(density((genomatrix %>% filter(group_class_breakpoints == '22q11.2_Duplication A-D'))$age_final, na.rm = TRUE), main = "age Density  across 22q11.2_Duplication A-D")
plot(density((genomatrix %>% filter(group_class_breakpoints == '16p11.2_Deletion BP4-5'))$age_final, na.rm = TRUE), main = "age Density across 16p11.2_Deletion BP4-5")
plot(density((genomatrix %>% filter(group_class_breakpoints == '16p11.2_Duplication BP4-5'))$age_final, na.rm = TRUE), main = "age Density across 16p11.2_Duplication BP4-5")
  
# raw IQ distribution across group classes
plot1 <- ggplot(genomatrix, aes(x = group_class_breakpoints, y = smry_iq_fsiq, fill = group_class_breakpoints)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Class / ID", y = "IQ", title = "IQ Distribution by 16p/22q CNV") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
ggsave('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/main_effects/iq_distributions.png', plot=plot1, width = 8, height = 6)

# train fsiq models
results_list <- list()
genomatrix1 <- genomatrix %>%
  drop_na(smry_iq_fsiq, has_secondary_cnv, age_final, age2, dem_sex) %>%
  filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5'))
for (cnv_class in c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5') ) {
  print(cnv_class)
  sample_subset <- genomatrix1 %>% filter(group_class_breakpoints == cnv_class)
  sample_subset <- sample_subset[complete.cases(sample_subset[, c("smry_iq_fsiq", "has_secondary_cnv", "age_final", "age2", "dem_sex", "redcap_data_access_group")]), ]
  model <- lm(smry_iq_fsiq ~ has_secondary_cnv + age_final + age2 + dem_sex + redcap_data_access_group, data = sample_subset)
  tidy_fit <- broom::tidy(model)
  eff <- tidy_fit[tidy_fit$term == "has_secondary_cnvTRUE", ]
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

ggplot(df_plot, aes(x = has_secondary_cnv, y = smry_iq_fsiq, fill = has_secondary_cnv)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6,
               position = position_dodge(width = 0.8)) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~ group_class_breakpoints, nrow = 1) +  # one panel per CNV group
  theme_minimal() +
  labs(
    x = "has secondary CNV",
    y = "FSIQ scores",
    title = "FSIQ scores across G2MH participants with and without secondary CNV",
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

# train viq models
results_list <- list()
genomatrix1 <- genomatrix %>%
  drop_na(smry_iq_viq, has_secondary_cnv, age_final, age2, dem_sex) %>%
  filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5'))
for (cnv_class in c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5') ) {
  print(cnv_class)
  sample_subset <- genomatrix1 %>% filter(group_class_breakpoints == cnv_class)
  sample_subset <- sample_subset[complete.cases(sample_subset[, c("smry_iq_viq", "has_secondary_cnv", "age_final", "age2", "dem_sex", 'redcap_data_access_group')]), ]
  model <- lm(smry_iq_viq ~ has_secondary_cnv + age_final + age2 + dem_sex + redcap_data_access_group, data = sample_subset)
  tidy_fit <- broom::tidy(model)
  eff <- tidy_fit[tidy_fit$term == "has_secondary_cnvTRUE", ]
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
ggplot(df_plot, aes(x = has_secondary_cnv, y = smry_iq_viq, fill = has_secondary_cnv)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6,
               position = position_dodge(width = 0.8)) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~ group_class_breakpoints, nrow = 1) +  # one panel per CNV group
  theme_minimal() +
  labs(
    x = "has secondary CNV",
    y = "VIQ scores",
    title = "VIQ scores across G2MH participants with and without secondary CNV",
    fill = "Secondary Mutation"
  ) +
  geom_text(
    data = results,
    aes(x = 1.5, y = max(genomatrix1$smry_iq_viq)*0.95,
        label = ifelse(is.na(p.value), "NA", paste0("p = ", signif(p.value, 2)))),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# train piq models
results_list <- list()
genomatrix1 <- genomatrix %>%
  drop_na(smry_iq_piq, has_secondary_cnv, age_final, age2, dem_sex) %>%
  filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5'))
for (cnv_class in c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5') ) {
  print(cnv_class)
  sample_subset <- genomatrix1 %>% filter(group_class_breakpoints == cnv_class)
  sample_subset <- sample_subset[complete.cases(sample_subset[, c("smry_iq_piq", "has_secondary_cnv", "age_final", "age2", "dem_sex", "redcap_data_access_group")]), ]
  model <- lm(smry_iq_piq ~ has_secondary_cnv + age_final + age2 + dem_sex + redcap_data_access_group, data = sample_subset)
  tidy_fit <- broom::tidy(model)
  eff <- tidy_fit[tidy_fit$term == "has_secondary_cnvTRUE", ]
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
ggplot(df_plot, aes(x = has_secondary_cnv, y = smry_iq_piq, fill = has_secondary_cnv)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6,
               position = position_dodge(width = 0.8)) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~ group_class_breakpoints, nrow = 1) +  # one panel per CNV group
  theme_minimal() +
  labs(
    x = "has secondary CNV",
    y = "PIQ scores",
    title = "PIQ scores across G2MH participants with and without secondary CNV",
    fill = "Secondary Mutation"
  ) +
  geom_text(
    data = results,
    aes(x = 1.5, y = max(genomatrix1$smry_iq_piq)*0.95,
        label = ifelse(is.na(p.value), "NA", paste0("p = ", signif(p.value, 2)))),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# now, test has_secondary_cnv with all the group classes, and control 
genomatrix1 <- genomatrix %>% filter(group_class_breakpoints %in% c('22q11.2_Deletion A-D', '22q11.2_Duplication A-D', '16p11.2_Deletion BP4-5'))
sample_subset <- genomatrix1[complete.cases(genomatrix1[, c("smry_iq_fsiq", "group_class_breakpoints", "has_secondary_cnv", "age_final", "age2", "dem_sex", "redcap_data_access_group")]), ]
genomatrix1$group_class_breakpoints <- as.factor(genomatrix1$group_class_breakpoints)
genomatrix1$group_class_breakpoints <- relevel(genomatrix1$group_class_breakpoints, ref='22q11.2_Deletion A-D')
genomatrix1$group_class_breakpoints <- droplevels(genomatrix1$group_class_breakpoints)
model <- lm(smry_iq_fsiq ~ has_secondary_cnv + group_class_breakpoints + age_final + age2 + dem_sex + redcap_data_access_group, data = sample_subset)
pval <- tidy(model) %>%
  filter(term == "has_secondary_cnvTRUE") %>%
  pull(p.value)
ggplot(genomatrix1, aes(x = has_secondary_cnv, y = smry_iq_fsiq, fill = has_secondary_cnv)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  theme_minimal() +
  labs(
    x = "Has Secondary Mutation",
    y = "FSIQ",
    title = "FSIQ scores across G2MH participants with and without secondary CNV"
  ) +
  geom_text(
    aes(x = 1.5, y = 100,
        label = paste0("p = ", signif(pval, 2))),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme(legend.position = "none")


# adding ancestry PCs for whoever we have it
genetics <- read.csv('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/data/genomatrix_20250214.csv') %>%
  select(GUID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
genomatrix_genetics_combined <- merge(genomatrix, genetics, by.x = 'guid', by.y = 'GUID', all.y=TRUE)








corepsych_vars <- c('corepsych_dep_dx', 'corepsych_bpd_dx', 'corepsych_psy_dx',
                    'corepsych_anx_dx', 'corepsych_ocd_dx', 'corepsych_adhd_dx',
                    'corepsych_odd_dx', 'corepsych_asd_dx')
symptom_vars <- c('corepsych_psy_symptoms___1', 'corepsych_psy_symptoms___2', 
                  'corepsych_adhd_symptoms___1', 'corepsych_adhd_symptoms___2',
                  'corepsych_odd_symptoms___1', 'corepsych_odd_symptoms___2',
                  'corepsych_asd_core_1', 'corepsych_asd_core_2', 
                  'corepsych_ld___0', 'corepsych_ld___1', 'corepsych_ld___2', 'corepsych_id', 'corepsych_sub_dx') # note that sub_dx is effectively binary...

for (corepsych_var in corepsych_vars) {
  
  # recode corepsych variable
  x <- genomatrix[[corepsych_var]]
  x[x == 1] <- 99   # temporary code
  x[x == 2] <- 1
  x[x == 99] <- 2
  genomatrix[[corepsych_var]] <- x
  
  genomatrix[[corepsych_var]] <- as.ordered(genomatrix[[corepsych_var]])
}
################################################################################
# group-wise associations
################################################################################
# get CNV list

# filter genomatrix to only inlcude groups with at least 8 people
genomatrix_filtered <- genomatrix %>%
  filter(!(group_class_breakpoints == 'no_cnv')) %>%
  group_by(group_class_breakpoints) %>%
  filter(n() >= 5) %>%
  ungroup()
CNVs = unique(genomatrix_filtered$group_class_breakpoints)
genomatrix_filtered[genomatrix_filtered == "NI"] <- '0'
# set reference to biggest group (22qA-D Del)
genomatrix_filtered$group_class_breakpoints <- relevel(factor(genomatrix_filtered$group_class_breakpoints),
                             ref = names(sort(table(genomatrix_filtered$group_class_breakpoints), decreasing = TRUE))[1])
genomatrix_filtered$group_class_breakpoints <- droplevels(genomatrix_filtered$group_class_breakpoints)

results_global <- data.frame()
# recode corepsych vars and run global model
for (corepsych_var in corepsych_vars) {
  print(corepsych_var)
  # setup model df
  genomatrix_mod <- genomatrix_filtered[, c(corepsych_var, "group_class_breakpoints", "dem_age", "dem_age2", "dem_sex")]
  genomatrix_mod <- genomatrix_mod[complete.cases(genomatrix_mod), ]
  genomatrix_mod$group_class_breakpoints <- droplevels(genomatrix_mod$group_class_breakpoints)
  no_variation = (length(unique(genomatrix_mod[[corepsych_var]])) < 2)
  
  # threshold age if the variable is corepsych_psy_dx
  if (corepsych_var == 'corepsych_psy_dx') {
    genomatrix_mod <- genomatrix_mod %>% filter(dem_age > 16)
  }
  
  # # make model
  model <- polr(as.formula(paste(corepsych_var, "~ group_class_breakpoints + dem_age + dem_age2 + dem_sex")), data = genomatrix_mod, Hess = TRUE)
  coef_summary <- summary(model)$coefficients
  coef_names <- rownames(coef_summary)
  for (grp in coef_names) {
    grp_name <- gsub("^group_class_breakpoints", "", grp)
    sub_grp <- genomatrix_mod[genomatrix_mod$group_class_breakpoints == grp_name, ]
    no_variation_flag <- length(unique(sub_grp[[corepsych_var]])) == 1
    if (grp %in% rownames(coef_summary)) {
      est <- coef_summary[grp, "Value"]
      OR <- exp(est)
      pval <- 2 * pnorm(abs(coef_summary[grp, "t value"]), lower.tail = FALSE)
    } else {
      est <- NA
      OR <- NA
      pval <- NA
    }

    # add to results
    results_global <- rbind(results_global, data.frame(
      phenotype = corepsych_var,
      group = grp_name,
      estimate = est,
      OR = OR,
      p_value = pval,
      flag_no_variation = no_variation_flag,
      stringsAsFactors = FALSE
    ))
  }
}

################################################################################
# pairwise associations
################################################################################

# main effects separately - tbd if we have enough power
results_pairwise_main_effects <- data.frame()
for (phen in corepsych_vars) {
  
  for (grp in CNVs) {
    
    # Create a binary indicator for the current group vs all others and filter for complete cases
    print(grp)
    genomatrix_model <- genomatrix_filtered[, c(phen, "group_class_breakpoints", "dem_age", 'dem_age2', "dem_sex")]
    genomatrix_model <- genomatrix_model[complete.cases(genomatrix_model), ]
    genomatrix_model$in_group <- ifelse(genomatrix_model$group_class_breakpoints == grp, 1, 0)
    
    # checking for sparsity
    counts <- table(genomatrix_model$in_group, genomatrix_model[[phen]])
    min_count <- min(counts)
    n_group <- sum(genomatrix_model$in_group == 1)
    n_control <- sum(genomatrix_model$in_group == 0)
    sparse_flag <- min_count < 1
    phen_vals_in_group <- unique(genomatrix_model[[phen]][genomatrix_model$in_group == 1])
    no_variation_flag <- length(phen_vals_in_group) == 1
    
    # Run ordinal regression
    model <- polr(as.formula(paste(phen, "~ in_group + dem_age + dem_age2 + dem_sex")), data = genomatrix_model, Hess = TRUE)
    
    # Extract summary stats
    coef_summary <- summary(model)$coefficients
    pval <- pnorm(abs(coef_summary[,"t value"]), lower.tail = FALSE) * 2
    results_pairwise_main_effects <- rbind(results_pairwise_main_effects, 
                     data.frame(
                       phenotype = phen,
                       group_class = grp,
                       estimate = coef_summary["in_group", "Value"],
                       p_value = pval["in_group"],
                       n_group = n_group,
                       n_control = n_control,
                       flag_sparse = sparse_flag,
                       no_variation = no_variation_flag
                     ))
  }
  #break
}



################################################################################
# plotting
################################################################################

# function to calibrate significance levels
p_to_asterisks <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return(NA)  # not significant
}

# global associations
plot_df <- results_global %>%
  mutate(
    estimate_plot = ifelse(flag_no_variation | is.na(estimate), NA, estimate),
    sig_label = ifelse(!flag_no_variation & !is.na(p_value) & p_value < 0.05,
                       sapply(p_value, p_to_asterisks),
                       NA)   # NA will be blank in geom_text
  )
ggplot(plot_df, aes(x = group, y = phenotype, fill = estimate_plot)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90",   # flagged / blank cells
    name = "Estimate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Group",
    y = "Phenotype",
    title = "Effect of Group on Each Phenotype",
    subtitle = "Blank cells = no variation within group or NA estimates"
  )

# pairwise associations
plot_df <- results_pairwise_main_effects %>%
  mutate(log_p = -log10(p_value)) %>%
  mutate(estimate_plot = ifelse(no_variation, NA, estimate)) %>%
  mutate(sig_label = ifelse(!no_variation & p_value < 0.05,
                            sapply(p_value, p_to_asterisks), NA))
plot_df$phenotype <- factor(plot_df$phenotype, levels = unique(plot_df$phenotype))
plot_df$group <- factor(plot_df$group, levels = unique(plot_df$group))

ggplot(plot_df, aes(x = group_class, y = phenotype, fill = estimate_plot)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90",    # grey for missing or flagged cells
    name = "Estimate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Group",
    y = "Phenotype",
    title = "Effect of Group on Each Phenotype",
    subtitle = "Blank cells = no variation within group"
  )

################################################################################
# secondary CNV effects (using CNV allle counts and controlling for group_class_breakpoints)
################################################################################

secondary_group_cols <- c("locus.7q11.23_WBS", "locus.2q13_NPHP1", "locus.2q13", "locus.2q21.1",
                          "locus.7q11.23_distal", "locus.10q11.21_11.23", "locus.13q12.12",
                          "locus.17p12", "locus.17q12", "locus.1q21.1_TAR_BP1_BP2",
                          "locus.1q21.1_TAR_BP2_BP3", "locus.1q21.1_TAR_BP1_BP3", "locus.15q11.2_BP1_BP2",
                          "locus.15q11.2.13.1_BP2_BP3", "locus.15q11.2.13.1_BP1_BP3", 
                          "locus.15q13.1.13.2_BP3_BP4", "locus.15q13.1.13.3_BP3_BP5", 
                          "locus.15q13.3_BP4_BP5", "locus.15q13.3_BP4_BP5_exclCHRNA7",
                          "locus.15q13.3_CHRNA7", "locus.16p12.2_BP1_BP2", "locus.16p12.2_BP2_BP3",
                          "locus.16p12.2_BP1_BP3", "locus.16p13.11_BP1_BP2",
                          "locus.16p13.11_BP2_BP3", "locus.16p13.11_BP1_BP3")

# run associations for all secondary group classes while controlling for group class breakpoints
results_secondary <- data.frame()

# loop through phenotypes and secondary variables
for (phen in corepsych_vars) {
  #print(phen)
  for (col in secondary_group_cols) {
    #print(col)
    
    genomatrix_secondary <- genomatrix_filtered[, c(phen, col, 'group_class_breakpoints', "dem_age", 'dem_age2', "dem_sex")]
    genomatrix_secondary <- genomatrix_secondary[complete.cases(genomatrix_secondary), ]
    
    # Skip if phenotype has no variation at all
    if(length(unique(genomatrix_secondary[[phen]])) < 2) next
    
    # model for DEL
    genomatrix_secondary$in_group_DEL <- ifelse(genomatrix_secondary[[col]] == 1, 1, 0)
    no_variation_DEL <- length(unique(genomatrix_secondary[[phen]][genomatrix_secondary$in_group_DEL == 1])) < 2
    n_group <- sum(genomatrix_secondary$in_group_DEL)
    if (!no_variation_DEL) { # & !(sum(genomatrix_secondary$in_group_DEL) <= 2)) {
      #print(paste0('phenotype: ', phen))
      #print(paste0('cnv: ', col))
      #print('')
      model_DEL <- polr(as.formula(paste(phen, "~ in_group_DEL + group_class_breakpoints + dem_age + dem_age2 + dem_sex")),
             data = genomatrix_secondary, Hess = TRUE)
      coef_SUM <- summary(model_DEL)$coefficients
      if ("in_group_DEL" %in% rownames(coef_SUM)) {
        pval <- 2 * pnorm(abs(coef_SUM["in_group_DEL","t value"]), lower.tail = FALSE)
        est <- coef_SUM["in_group_DEL","Value"]
       } else {
        pval <- NA
        est <- NA
        }
      } else {
        pval <- NA
        est <- NA
      }
    results_secondary <- rbind(results_secondary, data.frame(
      phenotype = phen,
      group = paste(col, "DEL", sep="_"),
      estimate = est,
      p_value = pval,
      flag_no_variation = no_variation_DEL,
      n = n_group
    ))
    
    # model for DUP
    genomatrix_secondary$in_group_DUP <- ifelse(genomatrix_secondary[[col]] == 3, 1, 0)
    no_variation_DUP <- length(unique(genomatrix_secondary[[phen]][genomatrix_secondary$in_group_DUP == 1])) < 2
    #print(paste0('num dup: ', sum(genomatrix_secondary$in_group_DUP)))
    n_group <- sum(genomatrix_secondary$in_group_DUP)
    if (!no_variation_DUP) { # & !(sum(genomatrix_secondary$in_group_DUP) <= 2)) {
      model_DUP <- polr(as.formula(paste(phen, "~ in_group_DUP + group_class_breakpoints + dem_age + dem_age2 + dem_sex")),
             data = genomatrix_secondary, Hess = TRUE)
      coef_SUM <- summary(model_DUP)$coefficients
      if ("in_group_DUP" %in% rownames(coef_SUM)) {
        pval <- 2 * pnorm(abs(coef_SUM["in_group_DUP","t value"]), lower.tail = FALSE)
        est <- coef_SUM["in_group_DUP","Value"]
      } else {
        pval <- NA
        est <- NA
      } 
      
      results_secondary <- rbind(results_secondary, data.frame(
      phenotype = phen,
      group = paste(col, "DUP", sep="_"),
      estimate = est,
      p_value = pval,
      flag_no_variation = no_variation_DUP,
      n = n_group
    ))
    }
  }
}

plot_df <- results_secondary %>%
  filter(!is.na(p_value)) %>%
  mutate(log_p = -log10(p_value)) %>%
  mutate(estimate_plot = ifelse(flag_no_variation, NA, estimate)) %>%
  mutate(sig_label = ifelse(!flag_no_variation & p_value < 0.05,
                            sapply(p_value, p_to_asterisks), NA))
plot_df$phenotype <- factor(plot_df$phenotype, levels = unique(plot_df$phenotype))
plot_df$group <- factor(plot_df$group, levels = unique(plot_df$group))
plot_df$group <- droplevels(plot_df$group)
ggplot(plot_df, aes(x = group, y = phenotype, fill = estimate_plot)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90",    # grey for missing or flagged cells
    name = "Estimate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Group",
    y = "Phenotype",
    title = "Effect of Group on Each Phenotype",
    subtitle = "Blank cells = no variation within group"
  )



