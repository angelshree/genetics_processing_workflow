
setwd('~/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByAnjali/analysis/post_data_freeze/genetics_processing/cnv_genotyping/cnv_calls/')

gsa_cnv_calls <- read.csv('gsa_cnv_calls.csv')

test <- gsa_cnv_calls %>% 
  #filter(!is.na(gen_group_class) | group_class != '') %>%
  select(identifier, gen_group_class, gen_group_22q_type, gen_group_16p_type, group_class_breakpoints,
         all_of(grep('16p11.2', names(gsa_cnv_calls))), 
         all_of(grep('22q11.2', names(gsa_cnv_calls))))

test <- test %>% select(-all_of(grep('oneway', names(test))))
