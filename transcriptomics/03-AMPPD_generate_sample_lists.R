library(tidyverse)

########## 1218

case.ctrl.1218 = read_csv("./metadata/2020_v2release_1218/amppd_case_control_1218.csv")
rna.samples.1218 = read_csv("./metadata/2020_v2release_1218/rna_sample_inventory_1218.csv")

rna.case.control.1218 = dplyr::inner_join(x=case.ctrl.1218, 
                                          y=rna.samples.1218, 
                                          by="participant_id") %>% 
  dplyr::mutate(cohort=stringr::str_extract(participant_id, "^.{2}")) %>% 
  dplyr::filter(cohort %in% c("BF", "PP", "PD")) %>% 
  dplyr::filter(visit_month %in% c(0.5, 0)) %>% 
  dplyr::filter(case_control_other_latest!="Other")

rna.case.control.1218 %>% 
  dplyr::filter(cohort=="BF") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/BF-visit_month=0-case_control_other_latest!=Other_unique_1218.csv")

rna.case.control.1218 %>% 
  dplyr::filter(cohort=="PP") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PP-visit_month=0-case_control_other_latest!=Other_unique_1218.csv")

rna.case.control.1218 %>% 
  dplyr::filter(cohort=="PD") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PD-visit_month=0-case_control_other_latest!=Other_unique_1218.csv")

########## 1015

case.ctrl.1015 = read_csv("./metadata/2019_v1release_1015/amppd_case_control_1015.csv")
rna.samples.1015 = read_csv("./metadata/2019_v1release_1015/rna_sample_inventory_1015.csv")

rna.case.control.1015 = dplyr::inner_join(x=case.ctrl.1015, 
                                          y=rna.samples.1015, 
                                          by="participant_id") %>% 
  dplyr::mutate(cohort=stringr::str_extract(participant_id, "^.{2}")) %>% 
  dplyr::filter(cohort %in% c("BF", "PP", "PD")) %>% 
  dplyr::filter(visit_month %in% c(0.5, 0)) %>% 
  dplyr::filter(case_control_other_latest!="Other")

# --- 

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="BF") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/BF-visit_month=0-case_control_other_latest!=Other_unique_1015.csv")

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="BF") %>% 
  dplyr::filter(case_control_other_latest=="Case") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/BF-visit_month=0-case_control_other_latest==Case_unique_1015.csv")

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="BF") %>% 
  dplyr::filter(case_control_other_latest=="Control") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/BF-visit_month=0-case_control_other_latest==Control_unique_1015.csv")

# --- 

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PP") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PP-visit_month=0-case_control_other_latest!=Other_unique_1015.csv")

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PP") %>% 
  dplyr::filter(case_control_other_latest=="Case") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PP-visit_month=0-case_control_other_latest==Case_unique_1015.csv")

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PP") %>% 
  dplyr::filter(case_control_other_latest=="Control") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PP-visit_month=0-case_control_other_latest==Control_unique_1015.csv")

# --- 

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PD") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PD-visit_month=0-case_control_other_latest!=Other_unique_1015.csv")

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PD") %>% 
  dplyr::filter(case_control_other_latest=="Case") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PD-visit_month=0-case_control_other_latest==Case_unique_1015.csv")

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PD") %>% 
  dplyr::filter(case_control_other_latest=="Control") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PD-visit_month=0-case_control_other_latest==Control_unique_1015.csv")

########## for the 1015 all timepoints analysis

dplyr::inner_join(x=case.ctrl.1015, 
                  y=rna.samples.1015, 
                  by="participant_id") %>% 
  dplyr::mutate(cohort=stringr::str_extract(participant_id, "^.{2}")) %>% 
  dplyr::filter(cohort %in% c("PD")) %>% 
  dplyr::filter(case_control_other_latest=="Case") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PD-visit_month=ALL-case_control_other_latest==Case_unique_1015.csv")

dplyr::inner_join(x=case.ctrl.1015, 
                  y=rna.samples.1015, 
                  by="participant_id") %>% 
  dplyr::mutate(cohort=stringr::str_extract(participant_id, "^.{2}")) %>% 
  dplyr::filter(cohort %in% c("PD")) %>% 
  dplyr::filter(case_control_other_latest=="Control") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PD-visit_month=ALL-case_control_other_latest==Control_unique_1015.csv")

dplyr::inner_join(x=case.ctrl.1015, 
                  y=rna.samples.1015, 
                  by="participant_id") %>% 
  dplyr::mutate(cohort=stringr::str_extract(participant_id, "^.{2}")) %>% 
  dplyr::filter(cohort %in% c("PP")) %>% 
  dplyr::filter(case_control_other_latest=="Case") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PP-visit_month=ALL-case_control_other_latest==Case_unique_1015.csv")

dplyr::inner_join(x=case.ctrl.1015, 
                  y=rna.samples.1015, 
                  by="participant_id") %>% 
  dplyr::mutate(cohort=stringr::str_extract(participant_id, "^.{2}")) %>% 
  dplyr::filter(cohort %in% c("PP")) %>% 
  dplyr::filter(case_control_other_latest=="Control") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  write.csv(x=., "./data/AMPPD_sample_lists/PP-visit_month=ALL-case_control_other_latest==Control_unique_1015.csv")



########## check

rna.case.control.1218 %>% 
  dplyr::filter(cohort=="BF") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  length()

rna.case.control.1218 %>% 
  dplyr::filter(cohort=="PP") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  length()

rna.case.control.1218 %>% 
  dplyr::filter(cohort=="PD") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  length()

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="BF") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  length()

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PP") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  length()

rna.case.control.1015 %>% 
  dplyr::filter(cohort=="PD") %>% 
  dplyr::pull(sample_id) %>% 
  unique() %>% 
  length()

