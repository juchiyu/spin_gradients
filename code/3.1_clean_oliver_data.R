###########################################
## Exclusion criteria for Lindsay's data ##
###########################################
library(tidyverse)
lol_spins_behav <- read_csv('data/spins_behav_data_full_03-03-2022.csv')
lol_spins_behav_old <- read_csv('data/spins_lolivers_subject_info_for_grads_2022-02-02.csv')
## People in early termination have already been removed

spins_RS_fd <- read.csv(file="data/SPINS_fd_by_run_12-14-2020.csv",stringsAsFactors=F)
early_term <- read.csv(file="data/spins_termination_info_11-20-2020.csv", header=T, stringsAsFactors=F)

select_behav <- lol_spins_behav[,c("record_id","scanner","diagnostic_group","demo_sex","demo_age_study_entry",
                    "scog_rmet_total","scog_er40_total","scog_mean_ea","scog_tasit1_total","scog_tasit2_total","scog_tasit2_sinc",
                    "scog_tasit2_simpsar","scog_tasit2_parsar","scog_tasit3_total","scog_tasit3_lie","scog_tasit3_sar",
                    "np_composite_tscore","np_domain_tscore_process_speed","np_domain_tscore_att_vigilance","np_domain_tscore_work_mem",
                    "np_domain_tscore_verbal_learning","np_domain_tscore_visual_learning","np_domain_tscore_reasoning_ps",
                    "bsfs_sec1_total", "bsfs_sec2_total", "bsfs_sec3_total", "bsfs_sec4_total", "bsfs_sec5_total", 
                    "bsfs_sec6_total", "bsfs_sec7_y_total_7a", "bsfs_sec7_n_total_7b", 
                    "qls20_empathy", "qls_factor_interpersonal", "qls_factor_instrumental_role", "qls_factor_intrapsychic", "qls_factor_comm_obj_activities", 
                    "bprs_factor_neg_symp", "bprs_factor_pos_symp", "bprs_factor_anxiety_depression", "bprs_factor_activation", "bprs_factor_hostility")]

# exclusion list based on imaging (fmriprep and ciftify) QC (participants.csv)
exclude_img <- c("SPN01_CMH_0104","SPN01_CMH_0136","SPN01_MRP_0143","SPN01_ZHH_0040","SPN01_ZHH_0047",
                 "SPN01_ZHH_0048","SPN01_ZHH_0052","SPN01_ZHP_0063","SPN01_ZHP_0093","SPN01_ZHP_0105")
select_behav$exclude_MRI <- (select_behav$record_id %in% exclude_img)

# read in mean FD info and exclude participants with mean FD > 0.5 for rest
spins_RS_hi_fd <- na.omit(spins_RS_fd[spins_RS_fd$fd_mean.rest_bold>0.5, c("record_id","fd_mean.rest_bold")])
spins_RS_hi_fd$record_id # 6 lost to motion
# add mean FD to the data table
spins_RS <- spins_RS_fd[,c("record_id","fd_mean.rest_bold")]
rownames(spins_RS) <- spins_RS$record_id
select_behav$fd_mean_rest <- spins_RS[select_behav$record_id,"fd_mean.rest_bold"]
# exclusion based on mean fd
select_behav$exclude_meanFD <- (select_behav$record_id %in% spins_RS_hi_fd$record_id)

# read in  early termination IDs - check to see who remains
# early_term[early_term$record_id %in% names(RS_ts),"record_id"]
# 9 participants remain
# CMP_0178 exclude due to positive UDS (only determined after completion)
# CMP_0183 exclude due to no cog data - unable to contact
# CMP_0202 exclude due to no cog data - did not want to continue
# ZHH_0019 exclude due to positive UDS
# ZHH_0020 exclude due to high hba1c
# ZHH_0045 exclude no longer eligible after baseline screening
# ZHP_0100 exclude due to no cog data - patient lost to follow up
# ZHP_0103 exclude due to no cog data - unable to contact
# ZHP_0162 exclude due to no cog data - did not want to continue
early_term_exc <- c("SPN01_CMP_0178","SPN01_CMP_0183","SPN01_CMP_0202","SPN01_ZHH_0019","SPN01_ZHH_0020","SPN01_ZHH_0045",
                    "SPN01_ZHP_0100", "SPN01_ZHP_0103", "SPN01_ZHP_0162")
select_behav$exclude_earlyTerm <- (select_behav$record_id %in% early_term_exc)

## save to .csv file
write.csv(select_behav, file = "data/spins_lolivers_subject_info_for_grads_2022-04-21(withcomposite).csv")
