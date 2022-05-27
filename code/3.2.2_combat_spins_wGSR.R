# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SummarizedExperiment")
# library(devtools)
# install_github("jfortin1/neuroCombatData")
# install_github("jfortin1/neuroCombat_Rpackage")

library(tidyverse)

# load correlation data
load("data/spins_RS_2mm_GSR_GlasserTian_zcor.rda")
rs_corrs_z <- rs_zcor %>% simplify2array %>% apply(3, function(x) x[upper.tri(x)])
##------------------------------From Lindsay, May-11-2022---------
# run ComBat on connectivity data to harmonize across scanners
library(neuroCombat)
# path to behavioural 
spins_behav_combat <- 
  read.csv("data/spins_lolivers_subject_info_for_grads_2022-04-21(withcomposite).csv") %>%
  filter(exclude_MRI==FALSE,
         exclude_meanFD==FALSE,
         exclude_earlyTerm==FALSE,
         exclude_nocog==FALSE,
         # exclude_noEA==FALSE,
         exclude_incmpltcog==FALSE)
spins_behav_combat_match <- spins_behav_combat[spins_behav_combat$record_id %in% colnames(rs_corrs_z),]
rownames(spins_behav_combat_match) <- spins_behav_combat_match$record_id
# dat is a data matrix of the data to harmonize - rows are features (connections) and columns are participants
rs_corrs_com <- as.matrix(rs_corrs_z[,colnames(rs_corrs_z) %in% spins_behav_combat$record_id])

## check complete cases --------------------------------------------------
data2impute <- spins_behav_combat_match %>% select(scanner, diagnostic_group, demo_sex, demo_age_study_entry,
                                               scog_rmet_total, scog_er40_total, #scog_mean_ea, 
                                               scog_tasit1_total, scog_tasit2_sinc, scog_tasit2_simpsar, 
                                               scog_tasit2_parsar, scog_tasit3_lie, scog_tasit3_sar,
                                               np_domain_tscore_process_speed, np_domain_tscore_att_vigilance, np_domain_tscore_work_mem,
                                               np_domain_tscore_verbal_learning, np_domain_tscore_visual_learning, np_domain_tscore_reasoning_ps)

## impute missing values
library(mice)
# impute missing values for combat variables
# explore rows with NAs - 9 people missing one beh value
data2impute[!complete.cases(data2impute),] %>% nrow
# impute remaining NA values using mice (1 each for 15 participants)
# ppm (predictive mean matching) is default method; m = 5 is the default number of imputations
# logged events are related to the categorical variables (not included in imputation)
spins_behav_combat_impt <- complete(mice(data2impute)) # n = 420
##-------------------------------------------------------------------------
# mod is a design matrix specifying biological covariates that should be protected - here diagnosis, age, sex, and cog variables
modcombat <- model.matrix(~ diagnostic_group + demo_sex + demo_age_study_entry +
                            scog_rmet_total + scog_er40_total + scog_tasit1_total + scog_tasit2_sinc +
                            scog_tasit2_simpsar + scog_tasit2_parsar + scog_tasit3_lie + scog_tasit3_sar +
                            np_domain_tscore_process_speed + np_domain_tscore_att_vigilance + np_domain_tscore_work_mem +
                            np_domain_tscore_verbal_learning + np_domain_tscore_visual_learning + np_domain_tscore_reasoning_ps, 
                          data=spins_behav_combat_impt)
# R run ComBat
# batch is a vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for
rs_combat <- neuroCombat(dat=rs_corrs_com, batch=c(spins_behav_combat_match$scanner), mod=modcombat)
# transpose the harmonized data matrix
rs_combat_data_z <- t(rs_combat$dat.combat)

## From z to r values again
source("code/functions/z2r.R")
rs_combat_data_r <- z2r(rs_combat_data_z)
## reformat the vectorized data as correlation matrices
source("code/functions/vec2sqmat.R")
rs_combat_cor_r <- list()
for (i in 1:nrow(rs_combat_data_r)){
  rs_combat_cor_r[[i]] <- vec2sqmat(rs_combat_data_r[i,]) 
  diag(rs_combat_cor_r[[i]]) <- 1
  write.table(rs_combat_cor_r[[i]], file = sprintf("data/spins_RS_2mm_GSR_GlasserTian_combat_cor/%s.txt", rownames(rs_combat_data_r)[i]))
}
names(rs_combat_cor_r) <- rownames(rs_combat_data_r)

## Imputed behavioral data for output
spins_behav_impt <- data.frame(record_id = rownames(spins_behav_combat_impt), spins_behav_combat_impt)
## save output
save(rs_combat_data_r, rs_combat_data_z, rs_combat_cor_r, spins_behav_impt, file = "data/spins_RS_2mm_GSR_GlasserTian_combated_imputed.RData")
##-----------------------------------------------------------------