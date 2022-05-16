# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SummarizedExperiment")
# # library(devtools)
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
spins_behav_combat <- read.csv("/projects/loliver/SPINS_PLS_Conn/data/processed/spins_behav_data_full_03-03-2022.csv")
spins_behav_combat_match <- spins_behav_combat[spins_behav_combat$record_id %in% colnames(rs_corrs_z),]
# dat is a data matrix of the data to harmonize - rows are features (connections) and columns are participants
rs_corrs_com <- as.matrix(rs_corrs_z[,colnames(rs_corrs_z) %in% spins_behav_combat$record_id])
# mod is a design matrix specifying biological covariates that should be protected - here diagnosis, age, sex, and cog variables
modcombat <- model.matrix(~ diagnostic_group + demo_sex + demo_age_study_entry +
                            scog_rmet_total + scog_er40_total + scog_mean_ea + scog_tasit1_total + scog_tasit2_sinc +
                            scog_tasit2_simpsar + scog_tasit2_parsar + scog_tasit3_lie + scog_tasit3_sar +
                            np_domain_tscore_process_speed + np_domain_tscore_att_vigilance + np_domain_tscore_work_mem +
                            np_domain_tscore_verbal_learning + np_domain_tscore_visual_learning + np_domain_tscore_reasoning_ps, 
                          data=spins_behav_combat_match)
# R run ComBat
# batch is a vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for
rs_combat <- neuroCombat(dat=rs_corrs_com, batch=c(spins_behav_combat_match$scanner), mod=modcombat)
# transpose the harmonized data matrix
rs_combat_data <- t(rs_combat$dat.combat)
##-----------------------------------------------------------------