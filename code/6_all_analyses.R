#___________________________________________
# Gather analysis/data for results plotting
#-------------------------------------------

# Libraries ----
library(tidyverse)
library(ggseg)
library(ggsegGlasser)
library(broom)
library(TExPosition)
library(PTCA4CATA)
library(plotly)
library(colorspace)
library(gridExtra)
library(MKinfer) # for bootstrap t
library(tableone)
library("RColorBrewer")

# Read data ----
## gradients
spins_grads <- read_csv("data/spins_gsr_RS_gradients.csv")

spins_grads_num <- data.frame(spins_grads[,c(2:5,7)])
spins_grads_num_full <- data.frame(spins_grads[,c(2:7)])
spins_grads_wide <- reshape(spins_grads_num, idvar = "Subject", timevar = "ROI", direction = "wide")

## subject data ----
load("data/spins_RS_2mm_GSR_GlasserTian_combated_imputed.RData")
rm(rs_combat_cor_r, rs_combat_data_r, rs_combat_data_z)

lol_spins_behav <- spins_behav_impt
names(lol_spins_behav)
lol_spins_behav$subject <- sub("SPN01_", "sub-", lol_spins_behav$record_id) %>% sub("_", "", .)

## add motion data from the original behavioral set ----
lol_original <- 
  read_csv('data/spins_lolivers_subject_info_for_grads_2022-05-30.csv') %>%
  filter(exclude_MRI==FALSE, 
         exclude_meanFD==FALSE, 
         exclude_earlyTerm==FALSE) %>% as.data.frame
lol_original$subject <- sub("SPN01_", "sub-", lol_original$record_id) %>% sub("_", "", .)
rownames(lol_original) <- lol_original$subject
lol_spins_behav$fd_mean_rest <- lol_original[lol_spins_behav$subject,"fd_mean_rest"]

## design matrix for subjects ----
spins_dx <- lol_spins_behav %>%
  select(subject,scanner,diagnostic_group,demo_sex,demo_age_study_entry)
spins_dx_org <- spins_dx[,-1] %>% data.frame

## demographics data ----
lol_demo <- 
  read_csv('data/spins_lolivers_subject_info_for_grads_2022-04-21(withcomposite).csv') %>%
  filter(exclude_MRI==FALSE, 
         exclude_meanFD==FALSE, 
         exclude_earlyTerm==FALSE) %>% as.data.frame
lol_demo$subject <- sub("SPN01_", "sub-", lol_demo$record_id) %>% sub("_", "", .)
rownames(lol_demo) <- lol_demo$record_id

## numeric data ----
spins_behav_num <- lol_spins_behav %>% 
  select(scog_rmet_total, scog_er40_total,
         scog_tasit1_total,
         scog_tasit2_parsar, scog_tasit2_simpsar, scog_tasit2_sinc,
         scog_tasit3_lie, scog_tasit3_sar, np_domain_tscore_att_vigilance,
         np_domain_tscore_process_speed, np_domain_tscore_work_mem,
         np_domain_tscore_verbal_learning, np_domain_tscore_visual_learning,
         np_domain_tscore_reasoning_ps, 
         fd_mean_rest
  ) %>% data.frame
rownames(spins_behav_num) <- lol_spins_behav$subject

## participants' variables of which the effects should be regressed out ----
var2regout <- lol_spins_behav %>%
  select(demo_sex, demo_age_study_entry, `fd_mean_rest`, scanner) %>% data.frame
rownames(var2regout) <- lol_spins_behav$record_id
var2regout$demo_sex_num <- as.numeric(as.factor(var2regout$demo_sex))-1
var2regout_num <- var2regout[,-1]

# Check subject overlap ----
grad.sub <- spins_grads_wide$Subject[order(spins_grads_wide$Subject)]
behav.sub <- lol_spins_behav$record_id[order(lol_spins_behav$record_id)]

kept.sub <- lol_spins_behav$record_id[complete.cases(lol_spins_behav)==TRUE] # 420
kept.sub <- kept.sub[kept.sub != "SPN01_ZHP_0150"] # non-eligible 

## grab the matching data ----

behav.dat <- lol_spins_behav[kept.sub,c(6:19)]
spins_grads_wide_org <- spins_grads_wide[,-1]
rownames(spins_grads_wide_org) <- spins_grads_wide$Subject
grad.dat <- spins_grads_wide_org[kept.sub,]
lol_demo_match <- lol_demo[kept.sub,]

## demo and numeric data ----
spins_demo <- lol_demo_match %>% 
  select(demo_sex, demo_age_study_entry, diagnostic_group, scog_rmet_total, scog_er40_total, #scog_mean_ea,
         scog_tasit1_total,
         scog_tasit2_sinc,
         scog_tasit2_simpsar,
         scog_tasit2_parsar,
         scog_tasit3_lie,
         scog_tasit3_sar, 
         np_domain_tscore_att_vigilance,
         np_domain_tscore_process_speed,
         np_domain_tscore_work_mem,
         np_domain_tscore_verbal_learning,
         np_domain_tscore_visual_learning,
         np_domain_tscore_reasoning_ps, 
  ) %>% data.frame
colnames(spins_demo)
rownames(spins_demo) <- lol_demo_match$subject

## variables to regress out ----
regout.dat <- var2regout_num[kept.sub,]

# Regress out the effects ----
behav.reg <- apply(behav.dat, 2, function(x) lm(x~regout.dat$demo_sex + regout.dat$demo_age_study_entry + regout.dat$fd_mean_rest)$residual)

grad.reg <- apply(grad.dat, 2, function(x) lm(x~regout.dat$demo_sex + regout.dat$demo_age_study_entry + regout.dat$fd_mean_rest)$residual)

grad.reg2plot <- apply(grad.dat, 2, function(x){
  model <- lm(x~regout.dat$demo_sex + regout.dat$demo_age_study_entry + regout.dat$fd_mean_rest)
  return(model$residual + model$coefficient[1])
} )

# Network colors ----
networks <- read_delim("networks.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(NETWORK, NETWORKKEY, RED, GREEN, BLUE, ALPHA) %>%
  distinct() %>%
  add_row(NETWORK = "Subcortical", NETWORKKEY = 13, RED = 0, GREEN=0, BLUE=0, ALPHA=255) %>%
  mutate(hex = rgb(RED, GREEN, BLUE, maxColorValue = 255)) %>%
  arrange(NETWORKKEY)

networks$hex <- darken(networks$hex, 0.2)

# Row and column designs ----
## match ROIs to networks ----
ROI.network.match <- cbind(spins_grads$ROI, spins_grads$Network) %>% unique
ROI.idx <- ROI.network.match[,2]
names(ROI.idx) <- ROI.network.match[,1]
### match networks with colors ----
net.col.idx <- networks$hex
names(net.col.idx) <- networks$NETWORK

## design matrix for subjects ----
sub.dx <- spins_dx_org[kept.sub,]

diagnostic.dx <- sub.dx$diagnostic_group %>% as.matrix
diagnostic.dx <- recode(diagnostic.dx, !!!c("case" = "SSD"))
diagnostic.col.idx <- c("SSD" = "darkorchid3",
                        "control" = "gray50")
diagnostic.col <- list()
diagnostic.col$oc <- recode(diagnostic.dx, !!!diagnostic.col.idx) %>% as.matrix()
diagnostic.col$gc <- diagnostic.col.idx %>% as.matrix

## design matrix for columns - behavioral ----
behav.dx <- matrix(nrow = ncol(behav.dat), ncol = 1, dimnames = list(colnames(behav.dat), "type")) %>% as.data.frame

behav.col <- c("scog" = "#CD661D",
               "np" = "#3F7538",#"#59A14F",
               "bsfs" = "#D37295")

behav.dx$type <- sub("(^[^_]+).*", "\\1", colnames(behav.dat))
behav.dx$type.col <- recode(behav.dx$type, !!!behav.col)

## design matrix for columns - gradient ----
grad.dx <- matrix(nrow = ncol(grad.dat), ncol = 4, dimnames = list(colnames(grad.dat), c("gradient", "ROI", "network", "network.col"))) %>% as.data.frame

grad.dx$gradient <- sub("(^[^.]+).*", "\\1", colnames(grad.dat))
grad.dx$ROI <- sub("^[^.]+.(*)", "\\1", colnames(grad.dat))
grad.dx$network <- recode(grad.dx$ROI, !!!ROI.idx)
grad.dx$network.col <- recode(grad.dx$network, !!!net.col.idx)

## get different alpha for gradients ----
grad.col.idx <- c("grad1" = "grey30",
                  "grad2" = "grey60",
                  "grad3" = "grey90")
grad.dx$gradient.col <- recode(grad.dx$gradient, !!!grad.col.idx)

# PLSC ----
pls.res <- tepPLS(behav.reg, grad.reg, DESIGN = sub.dx$diagnostic_group, make_design_nominal = TRUE, graphs = FALSE)

pls.boot <- data4PCCAR::Boot4PLSC(behav.reg, grad.reg, scale1 = "SS1", scale2 = "SS1", nIter = 1000, nf2keep = 5, eig = TRUE)

pls.boot$bootRatiosSignificant.j[abs(pls.boot$bootRatios.j) < 2.88] <- FALSE
pls.boot$bootRatiosSignificant.i[abs(pls.boot$bootRatios.i) < 2.88] <- FALSE

pls.inf <- data4PCCAR::perm4PLSC(behav.reg, grad.reg, scale1 = "SS1", scale2 = "SS1", nIter = 1000)

## swith direction for dimension 3 ----
pls.res$TExPosition.Data$fi[,1] <- pls.res$TExPosition.Data$fi[,1]*-1
pls.res$TExPosition.Data$fj[,1] <- pls.res$TExPosition.Data$fj[,1]*-1
pls.res$TExPosition.Data$pdq$p[,1] <- pls.res$TExPosition.Data$pdq$p[,1]*-1
pls.res$TExPosition.Data$pdq$q[,1] <- pls.res$TExPosition.Data$pdq$q[,1]*-1
pls.res$TExPosition.Data$lx[,1] <- pls.res$TExPosition.Data$lx[,1]*-1
pls.res$TExPosition.Data$ly[,1] <- pls.res$TExPosition.Data$ly[,1]*-1

## bootstrap for means ----
check.dim = 1
lxly <- cbind(pls.res$TExPosition.Data$lx[,check.dim], pls.res$TExPosition.Data$ly[,check.dim])
colnames(lxly) <- c(paste0("Dim", check.dim, c(".Behavioural", ".gradient")))

lxly.boot <- Boot4Mean(lxly, diagnostic.dx, niter = 1000)
colnames(lxly.boot$GroupMeans) <- colnames(lxly.boot$BootCube) <- c(paste0("Dim", check.dim, c(".Behavioural", ".gradient")))

## Parse q ----
q.grad1 <- pls.res$TExPosition.Data$fj[which(grad.dx$gradient == "grad1"),]
q.sig.grad1 <- pls.boot$bootRatiosSignificant.j[which(grad.dx$gradient == "grad1"),]
col.grad1 <- grad.dx$network.col[which(grad.dx$gradient == "grad1")]
dx.grad1 <- grad.dx$network[which(grad.dx$gradient == "grad1")]

q.grad2 <- pls.res$TExPosition.Data$fj[which(grad.dx$gradient == "grad2"),]
q.sig.grad2 <- pls.boot$bootRatiosSignificant.j[which(grad.dx$gradient == "grad2"),]
col.grad2 <- grad.dx$network.col[which(grad.dx$gradient == "grad2")]
dx.grad2 <- grad.dx$network[which(grad.dx$gradient == "grad2")]

q.grad3 <- pls.res$TExPosition.Data$fj[which(grad.dx$gradient == "grad3"),]
q.sig.grad3 <- pls.boot$bootRatiosSignificant.j[which(grad.dx$gradient == "grad3"),]
col.grad3 <- grad.dx$network.col[which(grad.dx$gradient == "grad3")]
dx.grad3 <- grad.dx$network[which(grad.dx$gradient == "grad3")]

# Data for brain plotting ----
## reorganize data frame ----
### factor scores ----
full.pq <- data.frame(grad.dx[,1:4], pls.res$TExPosition.Data$fj[,1:4], colMeans(grad.dat))
colnames(full.pq)[5:9] <- c(paste0("fj", c(1:4)), "raw")

### filtered by significance ----
full.pq.sig <- full.pq[,5:8]
full.pq.sig[pls.boot$bootRatiosSignificant.j[,1:4] == FALSE] <- 0
full.pq[,10:13] <- full.pq.sig
colnames(full.pq)[10:13] <- c(paste0("fj.sig", c(1:4)))

### contributions ----
full.cj <- data.frame(grad.dx[,1:4],
                      pls.res$TExPosition.Data$cj[,1:4],
                      colMeans(grad.dat))
colnames(full.cj)[5:9] <- c(paste0("cj", c(1:4)), "raw")
full.cj.sig <- full.cj[,5:8]
full.cj.sig[pls.boot$bootRatiosSignificant.j[,1:4] == FALSE] <- 0
full.cj[,10:13] <- full.cj.sig
colnames(full.cj)[10:13] <- c(paste0("cjsig", c(1:4)))

## filtered by importance (1/#J) ----
full.cj.imp <- full.cj[,5:8]
full.cj.imp[full.cj[,5:8] < 1/nrow(full.cj)] <- 0
full.cj[,14:17] <- full.cj.imp
colnames(full.cj)[14:17] <- c(paste0("cjimp", c(1:4)))

### sign ----
full.cj[,c("dir.1", "dir.2", "dir.3", "dir.4")] <- sign(pls.res$TExPosition.Data$fj[,1:4])

### full positive contributions ----
full.cj[,c("cj1_x_mean", "cj2_x_mean", "cj3_x_mean", "cj4_x_mean")] <- full.cj[,c("cj1", "cj2", "cj3", "cj4")]*full.cj$raw

full.cj[,c("cj1_x_mean.pos", "cj2_x_mean.pos", "cj3_x_mean.pos", "cj4_x_mean.pos")] <- full.cj[,c("cj1_x_mean.neg", "cj2_x_mean.neg", "cj3_x_mean.neg", "cj4_x_mean.neg")] <- full.cj[,c("cj1_x_mean", "cj2_x_mean", "cj3_x_mean", "cj4_x_mean")] 

### filtered by significance (bootstrap of fj) ----
full.cj[,c("cjsig1_x_mean", "cjsig2_x_mean", "cjsig3_x_mean", "cjsig4_x_mean")] <- full.cj[,c("cjsig1", "cjsig2", "cjsig3", "cjsig4")]*full.cj$raw

### filtered by importance (q^2) ----
full.cj[,c("cjimp1_x_mean", "cjimp2_x_mean", "cjimp3_x_mean", "cjimp4_x_mean")] <- full.cj[,c("cjimp1", "cjimp2", "cjimp3", "cjimp4")]*full.cj$raw

### separate by sign ----
full.cj <- full.cj %>% mutate(
  cjimp1.pos = if_else(dir.1 > 0,cjimp1,0),
  cjimp2.pos = if_else(dir.2 > 0,cjimp2,0),
  cjsig1.pos = if_else(dir.1 > 0,cjsig1,0),
  cjsig2.pos = if_else(dir.2 > 0,cjsig2,0),
  cj1_x_mean.pos = if_else(dir.1 > 0,cj1_x_mean,0),
  cj2_x_mean.pos = if_else(dir.2 > 0,cj2_x_mean,0),
  cj3_x_mean.pos = if_else(dir.3 > 0,cj3_x_mean,0),
  cj4_x_mean.pos = if_else(dir.4 > 0,cj4_x_mean,0),
  cjsig1_x_mean.pos = if_else(dir.1 > 0,cjsig1_x_mean,0),
  cjsig2_x_mean.pos = if_else(dir.2 > 0,cjsig2_x_mean,0),
  cjsig3_x_mean.pos = if_else(dir.3 > 0,cjsig3_x_mean,0),
  cjsig4_x_mean.pos = if_else(dir.4 > 0,cjsig4_x_mean,0),
  cjimp1_x_mean.pos = if_else(dir.1 > 0,cjimp1_x_mean,0),
  cjimp2_x_mean.pos = if_else(dir.2 > 0,cjimp2_x_mean,0),
  cjimp3_x_mean.pos = if_else(dir.3 > 0,cjimp3_x_mean,0),
  cjimp4_x_mean.pos = if_else(dir.4 > 0,cjimp4_x_mean,0),
  cjimp1.neg = if_else(dir.1 < 0,cjimp1,0),
  cjimp2.neg = if_else(dir.2 < 0,cjimp2,0),
  cjsig1.neg = if_else(dir.1 < 0,cjsig1,0),
  cjsig2.neg = if_else(dir.2 < 0,cjsig2,0),
  cj1_x_mean.neg = if_else(dir.1 < 0,cj1_x_mean,0),
  cj2_x_mean.neg = if_else(dir.2 < 0,cj2_x_mean,0),
  cj3_x_mean.neg = if_else(dir.3 < 0,cj3_x_mean,0),
  cj4_x_mean.neg = if_else(dir.4 < 0,cj4_x_mean,0),
  cjsig1_x_mean.neg = if_else(dir.1 < 0,cjsig1_x_mean,0),
  cjsig2_x_mean.neg = if_else(dir.2 < 0,cjsig2_x_mean,0),
  cjsig3_x_mean.neg = if_else(dir.3 < 0,cjsig3_x_mean,0),
  cjsig4_x_mean.neg = if_else(dir.4 < 0,cjsig4_x_mean,0),
  cjimp1_x_mean.neg = if_else(dir.1 < 0,cjimp1_x_mean,0),
  cjimp2_x_mean.neg = if_else(dir.2 < 0,cjimp2_x_mean,0),
  cjimp3_x_mean.neg = if_else(dir.3 < 0,cjimp3_x_mean,0),
  cjimp4_x_mean.neg = if_else(dir.4 < 0,cjimp4_x_mean,0))

### Organize data ----
pq_for_plot <- full.pq %>% 
  mutate(roi = str_remove(ROI, "_ROI"),
         label = case_when(str_starts(ROI, "L") ~ str_c("lh_", roi),
                           str_starts(ROI, "R_") ~ str_c("rh_", roi)))

pqxcj_for_plot <- full.cj %>% 
  mutate(roi = str_remove(ROI, "_ROI"),
         label = case_when(str_starts(ROI, "L") ~ str_c("lh_", roi),
                           str_starts(ROI, "R_") ~ str_c("rh_", roi)))


# Gradients: group differences ----
## Different formatting ----
grad.reg.long <- data.frame(t(grad.reg2plot), grad.dx[,1:3])
grad.reg.wide <- grad.reg.long %>% pivot_longer(cols = "SPN01_CMH_0001":"SPN01_ZHP_0171",
                                                names_to = "Subject",
                                                names_prefix = "",
                                                values_to = "value") 
grad.reg.wide4plot <- grad.reg.wide %>% pivot_wider(names_from = "gradient",
                                                    names_prefix = "",
                                                    values_from = "value",) %>% data.frame
grad.reg.wide4plot$Subject <- sub("sub.", "sub-", grad.reg.wide4plot$Subject)
names(grad.reg.wide4plot)[2] <- "Network"

## organize data - contributions ----
spins_cj <- list(grad1 = pqxcj_for_plot %>% filter(gradient == "grad1"),
                 grad2 = pqxcj_for_plot %>% filter(gradient == "grad2"),
                 grad3 = pqxcj_for_plot %>% filter(gradient == "grad3"))

spins_cj <- lapply(spins_cj, function(x) {
  rownames(x) <- spins_cj$grad1$ROI
  return(x)
})

## organize data - ROI average ----
spins_mean <- getMeans(grad.reg.wide4plot[,4:6], grad.reg.wide4plot$ROI)
grad2plot <- spins_mean
colnames(grad2plot) <- c(paste0("Gradient ", c(1:3)))
grad4grp <- unique(grad.reg.wide4plot[,c(1,2)]) %>% data.frame
rownames(grad4grp) <- grad4grp$ROI

## organize data - network mean ----
net_mean <- getMeans(spins_mean, grad4grp[rownames(spins_mean), "Network"])
colnames(net_mean) <- c(paste0("Gradient ", c(1:3)))

## organize data - ROI & group average ----
spins_grad_everything <- merge(grad.reg.wide4plot, sub.dx, by.x = "Subject", by.y = 0)

spins_mean_bygrp<- spins_grad_everything %>% 
  group_by(ROI, Network,diagnostic_group) %>% 
  dplyr::summarize(across(starts_with("grad"), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame %>%
  reshape(idvar = c("ROI", "Network"), timevar = "diagnostic_group", direction = "wide")

spins_sd_bygrp<- spins_grad_everything %>% 
  group_by(ROI, Network,diagnostic_group) %>% 
  dplyr::summarize(across(starts_with("grad"), ~ sd(.x, na.rm = TRUE))) %>%
  as.data.frame %>%
  reshape(idvar = c("ROI", "Network"), timevar = "diagnostic_group", direction = "wide")

grad2plot2 <- spins_mean_bygrp[,c(3:8)]
grad4grp2 <- unique(spins_mean_bygrp[,c(1:2)]) %>% as.data.frame
rownames(grad4grp2) <- grad4grp2$ROI
rownames(grad2plot2) <- grad4grp2$ROI

## organize data - network mean by groups ----
net_mean_bygrp <- spins_mean_bygrp %>%
  group_by(Network) %>%
  dplyr::summarize(across(starts_with("grad"), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame
rownames(net_mean_bygrp) <- net_mean_bygrp$Network

## get colors ----
roi.col <- list()
roi.col$oc <- recode(grad4grp$Network, !!!net.col.idx) %>% as.matrix
rownames(roi.col$oc) <- grad4grp$ROI
roi.col$gc <- as.matrix(net.col.idx)

## color for glasser atlas ----
net.col.idx.light <- lighten(net.col.idx)
names(net.col.idx.light) <- names(net.col.idx)

## data for arrows with fj ----
grad_fjxgrad <- data.frame(pls.res$TExPosition.Data$pdq$q[,1:3]*5, grad.dx[,c("ROI", "gradient")])
colnames(grad_fjxgrad)[1:3] <- paste0("Dim", c(1:3))
grad_fj_wide <- reshape(grad_fjxgrad, idvar = "ROI", timevar = "gradient", direction = "wide") %>% as.data.frame
grad_fj_wide$Network <- grad4grp[grad_fj_wide$ROI, "Network"]
rownames(grad_fj_wide) <- grad_fj_wide$ROI

net_mean.fj <- getMeans(grad_fj_wide[,2:10], grad_fj_wide$Network)

## color of arrows ----
arrow.col <- darken(roi.col$gc[rownames(net_mean),], 0.4)
names(arrow.col) <- rownames(net_mean)

## organize data - ROI average ----
spins_mean <- getMeans(grad.reg.wide4plot[,4:6], grad.reg.wide4plot$ROI)
grad2plot <- spins_mean
colnames(grad2plot) <- c(paste0("Gradient ", c(1:3)))
grad4grp <- unique(grad.reg.wide4plot[,c(1,2)]) %>% data.frame
rownames(grad4grp) <- grad4grp$ROI

## compute the difference with GLM ----
grad.reg.wide.diff <- grad.reg.wide
grad.reg.wide.diff$diagnostic_group <- sub.dx[grad.reg.wide$Subject,"diagnostic_group"]

all_roi_results <- grad.reg.wide.diff %>%
  select(-Subject) %>%
  ungroup() %>%
  group_by(network, gradient, ROI) %>%
  do(tidy(lm(value ~ diagnostic_group, data = .))) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(p_FDR = p.adjust(p.value, method = "fdr"))

## results from t-test filtered by p_FDR ----
all_roi_test <- all_roi_results %>%
  filter(term == "diagnostic_groupcontrol") %>%
  mutate(roi = str_remove(ROI, "_ROI"),
         label = case_when(str_starts(ROI, "L") ~ str_c("lh_", roi),
                           str_starts(ROI, "R_") ~ str_c("rh_", roi)),
         statistic.cor = ifelse(p_FDR < 0.05, statistic, 0)) %>%
  filter(ROI != "L_10pp_ROI") %>%
  as.data.frame() %>%
  group_by(gradient)

## results from t-test filtered by p_FDR (subcortical) ----
grad.reg.wide.sub <- grad.reg.wide.diff %>%
  select(-Subject) %>%
  filter(network == "Subcortical") %>%
  ### Mapping atlases ----
  mutate(Aseg = rep(rep(c(rep("Right-Hippocampus", 2), # right hippocampus
                          rep("Right-Amygdala", 2), # right amygdala
                          rep("Right-Thalamus-Proper", 4), # right thalamas
                          rep("Right-Nucleus-Accumbuns", 2), # right NAcc
                          rep("Right-Pallidum", 2), # right GP
                          rep("Right-Putamen", 2), # right putaman
                          rep("Right-Caudate", 2), # right caudate
                          rep("Left-Hippocampus", 2), # left hippocampus
                          rep("Left-Amygdala", 2), # left amygdata
                          rep("Left-Thalamus-Proper", 4), # left thalamas
                          rep("Left-Nucleus-Accumbuns", 2), # left Nacc
                          rep("Left-Pallidum", 2), # left GP
                          rep("Left-Putamen", 2), # left putaman
                          rep("Left-Caudate", 2)), # left caudate
                        each = 3), each = 419)) %>%
  ungroup() %>%
  group_by(network, gradient, Aseg) %>%
  do(tidy(lm(value ~ diagnostic_group, data = .))) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(p_FDR = p.adjust(p.value, method = "fdr"))

all_subroi_results <- grad.reg.wide.sub %>% 
  filter(term == "diagnostic_groupcontrol") %>%
  mutate(statistic.cor = ifelse(p_FDR < 0.05, statistic, 0)) %>%
  filter(!(Aseg %in% c("Left-Nucleus-Accumbuns", "Right-Nucleus-Accumbuns"))) %>%
  as.data.frame() %>%
  group_by(gradient)

# subT2plot$Aseg <- rep(c("Left-Caudate",
#                           "Right-Caudate",
#                           "Left-Pallidum",
#                           "Right-Pallidum",
#                           "Left-Hippocampus",
#                           "Right-Hippocampus",
#                           "Left-Putamen",
#                           "Right-Putamen",
#                           rep(c("Left-Amygdala",
#                           "Right-Amygdala"), 2),
#                           rep(c("Left-Nucleus-Accumbuns",
#                                 "Right-Nucleus-Accumbuns"), 2),
#                           "Left-Caudate",
#                           "Right-Caudate",
#                           "Left-Pallidum",
#                           "Right-Pallidum",
#                           "Left-Hippocampus",
#                           "Right-Hippocampus",
#                           "Left-Putamen",
#                           "Right-Putamen",
#                           rep(c("Left-Thalamus-Proper",
#                                 "Right-Thalamus-Proper"), 4))
#                         ,3)

### plot T in subcortical ----
for (grad in c(1:3)){
  gradient_raw_data_sub.G <- all_subroi_results %>%
    as.data.frame() %>%
    select(gradient, network, statistic.cor, label = Aseg) %>%
    group_by(label, gradient) %>%
    filter(gradient == paste0("grad",grad))
  
  gradient_raw_aseg_sub.G <- aseg %>%
    as.tibble %>% left_join(gradient_raw_data_sub.G)
  gradient_raw_aseg_sub.G$gradient <- paste0("grad",grad)
  
  if (grad == 1){
    gradient_raw_aseg_sub <- gradient_raw_aseg_sub.G
  }else{
    gradient_raw_aseg_sub <- rbind(gradient_raw_aseg_sub, gradient_raw_aseg_sub.G)
  }
}

gradAseg.test <- gradient_raw_aseg_sub %>% as_brain_atlas()


## brain data ----
brain_sig_grad1dat <- pqxcj_for_plot %>%
  as.data.frame() %>%
  group_by(gradient) %>%
  filter(gradient == "grad1")
cjsig1pick <- brain_sig_grad1dat$cjsig1 %>% setNames(brain_sig_grad1dat$ROI)

brain_sig_grad2dat <- pqxcj_for_plot %>%
  as.data.frame() %>%
  group_by(gradient) %>%
  filter(gradient == "grad2")
cjsig2pick <- brain_sig_grad2dat$cjsig1 %>% setNames(brain_sig_grad2dat$ROI)

brain_sig_grad3dat <- pqxcj_for_plot %>%
  as.data.frame() %>%
  group_by(gradient) %>%
  filter(gradient == "grad3")
cjsig3pick <- brain_sig_grad3dat$cjsig1 %>% setNames(brain_sig_grad3dat$ROI)

# Subcortical ----
subcor.col <- c("grey60", "grey60", "#a6cee3", "#1f78b4", 
                        "grey60", "#b2df8a", "grey60", "#B03060", 
                        "#ff7f00", "#CD6889", "grey60", "grey60", 
                        "#33a02c", "#EE9A49", "#B452CD", "#68228B", 
                        "#FFD700", "#D1AE22", rep("grey60", 11)) %>% setNames(aseg$data$label)
                        
names(subcor.col)[c(1,2,29)] <- NA
## by group (original gradients) ----
subroi2plot2 <- grad2plot2[grad4grp2$Network == "Subcortical",]

## Mapping Tien to aseg ----
subroi2plot <- pq_for_plot[pq_for_plot$network == "Subcortical",]
subroi2plot$Aseg <- rep(c(rep("Right-Hippocampus", 2), # right hippocampus
                          rep("Right-Amygdala", 2), # right amygdala
                          rep("Right-Thalamus-Proper", 4), # right thalamas
                          rep("Right-Nucleus-Accumbuns", 2), # right NAcc
                          rep("Right-Pallidum", 2), # right GP
                          rep("Right-Putamen", 2), # right putaman
                          rep("Right-Caudate", 2), # right caudate
                          rep("Left-Hippocampus", 2), # left hippocampus
                          rep("Left-Amygdala", 2), # left amygdata
                          rep("Left-Thalamus-Proper", 4), # left thalamas
                          rep("Left-Nucleus-Accumbuns", 2), # left Nacc
                          rep("Left-Pallidum", 2), # left GP
                          rep("Left-Putamen", 2), # left putaman
                          rep("Left-Caudate", 2)), # left caudate
                        each = 3)

## plot gradients in subcortical ----
for (grad in c(1:3)){
  gradient_raw_data.G <- subroi2plot %>%
    as.data.frame() %>%
    select(gradient, network, raw, label = Aseg) %>%
    group_by(label, gradient) %>%
    summarize(raw = mean(raw)) %>%
    filter(!(label %in% c("Left-Nucleus-Accumbuns", "Right-Nucleus-Accumbuns"))) %>%
    filter(gradient == paste0("grad",grad))
  
  gradient_raw_aseg.G <- aseg %>%
    as.tibble %>% left_join(gradient_raw_data.G)
  gradient_raw_aseg.G$gradient <- paste0("grad",grad)
  
  if (grad == 1){
    gradient_raw_aseg <- gradient_raw_aseg.G
  }else{
    gradient_raw_aseg <- rbind(gradient_raw_aseg, gradient_raw_aseg.G)
  }
}

gradAseg <- gradient_raw_aseg %>% as_brain_atlas()


## by gradients (cj) ----
spins_cj_subcor <- spins_cj
lapply(spins_cj_subcor, function(x){
  subcort.dat <- x[ROI.idx == "Subcortical",]
  return(subcort.dat)
})

## results from PLS filtered by bootstrap (subcortical) ----
pq_for_plot_sub <- pq_for_plot %>%
  filter(network == "Subcortical") %>%
  ### Mapping atlases ----
mutate(Aseg = rep(c(rep("Right-Hippocampus", 2), # right hippocampus
                        rep("Right-Amygdala", 2), # right amygdala
                        rep("Right-Thalamus-Proper", 4), # right thalamas
                        rep("Right-Nucleus-Accumbuns", 2), # right NAcc
                        rep("Right-Pallidum", 2), # right GP
                        rep("Right-Putamen", 2), # right putaman
                        rep("Right-Caudate", 2), # right caudate
                        rep("Left-Hippocampus", 2), # left hippocampus
                        rep("Left-Amygdala", 2), # left amygdata
                        rep("Left-Thalamus-Proper", 4), # left thalamas
                        rep("Left-Nucleus-Accumbuns", 2), # left Nacc
                        rep("Left-Pallidum", 2), # left GP
                        rep("Left-Putamen", 2), # left putaman
                        rep("Left-Caudate", 2)), # left caudate
                      each = 3)) %>%
  ungroup() %>%
  group_by(network, gradient, Aseg) %>%
  filter(!(Aseg %in% c("Left-Nucleus-Accumbuns", "Right-Nucleus-Accumbuns"))) %>%
  as.data.frame()

### plot T in subcortical ----
for (grad in c(1:3)){
  pq_sub <- pq_for_plot_sub %>%
    as.data.frame() %>%
    select(gradient, network, fj.sig1, label = Aseg) %>%
    group_by(label, gradient) %>%
    filter(gradient == paste0("grad",grad))
  
  pq_aseg_sub.G <- aseg %>%
    as.tibble %>% left_join(pq_sub)
  pq_aseg_sub.G$gradient <- paste0("grad",grad)
  
  if (grad == 1){
    pq_aseg_sub <- pq_aseg_sub.G
  }else{
    pq_aseg_sub <- rbind(pq_aseg_sub, pq_aseg_sub.G)
  }
}

pqAseg.test <- pq_aseg_sub %>% as_brain_atlas()

## grab color ----
subroi.col <- c(rep("#33a02c", 2), # right hippocampus
                rep("#ff7f00", 2), # right amygdala
                rep("#1f78b4", 4), # right thalamas
                rep("#b15928", 2), # right NAcc
                rep("#68228B", 2), # right GP
                rep("#B03060", 2), # right putaman
                rep("#D1AE22", 2), # right caudate
                rep("#b2df8a", 2), # left hippocampus
                rep("#EE9A49", 2), # left amygdata
                rep("#a6cee3", 4), # left thalamas
                rep("#B57957", 2), # left Nacc
                rep("#B452CD", 2), # left GP
                rep("#CD6889", 2), # left putaman
                rep("#FFD700", 2)) # left caudate


roi.subcor.idx <- ROI.idx[ROI.idx == "Subcortical"] %>% as.data.frame
colnames(roi.subcor.idx) <- "Network"
roi.subcor.idx$color <- subroi.col

## By contributions and significance to PLS ----
### G1G2 ----
G1G2.subcorpoints2plot.sigG1 <- cbind(subroi2plot2$grad1.control, subroi2plot2$grad2.control)

sigsubcorROI.G1 <- subcor.col %in% ifelse(spins_cj$grad1[rownames(subroi2plot2), "cjsig1"] != 0, roi.subcor.idx[rownames(subroi2plot2),"color"], NA)
subcor.col.G1 <- ifelse(sigsubcorROI.G1 == TRUE, subcor.col, "grey96") %>% setNames(names(subcor.col))

### G2G3 ----
G2G3.subcorpoints2plot.sigG2 <- cbind(subroi2plot2$grad2.control, subroi2plot2$grad3.control)

sigsubcorROI.G2 <- subcor.col %in% ifelse(spins_cj$grad2[rownames(subroi2plot2), "cjsig1"] != 0, roi.subcor.idx[rownames(subroi2plot2),"color"], NA)
subcor.col.G2 <- ifelse(sigsubcorROI.G2 == TRUE, subcor.col, "grey96") %>% setNames(names(subcor.col))

sigsubcorROI.G3 <- subcor.col %in% ifelse(spins_cj$grad3[rownames(subroi2plot2), "cjsig1"] != 0, roi.subcor.idx[rownames(subroi2plot2),"color"], NA)
subcor.col.G3 <- ifelse(sigsubcorROI.G3 == TRUE, subcor.col, "grey96") %>% setNames(names(subcor.col))

# Correlation to symptoms ----
behav_sympt <- read.csv("data/spins_behav_data_full_03-03-2022.csv") %>%
  select(record_id, bsfs_total, qls_total, bprs_factor_total, sans_total_sc, 
         bsfs_sec1_total, bsfs_sec2_total, bsfs_sec3_total, bsfs_sec4_total, bsfs_sec5_total, bsfs_sec6_total,
         qls20_empathy, qls_factor_interpersonal, qls_factor_instrumental_role,
         qls_factor_intrapsychic, qls_factor_comm_obj_activities, bprs_factor_neg_symp,
         bprs_factor_pos_symp, bprs_factor_anxiety_depression, bprs_factor_activation, bprs_factor_hostility, sans_sub_affective_flat_blunt, sans_sub_alogia, sans_sub_avolition_apathy, sans_sub_asocial_anhedonia)
rownames(behav_sympt) <- behav_sympt$record_id

lol_spins_behav_smp <- cbind(lol_spins_behav, behav_sympt[lol_spins_behav$record_id, ])
lol_spins_behav_ssd <- lol_spins_behav_smp[lol_spins_behav_smp$diagnostic_group == "case",]
# lol_spins_behav_ssd <- lol_spins_behav_ssd[complete.cases(lol_spins_behav_ssd),]

## numeric data ----
spins_symp_ssd_indiv <- lol_spins_behav_ssd %>% 
  select(
    bsfs_sec1_total, bsfs_sec2_total, bsfs_sec3_total, bsfs_sec4_total, bsfs_sec5_total, bsfs_sec6_total,
    qls_factor_interpersonal, qls_factor_instrumental_role,
    qls_factor_intrapsychic, qls_factor_comm_obj_activities, bprs_factor_neg_symp,
    bprs_factor_pos_symp, bprs_factor_anxiety_depression, bprs_factor_activation, bprs_factor_hostility, sans_sub_affective_flat_blunt, sans_sub_alogia, sans_sub_avolition_apathy, sans_sub_asocial_anhedonia
  ) %>% data.frame
rownames(spins_symp_ssd_indiv) <- lol_spins_behav_ssd$record_id

# Save everything
save.image(file = "data/res4plots.rda")

