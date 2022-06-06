##################################
## Save loadings to Nifti files ##
##################################

rm(list = ls())

library(dplyr)
library(tidyverse)
library(InPosition)
library(ggplot2)
library(readxl)
library(oro.nifti)
library(PTCA4CATA)
library(superheat)
library(oro.nifti)

# set up paths ----
mask.path <- "../../Data/Mask/"
pet.path <- "../../Data/PET/"
indiv.path <- "../../Data/PET/indiv/"
mni.path <- "../../Data/canonical/"
demo.path <- "../../Data/Demographics/"
func.dir <- "./functions"
home.dir <- getwd()
# source all functions==
setwd(func.dir)
sapply(list.files(), source)
setwd(home.dir)
#=======================

load("ca.hg.rda")

mniBrain <- readNIfTI(paste0(mni.path,"avg152T1.nii"))
LatentBrain <- readNIfTI(paste0(indiv.path,"sub-007_ses-01_SRTM_BPnd_MNI.nii.gz"))

str(LatentBrain@.Data)

# sort data ----
FactorScores <- sta.dxmat
FactorScores$`Dimension2` <- ca.hg$ExPosition.Data$fj[as.character(sta.dxmat$voxID),2]
FactorScores$`Dimension4` <- ca.hg$ExPosition.Data$fj[as.character(sta.dxmat$voxID),4]
head(FactorScores)

# map values back to brain array ----
## create empty array ----
loadmap <- array(0, dim = dim(LatentBrain@.Data))

# save fj of dimension 2 and 4
brain.lv2 <- loadmap
brain.lv4 <- loadmap
for (i in 1:nrow(FactorScores)){
  brain.lv2[FactorScores$rind[i], FactorScores$cind[i], FactorScores$tind[i]] <- FactorScores$Dimension2[i]
  brain.lv4[FactorScores$rind[i], FactorScores$cind[i], FactorScores$tind[i]] <- FactorScores$Dimension4[i]
}


# Back to nifti ----
Brain_lv2 <- Brain_lv4 <- LatentBrain

Brain_lv2@.Data <- brain.lv2
Brain_lv4@.Data <- brain.lv4

# write nifti ----
writeNIfTI(
  nim = Brain_lv2,
  filename = "./Results/0.4_2_PET_cahg_Brain_exWM_lv2")

writeNIfTI(
  nim = Brain_lv4,
  filename = "./Results/0.4_2_PET__cahg_Brain_exWM_lv4"
)
