##################################
## Save loadings to Nifti files ##
##################################

rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(oro.nifti)
library(PTCA4CATA)
library(superheat)
library(oro.nifti)

# set up paths ----
mni.path <- "data/canonical/"
#=======================

load("results/data4Plot_from_4.RData")

mniBrain <- readNIfTI(paste0(mni.path,"avg152T1.nii"))
LatentBrain <- readNIfTI(paste0(mni.path,"Tian_Subcortex_S2_3T.nii"))
label.dx <- read.table(paste0(mni.path, "Tian_Subcortex_S2_3T_label.txt"))
rownames(label.dx) <- 1:32

str(LatentBrain@.Data)

# sort data ----
PlotData <- list(grad1 = pqxcj_for_plot[pqxcj_for_plot$gradient == "grad1" & pqxcj_for_plot$network == "Subcortical",],
                 grad2 = pqxcj_for_plot[pqxcj_for_plot$gradient == "grad2" & pqxcj_for_plot$network == "Subcortical",],
                 grad3 = pqxcj_for_plot[pqxcj_for_plot$gradient == "grad3" & pqxcj_for_plot$network == "Subcortical",])

rownames(PlotData$grad1) <- PlotData$grad1$ROI
rownames(PlotData$grad2) <- PlotData$grad2$ROI
rownames(PlotData$grad3) <- PlotData$grad3$ROI
PlotData

# map values back to brain array ----
## create empty array ----
loadmap <- array(NA, dim = dim(LatentBrain@.Data))

# save fj of dimension 2 and 4
brain <- list()
brain$raw <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjimp1_x_mean.pos <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjimp1_x_mean.neg <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjimp2_x_mean.pos <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjimp2_x_mean.neg <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjsig1_x_mean.pos <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjsig1_x_mean.neg <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjsig2_x_mean.pos <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)
brain$cjsig2_x_mean.neg <- list(grad1 = loadmap, grad2 = loadmap, grad3 = loadmap)

# Replace values
gothrough <- c("raw", 
               "cjimp1_x_mean.pos", "cjimp1_x_mean.neg",
               "cjimp2_x_mean.pos", "cjimp2_x_mean.neg",
               "cjsig1_x_mean.pos", "cjsig1_x_mean.neg",
               "cjsig2_x_mean.pos", "cjsig2_x_mean.neg")


for (value in 1:length(gothrough)){
  for (grad in 1:3){
    for (roi in 1:32){
      brain[[gothrough[value]]][[grad]][LatentBrain@.Data == roi] <- as.numeric(PlotData[[grad]][label.dx[roi,], gothrough[value]])  
    }
    tmp <- LatentBrain
    tmp@.Data <- brain[[gothrough[value]]][[grad]]
    assign(sprintf("BrainFrom4.%s.grad%s", gothrough[value],grad), tmp)
    
    nifti.name <- sprintf("data/pls_subcort_nii/BrainFrom4.%s.grad%s", gothrough[value],grad)
    
    # write to nifti
    writeNIfTI(
      nim = get(sprintf("BrainFrom4.%s.grad%s", gothrough[value],grad)),
      filename = sprintf("data/pls_subcort_nii/BrainFrom4.%s.grad%s", gothrough[value],grad))
  }
}

# Back to nifti ----
# write nifti ----

