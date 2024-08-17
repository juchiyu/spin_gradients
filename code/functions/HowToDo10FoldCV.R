### 10-fold cross validation for PLSC ###
### Author: Ju-Chi Yu
### Date: Nov. 30th, 2023
##---------------------------------------

## Reading packages ----
library(TExPosition) # to run PLSC
library(purrr) # used in the kfold function

## Read needed functions ----
source("ProjectSupplementaryData4PLS.R")
source("PLS.kFoldCV.R")

## Read example data ----
data(beer.tasting.notes)
data1<-beer.tasting.notes$data[,1:8]
data2<-beer.tasting.notes$data[,9:16]

## Run PLSC ----
pls.res <- tepPLS(data1,data2, graphs = FALSE)

## Perform 10-fold validation
pls.cv <- PLS.kFoldCV(data1, data2, pls.res, k = 10)

## plot

plot(pls.res$TExPosition.Data$lx[,1], pls.cv$cross.validation.res$lx.hat[,1])
plot(pls.res$TExPosition.Data$ly[,1], pls.cv$cross.validation.res$ly.hat[,1])
