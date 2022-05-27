## Compute correlation
# set working dir
# setwd("/scratch/loliver/SPINS_PLS_Conn")
# load function to do Fisher's transform
source("code/functions/r2z.R")
##------------------------------From Lindsay, May-11-2022---------
# find RS time series files  # pattern glob2rx("*_RS_2mm_noGSR_glasser_meants.csv")
files_RS_ts <- list.files(path= "/scratch/loliver/SPINS_PLS_Conn", recursive=T, full.names=T, pattern="^.*_RS_2mm_GSR_glasser_tian_meants\\.csv$")
# confirm csvs aren't empty
files_RS_ts[file.size(files_RS_ts) == 0]
# create list of IDs
ptlist <- paste("SPN01", substring(files_RS_ts,32+5,32+7), substring(files_RS_ts,32+8,32+11), sep = "_")
# read in time series files
RS_ts <- lapply(files_RS_ts, read.csv, header=F)
# transpose dfs
RS_ts <- lapply(RS_ts, t)
# Name dfs with participant IDs
names(RS_ts) <- ptlist

# compute correlation matrix
rs_cor <- lapply(RS_ts, cor)
##-----------------------------------------------------------------
# Fisher's z-transformed
rs_zcor <- lapply(rs_cor, r2z)

save(rs_cor, file = "data/spins_RS_2mm_GSR_GlasserTian_cor.rda")
save(rs_zcor, file = "data/spins_RS_2mm_GSR_GlasserTian_zcor.rda")