Sure thing. The SPINS behavioural data is here: /projects/loliver/SPINS_PLS_Conn/data/processed/spins_behav_data_full_11-20-2020.csv

loliver  6:03 PM
As for the QC exclusions, here's a list based on the RS imaging (fmriprep and ciftify) QC
exclude_img <- c("SPN01_CMH_0104","SPN01_CMH_0136","SPN01_MRP_0143","SPN01_ZHH_0040","SPN01_ZHH_0047",
                  "SPN01_ZHH_0048","SPN01_ZHH_0052","SPN01_ZHP_0063","SPN01_ZHP_0093","SPN01_ZHP_0105")





6:04
I also have an exclusion list including EA task, motion (mean FD > .5), and imaging QC for the EA task data, but you're just using rest right?
6:08
I always check the termination info after other exclusions to see if anyone remains too, and then decide on a case by case basis if they should be included or not, depending on the data you're looking at etc. I extracted the termination info here: /projects/loliver/SPINS_PLS_Conn/data/processed/spins_termination_info_11-20-2020.csv