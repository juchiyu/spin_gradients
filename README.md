# spin_gradients
project for tigrlab

all output csv files can be found under the main directory /scratch/a/arisvoin/lbassman/spins_gradients/spin_gradients/

four csv files containing first 10 gradients for all subjects with networks:
  spins_concat_full has the full ea timeseries and no gsr applied,
  spins_concat_shortened has the shortened ea timeseries (209 timepoints to match rest) and no gsr applied,
  gsr_spins_concat_full has the full ea timeseries and gsr applied,
  gsr_spins_concat_shortened has the short ea and gsr applied,
  
four csv files containing the network centroids and average dispersion for every subject:
  network_averages_full has the full ea timeseries and no gsr applied, 
  network_averages_shortened has the shortened ea timeseries and no gsr applied,
  gsr_network_averages_full has the full ea timeseries with gsr applied,
  gsr_network_averages_shortened has the short ea timeseries and gsr applied,
  
the spins_gradients folder (directory above spin_gradients) has a folder for each subject in the study, plus others that were excluded, with pscalars for every gradient variation (gsr, task, concatenated, shortened) 
  
