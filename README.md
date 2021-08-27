# spin_gradients
project for tigrlab

by Lucy Bassman e:lbassman@uoguelph.ca

May 20 2021 - August 27 2021

Report Paper can be found at https://docs.google.com/document/d/1Qr25lIVcH24J8fKIM8v2vZNBFHagVYNCdLWrhlToU7w/edit

Presentation can be found at https://docs.google.com/presentation/d/1y-f1USF6cm8X_cUrxV9pGlbY1QvMKObm_CIh7cFLqWI/edit#slide=id.ge5e02efbea_0_18

Derived SPINS study data can be found at https://github.com/edickie/ssd_gradients , or locally on KIMEL server at /scratch/edickie/ssd_gradients

Requirements for the virtual environment for Jupyter can be found at nilearn_brainspace_requirements.txt

Software used for this project includes Jupyter notebook and RStudio. Brainspace packages were installed to the python environment to facilitate gradient creation.

all output csv files can be found under the main directory /scratch/a/arisvoin/lbassman/spins_gradients/spin_gradients/

- four csv files containing first 10 gradients for all subjects with networks:
  - spins_concat_full.csv has the full ea timeseries and no gsr applied,
  - spins_concat_shortened.csv has the shortened ea timeseries (209 timepoints to match rest) and no gsr applied,
  - gsr_spins_concat_full.csv has the full ea timeseries and gsr applied,
  - gsr_spins_concat_shortened.csv has the short ea and gsr applied,
  
- four csv files containing the network centroids and average dispersion for every subject:
  - network_averages_full.csv has the full ea timeseries and no gsr applied, 
  - network_averages_shortened.csv has the shortened ea timeseries and no gsr applied,
  - gsr_network_averages_full.csv has the full ea timeseries with gsr applied,
  - gsr_network_averages_shortened.csv has the short ea timeseries and gsr applied,
  
- the spins_gradients folder (directory above spin_gradients) has a folder for each subject in the study, with pscalars for every gradient variation (gsr, task, concatenated, shortened) 

- template gradients used for alignment (HCP and margulies):
- HCP gradients found at `/group-HCPS1200_atlas-GlasserTian_desc-subcorticalS2_conn.pconn.nii`
- margulies gradients found at `/tpl-fsLR_den-32k_atlas-Glasser2016Tian2019S2_desc-margulies2016_gradients.pscalar.nii`
