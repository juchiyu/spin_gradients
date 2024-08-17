This folder includes the code used in this study.

These scripts include some preliminary analyses and attempts. The main findings of the paper are from the following scripts:

- Extract gradients: `spins_gradient.ipynb`
- Rmarkdown for main PLSC results: `4_PLSC_grad_GSR_combat.rmd`
- Code to generate figures:
    + `6_all_analyses.R`: R code copied from `4_PLSC_grad_GSR_combat.rmd` to generate a single data file
    + `6_figures_for_paper_(from4).Rmd`: Takes the data generated above to create the 4 figures in the paper
- Code for analyses requested by reviewers: `7.1_PLSC_grad_GSR_combat_CV.rmd`
  > This includes cross-validation, top-half vs. bottom-half movers, separate PLSC for SSDs and Controls, the correlation between medication and group effects in gradients, and details in medication, wakefulness, and diagnoses.

Other useful functions can be found in the `functions` folder, including:
+ `z2r.R`: Function that transforms Fisher's Z-transformed matrix back to correlation matrix
+ `r2x.R`: Function that performs Fisher's Z-transformation on correlation matrix
+ `vec2sqmat.R`: Function that puts the vectorized correlation matrix back in to a square matrix
+ `PLS.kFoldCV.R`: Core function for PLSC k-fold cross-validation
+ `ProjectSupplementaryData4PLS.R`: Function called by `PLS.kFoldCV.R` to predict out-of-sample data from a PLSC model
+ `HowToDo10FoldCV.R`: Example code that demonstrates how to perform k-fold cross-validation of PLSC with the `PLS.kFoldCV` function
