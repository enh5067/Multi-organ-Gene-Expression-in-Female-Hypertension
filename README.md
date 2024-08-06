# Multi-organ-Gene-Expression-in-Female-Hypertension
This repository contains all the code and data needed to reproduce the results of the paper “Multi-organ gene expression analysis and network modeling reveal regulatory control cascades during the development of hypertension in female spontaneously hypertensive rat".

## This repository contains 13 files.

### Code-Folder containing all code necessary to reproduce the results of the paper.
  - Multiorgan_HMF.R-Calling script to load the female expression data and run the HMF method. Loads the multiorgan data to run the HMF method on the SHR and WKY data separately. Also loads and integrates brainstem RNAseq data with common genes to the multiorgan RNAseq data. This data is also run through the HMF method on the SHR and WKY data separately. 
  - Multiorgan_HMF_Male-Female.R-Calling script, highly similar to 'Multiorgan_HMF.R' except it integrates Male SHR and WKY data from Anderson et al., 2017 to allow for a male-female comparison.
  - Negddct_expression_heatmaps.R-Calling script to load the multiorgan negddct data and visualize the data in heatmaps with annotations shown. 

### Data-Folder containing all data necessary to reproduce the results of the paper.
  - Organ_Norm_AllData.txt-Negddct data from HT-qPCR across 5 organs in the female SHR and WKY.
  - Brainstem_Subsetted_Genes.txt-Brainstem RNAseq data for genes in common with HT-qPCR data from female SHR and WKY.
  - Male_Norm_AllData-Negddct data from HT-qPCR across 5 organs in the male SHR and WKY, from Anderson et al., 2017.

### Functions—Folder containing all functions necessary to run the HMF method on a given dataset. All functions in this folder are accessed from the calling script Multiorgan_HMF.
  - Hartley_functions.R-This function performs linear interpolation between data points and implements the continuous Hartley transform.
  - Master.R-This function performs system identification to infer the connectivity and dynamics of a linear network model (Anderson et al., 2017).
  - m-value_functions.R-This function selects a set of values for M = max(m).
  - plotting_functions.R-This function plots simulation data.
  - regression_functions.R-This function establishes the Y vector for a specific set of m-values (Equation 11, Anderson et al., 2017).
  - scaling_functions.R-This function normalizes input data to the range (0,1).
  - simulation_functions.R-This function extracts the initial conditions for the dynamic variables corresponding to measurements.
