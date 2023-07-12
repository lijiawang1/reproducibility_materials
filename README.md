

# Part 1: Data

## Abstract

The dataset used in the case study is integrate from 20 collection of scRNA-seq datasets from peripheral blood mononuclear cells (PBMCs).   We note that some of these datasets contain patients with longitudinal records; we take only one sample from these multiple measurements to ensure independence. A total of 864 patients are available in the dataset  with three severity levels marked as "Severe/Critical" (318 patients), "Mild/Moderate" (353 patients), and "Healthy" (193 patients). For each patient, PBMC scRNA-seq data is available in the form of a matrix recording the expression levels of genes in several cell types. More details of the integration process can be found in Supplementary Section A. 

The following datasets are available in the data folder:

1. "patient_meta.rds": This dataset contains sample metadata where each row represents a patient (donor). 
2. "aggExprs_meta.rds": This dataset represents pseudo-bulk metadata where each row corresponds to a pseudo-cell. 
3. "aggExprs_scMerge.rds": This dataset consists of a pseudo-bulk matrix where genes are represented as rows and pseudo-cells are represented as columns.

A detailed description of the datasets is available in "README_data" file in the "data" folder.


# Part 2: Code

## Abstract



The code includes all the functions necessary to generate the major results in Section 2.4 (simulation study) and Section 3 (case study). 

For the simulation study results in Section 2.4, the R code ("Figure 2 results.R", "Figure 3 results.R", "Figure 4 results.R", "Figure 5 results.R") can be found in the "code/simulation code" folder. As the some simulations have a high computational cost, we provide simpler examples (not the exact results in the manuscript) for checking the reproducibility of the process in the R files ("Figure 3 quick example (not exact).R", "Figure 4 quick example (not exact).R", "Figure 5 quick example (not exact).R") located at the same folder.

For the case study results in Section 3, the R code ("M.1 feature.R","M.2 feature.R","M.3 feature.R","M.4 feature.R") can be found in the "code/featurization code" folder. The manuscript summarizes these case study results in Table 1 and Figure 6 in Section 3.2. We provide the R code ("Table 1 results.R", "Figure 6 results.R") to produce the table and figure.  

A detailed description of the these R code files is available in "README_code" file in the "code" folder.


## Description

### Supporting software requirements

#### Version of primary software used

R version 4.2.2

#### Libraries and dependencies used by the code

R CRAN packages:
nnet 7.3-18;
randomForest 4.7-1.1;
e1071 1.7-12;
pROC 1.18.0;
dplyr 1.1.1;
ggplot2 3.4.0;
gridExtra 2.3;
readr 2.1.3;
resample 0.6;
stats 4.2.2;

These packages can be downloaded by \texttt{install.packages("...")} and loaded by \texttt{library("...")}.


### Instructions

<!--
Describe how to use the materials provided to reproduce analyses in the manuscript. Additional details can be provided in file(s) accompanying the reproducibility materials. If no workflow is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

Each R code file can be run individually as they contain all the necessary functions required for execution.

Before these running the R files ("M.1 feature.R","M.2 feature.R","M.3 feature.R","M.4 feature.R","Table 1 results.R","Figure 6 results.R"), please ensure that the user specifies the path to the entire folder, currently named "reproducibility_materials". This can be done by modifying the code "setwd("...")" to reflect the correct directory where the "reproducibility_materials" folder is located. Ensuring that the path is correctly specified will allow the R files to access the necessary files and directories for execution.

The output of each R code file and approximate run-time:

  i.) Located at "code/simulation code",

  * "Figure 2 results.R": reproduces the exact result in Figure 2 (<20min)
  * "Figure 3 results.R": reproduces the exact result in Figure 3 (~2hr)
  * "Figure 4 results.R": reproduces the exact result in Figure 4 (~2hr)
  * "Figure 5 results.R": reproduces the exact result in Figure 5 (~2hr)
  * "Figure 3 quick example (not exact).R": a quick example for Figure 3 (<10min)
  * "Figure 4 quick example (not exact).R": a quick example for Figure 4 (<10min)
  * "Figure 5 quick example (not exact).R": a quick example for Figure 5 (<10min)


  ii.) Located at "code/featurization code",
  
  * "M.1 feature.R": output classification performance results for M.1 feature. (<45min)
  * "M.2 feature.R": output classification performance results for M.2 feature. (<45min)
  * "M.3 feature.R": output classification performance results for M.3 feature. (<45min)
  * "M.4 feature.R": output classification performance results for M.2 feature. (<45min)
  
  These results are summarized in Table 1 and Figure 6 in the manuscript.

  iii.) Located at "code",
  
  * "Table 1 results.R": reproduces Table 1. (<3min)
  * "Figure 6 results.R": reproduces Figure 6. (<3min)
  
  In "Figure 6 results.R", by changing the code "feature\_method $= ``..."$", user can get the half violin plots for the corresponding featurization method. 
The input values that can be used are "m1", "m2", "m3", and "m4", representing the M.1, M.2, M.3, and M.4 featurization methods, respectively. 


# Notes 
More detailed descriptions of data and code are in README_data and README_code in the corresponding folders.


