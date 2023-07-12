The "code" folder contains all the necessary functions to generate the main results presented in Section 2.4 (simulation study) and Section 3 (case study). 

For the simulation study results in Section 2.4, the R code can be found in the "code/simulation code" folder. The following R files are available:

	* "Figure 2 results.R": reproduces the exact result in Figure 2 (<20min)
	* "Figure 3 results.R": reproduces the exact result in Figure 3 (~2hr)
	* "Figure 4 results.R": reproduces the exact result in Figure 4 (~2hr)
	* "Figure 5 results.R": reproduces the exact result in Figure 5 (~2hr)
As the last three simulations have a high computational cost, we provide simpler examples (not the exact results in the manuscript) for checking the reproducibility of the process as follows:

	* "Figure 3 quick example (not exact).R": provides a quick example for Figure 3 (<10min)
	* "Figure 4 quick example (not exact).R": provides a quick example for Figure 4 (<10min)
	* "Figure 5 quick example (not exact).R": provides a quick example for Figure 5 (<10min)

For the case study results in Section 3, the R code can be found in the "code/featurization code" folder. The following R files are available:

	* "M.1 feature.R": reproduces the classification performance results with respect to the M.1 featurization method. The outputs are stored in "logistic_m1.Rdata", "forest_m1.Rdata", "svm_m1.Rdata" files in the "output" folder. (<45min)

	* "M.2 feature.R": reproduces the classification performance results with respect to the M.2 featurization method. The outputs are stored in "logistic_m2.Rdata", "forest_m2.Rdata", "svm_m2.Rdata" files in the "output" folder. (<45min)
	
	* "M.3 feature.R": reproducing the classification performance results with respect to the M.3 featurization method. The outputs are stored in "logistic_m3.Rdata", "forest_m3.Rdata", "svm_m3.Rdata" files in the "output" folder. (<45min)

	* "M.4 feature.R": reproduces the classification performance results with respect to the M.1 featurization method. The outputs are stored in "logistic_m4.Rdata", "forest_m4.Rdata", "svm_m4.Rdata" files in the "output" folder. (<45min)

Section 3.1 of the manuscript provides detailed information about the featurization methods. The necessary datasets ("patient_meta.rds", "aggExprs_meta.rds", "aggExprs_scMerge.rds") required by these R files are located in the "data" folder. To facilitate understanding, we have included a "README_data" file in the "data" folder, which explains the format and meaning of these datasets in more detail. You can refer to the "README_data" file for further clarification on the dataset structure and content. Before these running the R files, please ensure that the user specifies the path to the entire folder, i.e., the "reproducibility_materials" folder (if you downloaded from github, the folder name could be "reproducibility_materials-main"), by modifying the code "setwd("~/Desktop/reproducibility_materials")".


The manuscript summarizes the case study results in Table 1 and Figure 6 in Section 3.2. The R code to produce the table and figure are available in the files located at the "code" folder: 

	* "Table 1 results.R": reproduces Table 1 as presented in the manuscript. (<3min)
	* "Figure 6 results.R": reproduces plots in Figure 6 as presented in the manuscript. By change the code 'feature_method = "m1"', user can get the half violin plots for the corresponding featurization method where "mi" stands for featurization method M.i proposed in section 3.1. e.g., "m1", "m2", "m3", "m4". (<3min)

These R code files take the results generated from previous steps, which are stored in the "output" folder, as inputs. By running these files, you can reproduce the results and generate Table 1 and Figure 6 based on the stored outputs. Again, before running these R files, please ensure that the user specifies the path to the entire folder, by modifying the code "setwd("~/Desktop/reproducibility_materials")".




