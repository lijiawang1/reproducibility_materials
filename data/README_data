
The dataset used in the case study is integrate from 20 collection of scRNA-seq datasets.   We note that some of these datasets contain patients with longitudinal records; we take only one sample from these multiple measurements to ensure independence. A total of 864 patients are available in the dataset  with three severity levels marked as "Severe/Critical" (318 patients), "Mild/Moderate" (353 patients), and "Healthy" (193 patients). For each patient, PBMC scRNA-seq data is available in the form of a matrix recording the expression levels of genes in several cell types. More details of the integration process can be found in Supplementary Section A. 


1. "patient_meta.rds": This dataset contains sample metadata where each row represents a patient(donor). It provides information about the original dataset, the condition of the patient (e.g., "Severe/Critical", "Mild/Moderate", "Healthy"), the outcome of Covid-19, and additional clinical variables such as age, gender, ethnicity, BMI, and previous medical history (e.g., hypertension, heart disease, etc.). The "rp_id" feature assigns a unique ID to each patient, facilitating matching across the three datasets.

2. "aggExprs_meta.rds": This dataset represents pseudo-bulk metadata where each row corresponds to a pseudo-cell. It includes information about the cell type, the sample it belongs to, the original dataset, the condition of the corresponding patient, the outcome of Covid-19, and additional clinical variables.

3. "aggExprs_scMerge.rds": This dataset consists of a pseudo-bulk matrix where genes are represented as rows and pseudo-cells are represented as columns. The matrix is obtained by averaging across the samples and cell types of the scMerge batch-corrected matrix. The column names follow the format "rp_id|cell-type", indicating the specific cell type and corresponding patient.

