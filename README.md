# CSE284_Group_16

TCGA Website Link: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Squamous%20Cell%20Carcinoma%20(LUSC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 

Gene data: “HTSeq - Counts” 

https://xenabrowser.net/datapages/?dataset=TCGA-LUSC.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 

Phenotype data: “phenotype” 

https://xenabrowser.net/datapages/?dataset=TCGA-LUSC.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 

Survival data: “survival” 

https://xenabrowser.net/datapages/?dataset=TCGA-LUSC.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 

Data dictionary: https://docs.gdc.cancer.gov/Data_Dictionary/viewer 

From the TCGA websites, download:
- TCGA-LUSC.GDC_phenotype.tsv
- TCGA-LUSC.htseq_counts.tsv
- TCGA-LUSC.survival.tsv
and unzip the first two gz files.

Code documentation:

- `workspace`
Reads in the raw downloaded files and prepares df_train and df_test RData files as well as the annotations Excel file.

- `workspace2`
Reads in the df_test and df_train RData files and generates uninformative, OLS, and LMM models and figures.

- `ryo_analysis`
Explores data and methods of specifying random effects in LMMs. Compares LMMs with differing covariates and error terms. Generates Manhattan plots for protein-coding genes.

- `formulas`
Formats Null, OLS, and LMM formulas for slides and report.

