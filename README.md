
######## Author names, contact details
Sarthak P. Malusare*, Giacomo Zilio*, Emanuel A. Fronhofer*
*ISEM, Université de Montpellier, CNRS, IRD, EPHE, Montpellier, France

Contacts:  sarthak-pravin.malusare@etu.umontpellier.fr
	   giacomo.zilio@umontpellier.fr	  
	   emanuel.fronhofer@umontpellier.fr

Responsible for collecting data: Sarthak P. Malusare
Responsible for statistical analysis and code annotation: Giacomo Zilio

######## Title
Evolution of thermal performance curves: a meta-analysis of selection experiments

######## Summary
Review and meta-analysis of thermal performance curves before and after evolution in laboratory-based studies.

######## Preprint
https://www.biorxiv.org/content/10.1101/2022.05.09.491229v2

######## Overview of folders

### The results of the study are based on 4 MAIN analyses, each with its own folder.
0-analysis_JEB:
-R script to reproduce the statistical analysis
-dataset (0-analysis_JEB.txt)
-README file describing the dataset columns

2-analysis_JEB:
-R script to reproduce the statistical analysis
-dataset (2-analysis_JEB.txt)
-README file describing the dataset columns

Multi-analysis_JEB:
-R script to reproduce the statistical analysis
-dataset (Multi_analysis_JEB.txt)
-README file describing the dataset columns

Hotter_JEB:
-R script to reproduce the statistical analysis
-dataset (Hotter_JEB.txt)
-README file describing the dataset columns



### Investigation of variance and heterogeneity in I^2_code folder and Posteriors_JEB folder.
I^2_code:
-The 4 datasets and README of each main analysis
-4 R scripts to calculate variance and heterogeneity

Posteriors_JEB:
-Posteriors for the best models of the 4 main analyses, for variance and heterogeneity.
To be used in the 4 R scripts of the I^2_code folder to calculate variance and heterogeneity



### SUPPLEMENTARY analysis
Funnel_JEB folder:
-R script to reproduce the supplementary analysis
-datasets of the 4 main analyses + complete dataset with all Curve_ID (All_data_JEB.txt)
-README file describing each dataset columns


### phylogenetic tree
The Phylogenetic_Tree_JEB folder contains the phylogenetic tree generated with the taxize package as described in the main text as a .nex file.
