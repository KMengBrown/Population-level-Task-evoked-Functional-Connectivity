# Population-level-Task-evoked-Functional-Connectivity

1. Description: This repository is devoted to the ptFCE and AMUSE-ptFCE algorithms. These algorithms come from the paper "Population-level Task-evoked Functional Connectivity" by Kun Meng and Ani Eloyan (hereafter ME) available at https://arxiv.org/pdf/2102.12039.pdf. Additionally, this repository contains the code script files for simulation studies and data analyses in ME.
2. The version of primary software used: R version 3.6.2.
3. Maintainer: Kun Meng <kun_meng@brown.edu>


The following file contains functions "ptFCE" and "AMUSE_ptFCE." These two functions implement Algorithms 1 and 2, respectively, of ME. Running the code in this file is needed before reproducing any simulations or data analyses in ME.
* ptFCE and AMUSE-ptFCE Algorithms.R

The following file illustrates the implementation of functions "ptFCE" and "AMUSE_ptFCE" using simulations.
* Illustrations.R


The paper ME has two data components: simulated data and the motor-task data set publicly available from the Human Connectome Project (HCP) database. We provide the complete code for replicating all simulation studies included in ME herein. To replicate the motor-task dataset analysis, we provide the full code that can be applied to data downloaded from HCP.


The following two files replicate the simulation studies presented in Tables 1 and 2 of ME.
* Simulation_Table_1.R
* Simulation_Table_2.R

The following file replicates the simulation study results presented in Figure 5 of ME.
* Influence_of_noise.R


The referred HCP data are available online at https://protocols.humanconnectome.org/HCP/3T/task-fMRI-protocol-details.html. These data are publicly available by request, following the process described here: They can be obtained from ConnectomeDB using the link provided here. An HCP user account and Aspera browser plugin are required for downloading the data. Since we used a subset of the subjects, we provide the list of the subjects used in our analysis. Once the fMRI data are downloaded from HCP, the code we provide herein can be used to replicate the full analysis presented in ME.


The following file presents the IDs of subjects in the HCP data set used in our data analyses.
* Subjects_IDs.csv

The following file replicates the procedure of computing region-specific time courses corresponding to the AAL atlas. Its output is a collection of 308 ".rda" files. Running this code is needed before replicating the data analyses presented in Figures 6 and 7 of ME.
* HCPdata.R

The following file replicates the data analyses presented in Figures 6 and 7 of ME. Before running this code, you need to put the "Subjects_IDs.csv" file and all the 308 ".rda" data files in the same folder and set the corresponding working directory, where the 308 ".rda" files are the output from running the "HCPdata.R" code.
* Data_analysis.R

Further details are provided in each script file.


Kun (Michael) Meng,

Ph.D. Candidate, 
Department of Biostatistics, 
Brown University
