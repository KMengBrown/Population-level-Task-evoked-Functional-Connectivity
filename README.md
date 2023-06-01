# Population-level-Task-evoked-Functional-Connectivity

1. Description: This repository is devoted to the ptFCE algorithm. This algorithm comes from the paper "Population-level Task-evoked Functional Connectivity" by Kun Meng and Ani Eloyan (hereafter ME) available at https://arxiv.org/pdf/2102.12039.pdf. Additionally, this repository contains the code script files for simulation studies and data analyses in ME.
2. The version of primary software used: R version 3.6.2.
3. Maintainer: Kun Meng <kun_meng@brown.edu>


The following file contains function "ptFCE." This function implements Algorithms 1 of ME. Running the code in this file is needed before reproducing any simulations or data analyses in ME.
* ptFCE Algorithm.R

The following file illustrates the implementation of function "ptFCE" using simulations.
* Illustrations.R


The paper ME has two data components: simulated data and the motor-task data set publicly available from the Human Connectome Project (HCP) database. We provide the complete code for replicating all simulation studies included in ME herein. To replicate the motor-task dataset analysis, we provide the full code that can be applied to data downloaded from HCP.


The following five files replicate the simulation studies for methods comparison presented in ME.
* Comparison_Mech_1.R
* Comparison_Mech_2.R
* MLM_based_Mechanism_1
* MLM_based_Mechanism_2
* 50_node_study

The following file replicates the simulation study for ptFCE estimation bias/variance in the Supporting Information of ME.
* Estimation_bias_and_variance_ptFCE.R

The following file replicates the simulation study for the influence of random noise on estimation bias/variance.
* Influence_of_noise.R


The referred HCP data are available online at https://protocols.humanconnectome.org/HCP/3T/task-fMRI-protocol-details.html. These data are publicly available by request, following the process described here: They can be obtained from ConnectomeDB using the link provided here. An HCP user account and Aspera browser plugin are required for downloading the data. Since we used a subset of the subjects, we provide the list of the subjects used in our analysis. Once the fMRI data are downloaded from HCP, the code we provide herein can be used to replicate the full analysis presented in ME.


The following file presents the IDs of subjects in the HCP data set used in our data analyses.
* Subjects_IDs.csv

The following file replicates the procedure of computing region-specific time courses corresponding to the AAL atlas. Its output is a collection of 308 ".rda" files. Running this code is needed before replicating the data analyses of ME.
* HCPdata.R

The following file replicates the data analyses presented in ME. Before running this code, you need to put the "Subjects_IDs.csv" file and all the 308 ".rda" data files in the same folder and set the corresponding working directory, where the 308 ".rda" files are the output from running the "HCPdata.R" code.
* Data_analysis.R

Further details are provided in each script file.


Kun Meng, Ph.D.
Prager Assistant Professor
Division of Applied Mathematics
Brown University
