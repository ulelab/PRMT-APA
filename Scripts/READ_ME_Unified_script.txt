
#################################
READ_ME before running Unfied_script.R
Questions: martina.2.pedna@kcl.ac.uk
#################################


INPUT FILES:

This script takes as input 2 different set of files:

-Time point files after PRMT inhibition, outputs of DESEQ analysis 
-siRNA inhibition of RBP files, outputs of DESEQ analysis  


TO RUN:

#need to give permission to run

Bash:
$chmod a+x Unfied_script.R

#required arguments:

-n_tp number of dataset PRMT inhibited
-n_RBP: number of dataset RBP

for each datset give path to dataset and the time point or RBP name:
-sample_tm_*: for timepoint ds name
-sample_RBP_*: for rbp ds name
-ds_tp_*: for the tp ds path
-ds_RBP_*: for the RBP ds path

N.B: * number from 1 to number of dataset

optional ones: 
-c1_prmt='value of first data change for Prmt category, default 0.2' 
-c2_prmt='value of other dataset change for prmt category default 0.15' 
-padj_prmt='default 0.05' 
-similarity_threshold_PRMT=default:0.7
-similarity_threshold_reg=default:0.7
-c1_reg= 'value for change for reg category, default:0.2'
-c2_reg = 'p-value for reg cat default=0.05'
-c1_ctrl='value for max change control site default 0.15'
-c2_ctrl='value for min use basal control site default=0.10'

#Mock command:

Rscript ./Unified_script.r -n_tp=5 -n_RBP=3 -ds_tp_1= $PATH/file_tp_24hrs.csv -sample_tp_1='24hrs' 
-ds_tp_2=$PATH/file_tp_48hrs.csv -sample_tp_2='48hrs' -ds_tp_3=$PATH/file_tp_72hrs.csv -sample_tp_3='72hrs' 
-ds_tp_4=$PATH/file_tp_96hrs.csv -sample_tp_4='96hrs'
-ds_RBP_1=$PATH/file_RBP_CFIM25.csv -sample_RBP_1='siCFIM25' -ds_RBP_2=$PATH/file_RBP_ELAVL.csv -sample_RBP_2='siELAVL1' 


WHAT IT DOES:

##reshapes the data into 1 unique dataframe where:
-ROWS: each row is a Poly Adenylation Site (PAS)
-Columns:The variables consist of a unique identifier for each PAS, the gene name, chromosome number, start and end in chromosome coordinates of the PAS, and the strand. Furthermore, a set of columns coded as ‘variable_name + dataset_name’ are present for each dataset. They contain information regarding: polyA usage in control sample and condition, mean difference change in usage of PAS inhibition compared to controls, and the adjusted p-value.

+saves reshaped file without filtering


##categorises sites into 3 different classes:
PRMT regulated:
-filters sites to have at least 20% change (default, can change using -c1_prmt='float:0.2') in 1 PRMT timepoint dataset + checks for significance: padjust <0.05(default, can change using -padj_prmt='float:0.05')
-checks that in all other timepoint ds the same PAS has at least 15% change (default, can change using -c2_prmt='float:0.15') in one other ds
-checks that the PRMT inhibition induced change is as strong as the siRNA effect with 70% threshold (default, can change using similarity_threshold_PRMT='float:0.7')
-selects one PAS per gene that has the biggest change and picks from all other PAS of that gene the one that has the strongest opposite change
-if PAS have a change in more than 1 RBP ds: keeps only the PAS pair that has the biggest change 
-writes file with just data from this category

Regulatable sites:
-filters PAS that have 20% change (default, can change using -c1_reg='float:0.2') in 1 of the RBP datasets and are significant:  padjust <0.05(default, can change using -c2_reg='float:0.05')
-checks that the siRNA effect is as strong as the PRMTi one: threshold 70% (default, can change using similarity_threshold_reg='float:0.7')
-picks 1 PAS per gene that has the biggest change and selects from all other PAS of that gene, the one with the biggest opposite effect
-removes genes previously found in the PRMT reg category
-writes a file with just this category

Control sites:
-removes all genes previously found in the other 2 categories
-picks pairs of PAS per gene that have a basal 10% usage in controls (default, can change using -c2_ctrl='float:0.1') and that have a change in usage less than 15% (default, can change using -c1_ctrl='float:0.15')
-writes a file with just data from this category


##Merges all categories together

##PAS Distal/Proximal + Isoform lengthening/shortening
-For each PAS checks the coordinates and strand and assigns if the PAS is proximal or distal compared to the other PAS of the pair
-Assigns if a gene has a lengtening or shortening isoform


##writes a final filtered file with all categories
 
