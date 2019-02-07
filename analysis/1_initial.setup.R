
# scRNAseq pipeline initial setup ______________________________________________________________________________________
# 
# This pipeline was designed, written and produced by Assaf Magen et al as part of CD4+ T tumor-infiltrating lymphocytes scRNAseq study
# Single-cell profiling of tumor-reactive CD4+ T-cells reveals unexpected transcriptomic diversity
# 
# Laboratory of Immune Cell Biology, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, Maryland
# Center for Bioinformatics and Computational Biology, University of Maryland, College Park, Maryland
# 
# Using any piece of this code or its modifications require the appropriate citations and agreement to code usage terms
# Default parameters should be adjusted to account for dataset-specific data features as needed and should not be taken as is
# 
# Contact: Assaf Magen, assaf.magen@nih.gov


# INITIAL SET UP (ONLY PRIOR TO FIRST RUN) ______________________________________________________________________________________
# 

# Set up session

# Replace USER_NAME and SERVER with your credentials

# Connect to biowulf
ssh USER_NAME@SERVER
mkdir data
mkdir scripts
mkdir analysis

# Set up SSH key and connect remote folders for browsing
scp /Users/USER_NAME/.ssh/id_rsa.pub USER_NAME@SERVER:~/tmp
ssh -X USER_NAME@SERVER
cat tmp >> ~/.ssh/authorized_keys
rm tmp
chmod 0600 ~/.ssh/authorized_keys

# Make data folder with dataset matrices and scripts folder with code files
scp -r /Users/USER_NAME/Desktop/TILsPd1 USER_NAME@SERVER:/home/USER_NAME/data/
scp -r /Users/USER_NAME/Desktop/LN1 USER_NAME@SERVER:/home/USER_NAME/data/
scp -r /Users/USER_NAME/Desktop/markers.csv USER_NAME@SERVER:/home/USER_NAME/data/
scp -r /Users/USER_NAME/Desktop/regev_lab_cell_cycle_genes.txt USER_NAME@SERVER:/home/USER_NAME/data/
scp -r /Users/USER_NAME/Desktop/sc.methods.pipeline.example.R USER_NAME@SERVER:/home/USER_NAME/scripts/

# Install packages
sinteractive
module purge;module load R/3.4;R
install.packages(pkgs = c('Matrix','limma','ggplot2','cccd','GGally','ggrepel','gplots','rslurm','dplyr','ggpubr','reshape2','RColorBrewer','xlsx','grid'))
if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("JinmiaoChenLab/Rphenograph")
