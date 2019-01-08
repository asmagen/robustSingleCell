
# scRNAseq pipeline initial setup ______________________________________________________________________________________

# This pipeline was designed, written and produced by Assaf Magen as part of CD4+ T tumor-infiltrating lymphocytes scRNAseq study
# Single-cell resolution profiling of tumor-reactive CD4+ T-cells reveals immune landscape alterations

# Laboratory of Immune Cell Biology, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, Maryland
# Center for Bioinformatics and Computational Biology, University of Maryland, College Park, Maryland

# Using any piece of this code or its modifications require the appropriate citations and agreement to code usage terms
# Default parameters should be adjusted to account for dataset-specific data features as needed and should not be taken as is

# Contact: Assaf Magen, assaf.magen@nih.gov


# INITIAL SET UP (ONLY PRIOR TO FIRST RUN) ______________________________________________________________________________________
# 

# Set up session

# Connect to biowulf
ssh magenay@biowulf.nih.gov
mkdir data
mkdir scripts
mkdir analysis
sinteractive
module purge;module load R/3.4;R

# Install packages
install.packages(pkgs = c('Matrix','limma','ggplot2','cccd','GGally','ggrepel','gplots','rslurm','dplyr','ggpubr','reshape2','RColorBrewer','xlsx','grid'))
if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("JinmiaoChenLab/Rphenograph")

# Make data folder with dataset matrices and scripts folder with code files
scp -r /Users/magenay/Desktop/TILsPd1 magenay@biowulf.nih.gov:/home/magenay/data/
scp -r /Users/magenay/Desktop/markers.csv magenay@biowulf.nih.gov:/home/magenay/data/
scp -r /Users/magenay/Desktop/regev_lab_cell_cycle_genes.txt magenay@biowulf.nih.gov:/home/magenay/data/
scp -r /Users/magenay/Desktop/sc.methods.pipeline.example.R magenay@biowulf.nih.gov:/home/magenay/scripts/

# Set up SSH key and connect remote folders for browsing
scp /Users/magenay/.ssh/id_rsa.pub magenay@biowulf.nih.gov:~/tmp
ssh -X magenay@biowulf.nih.gov
cat tmp >> ~/.ssh/authorized_keys
rm tmp
chmod 0600 ~/.ssh/authorized_keys
