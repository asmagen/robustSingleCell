
# TILs scRNAseq sample analysis ______________________________________________________________________________________
# 
# This pipeline was designed, written and produced by Assaf Magen as part of CD4+ T tumor-infiltrating lymphocytes scRNAseq study
# Single-cell resolution profiling of tumor-reactive CD4+ T-cells reveals immune landscape alterations
# 
# Laboratory of Immune Cell Biology, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, Maryland
# Center for Bioinformatics and Computational Biology, University of Maryland, College Park, Maryland
# 
# Using any piece of this code or its modifications require the appropriate citations and agreement to code usage terms
# Default parameters should be adjusted to account for dataset-specific data features as needed and should not be taken as is
# 
# Contact: Assaf Magen, assaf.magen@nih.gov

# SET UP ________________________________________________________________________________________________________________________
# 

# Mount directories locally (not on the server)
# create ~/projects/biowulf.home locally first
/usr/local/bin/sshfs USER_NAME@SERVER:/home/USER_NAME/ ~/projects/biowulf.home -o volname=biowulf.home

# Connect to biowulf
ssh USER_NAME@SERVER
screen -r -d TILs_example # run first time: screen -S TILs_example
sinteractive --mem=16g --time=18:00:00
module purge;module load R/3.4;R

# Load source code
options("width"=250)
source('/home/USER_NAME/scripts/Magen_scRNAseq.pipeline.R')

# Define gene signatures of interest for regression or visualization
library(data.table);markers = data.frame(fread('/home/USER_NAME/data/markers.csv'))[,1:5];marker.genes = unique(markers$symbol)
Exhaustion = c('Pdcd1','Cd244','Havcr2','Ctla4','Cd160','Lag3','Tigit','Cd96')

# Initialize project information and drive paths
environment = initialize.project (
	datasets = c('TILsPd1'), # List of dataset code names
	origins = c('TILs PD1+'), # List of dataset tissue origin/condition full name
	experiments = c('Exp 1'), # List of experiment design annotations (Replicate experiment annotation)
	data.path = '/home/USER_NAME/data', # Path to where the data is located
	work.path = '/home/USER_NAME/analysis', # path to where the analysis should create folders and store files
	marker.genes = marker.genes) # Set of genes of interest for visualization purposes

# PREPROCESSING _____________________________________________________________________________________________________________________
# 

environment = read.data () # Read dataset and perform quality filtering

environment = get.variable.genes () # Identify highly variable genes

environment = add.confounder.variables ( # Compute gene signatures activation for signatures of interest; used for regression and plotting
	ribosomal.score = ribosomal.score(),
	mitochondrial.score = mitochondrial.score(),
	cell.cycle.score = cell.cycle.score(),
	Exhaustion = controlled.mean.score(Exhaustion))

environment = PCA () # Perform PCA analysis with shuffled simulation to define the appropriate number of PCs

# CLUSTERING ________________________________________________________________________________________________________________________
# 

environment = cluster.analysis () # Perform clustering analysis scanning a range of clustering resolutions

# DIFFERENTIAL EXPRESSION AND VISUALIZATIONS ________________________________________________________________________________________
# 

summarize () # Perform differential expression analysis, compute tSNE and plot various visualizations

# Filter specific clusters then rerun analysis from first step (initialize.project)
filter.cluster.data (remove.clusters = 2) # run only once when filtration is needed

# CLUSTER NAME ASSIGNMENT ___________________________________________________________________________________________________________
# 

# Get provisional cluster names and set them in the dataset

types = rbind(
	data.frame(type='Tfh',gene=c('Tcf7','Cxcr5','Bcl6')),
	data.frame(type='Th1',gene=c('Cxcr6','Ifng','Tbx21')),
	data.frame(type='Isc',gene=c('Irf7','Stat1','Ifit3')),
	data.frame(type='Treg',gene=c('Foxp3','Il2ra')),
	data.frame(type='Ccr7',gene=c('Ccr7')),
	data.frame(type='Il7r',gene=c('Il7r')),
	data.frame(type='Cycle',gene=c('Mki67','Top2a','Birc5'))
)

cluster.names = get.cluster.names (types)
environment = set.cluster.names (names = cluster.names)
summarize () # Rerun summary to include new cluster names in figures
