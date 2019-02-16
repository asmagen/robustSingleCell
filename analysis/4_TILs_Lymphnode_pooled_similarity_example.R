
# Pooled TILs and Lymphnode scRNAseq samples analysis ______________________________________________________________________________________
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
screen -r -d TILs_LN_pooled_example # run first time: screen -S TILs_LN_pooled_example
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
	datasets = c('LN1','TILsPd1'), # List of dataset code names
	origins = c('Lymphnode','TILs PD1+'), # List of dataset tissue origin/condition full name
	experiments = c('Exp 1','Exp 1'), # List of experiment design annotations (Replicate experiment annotation)
	data.path = '/home/USER_NAME/data', # Path to where the data is located
	work.path = '/home/USER_NAME/analysis', # path to where the analysis should create folders and store files
	marker.genes = marker.genes) # Set of genes of interest for visualization purposes

# PREPROCESSING _____________________________________________________________________________________________________________________
# 

environment = read.preclustered.datasets () # Read datasets

environment = add.confounder.variables ( # Compute gene signatures activation for signatures of interest; used for regression and plotting
	ribosomal.score = ribosomal.score(),
	mitochondrial.score = mitochondrial.score(),
	cell.cycle.score = cell.cycle.score(),
	Exhaustion = controlled.mean.score(Exhaustion))

environment = PCA (clear.previously.calculated.clustering = F) # Perform PCA analysis with shuffled simulation to define the appropriate number of PCs

# DIFFERENTIAL EXPRESSION AND VISUALIZATIONS ________________________________________________________________________________________
# 

summarize (contrast = 'datasets') # Perform differential expression analysis, compute tSNE and plot various visualizations

# PERFORM CLUSTER SIMILARITY ANALYSES _______________________________________________________________________________________________
# 

cluster.similarity = compare.cluster.similarity (diff.exp.file = 'main.datasets.diff.exp.RData',cluster.similarity.function = pearson.correlation,label = 'pearson',rerun = F)
similarity  = cluster.similarity$similarity
map = cluster.similarity$map
visualize.cluster.cors.heatmaps (environment$work.path,similarity)

# Compute robust clusters - possible only if there are 2 replicate experiments per sample origin
filtered.similarity = get.robust.cluster.similarity(similarity,min.sd = qnorm(.9),max.q.val = 0.01,rerun = F)
robust.clusters = sort(unique(c(filtered.similarity$cluster1,filtered.similarity$cluster2)));length(robust.clusters)
head(filtered.similarity);dim(filtered.similarity)
summary(filtered.similarity)
visualize.cluster.cors.heatmaps (environment$work.path,filtered.similarity)
similarity = filtered.similarity


# Visualize cluster similarity metrics
plot.cluster.similarity.stats(environment,similarity)
