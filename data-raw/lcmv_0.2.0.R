library(robustSingleCell)
library(SingleCellExperiment)

download_LCMV()
LCMV1_sce <- read_10x_data(file.path(tempdir(), "LCMV/LCMV1"))
LCMV1_sce <- robustSingleCell::filter(LCMV1_sce, verbose = T)
LCMV1_sce <- log_normalize(LCMV1_sce, verbose = T)
LCMV1_sce <- get_variable_genes(LCMV1_sce, verbose = T)
LCMV1_sce <- ribosomal_score(LCMV1_sce)
LCMV1_sce <- mitochondrial_score(LCMV1_sce)
LCMV1_sce <- cell_cycle_score(LCMV1_sce, verbose = T)

exhaustion_markers <- c('Pdcd1', 'Cd244', 'Havcr2', 'Ctla4', 'Cd160', 'Lag3', 'Tigit', 'Cd96')
exhaustion_score <- controlled_mean_score(LCMV1_sce, exhaustion_markers)
LCMV1_sce <- add_cell_score(LCMV1_sce, "exhaustion", exhaustion_score)
colnames(colData(LCMV1_sce))
# TODO allow PCA to run on non-slurm environment
LCMV1_sce <- shuffled_PCA(LCMV1_sce, local = T, regress = NULL)
reducedDims(LCMV1_sce)
