context("PCA computation")
library(robustSingleCell)
library(rslurm)

SLURM = system('sinfo', ignore.stdout = TRUE, ignore.stderr = TRUE)
SLURM_MSG = 'Only test on Slurm head node.'

data.path <- system.file("extdata", package = "robustSingleCell")
raw_LCMV1 <- as.matrix(read.table(file.path(data.path, "LCMV1_small.txt"), check.names = F))
LCMV1 <- initialize.project(datasets = "LCMV1",
                            origins = "CD44+ cells",
                            experiments = "Rep1",
                            data.path = data.path,
                            work.path = file.path(tempdir(), "LCMV/LCMV_analysis"))

LCMV1 <- read.data(LCMV1,
                   raw.data.matrices = structure(list(raw_LCMV1), names = "LCMV1"),
                   min.genes.per.cell = 10,
                   max.genes.per.cell.quantile = 1,
                   max.UMIs.per.cell.quantile = 1,
                   min.cells.per.gene = 1, rerun = T)

LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
                            min.dispersion.scaled = 0.1, rerun = T)

test_that("The number of 'significant' PCs and number of clusters", {
  if (SLURM) skip(SLURM_MSG)
  LCMV1 <- PCA(LCMV1, local = T)
  expect_equal(dim(LCMV1$PCA), c(8, 27))
#  LCMV1 <- cluster.analysis(LCMV1, knn.ratios = c(0.1, 0.2), loadPreviousKnn = F,
#                            rerun = T, deleteCache = T, plot = F, local = T)
#  expect_equal(LCMV1$clustering$nclusters, 6)
# types = rbind(
#                 data.frame(type='Tfh',gene=c('Tcf7','Cxcr5','Bcl6')),
#                 data.frame(type='Th1',gene=c('Cxcr6','Ifng','Tbx21')),
#                 data.frame(type='Tcmp',gene=c('Ccr7','Bcl2','Tcf7')),
#                 data.frame(type='Treg',gene=c('Foxp3','Il2ra')),
#                 data.frame(type='Tmem',gene=c('Il7r','Ccr7')),
#                 data.frame(type='CD8',gene=c('Cd8a')),
#                 data.frame(type='CD4', gene = c("Cd4")),
#                 data.frame(type='Cycle',gene=c('Mki67','Top2a','Birc5'))
# )
# summarize(LCMV1)
# LCMV1_cluster_names <- get.cluster.names(LCMV1, types, min.fold = 1.0, max.Qval = 0.2)
# LCMV1 <- set.cluster.names(LCMV1, names = LCMV1_cluster_names)
# summarize(LCMV1)
# library(ggplot2)
# plot-contour_overlay_tSNE (LCMV1,genes = c('Cxcr5','Cxcr6','Cd4','Cd8a'),perplexity = 10,max_iter = 10000,width = 10, height = 10)
# plot_pair_scatter (LCMV1,gene1 = 'Cxcr5',gene2 = 'Tcf7',cluster_group1 = c('Unknown_2','Tcmp'),cluster_group2 = c('CD8','Cycle_CD8'),group1_label = 'Tfh_Tcmp',group2_label = 'CD8',width = 10, height = 10)
# library(ggplot2);library(ggrepel)
# diff.exp = get.robust.markers (environment = LCMV1,cluster_group1 = c('Unknown_2','Tcmp'),cluster_group2 = c('CD8','Cycle_CD8'),group1_label = 'Tfh_Tcmp',group2_label = 'CD8')
# head(diff.exp)
})

