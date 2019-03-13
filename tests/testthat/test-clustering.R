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
                            data.path = data.path)

LCMV1 <- read.data(LCMV1,
                   raw.data.matrices = list(LCMV1 = raw_LCMV1),
                   min.genes.per.cell = 100,
                   max.genes.per.cell.quantile = 1,
                   max.UMIs.per.cell.quantile = 1,
                   min.cells.per.gene = 1)

LCMV1 <- get.variable.genes(LCMV1)

test_that("The number of 'significant' PCs and number of clusters", {
  if (SLURM) skip(SLURM_MSG)
  LCMV1 <- PCA(LCMV1, local = T)
  expect_equal(dim(LCMV1$PCA), c(7, 97))
  LCMV1 <- cluster.analysis(LCMV1, knn.ratios = c(0.1, 0.2), loadPreviousKnn = F,
                            rerun = T, deleteCache = T, plot = F, local = T)
  expect_equal(LCMV1$clustering$nclusters, 6)
})

