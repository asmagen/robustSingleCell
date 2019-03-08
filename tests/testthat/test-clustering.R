context("Filtering and normalization")
library(robustSingleCell)

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


# Tests for PCA computation

context("PCA computation")

library(rslurm)

local_slurm_array <- function(slr_job) {
    olddir <- getwd()
    rscript_path <- file.path(R.home("bin"), "Rscript")
    setwd(paste0("_rslurm_", slr_job$jobname))
    tryCatch({
        writeLines(c(paste0("for (i in 1:", slr_job$nodes, " - 1) {"),
                     "Sys.setenv(SLURM_ARRAY_TASK_ID = i)",
                     "source('slurm_run.R')", "}"), "local_run.R")
        system(paste(rscript_path, "--vanilla local_run.R"))
    }, finally = setwd(olddir))
    return(slr_job)
}

LCMV1 <- PCA(LCMV1,local = T)

test_that("The number of 'significant' PCs ", {
  expect_equal(dim(LCMV1$PCA), c(7, 97))
})


# Tests for clustering computation

context("Clustering computation")

LCMV1 <- cluster.analysis(LCMV1,knn.ratios = c(0.1,0.2),loadPreviousKnn = F, rerun = T, deleteCache = T, plot = F, local = T)
0.1
LCMV1$clustering$
# Need to resolve resolution input parameter
test_that("The number of clusters ", {
  expect_equal(LCMV1$clustering$nclusters, 6)
})

