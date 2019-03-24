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
                   min.genes.per.cell = 10,
                   max.genes.per.cell.quantile = 1,
                   max.UMIs.per.cell.quantile = 1,
                   min.cells.per.gene = 1)

# Tests for filtering and normalization

test_that("The number of cells and genes left after filtering ", {
  expect_equal(dim(LCMV1$counts), c(1390, 27))
})

test_that("The normalization and scaling ", {
  expect_equal(LCMV1$normalized[1,2], 3.920285, tolerance = 1e-6)
  expect_equal(LCMV1$normalized[11,1], 5.009584, tolerance = 1e-6)
})


# Tests for getting variable genes

context("Get variable genes")

LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
                            min.dispersion.scaled = 0.1, rerun = T)
test_that("The number of variable genes ", {
  expect_equal(length(LCMV1$HVG), 641)
})

# Tests for score computation

context("Score computation")
exhaustion_markers <- c('Pdcd1', 'Cd244', 'Havcr2', 'Ctla4', 'Cd160', 'Lag3',
                        'Tigit', 'Cd96')
ribosomal.score <- ribosomal.score(LCMV1)
mitochondrial.score <- mitochondrial.score(LCMV1)
cell.cycle.score <- cell.cycle.score(LCMV1)
Exhaustion <- controlled.mean.score(LCMV1, exhaustion_markers)

test_that("The scoring function ", {
  expect_equal(unname(ribosomal.score[10]), 2.889561, tolerance = 1e-6)
  expect_equal(unname(mitochondrial.score[11]), 0.5267472, tolerance = 1e-6)
  expect_equal(unname(cell.cycle.score[12, 3]), -0.3081372, tolerance = 1e-6)
  expect_equal(unname(Exhaustion[13]), -0.1526918, tolerance = 1e-6)
})
