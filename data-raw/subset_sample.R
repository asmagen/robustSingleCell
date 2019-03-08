LCMV1 <- initialize.project(datasets = "LCMV1",
                            origins = "CD44+ cells",
                            experiments = "Rep1",
                            data.path = "~/LCMV/",
                            work.path = "~/LCMV/LCMV_analysis2")
LCMV1 <- read.data(LCMV1, subsample = 100)
LCMV1 <- get.variable.genes(LCMV1)
LCMV1_small <- LCMV1$counts[LCMV1$HVG,]
write.table(LCMV1_small, file = "inst/extdata/LCMV1_small.txt")


# randomly sample
LCMV2 <- initialize.project(datasets = "LCMV2",
                            origins = "CD44+ cells",
                            experiments = "Rep2",
                            data.path = "~/LCMV/",
                            work.path = "~/LCMV/LCMV_analysis2")
LCMV2 <- read.data(LCMV2, subsample = 100)
LCMV2 <- get.variable.genes(LCMV2)
LCMV2_small <- LCMV2$counts[LCMV2$HVG,]
write.table(LCMV2_small, file = "inst/extdata/LCMV2_small.txt")

# write examples
data.path <- system.file("extdata", package = "robustSingleCell")
raw_LCMV1 <- as.matrix(read.table(file = file.path(data.path, "LCMV1_small.txt"), check.names = F))
LCMV1_small <- initialize.project(datasets = "LCMV1",
                            origins = "CD44+ cells",
                            experiments = "Rep1",
                            data.path = data.path)
LCMV1_small <- read.data(LCMV1_small, raw.data.matrices = list(LCMV1 = raw_LCMV1),
min.genes.per.cell = 100, max.genes.per.cell.quantile = 1,
max.UMIs.per.cell.quantile = 1, min.cells.per.gene = 1)
LCMV1_small <- get.variable.genes(LCMV1_small)
usethis::use_data(LCMV1_small, overwrite = T)



LCMV1 <- setup_LCMV1_example()
# The data had been filtered already and parameters were chosen correspondingly
LCMV1 <- read.data(LCMV1,
raw.data.matrices = list(LCMV1 = raw_LCMV1),
min.genes.per.cell = 100,
max.genes.per.cell.quantile = 1,
max.UMIs.per.cell.quantile = 1,
min.cells.per.gene = 1)
