LCMV1 <- initialize.project(datasets = "LCMV1",
                            origins = "CD44+ cells",
                            experiments = "Rep1",
                            data.path = "~/LCMV/",
                            work.path = "~/LCMV/LCMV_analysis2")
LCMV1 <- read.data(LCMV1, subsample = 50)
mitogenes <- get.mito.genes(LCMV1$genes)
ribogenes <- get.ribo.genes(LCMV1$genes)
exhaustion_markers <- c('Pdcd1', 'Cd244', 'Havcr2', 'Ctla4', 'Cd160', 'Lag3',
                        'Tigit', 'Cd96')
LCMV1 <- get.variable.genes(LCMV1)
LCMV1_genes <- unique(c(ribogenes, mitogenes, LCMV1$HVG, intersect(exhaustion_markers, LCMV1$genes)))
LCMV1_small <- LCMV1$counts[LCMV1_genes,]
write.table(LCMV1_small, file = "inst/extdata/LCMV1_small.txt")

# only 205 genes that are variable
LCMV1 <- setup_LCMV_example()
LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
                            min.dispersion.scaled = 0.1, rerun = T)

# randomly sample
LCMV2 <- initialize.project(datasets = "LCMV2",
                            origins = "CD44+ cells",
                            experiments = "Rep2",
                            data.path = "~/LCMV/",
                            work.path = "~/LCMV/LCMV_analysis2")
LCMV2 <- read.data(LCMV2, subsample = 50)
LCMV2 <- get.variable.genes(LCMV2)
mitogenes <- get.mito.genes(LCMV2$genes)
ribogenes <- get.ribo.genes(LCMV2$genes)
LCMV2_genes <- unique(c(ribogenes, mitogenes, LCMV2$HVG, intersect(exhaustion_markers, LCMV2$genes)))
LCMV2_small <- LCMV2$counts[LCMV2_genes,]
write.table(LCMV2_small, file = "inst/extdata/LCMV2_small.txt")


# only 205 genes that are variable
LCMV2 <- setup_LCMV_example("LCMV2")
LCMV2 <- get.variable.genes(LCMV2, min.mean = 0.1, min.frac.cells = 0,
                            min.dispersion.scaled = 0.1, rerun = T)



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
