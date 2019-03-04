LCMV1 <- initialize.project(datasets = "LCMV1",
                            origins = "CD44+ cells",
                            experiments = "Rep1",
                            data.path = "~/LCMV/",
                            work.path = "~/LCMV/LCMV_analysis2")
LCMV1 <- read.data(LCMV1, subsample = 100)
LCMV1 <- get.variable.genes(LCMV1)


# randomly sample
LCMV2 <- initialize.project(datasets = "LCMV2",
                            origins = "CD44+ cells",
                            experiments = "Rep2",
                            data.path = "~/LCMV/",
                            work.path = "~/LCMV/LCMV_analysis2")
LCMV2 <- read.data(LCMV2, subsample = 100)
LCMV2 <- get.variable.genes(LCMV2)

subset_genes <- union(LCMV1$HVG, LCMV2$HVG)

LCMV1_normalized <- LCMV1$normalized[LCMV1$HVG,]
LCMV2_normalized <- LCMV2$normalized[LCMV2$HVG,]
devtools::use_data(LCMV1_normalized, robustSingleCell)
devtools::use_data(LCMV2_normalized, robustSingleCell)
