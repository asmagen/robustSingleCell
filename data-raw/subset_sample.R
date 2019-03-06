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

