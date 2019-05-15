library(robustSingleCell)

download_LCMV()
LCMV1_sce <- read_10x_data(file.path(tempdir(), "LCMV/LCMV1"), subsample = 500)
compute_summary_stats(LCMV1_sce)
