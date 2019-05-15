#' Parallelized PCA Analysis
#'
#' Run PCA analysis with a simulation analysis of shuffled data to determine the appropriate number of PCs.
#'
#' @param sce \code{SingleCellExperiment} object
#' @param regress gene signature activation scores to regress
#' @param groups experimental design annotation to guide dataset-specific regression
#' @param nShuffleRuns number of shuffled analyses
#' @param threshold FDR threshold
#' @param maxPCs maximum number of possible PCs
#' @param label optional analyses label folder
#' @param mem HPC memory
#' @param time HPC time
#' @param rerun whether to rerun the analysis rather than load from cache
#' @param clear.previously.calculated.clustering whether to clear previous clustering analysis
#' @param local whether to run jobs locally on slurm instead of submitting the job
#' @return \code{environment} parameter containing PC coordinates
#' @export
#' @import rslurm SingleCellExperiment
shuffled_PCA <- function(sce, regress = NULL, groups = NULL, nShuffleRuns = 10, threshold = 0.1,
                maxPCs = 100, label = NULL, mem = "2GB", time = "0:10:00", rerun = F,
                local = F) {
  browser()
  #check_not_slurm("PCA")
  if (length(regress) > 1 || !is.null(regress)) {
    config <- paste(colnames(regress), collapse = "+")
  } else {
    config <- "not.regressed"
  }

  if (length(groups) == 1 && is.null(groups))
    groups <- rep(1, ncol(sce))
  if (nShuffleRuns != 10)
    config <- paste(config, nShuffleRuns, sep = ".")
  if (threshold != 0.1)
    config <- paste(config, threshold, sep = ".")
  if (maxPCs != 100)
    config <- paste(config, maxPCs, sep = ".")
  if (!is.null(label))
    config <- paste(config, label, sep = ".")

  cache_dir <- file.path(tempdir(), "shuffled.PCA")
  cache <- file.path(cache_dir, paste(config, "PCA.rds", sep = "."))

  if (!rerun && file.exists(cache)) {
    print.message("Loading precomputed")
    precomputed <- readRDS(cache)
    PCA <- precomputed$PCA
    Rotation <- precomputed$Rotation
    rm(precomputed)
  } else {
    print.message("Computing")
    unlink(cache_dir, recursive = T, force = T)
    dir.create(cache_dir, showWarnings = F, recursive = T, mode = "700")

    data <- assay(sce, "normcounts")[rowData(sce)$is_HVG, ]
    print.message("Dim")
    print(dim(data))

    if (length(regress) > 1 || !is.na(regress)) {
      corrected <- regress.covariates(regress, data, groups, rerun, save = T)
      print.message("Regressed matrix")
      corner(corrected)
      data <- corrected
    }

    n <- min(maxPCs, ncol(data))
    ndf <- n - 1

    get.shuffled.var <- function(rep) {

      data.perm <- apply(data, 2, sample, replace = FALSE)
      data.perm <- t(data.perm)
      data.perm <- data.perm[,apply(data.perm, 2, var) > 0]
      pca.perm <- stats::prcomp(data.perm, retx = TRUE, center = T, scale = T)
      var.perm <- pca.perm$sdev[1:ndf]^2/sum(pca.perm$sdev[1:ndf]^2)
      saveRDS(list(pca.perm = pca.perm, var.perm = var.perm), file = file.path(cache_dir,
                                                                               paste("shuffled.PCA.rep", rep, "rds", sep = ".")))
      return(var.perm)
    }

    sopt <- list(mem = mem, time = time, share = TRUE)

    if (local) {
      sjob <- slurm_apply(get.shuffled.var, data.frame(rep = seq(nShuffleRuns)),
                          add_objects = c("shuffled.PCA.data.path", "data", "ndf"), pkgs = NULL,
                          nodes = nShuffleRuns, cpus_per_node = 1, submit = FALSE, slurm_options = sopt)
      local_slurm_array(sjob)
    } else {
      sjob <- slurm_apply(get.shuffled.var, data.frame(rep = seq(nShuffleRuns)),
                          add_objects = c("shuffled.PCA.data.path", "data", "ndf"), pkgs = NULL,
                          nodes = nShuffleRuns, cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
    }

    pc.time <- Sys.time()
    pca <- stats::prcomp(t(data), retx = TRUE, center = T, scale = T)
    print.message("Single PCA run time:")
    print(Sys.time() - pc.time)
    var <- pca$sdev[1:ndf]^2/sum(pca$sdev[1:ndf]^2)
    print.message("Real PCA Var")
    print(utils::head(var, 10))

    var.perm <- get_slurm_out(sjob)
    dim(var.perm)

    if (length(var.perm) == 0 || nrow(var.perm) < nShuffleRuns) {
      print.message("JOB ERROR: Not enough shuffled results:", nrow(var.perm),
                    "<", nShuffleRuns, "\nCHECK FOR FAILED JOBS\n\n\n")
      terminate <- readline(prompt = "Terminate? (y/n) ")
      if (terminate != "n") {
        tryCatch({
          cleanup_files(sjob)
          cleanup_files(sjob)
          cleanup_files(sjob)
        }, error = function(v) v)
        return(sce)
      }
    }

    tryCatch({
      cleanup_files(sjob)
      cleanup_files(sjob)
      cleanup_files(sjob)
    }, error = function(v) v)
    print.message("Shuffled PCA Var")
    corner(var.perm, 10)

    p <- rep(1, n)
    for (i in 1:ndf) {
      p[i] <- mean(var.perm[, i] >= var[i])
    }
    p
    for (i in 2:ndf) {
      p[i] <- max(p[(i - 1)], p[i])
    }
    nPCs <- sum(p <= threshold, na.rm = T)
    print.message("nPCs", nPCs)

    PCA <- t(pca$x)[seq(nPCs), ]
    print.message("PCA")
    corner(PCA)
    Rotation <- t(pca$rotation[, seq(nPCs)])
    print.message("Rotation")
    corner(Rotation)

    saveRDS(list(PCA = PCA, Rotation = Rotation), file = cache)
  }

  reducedDim(example_sce)$PCA <- AnnotatedDataFrame(PCA, varMetaData = Rotation)

  cat("# PCs = ", nrow(environment$PCA), "\n", sep = "")

  return(sce)
}
