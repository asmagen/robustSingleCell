check_not_slurm <- function(func_name) {
    not_slurm = suppressWarnings(system('sinfo', ignore.stdout = TRUE, ignore.stderr = TRUE))
    slurm_msg = paste0('SLURM not detected. Please run ',
                       func_name,
                       ' on a SLURM cluster.\n')
    if (not_slurm) cat(slurm_msg)
    return(not_slurm)
}

#' Parallelized PCA Analysis
#'
#' Run PCA analysis with a simulation analysis of shuffled data to determine the appropriate number of PCs.
#'
#' @param environment \code{environment} object
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
#' @import rslurm
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' }
PCA <- function(environment, regress = NA, groups = NA, nShuffleRuns = 10, threshold = 0.1,
    maxPCs = 100, label = NA, mem = "2GB", time = "0:10:00", rerun = F, clear.previously.calculated.clustering = T, local = F) {

    if (check_not_slurm("PCA")) {
        return(environment)
    }
    if (length(regress) > 1 || !is.na(regress)) {
        config <- paste(colnames(regress), collapse = "+")
    } else {
        config <- "not.regressed"
    }

    if (length(groups) == 1 && is.na(groups))
        groups <- rep(1, environment$nsamples)
    if (nShuffleRuns != 10)
        config <- paste(config, nShuffleRuns, sep = ".")
    if (threshold != 0.1)
        config <- paste(config, threshold, sep = ".")
    if (maxPCs != 100)
        config <- paste(config, maxPCs, sep = ".")
    if (!is.na(label))
        config <- paste(config, label, sep = ".")

    environment$work.path <- file.path(environment$baseline.work.path, config)
    environment$PCA <- environment$Rotation <- environment$PCA.path <- NA
    if (clear.previously.calculated.clustering)
        environment$clustering <- environment$seurat.cluster.association <- NA
    dir.create(file.path(environment$work.path, "tracking"), showWarnings = F, recursive = T,
        mode = "700")
    print.message("Transitioning to", config, "folder")

    environment$res.data.path <- file.path(environment$work.path, "data")
    shuffled.PCA.data.path <- file.path(environment$res.data.path, "shuffled.PCA")

    dir.create(environment$res.data.path, showWarnings = F, recursive = T, mode = "700")

    cache <- file.path(environment$res.data.path, paste(config, "PCA.rds", sep = "."))

    if (!rerun && file.exists(cache)) {
        print.message("Loading precomputed")
        precomputed <- readRDS(cache)
        PCA <- precomputed$PCA
        Rotation <- precomputed$Rotation
        rm(precomputed)
    } else {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"))
        on.exit(end(t))
        if (clear.previously.calculated.clustering)
            unlink(file.path(environment$res.data.path, "clustering"), recursive = T,
                force = T)

        unlink(shuffled.PCA.data.path, recursive = T, force = T)
        dir.create(shuffled.PCA.data.path, showWarnings = F, recursive = T, mode = "700")

        raw.data <- environment$counts
        data <- environment$normalized[environment$HVG, ]
        print.message("Dim")
        print(dim(data))

        if (length(regress) > 1 || !is.na(regress)) {
            corrected <- regress.covariates(environment = environment, regress = regress, data = data, groups = groups, rerun = rerun, save = T)
            print.message("Regressed matrix")
            corner(corrected)
            data <- corrected
        }

        n <- min(maxPCs, ncol(data))
        m <- nrow(data)
        ndf <- n - 1

        get.shuffled.var <- function(rep) {

            data.perm <- apply(data, 2, sample, replace = FALSE)
            data.perm <- t(data.perm)
            data.perm <- data.perm[,apply(data.perm, 2, var) > 0]
            pca.perm <- stats::prcomp(data.perm, retx = TRUE, center = T, scale = T)
            var.perm <- pca.perm$sdev[1:ndf]^2/sum(pca.perm$sdev[1:ndf]^2)
            saveRDS(list(pca.perm = pca.perm, var.perm = var.perm), file = file.path(shuffled.PCA.data.path,
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

        t <- start(file.path(environment$work.path, "tracking"), append = T, split = T)
        on.exit(end(t), add = T)

        if (length(var.perm) == 0 || nrow(var.perm) < nShuffleRuns) {
            print.message("JOB ERROR: Not enough shuffled results:", nrow(var.perm),
                "<", nShuffleRuns, "\nCHECK FOR FAILED JOBS\n\n\n")
            terminate <- readline(prompt = "Terminate? (y/n) ")
            if (terminate != "n") {
                #end(t)
                tryCatch({
                  cleanup_files(sjob)
                  cleanup_files(sjob)
                  cleanup_files(sjob)
                }, error = function(v) v)
                return(environment)
            }
        }


        t <- start(file.path(environment$work.path, "tracking"), append = T)
        on.exit(end(t), add = T)

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

    environment$PCA <- PCA
    environment$Rotation <- Rotation
    environment$PCA.path <- cache

    cat("# PCs = ", nrow(environment$PCA), "\n", sep = "")

    return(environment)
}

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
