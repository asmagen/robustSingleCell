#' Perform tSNE Analyses
#'
#' Run distributed tSNE analysis for multiple input hyperparameters
#'
#' @param environment \code{environment} object
#' @param perplexity perplexity parameter of tSNE
#' @param max_iter maximum number of iterations to run the tSNE
#' @param rerun whether to rerun or load from cache
#' @param local whether to run tSNE locally
#' @param mem Memory for each job; default 4 GB
#' @param time Time for each job; default 15 minutes
#' @return Distributed job identified object
run_tSNE <- function(environment, perplexity, max_iter, rerun, local = F, mem = "4GB", time = "0:15:00") {

    tSNEs.dir <- file.path(environment$res.data.path, "tSNEs")
    list.files(tSNEs.dir)
    sjob <- NA

    if (rerun || !dir.exists(tSNEs.dir)) {
        t <- start(file.path(environment$work.path, "tracking"), split = T)
        on.exit(end(t))
        print(tSNEs.dir)
        # unlink(tSNEs.dir,recursive=T,force=T)
        dir.create(tSNEs.dir, showWarnings = F)

        duplicated.indices <- duplicated(t(environment$PCA))
        if (sum(duplicated.indices) > 0) {
            print.message("Excluding", sum(duplicated.indices), "duplicated cells from tSNE analysis")
            print(table(environment$dataset.labels[duplicated.indices]))
        }

        data.path <- environment$PCA.path

        tSNE <- function(perplexity, max_iter) {

            precomputed <- readRDS(data.path)
            PCA <- precomputed$PCA
            rm(precomputed)

            duplicated.indices <- duplicated(t(PCA))
            tSNE <- Rtsne::Rtsne(t(PCA[, !duplicated.indices]), pca = F, initial_dims = nrow(PCA),
                perplexity = perplexity, max_iter = max_iter, verbose = T, whiten = F)$Y

            saveRDS(tSNE, file = file.path(tSNEs.dir, paste(perplexity, max_iter,
                "tSNE.rds", sep = ".")))
        }

        params <- data.frame(expand.grid(perplexity = perplexity, max_iter = max_iter))
        if (nrow(params) > 10) {
            print.message("Warning: possibly too many tSNE parameter combinations [",
                nrow(params), "combinations ]")
            terminate <- readline(prompt = "Terminate? (y/n) ")
            if (terminate != "n") {
                return()
            }
        }

        sopt <- list(mem = mem, time = time, share = TRUE)

        if (local) {
            sjob <- slurm_apply(tSNE, params, nodes = nrow(params),
                                add_objects = c("data.path", "tSNEs.dir"),
                                cpus_per_node = 1, submit = FALSE, slurm_options = sopt)
            local_slurm_array(sjob)
        } else {
            sjob <- slurm_apply(tSNE, params, nodes = nrow(params), cpus_per_node = 1,
                                add_objects = c("data.path", "tSNEs.dir"), submit = TRUE, slurm_options = sopt)
        }
    }

    return(sjob)
}
