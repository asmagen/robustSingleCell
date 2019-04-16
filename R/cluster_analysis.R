cluster <- function(knn.ratio, label, path, data.path, nPCs) {

    cat(paste("Params", "\nknn.ratio = ", knn.ratio, "\nlabel = ", label, "\npath = ",
        path, "\ndata.path = ", data.path, "\n", sep = ""))

    if (length(grep("shuffled", label)) > 0) {
        rep <- strsplit(as.character(label), split = "[.]")[[1]][2]
        file <- file.path(path, "shuffled.PCA", paste("shuffled.PCA.rep", rep, "rds",
            sep = "."))
        if (!file.exists(file))
            print.message("Shuffled PCA file", rep, "does not exist [Change PCA or clustering parameters]")
        precomputed <- readRDS(file)
        pca.perm <- precomputed$pca.perm
        PCA <- t(pca.perm$x)[seq(as.numeric(nPCs)), ]
    } else {
        precomputed <- readRDS(file.path(data.path))
        PCA <- precomputed$PCA
    }

    res <- Rphenograph(t(PCA), k = max(2,floor(knn.ratio * ncol(PCA))))
    community <- res[[2]]
    modularity <- igraph::modularity(community)
    membership <- igraph::membership(community)
    memberships <- community$memberships
    nclusters <- length(unique(membership))

    clustering <- list(modularity = modularity, memberships = memberships, membership = membership,
        knn.ratio = knn.ratio, nclusters = nclusters)

    saveRDS(clustering, file = file.path(path, "clustering", paste(label, knn.ratio,
        "rds", sep = ".")))
}

get.clustering.results <- function(clustering.dir, knn.ratio, shuffledKNN) {

    shuffled.result.files <- list.files(clustering.dir, pattern = paste("shuffled.*.*.rds",
        sep = ""))
    shuffled.result.files
    nshuffled <- length(shuffled.result.files)
    shuffled <- {
    }
    rep <- 1
    shuffled.membership <- list()
    for (rep in seq(nshuffled)) {
        clustering <- readRDS(file.path(clustering.dir, shuffled.result.files[rep]))
        shuffled <- rbind(shuffled, unlist(clustering[c(1, 4, 5)]))
        shuffled.membership[[rep]] <- as.vector(clustering[2])
    }
    shuffled <- as.data.frame(shuffled)

    real.clustering.file <- file.path(clustering.dir, paste("real", knn.ratio, "rds",
        sep = "."))
    if (!file.exists(real.clustering.file))
        cat("\nERROR: Couldn't find real clustering file knn.ratio =", knn.ratio,
            "\nCHECK FOR FAILED JOBS\n\n\n")
    clustering <- readRDS(real.clustering.file)
    membership <- clustering$membership
    nFinalClusters <- length(unique(membership))
    memberships <- clustering$memberships
    modularity <- clustering$modularity
    nclusters <- apply(memberships, 1, function(v) length(unique(v)))
    shuffled.indices <- order(abs(shuffled$nclusters - nFinalClusters))[seq(shuffledKNN)]
    shuffled.matches <- shuffled[shuffled.indices, ]
    shuffled.membership <- shuffled.membership[shuffled.indices]

    cat("knn.ratio = ", knn.ratio, " nclusters = ", nFinalClusters, " (nClustShuffled = ",
        ifelse(stats::var(shuffled.matches$nclusters) > 0, paste(min(shuffled.matches$nclusters),
            "~", max(shuffled.matches$nclusters), sep = ""), shuffled.matches$nclusters[1]),
        ")\n", sep = "")

    mean.shuffled <- mean(shuffled.matches$modularity)
    max.shuffled <- max(shuffled.matches$modularity)
    cat("mdlrty=", round(modularity, 2), " mean.shfl=", round(mean.shuffled, 2),
        " max.shfl=", round(max.shuffled, 2), " mdlrty/mean.shfl=", round(modularity/mean.shuffled,
            2), " mdlrty/max.shfl=", round(modularity/max.shuffled, 2), "\n", sep = "")

    return(list(membership = membership, memberships = memberships, modularity = modularity,
        shuffled = shuffled.matches, shuffled.membership = shuffled.membership, shuffled.modularity = shuffled.matches$modularity))
}

#' Distributed Clustering Analysis
#'
#' Perform clustering analysis for a range of hyperparameter (KNN Ratios) values and assess clustering quality relative to simulation analysis of shuffled data.
#'
#' @param environment \code{environment} object
#' @param knn.ratios range of KNN parameters to scan (corresponding to different resolutions)
#' @param nShuffleRuns number of shuffled clustering analyses to perform per KNN threshold
#' @param shuffledKNN number of closest KNN shuffled analyses to include in background clustering quality computation
#' @param loadPreviousKnn whether to load previous analysis results
#' @param rerun whether to rerun the analysis rather than load from cache
#' @param deleteCache whether to delete cache files
#' @param mem HPC memory
#' @param time HPC time
#' @param plot whether to plot the clustering qualities compared to shuffled
#' @param local whether to run jobs locally rather than using distributed slurm system
#' @return \code{environment} parameter containing clustering assignment and provisional cluster names
#' @import cccd
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' }
cluster.analysis <- function(environment, knn.ratios = c(0.01, 0.05, 0.1), nShuffleRuns = 10,
    shuffledKNN = 10, loadPreviousKnn = T, rerun = F, deleteCache = F, mem = "4GB",
    time = "0:15:00", plot = T, local = F) {

    if (check_not_slurm("cluster.analysis")) {
        return(environment)
    }

    clustering.dir <- file.path(environment$res.data.path, "clustering")
    nresults = length(list.files(clustering.dir, pattern = paste("*.rds", sep = "")))
    shuffled.result.files <- list.files(clustering.dir, pattern = paste("shuffled.*.*.rds",
        sep = ""))
    shuffled.result.files

    extended.knn.ratios <- unique(c(knn.ratios, seq(max(knn.ratios), max(knn.ratios) *
        1.25, max(knn.ratios) - knn.ratios[order(knn.ratios)][length(knn.ratios) -
        1])))
    path <- environment$res.data.path
    data.path <- environment$PCA.path
    params <- rbind(data.frame(knn.ratio = knn.ratios, label = "real", path, data.path,
        nPCs = nrow(environment$PCA)), data.frame(expand.grid(knn.ratio = extended.knn.ratios,
        label = paste("shuffled", seq(nShuffleRuns), sep = ".")), path, data.path,
        nPCs = nrow(environment$PCA)))
    cache <- file.path(environment$res.data.path, "clustering.rds")
    sjob <- NA

    if (rerun || !dir.exists(clustering.dir) || nresults == 0) {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"))
        on.exit(end(t))
        if (deleteCache)
            unlink(clustering.dir, recursive = T, force = T)
        dir.create(clustering.dir, showWarnings = F)

        print.message("knn.ratios:")
        print(knn.ratios)
        print.message("# params:", nrow(params))
        print.message("Head params:")
        print(utils::head(params))
        print.message("Tail params:")
        print(utils::tail(params))

        sopt <- list(mem = mem, time = time, share = TRUE)

        if (local) {
            sjob <- sjob <- slurm_apply(cluster, params, nodes = nrow(params), cpus_per_node = 1,
            submit = FALSE, slurm_options = sopt)
            local_slurm_array(sjob)
        } else {
            sjob <- slurm_apply(cluster, params, nodes = nrow(params), cpus_per_node = 1,
            submit = TRUE, slurm_options = sopt)
        }

        tryCatch({
            get_slurm_out(sjob)
            get_slurm_out(sjob)
            get_slurm_out(sjob)
        }, error = function(v) v)
        end(t)
    }

    if (loadPreviousKnn && file.exists(cache)) {
        print.message("Loading precomputed")
        clustering <- readRDS(cache)
    } else {
        t <- start(file.path(environment$work.path, "tracking"), name = "KNN.stats",
            split = T)
        on.exit(end(t), add = T)
        nResults <- length(list.files(clustering.dir, pattern = paste("*.rds", sep = "")))
        if (nResults < nrow(params)) {
            cat("\nERROR: Found just", nResults, "shuffled clusterings instead of",
                nrow(params), "\nCHECK FOR FAILED JOBS\n\n\n")
            terminate <- readline(prompt = "Terminate? (y/n) ")
            if (terminate != "n") {
                tryCatch({
                  cleanup_files(sjob)
                  cleanup_files(sjob)
                  cleanup_files(sjob)
                }, error = function(v) v)
                stop()
            }
        }
        tryCatch({
            cleanup_files(sjob)
            cleanup_files(sjob)
            cleanup_files(sjob)
        }, error = function(v) v)

        plot.stats <- list()
        clusters.aggregate <- {
        }

        for (knn.ratio in knn.ratios[length(knn.ratios):1]) {
            tryCatch({
                result <- get.clustering.results(clustering.dir, knn.ratio, shuffledKNN)
                nClusters <- paste(length(unique(result$membership)))
                stats <- plot.stats[[nClusters]]
                if (is.null(stats))
                  stats <- {
                  }
                stats <- rbind(stats, data.frame(knn.ratio = knn.ratio, Type = "Original data",
                  Modularity = result$modularity))
                stats <- rbind(stats, data.frame(knn.ratio = knn.ratio, Type = "Shuffled data",
                  Modularity = result$shuffled.modularity))
                plot.stats[[nClusters]] <- stats
                clusters.aggregate <- cbind(clusters.aggregate, result$membership)
            }, error = function(v) v)
        }

        if (plot) {
            ncluster <- names(plot.stats)[1]
            new.plot.stats <- {
            }
            fold.data <- {
            }
            for (ncluster in as.character(sort(as.numeric(names(plot.stats))))) {
                summary(plot.stats[[ncluster]])
                data <- plot.stats[[ncluster]]
                original <- data$Modularity[data$Type == "Original data"]
                shuffled <- data$Modularity[data$Type == "Shuffled data"]
                new.plot.stats <- rbind(new.plot.stats, data.frame(nClusters = ncluster,
                  Type = "Original Data", Modularity = mean(original), sd = stats::sd(original)))
                new.plot.stats <- rbind(new.plot.stats, data.frame(nClusters = ncluster,
                  Type = "Shuffled Data", Modularity = mean(shuffled), sd = stats::sd(shuffled)))
                fold.data <- rbind(fold.data, data.frame(nClusters = ncluster, Type = "Original Data",
                  Modularity = mean(original), Fold = mean(original)/mean(shuffled)))
            }

            grDevices::pdf(file.path(environment$work.path, "Clustering.modularity.pdf"),
                width = 10)
            print(ggplot(new.plot.stats, aes(x = nClusters, y = Modularity, fill = Type)) +
                geom_bar(stat = "identity", color = "black", position = position_dodge()) +
                geom_errorbar(aes(ymin = Modularity - sd, ymax = Modularity + sd),
                  width = 0.2, position = position_dodge(0.9)) + labs(title = "Clustering modularity analysis",
                x = "Number of Clusters", y = "Modularity") + theme_classic(base_size = 25) +
                scale_fill_manual(values = c("#999999", "#E69F00")) + geom_text(data = fold.data,
                aes(nClusters, max(Modularity) * 1.05, label = sprintf("%2.1f", Fold)),
                size = 7))
            print(ggplot(new.plot.stats, aes(x = nClusters, y = Modularity, fill = Type)) +
                geom_bar(stat = "identity", color = "black", position = position_dodge()) +
                geom_errorbar(aes(ymin = Modularity - sd, ymax = Modularity + sd),
                  width = 0.2, position = position_dodge(0.9)) + labs(title = "Clustering modularity analysis",
                x = "Number of Clusters", y = "Modularity") + theme_classic(base_size = 25) +
                scale_fill_manual(values = c("#999999", "#E69F00")) + geom_text(data = fold.data,
                aes(nClusters, max(Modularity) * 1.05, label = sprintf("%2.2f", Fold)),
                size = 7))
            grDevices::dev.off()
        }

        knn.choose <- as.numeric(readline(prompt = "Select KNN ratio: "))
        print.message("knn.choose =", knn.choose)
        clustering.results <- get.clustering.results(clustering.dir, knn.choose,
            nShuffleRuns)

        clustering <- {
        }
        clustering$knn.choose <- knn.choose
        clustering$membership <- as.vector(clustering.results$memberships[nrow(clustering.results$memberships), ])
        clustering$nclusters <- length(unique(clustering$membership))
        print.message("# clusters:", clustering$nclusters)
        print.message("Length:", length(clustering$membership))
        clustering$shuffled <- clustering.results$shuffled
        clustering$shuffled.membership <- clustering.results$shuffled.membership

        files <- c(list.files(environment$work.path, pattern = "*.pdf", full.names = T),
            list.files(environment$work.path, pattern = "*.csv", full.names = T),
            list.files(file.path(environment$work.path, "diff.exp", "main"), pattern = "*.csv",
                full.names = T), list.files(file.path(environment$work.path, "Seurat"),
                pattern = "*.pdf", full.names = T), list.files(file.path(environment$work.path,
                "tracking"), pattern = "*.txt", full.names = T))
        dir <- file.path(environment$work.path, format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"))
        dir.create(dir, showWarnings = F)
        if (length(files) > 0) {
            file.copy(files, dir)
            file.remove(files)
            file.remove(list.files(environment$res.data.path, pattern = "*.diff.exp.rds",
                full.names = T))
        }
        unlink(file.path(environment$work.path, "diff.exp"), recursive = T, force = T)
        saveRDS(clustering, file = cache)
    }

    environment$clustering <- clustering
    environment$cluster.names <- environment$clustering$membership

    print.message("Membership table:")
    print(table(clustering$membership))

    return(environment)
}
