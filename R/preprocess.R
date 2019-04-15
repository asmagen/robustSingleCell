#' Identify Highly Variable Genes
#'
#' Get highly variable genes by Heteroscedasticity controlled binning of gene expression measurements within each dataset separately.
#'
#' @param environment \code{environment} object
#' @param min.mean minimum mean expression per gene
#' @param min.frac.cells minimum fraction of cells expressing each gene
#' @param min.dispersion.scaled minimum dispersion value
#' @param rerun whether to rerun the analysis or load from cache
#' @return \code{environment} parameter containing highly variable genes selection
#' @export
#' @examples
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1)
get.variable.genes <- function(environment, min.mean = 0.05, min.frac.cells = 0,
    min.dispersion.scaled = 1, rerun = F) {

    cache <- file.path(environment$baseline.data.path, "HVG.rds")

    if (!rerun & file.exists(cache)) {
        print.message("Loading precomputed")
        HVG <- readRDS(cache)
    } else {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"))
        on.exit(end(t))
        normalized <- environment$normalized

        datasets <- environment$dataset_ids
        table(datasets)
        dataset <- sort(unique(datasets))[1]
        dataset
        merged <- {
        }
        for (dataset in sort(unique(datasets))) {
            print.message(dataset)
            filtered.normalized <- normalized[, datasets == dataset]
            means <- apply(filtered.normalized, 1, function(v) log(mean(exp(v) -
                1) + 1))
            frac.cells <- rowSums(filtered.normalized > 0)/ncol(filtered.normalized)
            vars <- apply(filtered.normalized, 1, function(v) log(stats::var(exp(v) -
                1) + 1))
            dispersion <- apply(filtered.normalized, 1, function(v) log(stats::var(exp(v) -
                1)/mean(exp(v) - 1)))
            dispersion[is.na(x = dispersion)] <- 0
            means[is.na(x = means)] <- 0

            grDevices::pdf(file.path(environment$baseline.work.path, paste(dataset,
                "VariableGenes.pdf", sep = ".")))
            plot.data <- data.frame(gene = names(means), means, dispersion)  #head(plot.data)
            print(ggplot(plot.data, aes(x = means, y = dispersion, label = gene)) +
                geom_text(check_overlap = TRUE, size = 2))
            graphics::smoothScatter(means, dispersion)
            grDevices::dev.off()

            num.bin <- 20
            bins <- cut(x = means, breaks = num.bin)
            names(x = bins) <- names(x = means)
            mean_y <- tapply(dispersion, bins, mean)
            sd_y <- tapply(dispersion, bins, stats::sd)
            dispersion.scaled <- (dispersion - mean_y[as.numeric(x = bins)])/sd_y[as.numeric(x = bins)]
            dispersion.scaled[is.na(x = dispersion.scaled)] <- 0
            names(x = dispersion.scaled) <- names(x = means)

            criterion <- means >= min.mean & frac.cells >= min.frac.cells & dispersion.scaled >=
                min.dispersion.scaled
            HVG <- names(means)[criterion]
            print.message("# qualifying genes", length(HVG))
            print.message("Qualifying markers")
            print(environment$marker.genes[environment$marker.genes %in% HVG])
            print.message("NOT qualifying markers")
            print(environment$marker.genes[!environment$marker.genes %in% HVG])
            genes <- c(get.ribo.genes(rownames(normalized)), get.mito.genes(rownames(normalized)))
            print.message("Qualifying Ribosomal & Mitochondrial")
            print(genes[genes %in% HVG])
            utils::write.csv(HVG, file = file.path(environment$baseline.work.path,
                paste(dataset, "VariableGenes.csv", sep = ".")))
            merged <- unique(c(merged, HVG))

        }
        HVG <- merged

        print.message("Overall # qualifying genes", length(HVG))
        print.message("Overall qualifying markers")
        print(environment$marker.genes[environment$marker.genes %in% HVG])
        print.message("Overall NOT qualifying markers")
        print(environment$marker.genes[!environment$marker.genes %in% HVG])
        genes <- c(get.ribo.genes(rownames(normalized)), get.mito.genes(rownames(normalized)))
        print.message("Overall qualifying Ribosomal & Mitochondrial")
        print(genes[genes %in% HVG])

        saveRDS(HVG, file = cache)
    }

    environment$HVG <- HVG

    cat("# highly variable genes = ", length(environment$HVG), "\n", sep = "")

    return(environment)
}

nUMIs <- function(environment) {
    return(colSums(environment$counts))
}

nGenes <- function() {
    return(colSums(environment$counts > 0))
}

#' Add Confounder Variables
#'
#' Add confounder variables' activation level per cell to \code{environment} object.
#'
#' @param environment \code{environment} object
#' @param ... vectors containing activation levels of confounding variables
#' @return \code{environment} parameter containing added confounder variable
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- add.confounder.variables(LCMV1, ribosomal.score = ribosomal.score(LCMV1))
#' }
add.confounder.variables <- function(environment, ...) {
    environment$confounders <- data.frame(environment$confounders, data.frame(...))
    print(utils::head(environment$confounders))
    return(environment)
}

#' Compute Ribosomal Score
#'
#' Compute the activation level of ribosomal genes.
#'
#' @param environment \code{environment} object
#' @param control whether to subtract the score defined by technically similar genes
#' @param knn number of nearest neighbor
#' @return a vector of ribosomal genes activation score
#' @export
#' @examples
#' LCMV1 <- setup_LCMV_example()
#' ribosomal.score <- ribosomal.score(LCMV1)
ribosomal.score <- function(environment, control = T, knn = 10) {
    t <- start(file.path(environment$work.path, "tracking"))
    on.exit(end(t))
    genes <- get.ribo.genes(environment$genes)
    print.message("Using genes:")
    print(genes)
    if (control) {
        score <- controlled.mean.score(environment, genes, knn)
    } else {
        score <- colMeans(environment$normalized[genes, ])
    }
    score
}

get.ribo.genes <- function(genes) {
    return(genes[c(grep("^Rpl", genes, ignore.case = T), grep("^Rps", genes, ignore.case = T))])
}

#' Compute Mitochondrial Score
#'
#' Compute the activation level of mitochondrial genes.
#'
#' @param environment \code{environment} object
#' @param control whether to subtract the score defined by technically similar genes
#' @param knn number of nearest neighbor
#' @return a vector of mitochondrial genes activation score
#' @export
#' @examples
#' LCMV1 <- setup_LCMV_example()
#' mitochondrial.score <- mitochondrial.score(LCMV1)
mitochondrial.score <- function(environment, control = F, knn = 10) {
    # browser()
    t <- start(file.path(environment$work.path, "tracking"))
    on.exit(end(t))
    genes <- get.mito.genes(environment$genes)
    print.message("Using genes:")
    print(genes)
    if (control) {
        score <- controlled.mean.score(environment, genes, knn)
    } else {
        score <- Matrix::colMeans(environment$normalized[genes %in% rownames(environment$normalized),
            ])
    }
    return(score)
}

get.mito.genes <- function(genes) {
    return(genes[grep("^Mt-", genes, ignore.case = T)])
}

#' Compute Cell Cycle Score
#'
#' Compute the activation of cell cycle genes defined in Kowalczyk, M. S. et al.
#'
#' @param environment \code{environment} object
#' @param knn number of nearest neighbor
#' @param cc.genes.path optional; path to defined cell cycle genes. Default uses gene sets defined in Kowalczyk, M. S. et al. Single-cell RNA-seq reveals changes in cell cycle and differentiation programs upon aging of hematopoietic stem cells. Genome Res 25, 1860-1872, doi:10.1101/gr.192237.115 (2015).
#' @return a matrix of cell cycle genes activation scores (S, G2M and aggregated S/G2M scores, separately)
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' cell.cycle.score <- cell.cycle.score(LCMV1)
#' }
cell.cycle.score <- function(environment, knn = 10, cc.genes.path = NA) {
    t <- start(file.path(environment$work.path, "tracking"))
    on.exit(end(t))

    if (is.na(cc.genes.path)) {
        cc.genes <- capwords(cell_cycle_genes)
    } else {
        cc.genes <- capwords(readLines(cc.genes.path))
    }

    s.genes <- cc.genes[1:43]
    s.genes <- s.genes[s.genes %in% capwords(environment$genes)]
    print(s.genes)
    g2m.genes <- cc.genes[44:98]
    g2m.genes <- g2m.genes[g2m.genes %in% capwords(environment$genes)]
    print(g2m.genes)

    s.score <- controlled.mean.score(environment, s.genes, knn)
    g2m.score <- controlled.mean.score(environment, g2m.genes, knn)
    cell.cycle.score <- controlled.mean.score(environment, c(s.genes, g2m.genes),
        knn)

    print.message("# s.score > 0:", sum(s.score > 0), "fraction", sum(s.score > 0)/length(s.score))
    print.message("# g2m.score > 0:", sum(g2m.score > 0), "fraction", sum(g2m.score >
        0)/length(g2m.score))

    return(data.frame(S.stage = s.score, G2M.stage = g2m.score, aggregate_S_G2M.stage = cell.cycle.score))
}

#' Compute Controlled Activation Score
#'
#' Compute mean gene signatures activation scores while controlling for technically similar genes.
#'
#' @param environment \code{environment} object
#' @param genes gene list upon which to calculate gene signature activate
#' @param knn number of nearest neighbors
#' @param exclude.missing.genes whether to exclude genes with missing values
#' @param constrain.cell.universe binary vector indicating in which subset of cells to calculate the gene signature activation. Default is all cells.
#' @return gene signature activation scores per cell
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' exhaustion_markers <- c('Pdcd1', 'Cd244', 'Havcr2', 'Ctla4', 'Cd160', 'Lag3',
#' 'Tigit', 'Cd96')
#' Exhaustion <- controlled.mean.score(LCMV1, exhaustion_markers)
#' }
controlled.mean.score <- function(environment, genes, knn = 10, exclude.missing.genes = T,
    constrain.cell.universe = NA) {
    # similarly to
    # http://science.sciencemag.org/content/sci/suppl/2016/04/07/352.6282.189.DC1/Tirosh.SM.pdf/Seurat
    # to reduce association with library size or other technical

    if (is.na(constrain.cell.universe))
        constrain.cell.universe <- rep(T, ncol(environment$normalized))
    if (knn > 0) {
        genes <- rownames(environment$normalized)[match(capwords(genes), capwords(rownames(environment$normalized)))]
        nExclude <- sum(is.na(genes))
        if (nExclude > 0) {
            if (exclude.missing.genes) {
                print.message("Excluding", nExclude, "out of", length(genes), "genes not found in the dataset")
                genes <- genes[genes %in% capwords(rownames(environment$normalized))]
            } else {
                print.message("Some signature genes are missing in dataset")
                return(NA)
            }
        }

        background.genes <- background.genes(environment, foreground.genes = genes,
            knn)

        return(Matrix::colMeans(environment$normalized[genes, constrain.cell.universe]) -
            Matrix::colMeans(environment$normalized[background.genes, constrain.cell.universe]))
    } else {
        return(Matrix::colMeans(environment$normalized[genes, constrain.cell.universe]))
    }
}

get.technically.similar.genes <- function(environment, knn = 10) {

    t <- start(file.path(environment$work.path, "tracking"))
    on.exit(t)

    cache <- file.path(environment$baseline.data.path, paste(knn, "technical.background.genes.distances.rds",
        sep = "."))

    if (file.exists(cache)) {
        print.message("Loading precomputed")
        precomputed <- readRDS(cache)
        knns <- precomputed$knns
        technical.variables <- precomputed$technical.variables
        rm(precomputed)
    } else {
        print.message("Computing")

        technical.variables <- data.frame(means = rowMeans(environment$normalized),
            vars = apply(environment$normalized, 1, stats::var))
        scaled.technical.variables <- apply(technical.variables, 2, scale)
        rownames(scaled.technical.variables) <- rownames(technical.variables)

        n <- nrow(scaled.technical.variables)
        dist_obj <- stats::dist(scaled.technical.variables)
        knns <- array("", c(length(environment$genes), knn))
        rownames(knns) <- environment$genes
        for (index in seq(length(environment$genes))) {
            gene.dist <- get_dist(dist_obj, index, n, environment$genes)
            knns[index, ] <- names(gene.dist[order(gene.dist)[2:(knn + 1)]])
        }
        saveRDS(list(knns = knns, technical.variables = technical.variables), file = cache)
    }

    return(list(knns = knns, technical.variables = technical.variables))
}

get_dist <- function(dist_obj, i, n, names) {
   stopifnot(i <= n)
   distance <- c(0)
   if (i < n) {
       idx <- c(seq(0, n - i - 1) + (2 * n - i) * (i - 1) / 2 + 1)
       distance <- c(distance, dist_obj[idx])
   }
   if (i > 1) {
       pre <- rev(seq(n - i + 1 , n - 2))
       pre <- Reduce(sum, pre, init = i - 1, accumulate = T)
       distance <- c(dist_obj[pre], distance)
   }
   names(distance) <- names
   return(distance)
}

background.genes <- function(environment, foreground.genes, knn) {

    t <- start(file.path(environment$work.path, "tracking"))
    on.exit(end(t))
    foreground.genes <- foreground.genes[foreground.genes %in% environment$genes]
    technically.similar.genes <- get.technically.similar.genes(environment, knn)
    knns <- technically.similar.genes$knns
    technical.variables <- technically.similar.genes$technical.variables

    background.genes <- unique(setdiff(as.vector(knns[foreground.genes, ]), foreground.genes))
    print.message("Head foreground.genes technical.variables")
    print(utils::head(technical.variables[foreground.genes, ]))
    print.message("Head background.genes technical.variables")
    print(utils::head(technical.variables[background.genes, ]))
    print.message("background.genes")
    print(background.genes)
    return(background.genes)
}

regress.covariates <- function(environment, regress, data, groups, rerun = F, save = F) {

    cache <- file.path(environment$res.data.path, paste(paste(colnames(regress),
        collapse = "+"), "HVG.regressed.covariates.rds", sep = "_"))

    if (!rerun && file.exists(cache)) {
        print.message("Loading precomputed")
        corrected <- readRDS(cache)
    } else {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"))
        on.exit(end(t))
        formula.str <- paste("gene", paste(colnames(regress), collapse = " + "),
            sep = " ~ ")
        formula <- stats::as.formula(formula.str)
        print.message("Regressing:", formula.str)
        print.message("Not regressed matrix")
        corner(data)
        corrected <- data

        for (group in unique(groups)) {
            group.indices <- groups == group
            for (gene in rownames(data)) {
                lm.data <- data.frame(gene = data[gene, group.indices], regress[group.indices,
                  ])
                colnames(lm.data)[2:ncol(lm.data)] <- colnames(regress)
                model <- stats::lm(formula, lm.data)
                corrected[gene, group.indices] <- model$residuals
            }
        }
        if (save)
            saveRDS(corrected, file = cache)
    }

    return(corrected)
}

#' Differential Expression and Figure Generation
#'
#' Summarize the clustering results by conducting differential expression analysis and plotting figures.
#'
#' @param environment \code{environment} object
#' @param perplexity perplexity parameters for tSNE analyses
#' @param max_iter maximum iterations for tSNE
#' @param rerun whether to rerun
#' @param order order in which to plot the clusters
#' @param contrast either 'all' indicating differential expression between one cluster against all others or 'datasets' indicating differential expression analysis comparing one cluster to all other within each dataset separately ('datasets' should be used in pooled analysis for optimal results)
#' @param min.fold minimum fold change for filtering final differentially expressed gene lists
#' @param quantile q-value cutoff for differential expression analysis
#' @param local Whether to run tSNE locally on SLURM
#' @param mem Memory for each job; default 4 GB
#' @param time Time for each job; default 15 minutes
#' @export
#' @examples
#' \donttest{
#' # after running cluster.analysis()
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' summarize(LCMV1)
#' }
summarize <- function(environment, perplexity = seq(10, 30, 10), max_iter = 10000,
    rerun = F, order = NA, contrast = "all", min.fold = 1.5, quantile = 0.95,
    local = F, mem = "4GB", time = "0:15:00") {

    cluster.size <- table(environment$cluster.names)
    if (length(order) == 1 && is.na(order))
        order <- names(cluster.size)[order(cluster.size, decreasing = T)]
    tSNE.job <- run_tSNE(environment, perplexity, max_iter, rerun, local = local,
                         mem = mem, time = time)
    plot_PCA(environment, quantile = 0.05, order)
    plot.cluster.stats(environment, membership = environment$cluster.names, order = order)
    if (length(environment$seurat.cluster.association) > 1)
        tryCatch({
            plot.cluster.stats(environment, membership = environment$seurat.cluster.association,
                label = "Seurat", order = order)
        }, error = function(v) v)

    final.diff <- run.diff.expression(environment, clustering = environment$clustering,
        min.fold, quantile, label = "main", rerun = rerun, contrast = contrast)

    order <- sort(unique(environment$cluster.names))
    plot.heatmaps(environment, diff.exp = final.diff, membership = environment$cluster.names,
        order = order)
    if (length(environment$seurat.cluster.association) > 1)
        tryCatch({
            plot.heatmaps(environment, diff.exp = final.diff, membership = environment$seurat.cluster.association,
                label = "Seurat")
        }, error = function(v) v)
    plot.tSNE(environment, tSNE.job, perplexity, max_iter)
}
