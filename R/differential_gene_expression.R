getDE.limma <- function(Y, group, filter = T) {
    design <- stats::model.matrix(~group)
    my.lm <- limma::lmFit(Y, design = design)
    my.lm <- limma::eBayes(my.lm)
    topTable <- limma::topTable(my.lm, number = Inf, coef = seq(2, ncol(design)))
    if (filter)
        topTable <- topTable[topTable$adj.P.Val < 0.1, ]
    topTable <- topTable[order(topTable$adj.P.Val, decreasing = F), ]
    colnames(topTable)[colnames(topTable) == "P.Value"] <- "PValue"
    colnames(topTable)[colnames(topTable) == "adj.P.Val"] <- "QValue"

    return(topTable)
}

run.diff.expression <- function(environment, clustering, min.fold, quantile,
    label, rerun = F, contrast = "all", contrast.groups = NA) {

    get.diff.exp.stats <- function(id) {
        t <- Sys.time()
        if (is.na(id)) {
            label <- "real"
            clustering <- as.vector(membership)
        } else {
            label <- "shuffled"
            clustering <- as.vector(unlist(shuffled.membership[[id]]))
        }

        fold <- function(v, group) {
            exp(mean(v[group == T]) - mean(v[group == F]))
        }

        cat(label, id, "\n")

        stats <- {
        }
        for (cluster in sort(unique(clustering))) {

            cat("cluster", cluster, "\n")
            tt <- Sys.time()
            group <- factor(clustering == cluster)

            stats <- rbind(stats, data.frame(cluster = cluster, gene = rownames(matrix),
                fold = apply(matrix, 1, fold, group)))

            print(Sys.time() - tt)
        }

        stats <- cbind(label, stats)

        print(Sys.time() - t)
        return(stats)
    }

    summarize.diff.exp.stats <- function(job.portion, min.fold, quantile) {

        real.stats <- stats[stats$label == "real", ]
        print(utils::head(real.stats))
        shuffled.stats <- stats[stats$label == "shuffled", ]
        print(utils::head(shuffled.stats))
        genes <- genes[job.portions == job.portion]
        results <- {
        }
        empirical.diff <- {
        }
        for (gene in genes) {
            real <- real.stats[real.stats$gene == gene, ]
            real <- real[real$fold >= min.fold, ]
            if (nrow(real) == 0)
                next
            shuffled <- shuffled.stats$fold[shuffled.stats$gene == gene]
            # quantile(shuffled,seq(0.01,1,0.01));quantile(shuffled,0.95);real
            significant <- real[real$fold >= quantile(shuffled, quantile), ]
            empirical.diff <- rbind(empirical.diff, significant)
        }
        rownames(empirical.diff) <- NULL
        empirical.diff <- empirical.diff[order(empirical.diff$fold, decreasing = T),
            ]
        print(utils::head(empirical.diff))
        return(empirical.diff)
    }

    cache <- file.path(environment$res.data.path, paste(label, contrast, "diff.exp.RData",
        sep = "."))

    if (!rerun && file.exists(cache)) {
        print.message("Loading precomputed")
        load(cache)
    } else {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"), split = T)
        membership <- as.vector(clustering$membership)

        empirical.diff <- NA

        diff.exp.dir <- file.path(environment$work.path, "diff.exp", label,
            contrast)
        unlink(diff.exp.dir, recursive = T, force = T)
        dir.create(diff.exp.dir, recursive = T)
        limma.diff <- {
        }

        if (contrast == "all") {
            contrast.groups <- rep(1, length(membership))
        } else if (contrast == "datasets") {
            contrast.groups <- environment$datasets
        }

        print.message("contrast =", contrast)

        print(table(contrast.groups))
        contrast.group <- unique(contrast.groups)[1]
        table(membership, contrast.groups)
        for (contrast.group in unique(contrast.groups)) {
            print.message("contrast.group =", contrast.group)
            contrast.group.indices <- contrast.groups == contrast.group
            sum(contrast.group.indices)
            contrast.group.clusters <- sort(unique(membership[contrast.group.indices]))
            cluster <- contrast.group.clusters[1]
            for (cluster in contrast.group.clusters) {
                print.message("cluster =", cluster)

                group <- membership[contrast.group.indices] == cluster
                group <- factor(group)

                diff.exp <- getDE.limma(Y = environment$normalized[, contrast.group.indices],
                  group = group, filter = F)
                diff.exp <- diff.exp[order(diff.exp$logFC, decreasing = T),
                  ]
                diff.exp <- data.frame(gene = rownames(diff.exp), logFC = diff.exp$logFC,
                  fold = exp(diff.exp$logFC), QValue = diff.exp$QValue, PValue = diff.exp$PValue,
                  AveExpr = diff.exp$AveExpr)
                utils::write.csv(diff.exp, file = file.path(diff.exp.dir, paste("cluster",
                  cluster, "csv", sep = ".")))
                markers.diff <- diff.exp[diff.exp$gene %in% environment$marker.genes,
                  ]
                utils::write.csv(markers.diff, file = file.path(diff.exp.dir, paste("markers.cluster",
                  cluster, "csv", sep = ".")))
                limma.diff <- rbind(limma.diff, data.frame(contrast.group = contrast.group,
                  cluster = cluster, diff.exp[, c(1, 3, 4)]))
            }
        }
        rownames(limma.diff) <- NULL
        limma.all <- limma.diff
        limma.diff <- limma.diff[limma.diff$QValue <= (1 - quantile) & limma.diff$fold >=
            min.fold, ]
        final.diff <- limma.diff <- limma.diff[order(limma.diff$fold, decreasing = T),
            ]
        print.message("head(limma.diff):")
        print(utils::head(limma.diff))
        utils::write.csv(limma.diff, file = file.path(diff.exp.dir, "limma.diff.csv"))
        utils::write.csv(final.diff, file = file.path(diff.exp.dir, "final.diff.csv"))

        save(final.diff, limma.diff, empirical.diff, limma.all, file = cache)
        end(t)
    }

    return(final.diff)
}
