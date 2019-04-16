#' Get Robust Cluster Similarity
#'
#' Use cross-replicate experiment cluster similarity to remove irreproducible clusters.
#'
#' @param environment \code{environment} object
#' @param similarity pearson correlation between clusters' FC vectors defined in assess.cluster.similarity
#' @param min.sd minimum standard deviation for cluster reproducibility assessment
#' @param max.q.val maximum q value for cluster correlation cutoff
#' @param rerun whether to rerun the analysis or load from cache
#' @return filtered cluster similarity matrix
#' @export
#' @import RColorBrewer
#' @import reshape2
#' @import gplots
#' @import dplyr
#' @import ggpubr
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' types = rbind(
#' data.frame(type='Tfh',gene=c('Tcf7','Cxcr5','Bcl6')),
#' data.frame(type='Th1',gene=c('Cxcr6','Ifng','Tbx21')),
#' data.frame(type='Tcmp',gene=c('Ccr7','Bcl2','Tcf7')),
#' data.frame(type='Treg',gene=c('Foxp3','Il2ra')),
#' data.frame(type='Tmem',gene=c('Il7r','Ccr7')),
#' data.frame(type='CD8',gene=c('Cd8a')),
#' data.frame(type='CD4', gene = c("Cd4")),
#' data.frame(type='Cycle',gene=c('Mki67','Top2a','Birc5'))
#' )
#' summarize(LCMV1)
#' cluster_names <- get.cluster.names(LCMV1, types, min.fold = 1.0, max.Qval = 0.01)
#' LCMV1 <- set.cluster.names(LCMV1, names = cluster_names)
#' LCMV2 <- setup_LCMV_example("LCMV2")
#' LCMV2 <- get.variable.genes(LCMV2, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV2 <- PCA(LCMV2)
#' LCMV2 <- cluster.analysis(LCMV2)
#' summarize(LCMV2)
#' cluster_names <- get.cluster.names(LCMV2, types, min.fold = 1.0, max.Qval = 0.01)
#' LCMV2 <- set.cluster.names(LCMV2, names = cluster_names)
#' pooled_env <- setup_pooled_env()
#' pooled_env <- read.preclustered.datasets(pooled_env)
#' pooled_env <- PCA(pooled_env, clear.previously.calculated.clustering = F)
#' summarize(pooled_env, contrast = "datasets")
#' cluster.similarity <- assess.cluster.similarity(pooled_env)
#' similarity <- cluster.similarity$similarity
#' map <- cluster.similarity$map
#' filtered.similarity <- get.robust.cluster.similarity(
#'    pooled_env, similarity, min.sd = qnorm(.9), max.q.val = 0.01, rerun = F)
#' }
get.robust.cluster.similarity <- function(environment, similarity, min.sd = stats::qnorm(0.95),
    max.q.val = 0.01, rerun = F) {

    cache <- file.path(environment$res.data.path, "filtered.cluster.similarity.rds")
    if (!rerun && file.exists(cache)) {
        print.message("Loading precomputed similarity")
        filtered.cluster.similarity <- readRDS(cache)
    } else {
        print.message("Computing")


        match.significance.stats <- {
        }
        origin <- unique(similarity$origin1)[1]
        origins <- sort(unique(c(as.vector(similarity$origin1), as.vector(similarity$origin2))))
        for (origin in origins) {
            filter <- similarity$origin1 == origin & similarity$origin2 == origin &
                similarity$experiment1 != similarity$experiment2
            filtered.similarity <- similarity[filter, ]
            cor.val <- filtered.similarity$similarity
            shapiro.test.p.value <- tryCatch({
                stats::shapiro.test(cor.val)$p.value
            }, error = function(v) return(NA))
            sd.dist <- (cor.val - mean(cor.val))/stats::sd(cor.val)
            match.significance.stats <- rbind(match.significance.stats, data.frame(origin,
                cor.val, sd.dist, shapiro.test.p.value, filtered.similarity))
        }

        name <- paste("get.robust.cluster.similarity", sep = "_")
        work.path <- file.path(environment$work.path, name)
        if (file.exists(work.path)) {
            new.dir <- file.path(environment$work.path, paste(name, format(Sys.time(),
                "%a_%b_%e_%Y__%H_%M_%S"), sep = "---"))
            file.rename(work.path, new.dir)
        }
        dir.create(work.path, showWarnings = T, recursive = T)

        plot.data <- data.frame(correlation = match.significance.stats$cor.val, SD = match.significance.stats$sd.dist,
            origin = match.significance.stats$origin)
        grDevices::pdf(file.path(work.path, "cluster.matching.similarity.histogram.pdf"))
        print(ggplot(plot.data, aes(x = SD, fill = origin)) + geom_density(alpha = 0.5) +
            geom_vline(xintercept = min.sd) + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black"), legend.position = "bottom") + ylab("Density") + theme_classic(base_size = 25))
        print(ggplot(plot.data, aes(x = correlation, fill = origin)) + geom_density(alpha = 0.5) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom") + ylab("Density") + theme_classic(base_size = 25))
        print(ggplot(plot.data, aes(x = SD)) + geom_density(alpha = 0.5) + geom_vline(xintercept = min.sd) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom") + ylab("Density") + theme_classic(base_size = 25))
        print(ggplot(plot.data, aes(x = correlation)) + geom_density(alpha = 0.5) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom") + ylab("Density") + theme_classic(base_size = 25))
        grDevices::dev.off()

        sum(match.significance.stats$sd.dist >= min.sd)/length(match.significance.stats$sd.dist)
        criteria <- match.significance.stats$sd.dist >= min.sd & match.significance.stats$Q.val.1 <=
            max.q.val
        print.message("Found", sum(criteria), "robust populations")
        robust.cluster.similarity <- match.significance.stats[criteria, ]
        robust.cluster.similarity <- robust.cluster.similarity[order(robust.cluster.similarity$origin,
            robust.cluster.similarity$similarity.1, decreasing = T), ]
        summary(robust.cluster.similarity)
        print(robust.cluster.similarity[1:min(nrow(robust.cluster.similarity), 5),
            ])
        utils::write.csv(robust.cluster.similarity, file = file.path(work.path, "robust.cluster.similarity.csv"))

        robust.clusters <- unique(c(robust.cluster.similarity$cluster1, robust.cluster.similarity$cluster2))
        nclusters.overall <- length(unique(c(similarity$cluster1, similarity$cluster2)))
        cat("found", length(robust.clusters), "/", nclusters.overall, "robust.clusters (",
            round(length(robust.clusters)/nclusters.overall, 2), "%)\n")
        filtered.cluster.similarity <- similarity[similarity$cluster1 %in% robust.clusters &
            similarity$cluster2 %in% robust.clusters, ]
        dim(filtered.cluster.similarity)
        summary(filtered.cluster.similarity)

        print(filtered.cluster.similarity[1:min(nrow(filtered.cluster.similarity),
            5), ])
        utils::write.csv(filtered.cluster.similarity, file = file.path(work.path,
            "filtered.cluster.similarity.csv"))

        similarity.summary.df <- data.frame(cluster1 = filtered.cluster.similarity$name1,
            cluster2 = filtered.cluster.similarity$name2, coef = filtered.cluster.similarity$similarity)
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.Cor.FC.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        similarity.summary.df <- data.frame(cluster1 = filtered.cluster.similarity$name1,
            cluster2 = filtered.cluster.similarity$name2, coef = filtered.cluster.similarity$similarity.1)
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.Cor.means.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        similarity.summary.df <- data.frame(cluster1 = filtered.cluster.similarity$name1,
            cluster2 = filtered.cluster.similarity$name2, coef = -scale(filtered.cluster.similarity$ocldist.FC))
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.ocldist.FC.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        similarity.summary.df <- data.frame(cluster1 = filtered.cluster.similarity$name1,
            cluster2 = filtered.cluster.similarity$name2, coef = -scale(filtered.cluster.similarity$ocldist))
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))

        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.ocldist.means.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        saveRDS(filtered.cluster.similarity, file = cache)
    }

    return(filtered.cluster.similarity)
}

pearson.correlation <- function(diff1, diff2) {

    cor <- stats::cor.test(diff1, diff2, method = "pearson")
    return(data.frame(similarity = cor$estimate, significance = cor$p.value))
}

#' Assess Cluster Similarity
#'
#' Assess similarity between pairs of clusters.
#'
#' @param environment \code{environment} object
#' @param diff.exp.file name of differential expression results file
#' @param cluster.similarity.function which similarity function to use (either 'pearson.correlation' or '?') Mamie - there was another similarity function using euclidean distance. Do you know where did it go to? Can you replace the '?' with the name of this other function?
#' @param label name of the similarity measure to use for the results folder
#' @param rerun whether to rerun the analysis or load from cache
#' @return pairwise cluster similarity measures
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' types = rbind(
#' data.frame(type='Tfh',gene=c('Tcf7','Cxcr5','Bcl6')),
#' data.frame(type='Th1',gene=c('Cxcr6','Ifng','Tbx21')),
#' data.frame(type='Tcmp',gene=c('Ccr7','Bcl2','Tcf7')),
#' data.frame(type='Treg',gene=c('Foxp3','Il2ra')),
#' data.frame(type='Tmem',gene=c('Il7r','Ccr7')),
#' data.frame(type='CD8',gene=c('Cd8a')),
#' data.frame(type='CD4', gene = c("Cd4")),
#' data.frame(type='Cycle',gene=c('Mki67','Top2a','Birc5'))
#' )
#' summarize(LCMV1)
#' cluster_names <- get.cluster.names(LCMV1, types, min.fold = 1.0, max.Qval = 0.01)
#' LCMV1 <- set.cluster.names(LCMV1, names = cluster_names)
#' LCMV2 <- setup_LCMV_example("LCMV2")
#' LCMV2 <- get.variable.genes(LCMV2, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV2 <- PCA(LCMV2)
#' LCMV2 <- cluster.analysis(LCMV2)
#' summarize(LCMV2)
#' cluster_names <- get.cluster.names(LCMV2, types, min.fold = 1.0, max.Qval = 0.01)
#' LCMV2 <- set.cluster.names(LCMV2, names = cluster_names)
#' pooled_env <- setup_pooled_env()
#' pooled_env <- read.preclustered.datasets(pooled_env)
#' pooled_env <- PCA(pooled_env, clear.previously.calculated.clustering = F)
#' summarize(pooled_env, contrast = "datasets")
#' cluster.similarity <- assess.cluster.similarity(pooled_env)
#' }
assess.cluster.similarity <- function(environment, diff.exp.file = "main.datasets.diff.exp.rds",
    cluster.similarity.function = pearson.correlation, label = "pearson", rerun = F) {

    # if (check_not_slurm("assess.cluster.similarity")) {
    #     return(environment)
    # }
    cache <- file.path(environment$res.data.path, paste(label, "cluster.similarity.rds",
        sep = "."))
    if (!rerun && file.exists(cache)) {
        print.message("Loading precomputed similarity")
        precomputed <- readRDS(cache)
        similarity <- precomputed$similarity
        map <- precomputed$map
        rm(precomputed)
    } else {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"), split = F)
        on.exit(end(t))
        name <- paste("assess.cluster.similarity", label, sep = "_")
        work.path <- file.path(environment$work.path, name)
        if (file.exists(work.path)) {
            new.dir <- file.path(environment$work.path, paste(name, format(Sys.time(),
                "%a_%b_%e_%Y__%H_%M_%S"), sep = "---"))
            file.rename(work.path, new.dir)
        }
        dir.create(work.path, showWarnings = T, recursive = T)

        precomputed <- readRDS(file.path(environment$work.path, "data", diff.exp.file))
        limma.all <- precomputed$limma.all
        final.diff <- limma.all
        utils::head(final.diff)
        configs <- apply(final.diff, 1, function(v) paste(c(v[1], v[2]), collapse = "_"))
        membership <- environment$clustering$membership
        original.membership <- environment$original.clustering
        nclusters <- length(unique(membership))
        clusters <- sort(unique(membership))
        map <- data.frame(t(sapply(clusters, function(c) c(unique(environment$cluster.names[membership ==
            c][1]), unique(environment$origins[membership == c][1]), unique(environment$experiments[membership ==
            c][1]), unique(environment$dataset.labels[membership == c][1]), unique(original.membership[membership ==
            c][1]), c))))
        colnames(map) <- c("name", "origin", "experiments", "sample", "original.membership",
            "membership")
        utils::head(map)

        cluster.descriptors <- {
        }
        for (cluster.index in seq(nclusters)) {
            cluster <- clusters[cluster.index]
            indices <- membership == cluster
            diff.exp.indices <- configs == paste(environment$dataset.labels[indices][1],
                membership[indices][1], sep = "_ ") | configs == paste(environment$dataset.labels[indices][1],
                membership[indices][1], sep = "_")
            diff <- final.diff[diff.exp.indices, ]
            rownames(diff) <- diff$gene
            if (is.null(cluster.descriptors)) {
                cluster.descriptors <- data.frame(diff$fold)
                rownames(cluster.descriptors) <- diff$gene
            } else {
                cluster.descriptors <- cbind(cluster.descriptors, diff$fold[match(rownames(cluster.descriptors),
                  diff$gene)])
            }
        }
        colnames(cluster.descriptors) <- map[match(seq(nclusters), map$membership),
            1]
        utils::head(cluster.descriptors)

        grDevices::pdf(file.path(environment$work.path, paste("hclust.cluster.descriptors.pdf",
            sep = "_")), width = 10, height = 10)
        hc.dist <- stats::hclust(stats::dist(t(cluster.descriptors)))
        plot(hc.dist)
        grDevices::dev.off()

        # Cluster similarity
        similarity <- {
        }
        for (cluster.index1 in seq(nclusters - 1)) {
            cluster1 <- clusters[cluster.index1]
            indices <- membership == cluster1
            diff.exp.indices1 <- configs == paste(environment$dataset.labels[indices][1],
                membership[indices][1], sep = "_ ") | configs == paste(environment$dataset.labels[indices][1],
                membership[indices][1], sep = "_")
            diff1 <- final.diff[diff.exp.indices1, ]
            rownames(diff1) <- diff1$gene
            measurements1 <- environment$normalized[, membership == cluster1]
            measurements1.means <- rowMeans(measurements1)
            for (cluster.index2 in (cluster.index1 + 1):nclusters) {
                cluster2 <- clusters[cluster.index2]
                indices <- membership == cluster2
                diff.exp.indices2 <- configs == paste(environment$dataset.labels[indices][1],
                  membership[indices][1], sep = "_ ") | configs == paste(environment$dataset.labels[indices][1],
                  membership[indices][1], sep = "_")
                diff2 <- final.diff[diff.exp.indices2, ]
                rownames(diff2) <- diff2$gene
                matches <- match(diff1$gene, diff2$gene)
                diff2 <- diff2[matches, ]
                if (any(rownames(diff1) != rownames(diff2)))
                  stop("\n\n\n\nGene mismatch\n\n\n\n")
                measurements2 <- environment$normalized[, membership == cluster2]
                measurements2.means <- rowMeans(measurements2)
                ocldist <- as.vector(stats::dist(rbind(measurements1.means, measurements2.means)))
                ocldist.FC <- as.vector(stats::dist(rbind(diff1$fold, diff2$fold)))
                over1 <- diff1[diff1$QValue <= 0.05 & diff1$fold >= 1, ]
                over2 <- diff2[diff2$QValue <= 0.05 & diff2$fold >= 1, ]
                intersect <- intersect(over1$gene, over2$gene)
                sort(intersect)
                union <- union(over1$gene, over2$gene)
                jaccard <- length(intersect)/length(union)
                over1 <- diff1[diff1$QValue <= 0.05 & diff1$fold >= 1.5, ]
                over2 <- diff2[diff2$QValue <= 0.05 & diff2$fold >= 1.5, ]
                intersect <- intersect(over1$gene, over2$gene)
                sort(intersect)
                union <- union(over1$gene, over2$gene)
                jaccard2 <- length(intersect)/length(union)
                diff1.compact <- diff1[diff1$fold > 1.05 | diff1$fold < (1/1.05),
                  ]
                diff2.compact <- diff2[diff2$fold > 1.05 | diff2$fold < (1/1.05),
                  ]
                diff.compact.genes <- union(diff1.compact$gene, diff2.compact$gene)
                diff1.compact.strict <- diff1[diff1$fold > 1.15 | diff1$fold < (1/1.15),
                  ]
                diff2.compact.strict <- diff2[diff2$fold > 1.15 | diff2$fold < (1/1.15),
                  ]
                diff.compact.strict.genes <- union(diff1.compact.strict$gene, diff2.compact.strict$gene)
                similarity <- rbind(similarity, data.frame(cluster1, cluster2, cluster.similarity.function(diff1$fold,
                  diff2$fold), cluster.similarity.function(measurements1.means, measurements2.means),
                  ocldist, ocldist.FC, jaccard, jaccard2))
            }
        }
        utils::head(similarity)
        dim(similarity)

        # Cluster mapping
        similarity$Q.val <- stats::p.adjust(similarity$significance, method = "BH")
        similarity$Q.val.1 <- stats::p.adjust(similarity$significance.1, method = "BH")

        similarity <- cbind(map[match(similarity$cluster1, map$membership), 1:5],
            map[match(similarity$cluster2, map$membership), 1:5], similarity)
        utils::head(similarity)
        colnames(similarity)[1:10] <- c("name1", "origin1", "experiment1", "sample1",
            "original.membership1", "name2", "origin2", "experiment2", "sample2",
            "original.membership2")
        similarity <- similarity[order(similarity$similarity, decreasing = T), ]
        utils::write.csv(similarity, row.names = F, file = file.path(work.path, paste(label,
            "cluster.similarity.csv", sep = "_")))

        similarity.summary.df <- data.frame(cluster1 = similarity$name1, cluster2 = similarity$name2,
            coef = similarity$similarity)
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.Cor.FC.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        similarity.summary.df <- data.frame(cluster1 = similarity$name1, cluster2 = similarity$name2,
            coef = similarity$similarity.1)
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.Cor.means.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        similarity.summary.df <- data.frame(cluster1 = similarity$name1, cluster2 = similarity$name2,
            coef = -scale(similarity$ocldist.FC))
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.ocldist.FC.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        similarity.summary.df <- data.frame(cluster1 = similarity$name1, cluster2 = similarity$name2,
            coef = -scale(similarity$ocldist))
        mirror <- similarity.summary.df
        cluster1 <- mirror$cluster1
        cluster2 <- mirror$cluster2
        mirror$cluster1 <- cluster2
        mirror$cluster2 <- cluster1
        similarity.summary.df <- rbind(similarity.summary.df, mirror)
        utils::head(similarity.summary.df)
        similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2, value.var = "coef")
        similarity.matrix[1:10, 1:10]
        hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
        grDevices::pdf(file.path(work.path, paste("filtered.cluster.similarity.heatmap.ocldist.means.pdf",
            sep = "_")), width = 20, height = 20)
        colors <- rev(brewer.pal(5, "PuOr"))
        color.palette <- grDevices::colorRampPalette(colors)
        print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
            cexCol = 1, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
            Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", cellnote = round(similarity.matrix,
                1), notecol = "white", main = "Cor", margins = c(10, 10)))
        grDevices::dev.off()

        saveRDS(list(similarity = similarity, map = map), file = cache)

    }

    return(list(map = map, similarity = similarity))
}

