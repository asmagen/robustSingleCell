globalVariables(c("Modularity", "Type", "sd", "Fold", "correlation", "SD", "gene",
    "clustering", "Dataset", "..count..", "tSNE.1", "tSNE.2", "Cluster", "Origin",
    "Experiment", "cell_type", "Background.ratio", "Foreground.ratio", "Gene",
    "activation"))

#' Plot PCA results
#'
#' Plot the results obtained from PCA analysis as cell embedding in 2D space and annotation of gene loadings
#'
#' @param environment \code{environment} object
#' @param quantile quantile of PCA loadings for which to define top genes driving PCs
#' @param order ordering by which to plot the in heatmap of top genes driving PCs
#' @import GGally
#' @import ggrepel
#' @importFrom graphics plot
#' @importFrom stats quantile rnorm
#' @importFrom utils head write.csv
plot_PCA <- function(environment, quantile, order) {

    work.path <- environment$work.path
    PCA <- environment$PCA
    Rotation <- environment$Rotation
    cluster.names <- environment$cluster.names
    dataset <- environment$dataset.labels
    marker.genes <- environment$marker.genes

    drivers <- apply(Rotation, 1, function(v) {
        q <- quantile(v, c(quantile, 1 - quantile))
        names(v[v <= q[1] | v >= q[2]])
    })

    nPCs <- 3

    drivers <- apply(Rotation, 1, function(v) {
        q <- quantile(v, c(quantile, 1 - quantile))
        names(v[v <= q[1] | v >= q[2]])
    })
    grDevices::pdf(file.path(work.path, "Rotation.PCA.pdf"))
    data <- data.frame(gene = colnames(Rotation), Rotation = t(Rotation))
    for (row in seq(nrow(Rotation) - 1)) {
        if (typeof(drivers) == "list") {
            gene.set <- as.vector(unlist(drivers[c(row, row + 1)]))
        } else {
            gene.set <- as.vector(drivers[, c(row, row + 1)])
        }
        plot.data <- data[data$gene %in% gene.set, c(1, row + 1, row + 2)]
        print(ggplot(plot.data, aes_string(x = colnames(plot.data)[2], y = colnames(plot.data)[3],
            label = "gene")) + geom_text(check_overlap = TRUE, size = 2) + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black")))
        print(ggplot(plot.data, aes_string(x = colnames(plot.data)[2], y = colnames(plot.data)[3],
            label = "gene")) + geom_text(check_overlap = TRUE, size = 3) + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black")))
        plot.data <- plot.data[plot.data$gene %in% marker.genes, ]
        print(ggplot(plot.data, aes_string(x = colnames(plot.data)[2], y = colnames(plot.data)[3],
            label = "gene")) + geom_point() + geom_text_repel() + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black")))
    }
    grDevices::dev.off()

    grDevices::pdf(file.path(work.path, "all.PCA.pdf"))
    PCA_t <- t(PCA)
    rownames(PCA_t) <- NULL
    data <- data.frame(PCA_t, Cluster = factor(cluster.names), Dataset = factor(dataset))
    utils::head(data)
    for (row in seq(nrow(PCA) - 1)) {
        print(ggplot(data, aes_string(x = rownames(PCA)[row], y = rownames(PCA)[row +
            1], color = "Cluster")) + geom_point() + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black")))
    }
    if (length(unique(dataset)) > 1) {
        for (row in seq(nrow(PCA) - 1)) {
            print(ggplot(data, aes_string(x = rownames(PCA)[row], y = rownames(PCA)[row +
                1], color = "Dataset")) + geom_point() + scale_colour_brewer(palette = "Set3") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")))
        }
    }
    grDevices::dev.off()

    data <- data.frame(PCA_t, Cluster = factor(cluster.names), Dataset = factor(dataset))
    utils::head(data)
    rm(PCA_t)
    grDevices::pdf(file.path(work.path, "PC.scores.histogram.pdf"), width = 10)
    for (row in seq(nrow(PCA) - 1)) {
        print(ggplot(data, aes_string(x = rownames(PCA)[row], fill = "Cluster")) +
            geom_density(alpha = 0.5) + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black")) + ylab("Density") + theme_classic(base_size = 25))
    }
    if (length(unique(dataset)) > 1) {
        for (row in seq(nrow(PCA) - 1)) {
            print(ggplot(data, aes_string(x = rownames(PCA)[row], fill = "Dataset")) +
                geom_density(alpha = 0.5) + scale_fill_brewer(palette = "Paired") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                ylab("Density") + theme_classic(base_size = 25))
        }
    }
    grDevices::dev.off()

    grDevices::pdf(file.path(work.path, "PC.scores.heatmap.pdf"), width = 5, height = 8)
    for (row in seq(nrow(PCA) - 1)) {
        pc <- environment$Rotation[row, ]
        high <- names(pc[pc > quantile(pc, 0.99)])
        high
        low <- names(pc[pc < quantile(pc, 0.01)])
        low
        gene.list <- list(markers = c(low, high))
        names(gene.list) <- paste("PC", row)
        plot.expression.heatmap.based.on.FC.marker(environment$normalized, environment$cluster.names,
            gene.list = gene.list, scale = "row", counts = F, order = order, filter.diff.exp = F,
            cellnote = F)
    }
    if (length(unique(dataset)) > 1) {
        for (row in seq(nrow(PCA) - 1)) {
            pc <- environment$Rotation[row, ]
            high <- names(pc[pc > quantile(pc, 0.99)])
            high
            low <- names(pc[pc < quantile(pc, 0.01)])
            low
            gene.list <- list(markers = c(low, high))
            names(gene.list) <- paste("PC", row)
            plot.expression.heatmap.based.on.FC.marker(environment$normalized, environment$dataset.labels,
                gene.list = gene.list, scale = "row", counts = F, order = NA, filter.diff.exp = F,
                cellnote = F)
        }
    }
    grDevices::dev.off()
}

plot.cluster.stats <- function(environment, membership, label = NA, order = NA) {
    file.name <- "cluster.size.pdf"
    work.path <- environment$work.path
    if (!is.na(label)) {
        work.path <- file.path(environment$work.path, label)
        dir.create(work.path, showWarnings = F)
        file.name <- paste(label, file.name, sep = ".")
    }
    cluster.size <- table(environment$cluster.names)
    if (length(order) == 1 && is.na(order))
        order <- names(cluster.size)[order(cluster.size, decreasing = T)]
    grDevices::pdf(file.path(work.path, file.name), width = 8, height = 5)
    data <- data.frame(clustering = factor(membership, levels = order), Dataset = factor(environment$dataset.labels))
    if (length(unique(environment$dataset.labels)) > 1) {
        print(ggplot(data, aes(clustering, fill = Dataset)) + geom_bar() + scale_fill_brewer(palette = "Set3") +
            xlab("Cluster ID") + ylab("Number of cells") + theme_classic(base_size = 15) +
            theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 25,
            hjust = 1)) + guides(fill = guide_legend(nrow = max(1, floor(length(unique(environment$dataset.labels))/2)),
            byrow = F)))
        print(ggplot(data, aes(clustering, fill = Dataset)) + geom_bar(aes(y = (..count..)/sum(..count..))) +
            scale_y_continuous(labels = scales::percent) + ylab("relative frequencies") +
            scale_fill_brewer(palette = "Set3") + xlab("Cluster ID") + theme_classic(base_size = 15) +
            theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 25,
            hjust = 1)) + guides(fill = guide_legend(nrow = max(1, floor(length(unique(environment$dataset.labels))/2)),
            byrow = F)))
        print(ggplot(data, aes(Dataset, fill = clustering)) + geom_bar() + xlab("Cluster ID") +
            ylab("Number of cells") + theme_classic(base_size = 15) + theme(axis.text.x = element_text(angle = 25,
            hjust = 1)))
    }
    for (dataset in unique(environment$dataset.labels)) {
        print(ggplot(data[data$Dataset == dataset, ], aes(clustering)) + geom_bar(aes(y = (..count..)/sum(..count..))) +
            scale_y_continuous(labels = scales::percent) + ylab("Relative frequencies") + theme_classic(base_size = 15) +
            theme(axis.text.x = element_text(angle = 25, hjust = 1)) + xlab("Cluster ID") + ggtitle(dataset))
    }
    print(ggplot(data, aes(clustering)) + geom_bar() + theme_classic(base_size = 15) + theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
            ylab("Number of cells") + xlab("Cluster ID"))
    print(ggplot(data, aes(clustering)) + geom_bar(aes(y = (..count..)/sum(..count..))) +
        scale_y_continuous(labels = scales::percent) + ylab("Relative frequencies") + xlab("Cluster ID") + theme_classic(base_size = 15) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1)))

    grDevices::dev.off()

    cell.confusion <- table(membership, environment$dataset.labels)
    utils::write.csv(cell.confusion, file = file.path(environment$work.path, "cell.confusion.csv"))

    confounders <- environment$confounders
    data <- data.frame(clustering = factor(as.vector(membership)), dataset = factor(environment$dataset.labels),
        confounders)
    utils::head(data)
    file.name <- "confounder.stats.violin.pdf"
    if (!is.na(label))
        file.name <- paste(label, file.name, sep = ".")
    grDevices::pdf(file.path(work.path, file.name))
    for (variable in colnames(confounders)) {
        print(ggplot(data, aes_string("clustering", variable)) + geom_violin(draw_quantiles = c(0.25,
            0.5, 0.75), scale = "width") + theme(axis.text.x = element_text(angle = 25,
            hjust = 1)) + ggtitle(variable))
    }
    # if (length(unique(environment$dataset.labels))>1) {
    for (variable in colnames(confounders)) {
        print(ggplot(data, aes_string("dataset", variable)) + geom_violin(draw_quantiles = c(0.25,
            0.5, 0.75), scale = "width") + theme(axis.text.x = element_text(angle = 25,
            hjust = 1)) + ggtitle(variable))
    }
    # }
    grDevices::dev.off()
}

plot.tSNE <- function(environment, tSNE.job, perplexity, max_iter, membership = NA) {

    if (length(tSNE.job) == 1 && !is.na(tSNE.job)) {
        tSNE <- get_slurm_out(tSNE.job)
        tryCatch({
            cleanup_files(tSNE.job)
            cleanup_files(tSNE.job)
            cleanup_files(tSNE.job)
        }, error = function(v) v)
    }

    if (is.na(membership))
        membership <- environment$cluster.names

    configs <- expand.grid(perplexity, max_iter)

    for (row in seq(nrow(configs))) {
        perplexity <- configs[row, 1]
        max_iter <- configs[row, 2]
        tryCatch({
            tSNE <- readRDS(file.path(environment$res.data.path, "tSNEs", paste(perplexity,
                max_iter, "tSNE.rds", sep = ".")))

            duplicated.indices <- duplicated(t(environment$PCA))
            data <- data.frame(Cluster = factor(membership[!duplicated.indices]),
                Origin = factor(environment$origins[!duplicated.indices]), Experiment = factor(environment$experiments[!duplicated.indices]),
                tSNE = tSNE)
            utils::head(data)

            grDevices::pdf(file.path(environment$work.path, paste("tSNE_perplexity",
                perplexity, "max_iter", max_iter, "pdf", sep = ".")), width = 15,
                height = 10)
            print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = Cluster, shape = Origin)) +
                geom_point(data = data, size = 4, alpha = 0.6) + scale_shape(solid = T) +
                xlab("tSNE 1") + ylab("tSNE 2") + theme_classic(base_size = 20) +
                theme(legend.position = "bottom"))
            print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = Cluster, shape = Experiment)) +
                geom_point(data = data, size = 7, alpha = 0.4) + scale_shape(solid = T) +
                xlab("tSNE 1") + ylab("tSNE 2") + theme_classic(base_size = 20) +
                theme(legend.position = "bottom"))
            print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = Experiment, shape = Origin)) +
                geom_point(data = data, size = 7, alpha = 0.4) + scale_shape(solid = T) +
                xlab("tSNE 1") + ylab("tSNE 2") + theme_classic(base_size = 20) +
                theme(legend.position = "bottom") + scale_color_brewer(palette = "Set2"))
            print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = Origin, shape = Experiment)) +
                geom_point(data = data, size = 7, alpha = 0.4) + scale_shape(solid = T) +
                xlab("tSNE 1") + ylab("tSNE 2") + theme_classic(base_size = 20) +
                theme(legend.position = "bottom") + scale_color_brewer(palette = "Set2"))
            print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = Origin, shape = Origin)) +
                geom_point(data = data, size = 7, alpha = 0.4) + scale_shape(solid = T) +
                xlab("tSNE 1") + ylab("tSNE 2") + theme_classic(base_size = 20) +
                theme(legend.position = "bottom") + scale_color_brewer(palette = "Set2"))
            print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = Experiment, shape = Experiment)) +
                geom_point(data = data, size = 7, alpha = 0.4) + scale_shape(solid = T) +
                xlab("tSNE 1") + ylab("tSNE 2") + theme_classic(base_size = 20) +
                theme(legend.position = "bottom") + scale_color_brewer(palette = "Set2"))

            if (ncol(environment$confounders) > 0) {
                data <- data.frame(Cluster = factor(membership[!duplicated.indices]),
                  Origin = factor(environment$origins[!duplicated.indices]), Experiment = factor(environment$experiments[!duplicated.indices]),
                  environment$confounders, tSNE = tSNE)
                utils::head(data)
                # name = colnames(environment$confounders)[2]
                for (name in colnames(environment$confounders)) {
                  signature.activation <- data[[name]]
                  if (is.numeric(signature.activation)) {
                    print(ggplot(data, aes(x = tSNE.1, y = tSNE.2, color = signature.activation,
                      shape = Origin)) + geom_point(data = data, size = 4, alpha = 0.6) +
                      scale_shape(solid = T) + xlab("tSNE 1") + ylab("tSNE 2") +
                      theme_classic(base_size = 20) + theme(legend.position = "bottom") +
                      scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0) + guides(color = FALSE) + ggtitle(name))
                  }
                }
            }
        }, error = function(e) e)
        grDevices::dev.off()
    }
}

plot.expression.heatmap.based.on.FC.marker <- function(measurements, clustering,
    gene.list, counts = F, order = NA, RowSideColors = NA, scale = "row", save = NA,
    filter.diff.exp = T, cellnote = T, exponent = F, doMeans = T, srtCol = 45, multiplication = 100,
    rounding = 0, breaks = 50, key = T, sort.rows = T, sort.cols = T, Rowv = F, Colv = F,
    dendrogram = "none", main.remove = T) {

    colors <- c("#d60c0c", "#ffb1ad", "#f4f4f4", "#aed4fc", "#0050ba")  #c('#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4')
    colors <- colors[length(colors):1]  #red should correspond to high
    color.palette <- grDevices::colorRampPalette(colors)
    list.row <- 1
    titles <- names(gene.list)
    plots <- list()
    nclusters <- length(unique(clustering))

    list.row <- 1

    for (list.row in seq(length(gene.list))) {
        tryCatch({
            genes <- gene.list[[list.row]]

            if (length(genes) == 0)
                next
            if (length(genes) == 1)
                genes <- rep(genes, 2)

            genes <- genes[genes %in% rownames(measurements)]

            measurements.diff.exp <- measurements[match(genes, rownames(measurements)),
                ]
            if (doMeans) {
                names <- rownames(measurements.diff.exp)
                mat <- {
                }
                row <- 1
                if (counts)
                  measurements.diff.exp <- measurements.diff.exp > 0
                for (row in seq(length(genes))) {
                  averages <- apply.by.group(clustering, measurements.diff.exp[row,
                    ], mean)
                  mat <- rbind(mat, averages)
                }
                rownames(mat) <- names
                if (sort.cols) {
                  mat <- mat[, order(colnames(mat))]
                }
            } else {
                mat <- measurements.diff.exp
            }

            if (exponent)
                mat <- exp(mat)

            if (length(order) == 1 && is.na(order))
                order <- seq(nclusters)
            others <- mat
            others <- others[, order]
            if (sort.rows) {
                row.order <- order(unlist(apply(others, 1, which.max)), decreasing = T)
            } else {
                row.order <- seq(nrow(others))
            }
            others <- others[row.order, ]
            RowSideColors <- RowSideColors[row.order]
            if (filter.diff.exp) {
                significant.diff <- getDE.limma(Y = measurements[rownames(measurements) %in%
                  rownames(others), ], group = factor(clustering))
                keep <- rownames(others) %in% rownames(significant.diff)
                others <- others[keep, ]
                RowSideColors <- RowSideColors[keep]
            }

            if (length(unique(RowSideColors)) > 1 & sort.rows) {
                blocks <- as.vector(unlist(apply(others, 1, which.max)))
                index <- 1
                for (index in unique(blocks)) {
                  indices <- which(blocks == index)
                  if (length(indices) < 2)
                    next
                  order <- order(RowSideColors[indices], as.vector(others[indices,
                    index]))
                  RowSideColors[indices] <- RowSideColors[indices][order]
                  others[indices, ] <- others[indices, ][order, ]
                  rownames(others)[indices] <- rownames(others)[indices][order]
                }
            }

            RowSideColors[is.na(RowSideColors)] <- "white"
            RowSideColors[RowSideColors == "-1"] <- "blue"
            RowSideColors[RowSideColors == "1"] <- "red"

            if (cellnote == T) {
                cellnote <- round(multiplication * others, rounding)
                plots[[titles[list.row]]] <- heatmap.2(as.matrix(others), col = color.palette,
                  scale = scale, key = key, cexRow = 1, cexCol = 1, density.info = "none",
                  trace = "none", srtCol = srtCol, adjCol = c(1, 1), Rowv = Rowv,
                  Colv = Colv, margins = c(8, 5), dendrogram = dendrogram, main = titles[list.row],
                  RowSideColors = RowSideColors, cellnote = cellnote, notecol = "white")
            } else {
                plots[[titles[list.row]]] <- heatmap.2(as.matrix(others), col = color.palette,
                  scale = scale, key = key, cexRow = 1, cexCol = 1, density.info = "none",
                  trace = "none", srtCol = srtCol, adjCol = c(1, 1), Rowv = Rowv,
                  Colv = Colv, margins = c(8, 5), dendrogram = dendrogram, main = titles[list.row],
                  RowSideColors = RowSideColors)
            }

            if (!is.na(save)) {
                utils::write.csv(as.matrix(others), file = save)
            }
        }, error = function(e) print(e))
    }

    rownames(others)
}

#' Plot heatmap
#'
#' Plot heatmap given a set of markers.
#'
#' @param environment The \code{environment} object
#' @param name The file name of the figure
#' @param markers The markers to be plotted
#' @param path The path where the plot is saved; by default in TMPDIR
#' @param membership The cluster membership
#' @param normalized The normalized data matrix
#' @param order The ordering of markers
#' @param width The width of the pdf figure
#' @param height The height of the pdf figure
#' @param save The path where the plot is saved
#' @param counts Plot count matrix or not
#' @param filter.diff.exp Whether to filter for differentially expressed genes
#' @param sort.rows Whether to sort rows
#' @param sort.cols Whether to sort columns
#' @inheritParams gplots::heatmap.2
#' @export
plot_simple_heatmap <- function(environment, name, markers, path = NA, membership = NA,
    normalized = NA, order = NA, width = 5, height = 5, scale = "row", RowSideColors = NA,
    counts = F, filter.diff.exp = F, cellnote = F, key = F, save = NA, sort.rows = T,
    sort.cols = T, Colv = F, Rowv = F, dendrogram = "none", main = NA) {
    if (check_not_slurm("plot_simple_heatmap")) {
        return(NULL)
    }
    if (is.na(path))
        path <- environment$work.path
    grDevices::pdf(file = file.path(path, paste(name, "heatmap.pdf", sep = ".")),
        width = width, height = height)

    if (is.na(membership))
        membership <- environment$cluster.names
    if (is.na(normalized))
        normalized <- environment$normalized

    cluster.size <- table(membership)
    if (length(order) == 1 && is.na(order))
        order <- names(cluster.size)[order(cluster.size, decreasing = T)]

    gene.list <- list(markers = markers)
    if (!is.null(main))
        names(gene.list) <- main
    plot.expression.heatmap.based.on.FC.marker(normalized, membership, gene.list = gene.list,
        scale = scale, RowSideColors = RowSideColors, counts = counts, order = order,
        filter.diff.exp = filter.diff.exp, cellnote = cellnote, doMeans = T, exponent = ifelse(counts,
            F, T), multiplication = ifelse(counts, 100, 10), rounding = 0, key = key,
        save = save, sort.rows = sort.rows, sort.cols = sort.cols, Rowv = Rowv, Colv = Colv,
        dendrogram = dendrogram)

    if (key & length(RowSideColors) > 1) {
        graphics::legend("topright", legend = unique(RowSideColors), col = unique(as.numeric(RowSideColors)),
            lty = 1, lwd = 5, cex = 0.7)
    }

    grDevices::dev.off()
}

plot.violin <- function(environment, genes, types, fore1exp1, fore2exp1, fore1exp2,
    fore2exp2, back1exp1, back2exp1, back1exp2, back2exp2, path, height = 5, width = 5,
    scale = T, palette = "Greys", separate.background = F) {
    # Pastel2

    violin.plot.data <- {
    }
    cluster.indices <- environment$cluster.names %in% c(fore1exp1, fore2exp1, fore1exp2,
        fore2exp2, back1exp1, back2exp1, back1exp2, back2exp2)
    for (gene in genes) {
        if (scale) {
            expression <- environment$normalized[gene, cluster.indices]
            expression <- (expression - min(expression))/(max(expression) - min(expression))
        } else {
            expression <- environment$normalized[gene, cluster.indices]
        }
        violin.plot.data <- rbind(violin.plot.data, data.frame(cluster = environment$cluster.names[cluster.indices],
            gene = gene, expression = expression))
    }

    violin.plot.data$cell_type <- NA
    violin.plot.data$cell_type[violin.plot.data$cluster %in% c(fore1exp1, fore1exp2)] <- types[1]
    violin.plot.data$cell_type[violin.plot.data$cluster %in% c(fore2exp1, fore2exp2)] <- types[2]

    if (separate.background) {
        violin.plot.data$cell_type[violin.plot.data$cluster %in% c(back1exp1, back1exp2)] <- "Background 1"
        violin.plot.data$cell_type[violin.plot.data$cluster %in% c(back2exp1, back2exp2)] <- "Background 2"
        violin.plot.data$cell_type <- factor(violin.plot.data$cell_type, levels = c("Background 1",
            types, "Background 2"))
    } else {
        violin.plot.data$cell_type[violin.plot.data$cluster %in% c(back1exp1, back2exp1,
            back1exp2, back2exp2)] <- "Others"
        violin.plot.data$cell_type <- factor(violin.plot.data$cell_type, levels = c(types,
            "Others"))
    }

    tryCatch({
        grDevices::pdf(file.path(path, paste(paste(types, collapse = "_"), "violin.pdf",
            sep = ".")), height = height, width = width)
        print(ggplot(violin.plot.data, aes(x = gene, y = expression, fill = cell_type)) +
            geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) + theme_classic(base_size = 15) +
            theme(axis.text.x = element_text(angle = 25, hjust = 1)) + scale_fill_brewer(palette = palette))
        grDevices::dev.off()
    }, error = function(v) {
        grDevices::pdf(file.path(path, paste(paste(types, collapse = "_"), "violin.pdf",
            sep = ".")), height = height, width = width)
        print(ggplot(violin.plot.data, aes(x = gene, y = expression, fill = cell_type)) +
            geom_violin(scale = "width", draw_quantiles = c(0.5)) + theme_classic(base_size = 15) +
            theme(axis.text.x = element_text(angle = 25, hjust = 1)) + scale_fill_brewer(palette = palette))
        grDevices::dev.off()
    })
}


plot.heatmaps <- function(environment, diff.exp, membership, order = NA, nTopRanked = 10,
    label = NA) {

    symbols <- c()
    for (cluster in unique(diff.exp$cluster)) {
        top.ranked <- diff.exp[diff.exp$cluster == cluster, ]
        top.ranked <- top.ranked[order(top.ranked$fold, decreasing = T), ]
        symbols <- c(symbols, utils::head(as.vector(top.ranked$gene), nTopRanked))
    }
    symbols <- unique(symbols)
    symbols

    clustering <- as.vector(membership)

    file.name <- "diff.genes.pdf"
    work.path <- environment$work.path
    if (length(label) == 1 && !is.na(label)) {
        work.path <- file.path(environment$work.path, label)
        dir.create(work.path, showWarnings = F)
        file.name <- paste(label, file.name, sep = ".")
    }


    grDevices::pdf(file = file.path(work.path, file.name), width = length(unique(clustering))/2,
        height = length(symbols)^(1/1.5))
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = symbols), order = order, filter.diff.exp = T,
        cellnote = T, doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = symbols), order = order, filter.diff.exp = F,
        cellnote = T, doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = symbols), counts = T, order = order, filter.diff.exp = F,
        cellnote = T, doMeans = T, exponent = F, multiplication = 100, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = symbols), counts = T, order = order, filter.diff.exp = T,
        cellnote = T, doMeans = T, exponent = F, multiplication = 100, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = symbols), order = order, filter.diff.exp = F,
        scale = "col", cellnote = T, doMeans = T, exponent = T, multiplication = 10,
        rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = symbols), order = order, filter.diff.exp = F,
        scale = "none", cellnote = T, doMeans = T, exponent = T, multiplication = 10,
        rounding = 0)
    grDevices::dev.off()

    extended.genes <- environment$marker.genes
    file.name <- "all.markers.pdf"
    if (length(label) == 1 && !is.na(label))
        file.name <- paste(label, file.name, sep = ".")
    grDevices::pdf(file = file.path(work.path, file.name), width = length(unique(clustering))/2,
        height = length(extended.genes)^(1/1.5))
    plot.expression.heatmap.based.on.FC.marker(measurements = environment$normalized,
        clustering, gene.list = list(markers = extended.genes), order = order, filter.diff.exp = T,
        doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(measurements = environment$normalized,
        clustering, gene.list = list(markers = extended.genes), counts = T, order = order,
        filter.diff.exp = T, doMeans = T, exponent = F, multiplication = 100, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = extended.genes), order = order, filter.diff.exp = T,
        scale = "col", doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = extended.genes), order = order, filter.diff.exp = F,
        scale = "col", doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = extended.genes), order = order, filter.diff.exp = T,
        scale = "none", doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    plot.expression.heatmap.based.on.FC.marker(environment$normalized, clustering,
        gene.list = list(markers = extended.genes), order = order, filter.diff.exp = F,
        scale = "none", doMeans = T, exponent = T, multiplication = 10, rounding = 0)
    grDevices::dev.off()
}

#' Plot Gene Expression on tSNE
#'
#' Visualize normalized expression of selected genes on tSNE plot with color-code and contour annotation.
#'
#' @param environment \code{environment} object
#' @param genes selected genes to visualize
#' @param perplexity tSNE perplexity parameter
#' @param max_iter tSNE max_iter parameter
#' @param width pdf file canvas width
#' @param height pdf file canvas height
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' plot_contour_overlay_tSNE(LCMV1,genes = c('Cd4','Cd8a'))
#' }
plot_contour_overlay_tSNE <- function (environment,genes,perplexity = 30,max_iter = 10000,width = 10, height = 10) {

    if (any(!genes %in% environment$genes)) {
        cat('Removing genes not found in dataset:')
        print(genes[!genes %in% environment$genes])
        genes = genes[genes %in% environment$genes]
    }

    tSNE <- readRDS(file.path(environment$res.data.path, "tSNEs", paste(perplexity,
                max_iter, "tSNE.rds", sep = ".")))

    duplicated.indices <- duplicated(t(environment$PCA))

    grDevices::pdf(file = file.path(environment$work.path, "contour.overlay.tSNE.pdf"),
        width = width, height = height)
    for (gene in genes) {
        data <- data.frame(tSNE = tSNE,activation = scale(environment$normalized[gene,]));head(data)
        data.filtered <- data[data$activation > quantile(data$activation,0.9),]
        print(ggplot() + geom_point(data = data,aes(x = tSNE.1, y = tSNE.2,color=activation),size = 4,alpha = 0.6) + scale_shape(solid = T) + geom_density_2d(data = data.filtered, aes(x = tSNE.1, y = tSNE.2),alpha = 1,colour = "gray65") + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + scale_colour_gradient(low = "#ededed", high = "#a50303") + theme(legend.position="bottom") + ggtitle(gene))
    }
    grDevices::dev.off()
}

#' Plot Pairwise Gene Scatter Plot
#'
#' Visualize normalized expression contours of a selected gene pair across selected cluster groups.
#'
#' @param environment \code{environment} object
#' @param gene1 selected gene number 1
#' @param gene2 selected gene number 2
#' @param cluster_group1 cluster group 1 to be visualized (one or more clusters)
#' @param cluster_group2 cluster group 2 to be visualized (one or more clusters)
#' @param group1_label label for group 1 legend and file name
#' @param group2_label label for group 2 legend and file name
#' @param width pdf file canvas width
#' @param height pdf file canvas height
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' cluster_names <- get.cluster.names(LCMV1, types, min.fold = 1.0, max.Qval = 0.01)
#' LCMV1 <- set.cluster.names(LCMV1, names = cluster_names)
#' plot_pair_scatter(LCMV1,
#' gene1 = 'Cd4',
#' gene2 = 'Cd8',
#' cluster_group1 = cluster_names[1:2],
#' cluster_group2 = cluster_names[3:4],
#' group1_label = 'CD4 T Cells',
#' group2_label = 'CD8 T Cells')
#' }
plot_pair_scatter <- function (environment,gene1,gene2,cluster_group1,cluster_group2,group1_label,group2_label,width = 10, height = 10) {
    if (check_not_slurm("plot_pair_scatter")) {
        return()
    }
    clusters <- c(cluster_group1, cluster_group2)
    plot.data <- data.frame(gene1 = environment$normalized[gene1,environment$cluster.name %in% clusters],gene2 = environment$normalized[gene2,environment$cluster.name %in% clusters]);head(plot.data)

    plot.data[,1][plot.data[,1]==0] <- rnorm(sum(plot.data[,1]==0),sd = min(plot.data[,1][plot.data[,1]!=0])/5)
    plot.data[,2][plot.data[,2]==0] <- rnorm(sum(plot.data[,2]==0),sd = min(plot.data[,2][plot.data[,2]!=0])/5)

    cluster <- environment$cluster.name[environment$cluster.name %in% clusters]

    cluster[cluster %in% c(cluster_group1)] <- group1_label
    cluster[cluster %in% c(cluster_group2)] <- group2_label

    plot.data <- cbind(plot.data,cluster = cluster)
    colnames(plot.data)[1:2] <- c(gene1,gene2)

    grDevices::pdf(file = file.path(environment$work.path, paste(group1_label,group2_label,gene1,gene2,'pdf',sep='.')),
        width = width, height = height)
    print(ggplot(plot.data, aes_string(x = gene1, y = gene2, colour = 'cluster')) + geom_point(aes(alpha=0.5,size=3)) + geom_density_2d() + theme_classic(base_size=15) + ggtitle(paste(gene1,'vs',gene2,'expression in',group1_label,'vs',group2_label,'clusters',sep=' ')))
    grDevices::dev.off()
}

#' Visualize Correlation
#'
#' Plot correlation heatmaps for each pair of datasets.
#'
#' @param environment \code{environment} object
#' @param work.path where to locate the figures
#' @param similarity similarity matrix defined in compare.cluster.similarity or get.robust.cluster.similarity
#' @param margins The margins to the plot
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
#' similarity <- cluster.similarity$similarity
#' map <- cluster.similarity$map
#' filtered.similarity <- get.robust.cluster.similarity(
#'    pooled_env, similarity, min.sd = qnorm(.9), max.q.val = 0.01, rerun = F)
#' robust.clusters <- sort(unique(c(filtered.similarity$cluster1,
#' filtered.similarity$cluster2)))
#' visualize.cluster.cors.heatmaps(pooled_env, pooled_env$work.path,
#'                                filtered.similarity)
#' }
visualize.cluster.cors.heatmaps <- function(environment, work.path, similarity, margins = c(17, 17)) {
    if (check_not_slurm("visualize.cluster.cors.heatmaps")) {
        return()
    }
    name <- "cross.sample.similarities"
    work.path <- file.path(environment$work.path, name)
    if (file.exists(work.path)) {
        new.dir <- file.path(environment$work.path, paste(name, format(Sys.time(),
            "%a_%b_%e_%Y__%H_%M_%S"), sep = "---"))
        file.rename(work.path, new.dir)
    }
    dir.create(work.path, showWarnings = T, recursive = T)


    samples.names <- unique(c(as.vector(similarity$sample1), as.vector(similarity$sample2)))
    samples.names <- sort(samples.names)
    sample.name1 <- samples.names[3]
    sample.name2 <- samples.names[4]
    for (sample1 in seq(length(samples.names) - 1)) {
        sample.name1 <- samples.names[sample1]
        for (sample2 in (sample1 + 1):length(samples.names)) {
            sample.name2 <- samples.names[sample2]
            sub.work.path <- file.path(work.path, sample.name1, sample.name2)
            dir.create(sub.work.path, showWarnings = T, recursive = T)
            file.name <- paste(paste(strsplit(sample.name1, split = " ")[[1]], collapse = "."),
                "to", paste(strsplit(sample.name2, split = " ")[[1]], collapse = "."),
                sep = "_")
            samples <- c(sample.name1, sample.name2)  #unique(c(as.vector(similarity$sample1[similarity$sample1 == sample.name1]),as.vector(similarity$sample2[similarity$sample2 == sample.name2])))
            sample <- similarity[(similarity$sample1 %in% samples[1] & similarity$sample2 %in%
                samples[2]) | (similarity$sample1 %in% samples[2] & similarity$sample2 %in%
                samples[1]), ]

            similarity.summary.df <- data.frame(cluster1 = sample$name1, cluster2 = sample$name2,
                coef = sample$similarity)
            similarity.matrix <- acast(similarity.summary.df, cluster1 ~ cluster2,
                value.var = "coef")
            utils::write.csv(sample, file = file.path(sub.work.path, paste(file.name,
                "similarity.csv", sep = "_")))
            similarity.matrix <- t(similarity.matrix)

            ocl.similarity.summary.df <- data.frame(cluster1 = sample$name1, cluster2 = sample$name2,
                coef = -scale(sample$ocldist))
            ocl.similarity.matrix <- acast(ocl.similarity.summary.df, cluster1 ~
                cluster2, value.var = "coef")
            ocl.similarity.matrix <- t(ocl.similarity.matrix)

            grDevices::pdf(file.path(sub.work.path, paste(file.name, "similarity.heatmap.pdf",
                sep = "_")), width = 10, height = 10)
            colors <- rev(brewer.pal(5, "PuOr"))
            color.palette <- grDevices::colorRampPalette(colors)
            print(heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 1,
                cexCol = 1, srtCol = 45, scale = "none", density.info = "none", trace = "none",
                Rowv = T, Colv = T, dendrogram = "both", margins = margins, cellnote = round(similarity.matrix,
                  1), notecol = "white", main = "Pearson Correlation Between Cluster FC"))
            print(heatmap.2(ocl.similarity.matrix, col = color.palette, key = T,
                cexRow = 1, cexCol = 1, srtCol = 45, scale = "none", density.info = "none",
                trace = "none", Rowv = F, Colv = F, dendrogram = "none", margins = margins,
                cellnote = round(ocl.similarity.matrix, 1), notecol = "white",
                main = "Cluster mean Euclidean similarity"))
            grDevices::dev.off()
        }
    }
}

#' Plot Similarity Results
#'
#' Perform hierarchical clustering and plot cluster similarities according to dendrogram.
#'
#' @param environment \code{environment} object
#' @param similarity similarity matrix defined in compare.cluster.similarity or get.robust.cluster.similarity
#' @param hclust.resolution clustering resolution to impose on hclust cutree function
#' @param margins The margins to the plot
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
#' similarity <- cluster.similarity$similarity
#' map <- cluster.similarity$map
#' filtered.similarity <- get.robust.cluster.similarity(
#'    pooled_env, similarity, min.sd = qnorm(.9), max.q.val = 0.01, rerun = F)
#' visualize.cluster.similarity.stats(pooled_env, filtered_similarity)
#' }
visualize.cluster.similarity.stats <- function(environment, similarity,hclust.resolution = 8, margins = c(40, 40)) {
    net <- igraph::graph_from_data_frame(d = similarity[similarity$similarity > 0.1,
        c("name1", "name2")], directed = F)
    deg <- igraph::degree(net, mode = "all")
    igraph::V(net)$size <- deg
    l <- igraph::layout_with_kk(net)

    similarity.summary.df <- data.frame(cluster1 = similarity$name1, cluster2 = similarity$name2,
        coef = similarity$similarity)
    mirror <- similarity.summary.df
    cluster1 <- mirror$cluster1
    cluster2 <- mirror$cluster2
    mirror$cluster1 <- cluster2
    mirror$cluster2 <- cluster1
    similarity.summary.df <- rbind(similarity.summary.df, mirror)
    similarity.matrix <- reshape2::acast(similarity.summary.df, cluster1 ~ cluster2,
        value.var = "coef")

    grDevices::pdf(file.path(environment$work.path, paste("hclust.dist.Cor.FC.pdf",
        sep = "_")), width = 30, height = 30)
    similarity.matrix[1:5, 1:5]
    stats::dist(similarity.matrix)
    hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
    plot(hc.dist)
    colors <- rev(RColorBrewer::brewer.pal(5, "PuOr"))
    color.palette <- grDevices::colorRampPalette(colors)
    hc.dist <- stats::hclust(stats::as.dist(1 - similarity.matrix))
    clusters <- stats::cutree(hc.dist, k = hclust.resolution)

    clusters.ordered <- clusters[match(names(igraph::V(net)), names(clusters))]
    mark.groups <- lapply(unique(clusters.ordered), function(c) as.vector(which(clusters.ordered ==
        c)))
    lapply(unique(clusters.ordered), function(c) names(clusters.ordered[clusters.ordered ==
        c]))
    plot(net, mark.groups = mark.groups, layout = l)
    print(gplots::heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 3.5,
        cexCol = 3.5, scale = "none", density.info = "none", trace = "none", Rowv = stats::as.dendrogram(hc.dist),
        Colv = stats::as.dendrogram(hc.dist), dendrogram = "both", notecol = "white",
        main = "", keysize = 1, margins = margins, cellnote = round(similarity.matrix, 1), notecex = 2))  # 
    # print(gplots::heatmap.2(similarity.matrix, col = color.palette, key = T, cexRow = 3.5,
    #     cexCol = 3.5, scale = "none", density.info = "none", trace = "none", keysize = 1,
    #     Rowv = T, Colv = T, dendrogram = "both", notecol = "white", main = "", margins = margins))  # cellnote = round(similarity.matrix, 1), notecex = 2, main = 'Pearson Correlation Between Cluster FC'
    grDevices::dev.off()
}
