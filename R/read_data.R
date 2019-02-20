#' Read 10X Data
#'
#' Load sparse data matrices from 10X genomics.
#'
#' @param path Path to directory containing matrix.mtx, genes.tsv, and barcodes.tsv
#' @return a matrix of genes by cells
#' @export
read.10x.data <- function(path) {
    barcode.path <- list.files(path, pattern = "barcodes.tsv", full = T)
    features.path <- list.files(path, pattern = "genes.tsv", full = T)
    matrix.path <- list.files(path, pattern = "matrix.mtx", full = T)
    mat <- Matrix::readMM(file = matrix.path)
    feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(mat) <- barcode.names$V1
    rownames(mat) <- feature.names$V2
    return(mat)
}


#' Read and Preprocess Data
#'
#' Read 10X data files or a raw data matrix and perform normalization, QC filtering and duplicates removal.
#'
#' @param environment \code{environment} object
#' @param genome genome annotation
#' @param min.genes.per.cell minimum required number of genes per cell
#' @param max.genes.per.cell.quantile upper quantile for number of genes per cell
#' @param max.UMIs.per.cell.quantile upper quantile for number of UMIs per cell
#' @param min.cells.per.gene minimum required number of cells per gene
#' @param max.mitochondrial.frac maximum fraction of reads mapped to mitochondrial
#' genes per cell
#' @param max.ribosomal.frac maximum fraction of reads mapped to ribosomal genes per cell
#' @param cell.filters filtering option for cells based on marker genes
#' @param raw.data.matrices logical indicating if data matrices is provided instead of 10X dataset
#' @param rerun whether to rerun loading the dataset or load from cache
#' @param subsample number of cells to subsample
#' @param seed seed for subsampling of cells
#' @export
#' @import ggplot2
read.data <- function(environment, genome = "mm10", min.genes.per.cell = 500,
    max.genes.per.cell.quantile = 0.98, max.UMIs.per.cell.quantile = 0.98, min.cells.per.gene = 1,
    max.mitochondrial.frac = 0.1, max.ribosomal.frac = NA, cell.filters = NA,
    raw.data.matrices = NA, rerun = F, subsample = NULL, seed = 0) {
    # browser()
    cache <- file.path(environment$baseline.data.path, "data.RData")

    if (!rerun & file.exists(cache)) {
        print.message("Loading precomputed")
        load(cache)
    } else {

        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"))

        merged <- NA
        dataset <- environment$datasets[1]
        merged.dataset.labels <- merged.origins <- merged.experiments <- merged.criteria <- {
        }
        set.seed(seed)
        for (dataset in environment$datasets) {
            if (is.na(raw.data.matrices)) {
                print.message("Loading", dataset)
                measurements <- read.10x.data(path = file.path(environment$data.path,
                  dataset))
                measurements <- as.matrix(measurements[, Matrix::colSums(measurements >
                  0) >= min.genes.per.cell])
            } else {
                print.message("Using input", dataset)
                measurements <- raw.data.matrices[[dataset]]
            }
            if (!is.null(subsample) & subsample < ncol(measurements)) {
                measurements <- measurements[, sample(seq(ncol(measurements)),
                  subsample)]
            }
            colnames(measurements) <- rep(environment$datasets[environment$datasets ==
                dataset], ncol(measurements))
            dataset.labels <- rep(paste(environment$origins[environment$datasets ==
                dataset], " (", environment$experiments[environment$datasets ==
                dataset], ")", sep = ""), ncol(measurements))
            origins <- rep(environment$origins[environment$datasets == dataset],
                ncol(measurements))
            experiments <- rep(environment$experiments[environment$datasets ==
                dataset], ncol(measurements))
            cat(dim(measurements))
            # corner(measurements) measurements = measurements[,measurements['Cd4',]>0]

            data <- data.frame(nUMI = colSums(measurements), nGenes = colSums(measurements >
                0), Ribosomal = colMeans(measurements[get.ribo.genes(rownames(measurements)),
                ]), Mitochondrial = colMeans(measurements[get.mito.genes(rownames(measurements)),
                ]), Ribosomal.frac = colSums(measurements[get.ribo.genes(rownames(measurements)),
                ])/colSums(measurements), Mitochondrial.frac = colSums(measurements[get.mito.genes(rownames(measurements)),
                ])/colSums(measurements))
            print(head(data))

            pdf(file.path(environment$baseline.work.path, paste(dataset, "pre.filter.dataset.stats.pdf",
                sep = ".")))
            for (ind1 in seq(length(colnames(data)) - 1)) {
                for (ind2 in (ind1 + 1):(length(colnames(data)))) {
                  v1 <- colnames(data)[ind1]
                  v2 <- colnames(data)[ind2]
                  print(ggplot(data, aes_string(v1, v2)) + geom_point())
                }
            }
            dev.off()

            print.message("Original dimensions")
            print(dim(measurements))

            genes.per.cell <- colSums(measurements > 0)
            print.message("Original genes.per.cell")
            print(summary(genes.per.cell))
            print(quantile(genes.per.cell, seq(0.01, 1, 0.01)))
            UMIs.per.cell <- colSums(measurements)
            print.message("Original UMIs.per.cell")
            print(summary(UMIs.per.cell))
            print(quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))
            max.genes.per.cell <- quantile(genes.per.cell, max.genes.per.cell.quantile)
            print.message("max.genes.per.cell", max.genes.per.cell)
            max.UMIs.per.cell <- quantile(UMIs.per.cell, max.UMIs.per.cell.quantile)
            print.message("max.UMIs.per.cell", max.UMIs.per.cell)
            print.message("max.mitochondrial.frac", max.mitochondrial.frac)

            print.message("# not qualify min.genes.per.cell", sum(genes.per.cell <
                min.genes.per.cell))
            print.message("# not qualify max.genes.per.cell", sum(genes.per.cell >
                max.genes.per.cell))
            print.message("# not qualify max.UMIs.per.cell", sum(UMIs.per.cell >
                max.UMIs.per.cell))
            print.message("# not qualify max.mitochondrial.frac", sum(data$Mitochondrial.frac >
                max.mitochondrial.frac))

            final.gene.filter.passed.cells <- rep(T, ncol(measurements))
            if (!is.na(cell.filters) && dataset %in% as.vector(cell.filters$dataset)) {
                dataset.filters <- cell.filters[as.vector(cell.filters$dataset) ==
                  dataset, ]
                filter.genes <- as.vector(dataset.filters$gene)
                if (sum(filter.genes %in% rownames(measurements)) == 0)
                  print.message("All filter genes do not exist in dataset - not filtering anything")
                not.exist.genes <- sum(!filter.genes %in% rownames(measurements))
                if (not.exist.genes > 0)
                  print.message("Not filtering based on", not.exist.genes, "filter genes do not exist in dataset")
                filter.genes <- filter.genes[filter.genes %in% rownames(measurements)]
                if (length(filter.genes) == 1)
                  filter.genes <- c(filter.genes, filter.genes)
                gene <- filter.genes[1]
                matched.conditions <- rep(T, ncol(measurements))
                for (gene in filter.genes) {
                  matched.conditions <- matched.conditions & ((measurements[gene,
                    ] > 0) == dataset.filters$expressed[dataset.filters$gene ==
                    gene])  #check if matching condition
                }
                sum.not.matching <- sum(matched.conditions)
                print.message("# not qualify dataset.filters", sum.not.matching,
                  "out of", ncol(measurements), "[", round(sum.not.matching/ncol(measurements) *
                    100, 1), "% ]")
                final.gene.filter.passed.cells <- matched.conditions
            }

            criteria <- genes.per.cell >= min.genes.per.cell & genes.per.cell <=
                max.genes.per.cell & UMIs.per.cell <= max.UMIs.per.cell & data$Mitochondrial.frac <=
                max.mitochondrial.frac & final.gene.filter.passed.cells
            if (!is.na(max.ribosomal.frac))
                criteria <- criteria & data$Ribosomal.frac <= max.ribosomal.frac
            print.message("# qualifying cells", sum(criteria), "# not qualifying cells",
                sum(!criteria))
            measurements <- measurements[, criteria]
            dataset.labels <- dataset.labels[criteria]
            origins <- origins[criteria]
            experiments <- experiments[criteria]
            genes.per.cell <- colSums(measurements > 0)
            print.message("Filtered genes.per.cell")
            print(summary(genes.per.cell))
            print(quantile(genes.per.cell, seq(0.01, 1, 0.01)))
            UMIs.per.cell <- colSums(measurements)
            print.message("Filtered UMIs.per.cell")
            print(summary(UMIs.per.cell))
            print(quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))

            data <- data.frame(nUMI = colSums(measurements), nGenes = colSums(measurements >
                0), Ribosomal = colMeans(measurements[get.ribo.genes(rownames(measurements)),
                ]), Mitochondrial = colMeans(measurements[get.mito.genes(rownames(measurements)),
                ]), Ribosomal.frac = colSums(measurements[get.ribo.genes(rownames(measurements)),
                ])/colSums(measurements), Mitochondrial.frac = colSums(measurements[get.mito.genes(rownames(measurements)),
                ])/colSums(measurements))

            print(head(data))
            pdf(file.path(environment$baseline.work.path, paste(dataset, "post.filter.dataset.stats.pdf",
                sep = ".")))
            for (ind1 in seq(length(colnames(data)) - 1)) {
                for (ind2 in (ind1 + 1):(length(colnames(data)))) {
                  v1 <- colnames(data)[ind1]
                  v2 <- colnames(data)[ind2]
                  print(ggplot(data, aes_string(v1, v2)) + geom_point())
                }
            }
            dev.off()

            if (length(merged) == 1 && is.na(merged)) {
                merged <- measurements
                merged.dataset.labels <- dataset.labels
                merged.origins <- origins
                merged.experiments <- experiments
                merged.criteria <- criteria
            } else {
                if (sum(rownames(merged) != rownames(measurements)) > 0)
                  stop("Feature genes mismatch - need to correct dataset binding matching")
                merged <- cbind(merged, measurements)
                merged.dataset.labels <- c(merged.dataset.labels, dataset.labels)
                merged.origins <- c(merged.origins, origins)
                merged.experiments <- c(merged.experiments, experiments)
                merged.criteria <- c(merged.criteria, criteria)
            }
            print(table(colnames(merged)))
        }

        counts <- merged
        dataset.labels <- merged.dataset.labels
        origins <- merged.origins
        experiments <- merged.experiments
        criteria <- merged.criteria
        rm(merged, measurements, merged.dataset.labels, merged.origins, merged.experiments,
            merged.criteria)

        cells.per.gene <- rowSums(counts > 0)
        print.message("Filtered cells.per.gene")
        print(summary(cells.per.gene))
        genes.filter <- cells.per.gene >= min.cells.per.gene
        print.message("# qualifying genes", sum(genes.filter), "# not qualifying genes",
            sum(!genes.filter))

        print.message("Aggregated dataset dim")
        print(dim(counts))
        genes.per.cell <- colSums(counts[genes.filter, ] > 0)
        print.message("Aggregated dataset genes.per.cell")
        print(summary(genes.per.cell))
        print(quantile(genes.per.cell, seq(0.01, 1, 0.01)))
        UMIs.per.cell <- colSums(counts[genes.filter, ])
        print.message("Aggregated dataset UMIs.per.cell")
        print(summary(UMIs.per.cell))
        print(quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))

        rownames <- data.frame(old = rownames(counts))
        if (sum(duplicated(rownames(counts))) > 0) {
            gene <- unique(rownames(counts)[duplicated(rownames(counts))])[2]
            for (gene in unique(rownames(counts)[duplicated(rownames(counts))])) {
                print(gene)
                indices <- which(rownames(counts) == gene)
                print(rownames(counts)[indices])
                rownames(counts)[indices] <- paste(rownames(counts)[indices],
                  seq(length(indices)), sep = ".")
                print(rownames(counts)[indices])
            }
        }
        rownames <- data.frame(rownames, new = rownames(counts))
        normalized <- sweep(counts, MARGIN = 2, FUN = "/", STATS = colSums(counts))
        # sum((counts[,1]/sum(counts[,1]))!=normalized[,1])
        normalized <- normalized * 10000
        normalized <- log(normalized + 1)
        print.message("Normalized")
        corner(normalized)

        t <- start(file.path(environment$work.path, "tracking"), append = T,
            split = T)
        duplicated.indices <- duplicated(t(counts[genes.filter, ])) | duplicated(t(counts[genes.filter,
            ]), fromLast = T)
        if (sum(duplicated.indices) > 0) {
            print.message("\nRemoving all", sum(duplicated.indices), "duplicated cells\n")
            print(table(colnames(counts[, duplicated.indices])))
            counts <- counts[, !duplicated.indices]
            normalized <- normalized[, !duplicated.indices]
            dataset.labels <- dataset.labels[!duplicated.indices]
            origins <- origins[!duplicated.indices]
            experiments <- experiments[!duplicated.indices]
            criteria[criteria == T][duplicated.indices] <- FALSE
        }
        end(t)

        save(genes.filter, counts, normalized, dataset.labels, origins, experiments,
            criteria, file = cache)

        end(t)
    }

    counts <- counts[genes.filter, ]
    normalized <- normalized[genes.filter, ]

    environment$genes.filter <- genes.filter
    environment$counts <- counts
    environment$normalized <- normalized
    environment$genes <- rownames(normalized)
    environment$datasets <- colnames(normalized)
    environment$dataset.labels <- dataset.labels
    environment$origins <- origins
    environment$experiments <- experiments
    if (!exists("criteria"))
        criteria <- NA
    environment$criteria <- criteria
    environment$confounders <- data.frame(nUMI = colSums(counts), nGenes = colSums(counts >
        0))
    environment$nsamples <- ncol(counts)

    return(environment)
}

#' Read Preclustered Datasets
#'
#' Read previous analysis of multiple datasets to perform integrated analysis.
#'
#' @param environment \code{environment} object
#' @param path search path for previous projects
#' @param recursive recursive path search
#' @param rerun whether to rerun the reading process or load from cache
#' @export
read.preclustered.datasets <- function(environment, path = NA, recursive = T,
    rerun = F) {

    cache <- file.path(environment$baseline.data.path, "preclustered.datasets.RData")

    if (!rerun & file.exists(cache)) {
        print.message("Loading precomputed")
        load(cache)
    } else {

        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"), split = T)

        if (is.na(path))
            path <- dirname(environment$baseline.work.path)
        dataset <- environment$datasets[1]
        merged.clustering <- {
        }
        merged.original.clustering <- {
        }
        merged.diff.exp <- {
        }
        merged.HVG <- {
        }
        merged.counts <- {
        }
        merged.normalized <- {
        }
        merged.dataset.labels <- {
        }
        merged.origins <- {
        }
        merged.experiments <- {
        }
        merged.cluster.names <- {
        }
        union.genes.filter <- {
        }
        dataset.genes <- NA
        sample.index <- 1
        for (sample.index in seq(length(environment$datasets))) {
            dataset <- environment$datasets[sample.index]
            origin <- environment$origins[sample.index]
            experiment <- environment$experiments[sample.index]
            data.files <- list.files(path = file.path(path, dataset), pattern = "clustering.RData",
                full.names = T, recursive = recursive)
            file.index <- 1
            if (length(data.files) > 1) {
                for (row in seq(length(data.files))) print.message(row, dirname(dirname(data.files[row])))
                file.index <- as.numeric(readline(prompt = "Select clustering: "))
            }
            print.message("Loading", dirname(dirname(data.files[file.index])),
                "\n")
            load(data.files[file.index])
            if (length(merged.clustering) == 0) {
                min <- 0
            } else {
                min <- max(merged.clustering)
            }

            # names(clustering) str(clustering$shuffled.membership)
            # str(clustering$shuffled.membership[[1]])
            # apply(clustering$shuffled.membership[[1]]$memberships,1,table) membership
            # = clustering$shuffled.membership[[2]]$memberships[2,] limma.diff = {} for(
            # cluster in seq(length(unique(membership))) ) { print.message('cluster
            # =',cluster)

            # group = membership == cluster group = factor(group)

            # diff.exp = getDE.limma( Y = normalized, group = group, filter = F )
            # diff.exp = diff.exp[order(diff.exp$logFC,decreasing=T),] diff.exp =
            # data.frame(gene=rownames(diff.exp),logFC=diff.exp$logFC,fold=exp(diff.exp$logFC),QValue=diff.exp$QValue,PValue=diff.exp$PValue,AveExpr=diff.exp$AveExpr)
            # limma.diff =
            # rbind(limma.diff,data.frame(cluster=cluster,diff.exp[,c(1,3,4)])) }
            # limma.diff = limma.diff[order(limma.diff$fold,decreasing=T),]
            # head(limma.diff) lcdif[lcdif$gene=='Foxp3',] lndif[lndif$gene=='Foxp3',]
            # res={} for(i in seq(length(unique(lcdif$cluster)))){ for(j in
            # seq(length(unique(lndif$cluster)))){ lc=lcdif[lcdif$cluster == i,]
            # ln=lndif[lndif$cluster == j,] match = match(ln$gene,lc$gene)
            # head(lc[match,]) head(ln)
            # ccc=cor.test(lc$fold[match],ln$fold,method='spearman') res =
            # rbind(res,data.frame(ccc$estimate,ccc$p.value)) } } head(res)
            # res[order(res$ccc.estimate),]

            merged.clustering <- c(merged.clustering, clustering$membership +
                min)
            merged.original.clustering <- c(merged.original.clustering, clustering$membership)

            # index = 1 merged.shuffled.clustering[[]] for (index in
            # seq(length(clustering$shuffled.membership))) { v =
            # clustering$shuffled.membership[[index]] shuffled.membership =
            # v$memberships[2,] if (length(shuffled.membership) ==
            # length(clustering$membership)) }

            load(file.path(dirname(data.files[file.index]), "main.all.diff.exp.RData"))
            if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] ==
                T) {
                print.message("Converting from human to mouse genes")
                human.genes <- toupper(limma.all$gene)
                mouse.gene.names <- convertHumanGeneList(human.genes)
                if (length(dataset.genes) > 1 || !is.na(dataset.genes))
                  mouse.gene.names <- mouse.gene.names[mouse.gene.names[, 2] %in%
                    dataset.genes, ]
                limma.all$gene <- mouse.gene.names[match(human.genes, mouse.gene.names[,
                  1]), 2]
                limma.all <- limma.all[!is.na(limma.all$gene), ]
                intersect.genes <- intersect(unique(limma.all$gene), unique(merged.diff.exp$gene))
                limma.all <- limma.all[limma.all$gene %in% intersect.genes,
                  ]
                merged.diff.exp <- merged.diff.exp[merged.diff.exp$gene %in%
                  intersect.genes, ]
            }
            limma.all <- cbind(dataset, origin, experiment, limma.all)
            merged.diff.exp <- rbind(merged.diff.exp, limma.all)
            base.path <- dirname(dirname(dirname(data.files[file.index])))
            load(file.path(base.path, "data/HVG.RData"))
            if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] ==
                T) {
                human.genes <- toupper(HVG)
                mouse.gene.names <- convertHumanGeneList(human.genes)
                if (length(dataset.genes) > 1 || !is.na(dataset.genes))
                  mouse.gene.names <- mouse.gene.names[mouse.gene.names[, 2] %in%
                    dataset.genes, ]
                HVG <- unique(mouse.gene.names[, 2])
            }
            merged.HVG <- unique(c(merged.HVG, HVG))
            load(file.path(base.path, "data/data.RData"))
            colnames(counts) <- colnames(normalized) <- rep(dataset, ncol(normalized))
            if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] ==
                T) {
                human.genes <- toupper(rownames(counts))
                mouse.gene.names <- convertHumanGeneList(human.genes)
                if (length(dataset.genes) > 1 || !is.na(dataset.genes))
                  mouse.gene.names <- mouse.gene.names[mouse.gene.names[, 2] %in%
                    dataset.genes, ]

                rownames(counts) <- mouse.gene.names[match(human.genes, mouse.gene.names[,
                  1]), 2]
                counts <- counts[!is.na(rownames(counts)), ]
                intersect.genes <- intersect(rownames(counts), rownames(merged.counts))
                counts <- counts[rownames(counts) %in% intersect.genes, ]
                merged.counts <- merged.counts[rownames(merged.counts) %in%
                  intersect.genes, ]
                counts <- counts[match(rownames(merged.counts), rownames(counts)),
                  ]

                rownames(normalized) <- mouse.gene.names[match(human.genes,
                  mouse.gene.names[, 1]), 2]
                normalized <- normalized[!is.na(rownames(normalized)), ]
                intersect.genes <- intersect(rownames(normalized), rownames(merged.normalized))
                normalized <- normalized[rownames(normalized) %in% intersect.genes,
                  ]
                merged.normalized <- merged.normalized[rownames(merged.normalized) %in%
                  intersect.genes, ]
                normalized <- normalized[match(rownames(merged.normalized),
                  rownames(normalized)), ]
            }
            if (length(merged.counts) > 1 && sum(rownames(merged.counts) !=
                rownames(counts)) > 0)
                stop("Feature genes mismatch - need to correct dataset binding matching")
            merged.counts <- cbind(merged.counts, counts)
            merged.normalized <- cbind(merged.normalized, normalized)
            merged.dataset.labels <- c(merged.dataset.labels, rep(dataset, ncol(normalized)))
            merged.origins <- c(merged.origins, rep(origin, ncol(normalized)))
            merged.experiments <- c(merged.experiments, rep(experiment, ncol(normalized)))
            load(file.path(dirname(data.files[file.index]), "cluster.names.RData"))
            merged.cluster.names <- c(merged.cluster.names, cluster.names)
            if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] ==
                T) {
                human.genes <- names(genes.filter)
                mouse.gene.names <- convertHumanGeneList(human.genes)
                if (length(dataset.genes) > 1 || !is.na(dataset.genes))
                  mouse.gene.names <- mouse.gene.names[mouse.gene.names[, 2] %in%
                    dataset.genes, ]
                names(genes.filter) <- mouse.gene.names[match(human.genes, mouse.gene.names[,
                  1]), 2]
                genes.filter <- genes.filter[!is.na(names(genes.filter))]
                intersect.genes <- intersect(names(genes.filter), names(union.genes.filter))
                genes.filter <- genes.filter[names(genes.filter) %in% intersect.genes]
                union.genes.filter <- union.genes.filter[names(union.genes.filter) %in%
                  intersect.genes]
                genes.filter <- genes.filter[match(names(union.genes.filter),
                  names(genes.filter))]
            }
            if (length(union.genes.filter) == 0) {
                union.genes.filter <- genes.filter
            } else {
                union.genes.filter <- union.genes.filter | genes.filter
            }

            dataset.genes <- rownames(merged.counts)[apply(merged.counts, 1,
                sd) > 0]
        }

        # merged.HVG[!merged.HVG%in%rownames(normalized[genes.filter,])]
        # sum(genes.filter)
        merged.HVG <- merged.HVG[merged.HVG %in% rownames(counts)]
        union.genes.filter <- union.genes.filter[names(union.genes.filter) %in%
            rownames(counts)]

        print.message("dim(merged.normalized)")
        print(dim(merged.normalized))
        print(table(merged.dataset.labels))
        print(table(merged.origins))
        print(table(merged.experiments))
        print(table(merged.clustering))

        genes.filter <- union.genes.filter
        counts <- merged.counts
        normalized <- merged.normalized
        dataset.labels <- merged.dataset.labels
        origins <- merged.origins
        experiments <- merged.experiments
        HVG <- merged.HVG
        clustering <- merged.clustering
        cluster.names <- merged.cluster.names
        print.message("Saving")
        save(genes.filter, counts, normalized, dataset.labels, origins, experiments,
            HVG, clustering, merged.diff.exp, merged.original.clustering, cluster.names,
            file = cache)

        end(t)
    }
    counts <- counts[genes.filter, ]
    normalized <- normalized[genes.filter, ]
    HVG <- HVG[HVG %in% rownames(normalized)]

    environment$counts <- counts
    environment$normalized <- normalized
    environment$genes <- rownames(normalized)
    environment$datasets <- colnames(environment$normalized)
    environment$dataset.labels <- dataset.labels
    environment$origins <- origins
    environment$experiments <- experiments
    environment$confounders <- data.frame(nUMI = colSums(counts), nGenes = colSums(counts >
        0))
    environment$nsamples <- ncol(counts)
    environment$HVG <- HVG
    environment$clustering <- {
    }
    environment$clustering$membership <- clustering
    environment$original.clustering <- merged.original.clustering
    environment$diff.exp <- merged.diff.exp
    environment$cluster.names <- apply(cbind(environment$datasets, cluster.names),
        1, function(v) paste(v, collapse = "_"))

    return(environment)
}
