#' Read 10X Data
#'
#' Load sparse data matrices from 10X genomics.
#'
#' @param path Path to directory containing matrix.mtx, genes.tsv, and barcodes.tsv
#' @return a matrix of genes by cells
read_10x_data <- function(path) {
    barcode.path <- list.files(path, pattern = "barcodes.tsv", full.names = T)
    features.path <- list.files(path, pattern = "genes.tsv", full.names = T)
    matrix.path <- list.files(path, pattern = "matrix.mtx", full.names = T)
    mat <- Matrix::readMM(file = matrix.path)
    feature.names <- utils::read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names <- utils::read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
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
#' @examples
#' LCMV1 <- setup_LCMV_example()
#' data.path <- system.file("extdata/LCMV1_small.txt", package = "robustSingleCell")
#' # name of list should be the same as LCMV1$datasets
#' raw_LCMV1 <- as.matrix(read.table(data.path, check.names = FALSE))
#' LCMV1 <- read.data(LCMV1,
#' raw.data.matrices = list(LCMV1 = raw_LCMV1),
#' min.genes.per.cell = 100,
#' max.genes.per.cell.quantile = 1,
#' max.UMIs.per.cell.quantile = 1,
#' min.cells.per.gene = 1)
read.data <- function(environment, genome = "mm10", min.genes.per.cell = 500, max.genes.per.cell.quantile = 0.98,
    max.UMIs.per.cell.quantile = 0.98, min.cells.per.gene = 1, max.mitochondrial.frac = 0.1,
    max.ribosomal.frac = NA, cell.filters = NA, raw.data.matrices = NA, rerun = F,
    subsample = NULL, seed = 0) {

    cache <- file.path(environment$baseline.data.path, "data.rds")

    if (!rerun & file.exists(cache)) {
        print.message("Loading precomputed")
        precomputed <- readRDS(cache)
        genes.filter <- precomputed$genes.filter
        counts <- precomputed$counts
        normalized <- precomputed$normalized
        dataset.labels <- precomputed$dataset.labels
        dataset_ids <- precomputed$dataset_ids
        origins <- precomputed$origins
        experiments <- precomputed$experiments
        criteria <- precomputed$criteria
        rm(precomputed)
    } else {

        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"))
        on.exit(end(t))

        merged <- NA
        dataset <- environment$datasets[1]
        merged.dataset.labels <- merged.origins <- merged.experiments <- merged.criteria <- {
        }
        set.seed(seed)
        for (dataset in environment$datasets) {
            if (is.na(raw.data.matrices)) {
                print.message("Loading", dataset)
                measurements <- read_10x_data(path = file.path(environment$data.path,
                  dataset))

                measurements <- as.matrix(measurements[, Matrix::colSums(measurements >
                  0) >= min.genes.per.cell])
                if (!is.null(subsample) && subsample < ncol(measurements)) {
                    measurements <- measurements[, sample(seq(ncol(measurements)), subsample)]
                }
            } else {
                print.message("Using input", dataset)
                measurements <- raw.data.matrices[[dataset]]
                stopifnot(is.matrix(measurements))
            }


            colnames(measurements) <- rep(environment$datasets[environment$datasets ==
                dataset], ncol(measurements))
            dataset.labels <- rep(paste(unique(environment$origins)[environment$datasets ==
                dataset], " (", unique(environment$experiments)[environment$datasets == dataset],
                ")", sep = ""), ncol(measurements))
            origins <- rep(unique(environment$origins)[environment$datasets == dataset],
                ncol(measurements))
            experiments <- rep(unique(environment$experiments)[environment$datasets == dataset],
                ncol(measurements))
            cat(dim(measurements))
            # corner(measurements) measurements = measurements[,measurements['Cd4',]>0]

            data <- data.frame(nUMI = colSums(measurements), nGenes = colSums(measurements >
                0), Ribosomal = colMeans(measurements[get.ribo.genes(rownames(measurements)),
                ]), Mitochondrial = colMeans(measurements[get.mito.genes(rownames(measurements)),
                ]), Ribosomal.frac = colSums(measurements[get.ribo.genes(rownames(measurements)),
                ])/colSums(measurements), Mitochondrial.frac = colSums(measurements[get.mito.genes(rownames(measurements)),
                ])/colSums(measurements))
            print(utils::head(data))

            grDevices::pdf(file.path(environment$baseline.work.path, paste(dataset,
                "pre.filter.dataset.stats.pdf", sep = ".")))
            for (ind1 in seq(length(colnames(data)) - 1)) {
                for (ind2 in (ind1 + 1):(length(colnames(data)))) {
                  v1 <- colnames(data)[ind1]
                  v2 <- colnames(data)[ind2]
                  print(ggplot(data, aes_string(v1, v2)) + geom_point())
                }
            }
            grDevices::dev.off()

            print.message("Original dimensions")
            print(dim(measurements))

            genes.per.cell <- colSums(measurements > 0)
            print.message("Original genes.per.cell")
            print(summary(genes.per.cell))
            print(stats::quantile(genes.per.cell, seq(0.01, 1, 0.01)))
            UMIs.per.cell <- colSums(measurements)
            print.message("Original UMIs.per.cell")
            print(summary(UMIs.per.cell))
            print(stats::quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))
            max.genes.per.cell <- stats::quantile(genes.per.cell, max.genes.per.cell.quantile)
            print.message("max.genes.per.cell", max.genes.per.cell)
            max.UMIs.per.cell <- stats::quantile(UMIs.per.cell, max.UMIs.per.cell.quantile)
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
                    ] > 0) == dataset.filters$expressed[dataset.filters$gene == gene])  #check if matching condition
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
            measurements <- measurements[, criteria] ###!!!
            dataset.labels <- dataset.labels[criteria]
            origins <- origins[criteria]
            experiments <- experiments[criteria]
            genes.per.cell <- colSums(measurements > 0)
            print.message("Filtered genes.per.cell")
            print(summary(genes.per.cell))
            print(stats::quantile(genes.per.cell, seq(0.01, 1, 0.01)))
            UMIs.per.cell <- colSums(measurements)
            print.message("Filtered UMIs.per.cell")
            print(summary(UMIs.per.cell))
            print(stats::quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))

            data <- data.frame(nUMI = colSums(measurements), nGenes = colSums(measurements >
                0), Ribosomal = colMeans(measurements[get.ribo.genes(rownames(measurements)),
                ]), Mitochondrial = colMeans(measurements[get.mito.genes(rownames(measurements)),
                ]), Ribosomal.frac = colSums(measurements[get.ribo.genes(rownames(measurements)),
                ])/colSums(measurements), Mitochondrial.frac = colSums(measurements[get.mito.genes(rownames(measurements)),
                ])/colSums(measurements))

            print(utils::head(data))
            grDevices::pdf(file.path(environment$baseline.work.path, paste(dataset,
                "post.filter.dataset.stats.pdf", sep = ".")))
            for (ind1 in seq(length(colnames(data)) - 1)) {
                for (ind2 in (ind1 + 1):(length(colnames(data)))) {
                  v1 <- colnames(data)[ind1]
                  v2 <- colnames(data)[ind2]
                  print(ggplot(data, aes_string(v1, v2)) + geom_point())
                }
            }
            grDevices::dev.off()

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
        print(stats::quantile(genes.per.cell, seq(0.01, 1, 0.01)))
        UMIs.per.cell <- colSums(counts[genes.filter, ])
        print.message("Aggregated dataset UMIs.per.cell")
        print(summary(UMIs.per.cell))
        print(stats::quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))

        rownames <- data.frame(old = rownames(counts))
        if (sum(duplicated(rownames(counts))) > 0) {
            gene <- unique(rownames(counts)[duplicated(rownames(counts))])[2]
            for (gene in unique(rownames(counts)[duplicated(rownames(counts))])) {
                print(gene)
                indices <- which(rownames(counts) == gene)
                print(rownames(counts)[indices])
                rownames(counts)[indices] <- paste(rownames(counts)[indices], seq(length(indices)),
                  sep = ".")
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

        t <- start(file.path(environment$work.path, "tracking"), append = T, split = T)
        on.exit(end(t), add = T)
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

        saveRDS(list(genes.filter = genes.filter, counts = counts, normalized = normalized,
            dataset.labels = dataset.labels, dataset_ids = colnames(normalized), origins = origins, experiments = experiments,
            criteria = criteria), file = cache)

    }

    counts <- counts[genes.filter, ]
    normalized <- normalized[genes.filter, ]

    environment$genes.filter <- genes.filter
    environment$counts <- counts
    environment$normalized <- normalized
    environment$genes <- rownames(normalized)
    environment$dataset_ids <- colnames(normalized)
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
#' }
read.preclustered.datasets <- function(environment, path = NA, recursive = T, rerun = F) {

    cache <- file.path(environment$baseline.data.path, "preclustered.datasets.rds")

    if (!rerun & file.exists(cache)) {
        print.message("Loading precomputed")
        precomputed <- readRDS(cache)
        counts <- precomputed$counts
        normalized <- precomputed$normalized
        HVG <- precomputed$HVG
        genes.filter <- precomputed$genes.filter
        dataset.labels <- precomputed$dataset.labels
        origins <- precomputed$origins
        origin_id <- precomputed$origin_id
        experiments <- precomputed$experiments
        clustering <- precomputed$clustering
        merged.original.clustering <- precomputed$merged.original.clustering
        merged.diff.exp <- precomputed$merged.diff.exp
        cluster.names <- precomputed$cluster.names
        dataset_ids <- precomputed$dataset_ids
        rm(precomputed)
    } else {

        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"), split = T)
        on.exit(end(t))
        if (is.na(path))
            path <- dirname(environment$baseline.work.path)
        environment$datasets <- environment$datasets
        environment$experiments <- environment$experiments
        environment$origin_id <- environment$origins 
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
        merged.dataset.ids <- {
        }
        dataset.genes <- NA
        sample.index <- 1
        for (sample.index in seq(length(environment$datasets))) {
            dataset <- environment$datasets[sample.index]
            origin <- environment$origin_id[sample.index]
            experiment <- environment$experiments[sample.index]
            data.files <- list.files(path = file.path(path, dataset), pattern = "clustering.rds",
                full.names = T, recursive = recursive)
            file.index <- 1
            if (length(data.files) > 1) {
                for (row in seq(length(data.files))) print.message(row, dirname(dirname(data.files[row])))
                file.index <- as.numeric(readline(prompt = "Select clustering: "))
            }
            print.message("Loading", dirname(dirname(data.files[file.index])), "\n")
            clustering <- readRDS(data.files[file.index])

            if (length(merged.clustering) == 0) {
                min <- 0
            } else {
                min <- max(merged.clustering)
            }

            merged.clustering <- c(merged.clustering, clustering$membership + min)
            merged.original.clustering <- c(merged.original.clustering, clustering$membership)

            precomputed <- readRDS(file.path(dirname(data.files[file.index]), "main.all.diff.exp.rds"))
            limma.all <- precomputed$limma.all
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
                limma.all <- limma.all[limma.all$gene %in% intersect.genes, ]
                merged.diff.exp <- merged.diff.exp[merged.diff.exp$gene %in% intersect.genes,
                  ]
            }
            limma.all <- cbind(dataset, origin, experiment, limma.all)
            merged.diff.exp <- rbind(merged.diff.exp, limma.all)
            base.path <- dirname(dirname(dirname(data.files[file.index])))
            HVG <- readRDS(file.path(base.path, "data/HVG.rds"))
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
            precomputed <- readRDS(file.path(base.path, "data/data.rds"))
            counts <- precomputed$counts
            normalized <- precomputed$normalized
            genes.filter <- precomputed$genes.filter
            dataset.labels <- precomputed$dataset.labels
            dataset_ids <- precomputed$dataset_ids
            origin_id <- precomputed$origin_id
            rm(precomputed)

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
                merged.counts <- merged.counts[rownames(merged.counts) %in% intersect.genes,
                  ]
                counts <- counts[match(rownames(merged.counts), rownames(counts)),
                  ]

                rownames(normalized) <- mouse.gene.names[match(human.genes, mouse.gene.names[,
                  1]), 2]
                normalized <- normalized[!is.na(rownames(normalized)), ]
                intersect.genes <- intersect(rownames(normalized), rownames(merged.normalized))
                normalized <- normalized[rownames(normalized) %in% intersect.genes,
                  ]
                merged.normalized <- merged.normalized[rownames(merged.normalized) %in%
                  intersect.genes, ]
                normalized <- normalized[match(rownames(merged.normalized), rownames(normalized)),
                  ]
            }
            if (length(merged.counts) > 1 && sum(rownames(merged.counts) != rownames(counts)) >
                0)
                stop("Feature genes mismatch - need to correct dataset binding matching")
            merged.counts <- cbind(merged.counts, counts)
            merged.normalized <- cbind(merged.normalized, normalized)
            merged.dataset.ids <- c(merged.dataset.ids, dataset_ids)
            merged.dataset.labels <- c(merged.dataset.labels, dataset.labels)
            merged.origins <- c(merged.origins, rep(origin, ncol(normalized)))
            merged.experiments <- c(merged.experiments, rep(experiment, ncol(normalized)))
            precomputed <- readRDS(file.path(dirname(data.files[file.index]), "cluster.names.rds"))
            cluster.names <- precomputed$cluster.names
            rm(precomputed)

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
                genes.filter <- genes.filter[match(names(union.genes.filter), names(genes.filter))]
            }
            if (length(union.genes.filter) == 0) {
                union.genes.filter <- genes.filter
            } else {
                union.genes.filter <- union.genes.filter | genes.filter
            }

            dataset.genes <- rownames(merged.counts)[apply(merged.counts, 1, stats::sd) >
                0]
        }

        # merged.HVG[!merged.HVG%in%rownames(normalized[genes.filter,])]
        # sum(genes.filter)
        merged.HVG <- merged.HVG[merged.HVG %in% rownames(counts)]
        union.genes.filter <- union.genes.filter[names(union.genes.filter) %in% rownames(counts)]

        print.message("dim(merged.normalized)")
        print(dim(merged.normalized))
        print(table(merged.dataset.labels))
        print(table(merged.origins))
        print(table(merged.experiments))
        print(table(merged.clustering))
        print(table(merged.dataset.ids))

        genes.filter <- union.genes.filter
        counts <- merged.counts
        normalized <- merged.normalized
        dataset.labels <- merged.dataset.labels
        dataset_ids <- merged.dataset.ids
        origins <- merged.origins
        experiments <- merged.experiments
        HVG <- merged.HVG
        clustering <- merged.clustering
        cluster.names <- merged.cluster.names
        print.message("Saving")
        saveRDS(list(genes.filter = genes.filter, counts = counts, normalized = normalized,
            dataset.labels = dataset.labels, dataset_ids = dataset_ids, origins = origins, origin_id = origin_id, experiments = experiments,
            HVG = HVG, clustering = clustering, merged.diff.exp = merged.diff.exp,
            merged.original.clustering = merged.original.clustering, cluster.names = cluster.names),
            file = cache)

    }
    counts <- counts[names(genes.filter)[genes.filter], ]
    normalized <- normalized[names(genes.filter)[genes.filter], ]
    HVG <- HVG[HVG %in% rownames(normalized)]

    environment$counts <- counts
    environment$normalized <- normalized
    environment$genes <- rownames(normalized)
    environment$genes.filter <- genes.filter
    environment$datasets <- colnames(environment$normalized)
    environment$dataset.labels <- dataset.labels
    environment$dataset_ids <- dataset_ids
    environment$origins <- origins
    environment$origin_id <- origin_id
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
