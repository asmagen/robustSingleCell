
#' Get/Set Cluster Names by Marker Gene Expression
#'
#' \code{get.cluster.names} uses predefined marker genes to assign clusters with
#' putative cell type or state labels.
#'
#' @param environment \code{environment} object
#' @param types data frame associating cell type or state with marker genes
#' @param min.fold minimum fold change to consider a marker as overexpressed
#' @param max.Qval maximum FDR q value to consider a marker as overexpressed
#' @param print whether to print output calculations
#' @return \code{get.cluster.names} returns a vector containing assigned cluster
#' name labels
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
#' }
get.cluster.names <- function(environment, types, min.fold = 1.25, max.Qval = 0.1,
    print = T) {

    precomputed <- readRDS(file.path(environment$res.data.path, paste("main", "all",
        "diff.exp.rds", sep = ".")))
    limma.all <- precomputed$limma.all
    if (print)
        print(summary(limma.all))
    diff.exp <- limma.all[limma.all$fold > min.fold & limma.all$QValue < max.Qval,
        ]
    if (print)
        print(summary(diff.exp))

    cluster <- 1
    cluster.names <- array("Unknown", environment$clustering$nclusters)
    for (cluster in seq(environment$clustering$nclusters)) {
        cluster.diff <- diff.exp[diff.exp$cluster == cluster, ]
        cluster.name <- get.cluster.names.with.diff(cluster.diff, types, print)
        if (print)
            print(cluster.name)
        if (!(length(cluster.names) == 1 && is.na(cluster.name)))
            cluster.names[cluster] <- paste(cluster.name, collapse = "_")
    }
    cluster.names

    for (name in unique(cluster.names)) {
        match <- cluster.names == name
        if (sum(match) > 1)
            cluster.names[match] <- paste(name, seq(sum(match)), sep = "_")
    }

    return(cluster.names)
}

get.cluster.names.with.diff <- function(cluster.diff, types, print) {

    types$gene <- as.vector(types$gene)
    minimum.genes.to.qualify <- table(types$type)/2

    expression <- cbind(types, cluster.diff[match(types$gene, cluster.diff$gene),
        ])
    if (print)
        print(expression[!is.na(expression$fold), ])
    table <- sort(table(expression$type[!is.na(expression$fold)]) - minimum.genes.to.qualify,
        decreasing = T)
    if (print)
        print(table)
    if (sum(table > 0) == 0)
        return("Unknown")
    table <- table[table > 0]
    cluster.name <- names(table)

    return(cluster.name)
}

#' Set Cluster Names in Environment
#'
#' \code{set.cluster.names} saves the cluster names in storage and in the \code{environment} object
#'
#' @param names cluster names defined in get.cluster.names
#' @return \code{set.cluster.names} returns an \code{environment} object coded
#' with cluster names
#' @export
#' @describeIn get.cluster.names set annotations to clusters
set.cluster.names <- function(environment, names) {

    cluster.name.map <- data.frame(id = seq(length(names)), name = names)
    environment$cluster.names <- cluster.names <- names[environment$clustering$membership]
    saveRDS(list(cluster.names = cluster.names, cluster.name.map = cluster.name.map),
        file = file.path(environment$res.data.path, "cluster.names.rds"))
    utils::write.csv(cluster.name.map, file = file.path(environment$work.path, "cluster.name.map.csv"))
    print(table(environment$cluster.names))
    return(environment)
}

load.cluster.names <- function(environment) {
    precomputed <- readRDS(file.path(environment$res.data.path, "cluster.names.rds"))
    environment$cluster.names <- precomputed$cluster.names
    print(table(environment$cluster.names))
    return(environment)
}

remove.cluster.names <- function(environment) {

    environment$cluster.names <- environment$clustering$membership

    return(environment)
}

#' Remove selected clusters
#'
#' Remove selected clusters from the environment object.
#'
#' @param environment The \code{environment} object
#' @param remove.clusters A character vector of the clusters to be removed
#' @return An environment object with selected clusters removed
#' @export
#' @examples
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- filter_cluster_data(LCMV1, "1")
filter_cluster_data <- function(environment, remove.clusters) {
    membership <- as.vector(environment$clustering$membership)
    keep <- !membership %in% remove.clusters
    filter.data(environment, keep)
}

filter.data <- function(environment, keep) {
    data.file <- file.path(environment$baseline.data.path, "data.rds")
    precomputed <- readRDS(data.file)
    genes.filter <- precomputed$genes.filter
    counts <- precomputed$counts
    normalized <- precomputed$normalized
    dataset.labels <- precomputed$dataset.labels
    origins <- precomputed$origins
    experiments <- precomputed$experiments
    rm(precomputed)

    file.rename(data.file, paste(data.file, format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),
        sep = "---"))
    counts <- counts[, keep]
    genes.filter <- genes.filter & apply(counts, 1, stats::var) > 0
    normalized <- normalized[, keep]
    dataset.labels <- dataset.labels[keep]
    origins <- origins[keep]
    experiments <- experiments[keep]
    unlink(environment$baseline.data.path, recursive = T, force = T)
    dir.create(environment$baseline.data.path)
    cache <- file.path(environment$baseline.data.path, "data.rds")
    saveRDS(list(genes.filter = genes.filter, counts = counts, normalized = normalized,
        dataset.labels = dataset.labels, origins = origins, experiments = experiments),
        file = cache)
    file.rename(environment$work.path, paste(environment$work.path, "pre.filter",
        format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"), sep = "_"))
}

filter.robust.clusters <- function(environment, robust.clusters) {

    precomputed <- readRDS(file.path(environment$baseline.data.path, "preclustered.datasets.rds"))
    genes.filter <- precomputed$genes.filter
    counts <- precomputed$counts
    normalized <- precomputed$normalized
    dataset.labels <- precomputed$dataset.labels
    origins <- precomputed$origins
    experiments <- precomputed$experiments
    HVG <- precomputed$HVG
    clustering <- precomputed$clustering
    merged.diff.exp <- precomputed$merged.diff.exp
    merged.original.clustering <- precomputed$merged.original.clustering
    rm(precomputed)

    membership <- as.vector(environment$clustering$membership)
    keep <- membership %in% robust.clusters

    counts <- counts[, keep]
    normalized <- normalized[, keep]
    genes.filter <- genes.filter & apply(counts, 1, stats::var) > 0
    dataset.labels <- dataset.labels[keep]
    origins <- origins[keep]
    experiments <- experiments[keep]
    HVG <- NA
    clustering <- clustering[keep]
    merged.original.clustering <- merged.original.clustering[keep]
    merged.diff.exp <- NA

    dir <- dirname(environment$work.path)
    new.dir <- file.path(dirname(dir), paste("filtered", basename(dir), sep = "_"),
        "data")
    dir.create(new.dir, recursive = T)

    saveRDS(list(genes.filter = genes.filter, counts = counts, normalized = normalized,
        dataset.labels = dataset.labels, origins = origins, experiments = experiments,
        HVG = HVG, clustering = clustering, merged.diff.exp = merged.diff.exp, merged.original.clustering = merged.original.clustering),
        file = file.path(new.dir, "preclustered.datasets.rds"))
}
