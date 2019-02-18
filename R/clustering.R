#' Get Cluster Names by Marker Gene Expression
#' 
#' Use predefined marker genes to assign clusters with putative cell type or state labels.
#'
#' @param environment \code{environment} object
#' @param types data frame associating cell type or state with marker genes
#' @param min.fold minimum fold change to consider a marker as overexpressed
#' @param max.Qval maximum FDR q value to consider a marker as overexpressed
#' @param print whether to print output calculations
#' @return vector containing assigned cluster name labels
#' @export
get.cluster.names <- function(environment, types, min.fold = 1.25, max.Qval = 0.1, 
    print = T) {
    
    load(file.path(environment$res.data.path, paste("main", "all", "diff.exp.RData", 
        sep = ".")))
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
#' Save the cluster names in storage and in the \code{environment} object
#' 
#' @param environment \code{environment} object
#' @param names cluster names defined in get.cluster.names
#' @return \code{environment} object coded with cluster names
#' @export
set.cluster.names <- function(environment, names) {
    
    cluster.name.map <- data.frame(id = seq(length(names)), name = names)
    environment$cluster.names <- cluster.names <- names[environment$clustering$membership]
    save(cluster.names, cluster.name.map, file = file.path(environment$res.data.path, 
        "cluster.names.RData"))
    write.csv(cluster.name.map, file = file.path(environment$work.path, "cluster.name.map.csv"))
    print(table(environment$cluster.names))
    return(environment)
}

load.cluster.names <- function(environment) {
    
    load(file.path(environment$res.data.path, "cluster.names.RData"))
    environment$cluster.names <- cluster.names
    print(table(environment$cluster.names))
    return(environment)
}

remove.cluster.names <- function(environment) {
    
    environment$cluster.names <- environment$clustering$membership
    
    return(environment)
}

filter.cluster.data <- function(environment, remove.clusters) {
    membership <- as.vector(environment$clustering$membership)
    keep <- !membership %in% remove.clusters
    filter.data(keep)
}

filter.data <- function(environment, keep) {
    data.file <- file.path(environment$baseline.data.path, "data.RData")
    load(data.file)
    file.rename(data.file, paste(data.file, format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"), 
        sep = "---"))
    counts <- counts[, keep]
    genes.filter <- genes.filter & apply(counts, 1, var) > 0
    normalized <- normalized[, keep]
    dataset.labels <- dataset.labels[keep]
    origins <- origins[keep]
    experiments <- experiments[keep]
    unlink(environment$baseline.data.path, recursive = T, force = T)
    dir.create(environment$baseline.data.path)
    cache <- file.path(environment$baseline.data.path, "data.RData")
    save(genes.filter, counts, normalized, dataset.labels, origins, experiments, 
        file = cache)
    file.rename(environment$work.path, paste(environment$work.path, "pre.filter", 
        format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"), sep = "_"))
}

filter.robust.clusters <- function(environment, robust.clusters) {
    
    load(file.path(environment$baseline.data.path, "preclustered.datasets.RData"))
    
    membership <- as.vector(environment$clustering$membership)
    keep <- membership %in% robust.clusters
    
    counts <- counts[, keep]
    normalized <- normalized[, keep]
    genes.filter <- genes.filter & apply(counts, 1, var) > 0
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
    
    save(genes.filter, counts, normalized, dataset.labels, origins, experiments, 
        HVG, clustering, merged.diff.exp, merged.original.clustering, file = file.path(new.dir, 
            "preclustered.datasets.RData"))
}
