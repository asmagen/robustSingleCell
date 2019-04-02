#' Intialize the Project Environment
#'
#' Set up a project environment variable mapped to project results folder.
#'
#' @param datasets list of dataset code names
#' @param origins list of dataset tissue origin/condition full name
#' @param experiments list of experiment design annotations
#' @param data.path path to where the data is located
#' @param work.path path to where the analysis results are stored; optional, by
#' default, a temporary directory
#' @param marker.genes set of genes of interest for visualization purposes
#' @param clear.history whether you would like to remove any previous project by this name
#' @param analysis.label whether you would like to add a specific label to the analysis folder
#' @param convert.to.mouse.gene.symbols whether you are using human gene symbols and would like to convert them to mouse gene symbols
#' @return \code{environment} parameter containing file paths and experiment parameters
#' @export
#' @examples
#' data.path <- system.file("extdata", package = "robustSingleCell")
#' LCMV1_proj <- initialize.project(datasets = "LCMV1",
#'                             origins = "CD44+ cells",
#'                             experiments = "Rep1",
#'                             data.path = data.path,
#'                             work.path = tempdir())
initialize.project <- function(datasets, origins, experiments, data.path, work.path = NULL,
    marker.genes = NULL, clear.history = F, analysis.label = NULL, convert.to.mouse.gene.symbols = NULL) {

    options(width = 205)
    environment <- list()
    environment$datasets <- datasets
    environment$origins <- origins
    environment$experiments <- experiments
    environment$labels <- labels
    environment$data.path <- data.path
    environment$project <- ifelse(length(datasets) > 1, paste("merged", paste(datasets,
        collapse = "."), sep = "."), datasets)
    if(is.null(work.path)) {
        work.path <- tempdir()
        cat("Results will be saved at ", work.path, ".\n", sep = "")
    }
    if (!is.null(analysis.label))
        environment$project <- paste(analysis.label, environment$project, sep = "_")
    if (length(convert.to.mouse.gene.symbols) > 1 || !is.null(convert.to.mouse.gene.symbols))
        environment$convert.to.mouse.gene.symbols <- convert.to.mouse.gene.symbols
    environment$baseline.work.path <- environment$work.path <- file.path(work.path,
        environment$project)
    if (clear.history)
        unlink(environment$work.path, recursive = T, force = T)
    dir.create(environment$work.path, showWarnings = F, recursive = T, mode = "700")
    dir.create(file.path(environment$work.path, "tracking"), showWarnings = F, recursive = T,
        mode = "700")
    environment$baseline.data.path <- environment$res.data.path <- file.path(environment$work.path,
        "data")
    dir.create(environment$res.data.path, showWarnings = F, recursive = T, mode = "700")
    if (is.null(marker.genes)) {
        marker.genes <- unique(marker_genes$symbol)
    }
    environment$marker.genes <- marker.genes
    t <- start(file.path(environment$work.path, "tracking"))
    on.exit(end(t))
    print(environment)


    return(environment)
}
