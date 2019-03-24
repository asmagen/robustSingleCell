#' Set up LCMV1 example
#'
#' Set the path in the LCMV1_small object
#' @return An environment object containing the LCMV1_small data
#' @export
#' @examples
#' setup_LCMV1_example()
setup_LCMV1_example <- function() {
    data.path <- system.file("extdata", package = "robustSingleCell")
    LCMV1_small <- initialize.project(datasets = "LCMV1",
                                origins = "CD44+ cells",
                                experiments = "Rep1",
                                data.path = data.path,
                                work.path = tempdir())
    data.path <- system.file("extdata/LCMV1_small.txt", package = "robustSingleCell")
    raw_LCMV1 <- as.matrix(read.table(data.path, check.names = FALSE))
    LCMV1_small <- read.data(LCMV1_small,
                       raw.data.matrices = list(LCMV1 = raw_LCMV1),
                       min.genes.per.cell = 100,
                       max.genes.per.cell.quantile = 1,
                       max.UMIs.per.cell.quantile = 1,
                       min.cells.per.gene = 1)
    LCMV1_small
}


#' Download Example Dataset
#'
#' Download two replicates of CD44+ T cell 10X scRNAseq data sets (Ciucci 2018).
#'
#' @param base_dir Full path to a directory where data and analysis will be stored
#' @export
#' @return 1 if download fails and 0 if succeeds
#' @examples
#' download_LCMV(tempdir())
download_LCMV <- function(base_dir) {
    base_dir_path <- file.path(base_dir, "LCMV")
    dir.create(base_dir_path, showWarnings = F)
    check_1_avail <- GEOquery::getGEOSuppFiles("GSM3423794")
    check_2_avail <- GEOquery::getGEOSuppFiles("GSM3423795")
    if (is.null(check_1_avail) | is.null(check_2_avail)) {
        cat("Example files cannot be downloaded. \nPlease check your network connection.\n")
        return(1)
    } else {
        GEOquery::getGEOSuppFiles("GSM3423794", baseDir = base_dir_path)
        GEOquery::getGEOSuppFiles("GSM3423795", baseDir = base_dir_path)
        file.rename(file.path(base_dir_path, "GSM3423795"),
                    file.path(base_dir_path, "LCMV1"))
        file.rename(file.path(base_dir_path, "GSM3423794"),
                    file.path(base_dir_path, "LCMV2"))
        cat(paste0("Data saved at ", base_dir_path, "\n"))
        return(0)
    }
}

capwords <- function(s, strict = FALSE) {
    s <- tolower(s)
    cap <- function(s) paste(toupper(substring(s, 1, 1)), {
        s <- substring(s, 2)
        if (strict)
            tolower(s) else s
    }, sep = "", collapse = " ")
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

corner <- function(matrix, n = 5, m = 5) {
    print(matrix[seq(min(n, nrow(matrix))), seq(min(m, nrow(matrix)))])
}

get_slurm_out <- function(sjob) {

    Sys.sleep(1)
    queued <- length(system(paste("squeue -hn", sjob$jobname), intern = T)) > 0
    while (length(system(paste("squeue -hn", sjob$jobname), intern = T)) > 0) {
        Sys.sleep(1)
    }

    res_files <- paste0("results_", 0:(sjob$nodes - 1), ".RDS")
    tmpdir <- paste0("_rslurm_", sjob$jobname)
    missing_files <- setdiff(res_files, dir(path = tmpdir))

    if (length(missing_files) > 0) {
        missing_list <- paste(missing_files, collapse = ", ")
        warning(paste("The following files are missing:", missing_list))
    }

    res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
    if (length(res_files) == 0)
        return(NA)

    slurm_out <- lapply(res_files, readRDS)
    slurm_out <- do.call(c, slurm_out)
    slurm_out <- as.data.frame(do.call(rbind, slurm_out))

    return(slurm_out)
}

start <- function(track_dir_path, name = NA, append = F, split = F, print = T) {
    time <- get.time()
    label <- as.character(sys.calls()[[sys.nframe() - 1]])[1]
    if (!is.na(name))
        label <- paste(label, name, sep = ".")
    file <- paste(label, "txt", sep = ".")
    sink(file.path(track_dir_path, file), type = "output", append, split)
    if (!split && print)
        print(as.list(sys.calls())[seq(sys.nframe() - 1)])
    return(time)
}

end <- function(time = NA) {
    if (!is.na(time))
        elapsed.time(time)
    sink()
}

print.message <- function(...) {
    cat(cat(..., sep = " "), "\n", sep = "")
}

get.time <- function() {
    return(Sys.time())
}

elapsed.time <- function(time) {
    print(Sys.time() - time)
    cat("\n")
}

apply.by.group <- function(groups, values, ...) {
    group <- group.by(groups, values)
    for (fun in list(...)) {
        group <- sapply(group, fun)
    }
    return(group)
}

group.by <- function(groups, values) {
    return(split(x = as.vector(values), f = as.vector(groups)))
}

convertHumanGeneList <- function(x) {

    human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    genesV2 <- biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
        values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)

    humanx <- genesV2

    return(humanx)
}
