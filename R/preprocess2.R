#' Identify Highly Variable Genes
#'
#' Get highly variable genes by Heteroscedasticity controlled binning of gene expression measurements within each dataset separately.
#'
#' @param sce A SingleCellExperiment object
#' @param min.mean minimum mean expression per gene
#' @param min.frac.cells minimum fraction of cells expressing each gene
#' @param min.dispersion.scaled minimum dispersion value
#' @param plot whether to plot mean vs dispersion plot
#' @param verbose whether to print diagnostic message
#' @return \code{SingleCellExperiment} parameter containing highly variable genes selection
#' @export
#' @importFrom Matrix rowSums rowMeans
get_variable_genes <- function(sce, min.mean = 0.05, min.frac.cells = 0,
                               min.dispersion.scaled = 1, plot = F, verbose = F) {
  normalized <- assay(sce, "normcounts")
  means <- apply(normalized, 1, function(v) log(mean(exp(v) - 1) + 1))
  frac.cells <- rowSums(normalized > 0)/ncol(normalized)
  vars <- apply(normalized, 1, function(v) log(stats::var(exp(v) - 1) + 1))
  dispersion <- apply(normalized, 1, function(v) log(stats::var(exp(v) - 1)/mean(exp(v) - 1)))
  dispersion[is.na(x = dispersion)] <- 0
  means[is.na(x = means)] <- 0

  if (plot) {
    plot.data <- data.frame(gene = names(means), means, dispersion)
    print(ggplot(plot.data, aes(x = means, y = dispersion, label = gene)) +
            geom_text(check_overlap = TRUE, size = 2))
    graphics::smoothScatter(means, dispersion)
    grDevices::dev.off()
  }

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

  if (verbose) {
    print.message("# qualifying genes", length(HVG))
    genes <- rownames(sce)[rowData(sce)$is_ribo | rowData(sce)$is_mito]
    print.message("Qualifying Ribosomal & Mitochondrial")
    print(genes[genes %in% HVG])
  }

  cat("# highly variable genes = ", length(HVG), "\n", sep = "")
  rowData(sce)$is_HVG <- rownames(sce) %in% HVG
  return(sce)
}


#' Compute Ribosomal Score
#'
#' Compute the activation level of ribosomal genes.
#'
#' @param sce A SingleCellExperiment object
#' @param control whether to subtract the score defined by technically similar genes
#' @param knn number of nearest neighbor
#' @param verbose Whether to print diagnostic messages
#' @param force Whether to rerun and overwrite already computed ribosomal score
#' @return a vector of ribosomal genes activation score
#' @export
ribosomal_score <- function(sce, control = T, knn = 10, verbose = F, force = F) {
  if (!is.null(colData(sce)$ribo_score) & !force) {
    print.message("Ribosomal score already computed. To rerun, use forece = T.")
    return(sce)
  }

  genes <- rowData(sce)$gene_name[rowData(sce)$is_ribo]
  if (verbose) {
    print.message("Using genes:")
    print(genes)
  }

  if (control) {
    score <- controlled_mean_score(sce, genes, knn)
  } else {
    score <- colMeans(assay(sce, "normcounts")[rowData(sce)$is_ribo, ])
  }
  colData(sce)$ribo_score <- score
  return(sce)
}

get_ribo_genes <- function(genes) {
  return(genes[c(grep("^Rpl", genes, ignore.case = T), grep("^Rps", genes, ignore.case = T))])
}

is_ribo <- function(genes) {
  grepl("^Rp[ls]", genes, ignore.case = T)
}


#' Compute Mitochondrial Score
#'
#' Compute the activation level of mitochondrial genes.
#'
#' @param sce A SingleCellExperiment object
#' @param control whether to subtract the score defined by technically similar genes
#' @param knn number of nearest neighbor
#' @param verbose Whether to print diagnostic messages
#' @param force Whether to rerun and overwrite already computed mitochondrial score
#' @return a vector of mitochondrial genes activation score
#' @export
#' @examples
#' LCMV1 <- setup_LCMV_example()
#' mitochondrial.score <- mitochondrial.score(LCMV1)
mitochondrial_score <- function(sce, control = F, knn = 10, force = F, verbose = F) {
  if (!is.null(colData(sce)$mito_score) & !force) {
    print.message("Mitochondrial score already computed. To rerun, use forece = T.")
    return(sce)
  }

  genes <- rowData(sce)$gene_name[rowData(sce)$is_mito]
  if (verbose) {
    print.message("Using genes:")
    print(genes)
  }

  if (control) {
    score <- controlled_mean_score(sce, genes, knn, verbose = verbose)
  } else {
    score <- Matrix::colMeans(assay(sce, "normcounts")[rowData(sce)$is_mito, ])
  }
  colData(sce)$mito_score <- score
  return(sce)
}

get_mito_genes <- function(genes) {
  return(genes[grep("^Mt-", genes, ignore.case = T)])
}

is_mito <- function(genes) {
  grepl("^Mt-", genes, ignore.case = T)
}

#' Compute Controlled Activation Score
#'
#' Compute mean gene signatures activation scores while controlling for technically similar genes.
#'
#' @param rbs A RobustSingleCell object
#' @param genes gene list upon which to calculate gene signature activate
#' @param knn number of nearest neighbors
#' @param exclude.missing.genes whether to exclude genes with missing values
#' @param constrain.cell.universe binary vector indicating in which subset of cells to calculate the gene signature activation. Default is all cells.
#' @return gene signature activation scores per cell
controlled_mean_score <- function(rbs, genes, knn = 10, exclude.missing.genes = T,
                                  constrain.cell.universe = NULL, verbose = F) {
  # similarly to
  # http://science.sciencemag.org/content/sci/suppl/2016/04/07/352.6282.189.DC1/Tirosh.SM.pdf/Seurat
  # to reduce association with library size or other technical

  if (is.null(constrain.cell.universe))
    constrain.cell.universe <- rep(T, ncol(rbs))
  if (knn > 0) {
    count_genes <- rowData(rbs)$gene_name
    normed_counts <- assay(rbs, "normcounts")
    genes <- rowData(rbs)$gene_name[match(capwords(genes), capwords(count_genes))]
    nExclude <- sum(is.na(genes))
    if (nExclude > 0) {
      if (exclude.missing.genes) {
        print.message("Excluding", nExclude, "out of", length(genes), "genes not found in the dataset")
        genes <- genes[genes %in% capwords(count_genes)]
      } else {
        print.message("Some signature genes are missing in dataset")
        return(NA)
      }
    }

    background.genes <- background_genes(rbs, foreground.genes = genes, knn, verbose = verbose)

    return(Matrix::colMeans(normed_counts[count_genes %in% genes, constrain.cell.universe]) -
             Matrix::colMeans(normed_counts[count_genes %in% background.genes, constrain.cell.universe]))
  } else {
    return(Matrix::colMeans(normed_counts[count_genes %in% genes, constrain.cell.universe]))
  }
}

background_genes <- function(rbs, foreground.genes, knn, verbose = F) {

  foreground.genes <- foreground.genes[foreground.genes %in% rowData(rbs)$gene_name]
  technically.similar.genes <- get_technically_similar_genes(rbs, knn)
  knns <- technically.similar.genes$knns
  technical.variables <- technically.similar.genes$technical.variables

  background.genes <- unique(setdiff(as.vector(knns[foreground.genes, ]), foreground.genes))
  if (verbose) {
    print.message("Head foreground.genes technical.variables")
    print(utils::head(technical.variables[rowData(rbs)$gene_name %in% foreground.genes, ]))
    print.message("Head background.genes technical.variables")
    print(utils::head(technical.variables[rowData(rbs)$gene_name %in% background.genes, ]))
    print.message("background.genes")
    print(background.genes)
  }

  return(background.genes)
}

get_technically_similar_genes <- function(rbs, knn = 10) {
  cache_path <- file.path(cache(rbs), "technically_similar_genes.rds")
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  print.message("Memory intensive step. Please allocate adequate memory for running.")
  normed_counts <- assay(rbs, "normcounts")
  row_genes <- rowData(rbs)$gene_name
  technical.variables <- data.frame(means = rowMeans(normed_counts),
                                    vars = apply(normed_counts, 1, stats::var))
  scaled.technical.variables <- apply(technical.variables, 2, scale)
  rownames(scaled.technical.variables) <- rownames(technical.variables)

  n <- nrow(scaled.technical.variables)
  dist_obj <- stats::dist(scaled.technical.variables)
  knns <- array("", c(length(row_genes), knn))
  rownames(knns) <- row_genes
  for (index in seq(length(row_genes))) {
    gene.dist <- get_dist(dist_obj, index, n, row_genes)
    knns[index, ] <- names(gene.dist[order(gene.dist)[2:(knn + 1)]])
  }

  ret <- list(knns = knns, technical.variables = technical.variables)
  saveRDS(ret, file = cache_path)
  print.message(paste0("Cached at ", cache_path, "\n"))
  return(ret)
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

