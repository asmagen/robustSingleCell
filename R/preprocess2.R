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

#' Compute Cell Cycle Score
#'
#' Compute the activation of cell cycle genes defined in Kowalczyk, M. S. et al.
#'
#' @param sce A SingleCellExperiment object
#' @param knn number of nearest neighbor
#' @param cc_genes_genes cell cycle gene set. Default uses gene sets defined in Kowalczyk, M. S. et al. Single-cell RNA-seq reveals changes in cell cycle and differentiation programs upon aging of hematopoietic stem cells. Genome Res 25, 1860-1872, doi:10.1101/gr.192237.115 (2015).
#' @param verbose Whether to print diagnostic messages
#' @return a matrix of cell cycle genes activation scores (S, G2M and aggregated S/G2M scores, separately)
#' @export
cell_cycle_score <- function(sce, knn = 10, cc_cycle_genes = NULL, verbose = F) {
  cc.genes <- capwords(cell_cycle_genes)
  s.genes <- cc.genes[1:43]
  row_genes <- rowData(sce)$gene_name
  s.genes <- s.genes[s.genes %in% capwords(row_genes)]
  if (verbose) {
    print.message("S phase siganture genes used are ")
    print.message(s.genes)
  }
  g2m.genes <- cc.genes[44:98]
  g2m.genes <- g2m.genes[g2m.genes %in% capwords(row_genes)]
  if (verbose) {
    print.message("G2/M phase siganture genes used are ")
    print.message(g2m.genes)
  }

  s.score <- controlled_mean_score(sce, s.genes, knn)
  g2m.score <- controlled_mean_score(sce, g2m.genes, knn)
  cell.cycle.score <- controlled_mean_score(sce, c(s.genes, g2m.genes),
                                            knn)
  if (verbose) {
    print.message("# s.score > 0:", sum(s.score > 0), "fraction", sum(s.score > 0)/length(s.score))
    print.message("# g2m.score > 0:", sum(g2m.score > 0), "fraction", sum(g2m.score >
                                                                            0)/length(g2m.score))
  }


  colData(sce)$S_stage <- s.score
  colData(sce)$G2M_stage <- g2m.score
  colData(sce)$aggregate_S_G2M_stage <- cell.cycle.score
  return(sce)
}

#' Add cell specific score to metadata
#'
#'
#' @param sce A SingleCellExperiment object
#' @param name The name of the variable
#' @param values The value of the variables
#' @export
add_cell_score <- function(sce, name, values) {
  colData(sce)[[name]] <- values
  return(sce)
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
#' @export
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

  print.message("Memory intensive step. Please allocate adequate memory if the function aborts.")
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
  p <- dplyr::progress_estimated(n)
  for (index in seq(n)) {
    gene.dist <- get_dist(dist_obj, index, row_genes)
    knns[index, ] <- names(gene.dist[order(gene.dist)[2:(knn + 1)]])
    p$pause(0.1)$tick()$print()
  }

  ret <- list(knns = knns, technical.variables = technical.variables)
  saveRDS(ret, file = cache_path)
  print.message(paste0("Cached at ", cache_path, "\n"))
  return(ret)
}

get_dist <- function(dist_obj, i, names) {
  n <- nrow(dist_obj)
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

#' Regress covariate
#' @param rbs A RobustSingleCell object
#' @param regress Character of variables to regress
#' @param data Data matrix of normalized count
#' @param groups Group of cells
#' @param verbose Whether to print diagnostic messages
regress_covariates <- function(rbs, regress, data, groups, verbose = F) {
  cache <- file.path(rbs@cache, paste0("HVG.regressed.", paste(sort(colnames(regress)), collapse = "."), ".rds"))

  if (file.exists(cache)) {
    if (verbose) print.message("Loading precomputed")
    corrected <- readRDS(cache)
  } else {
    if (verbose) print.message("Computing")
    formula.str <- paste("gene", paste(colnames(regress), collapse = " + "),
                         sep = " ~ ")
    formula <- stats::as.formula(formula.str)
    if (verbose) {
      print.message("Regressing:", formula.str)
      print.message("Not regressed matrix")
      corner(data)
    }
    corrected <- data

    for (group in unique(groups)) {
      group.indices <- groups == group
      p <- dplyr::progress_estimated(nrow(data))
      for (gene in rownames(data)) {
        lm.data <- data.frame(gene = data[gene, group.indices], regress[group.indices,])
        colnames(lm.data)[2:ncol(lm.data)] <- colnames(regress)
        model <- stats::lm(formula, lm.data)
        corrected[gene, group.indices] <- model$residuals
        p$pause(0.1)$tick()$print()
      }
    }
    saveRDS(corrected, file = cache)
  }

  return(corrected)
}