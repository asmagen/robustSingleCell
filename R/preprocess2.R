#' Identify Highly Variable Genes
#'
#' Get highly variable genes by Heteroscedasticity controlled binning of gene expression measurements within each dataset separately.
#'
#' @param environment A SingleCellExperiment object
#' @param min.mean minimum mean expression per gene
#' @param min.frac.cells minimum fraction of cells expressing each gene
#' @param min.dispersion.scaled minimum dispersion value
#' @param plot whether to plot mean vs dispersion plot
#' @param verbose whether to print diagnostic message
#' @return \code{environment} parameter containing highly variable genes selection
#' @export
#' @importFrom Matrix rowSums
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