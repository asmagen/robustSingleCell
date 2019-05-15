#' Read 10X Data
#'
#' Load sparse data matrices from 10X genomics into a SingleCellExperiment object.
#'
#' @param path Path to directory containing matrix.mtx, genes.tsv, and barcodes.tsv
#' @param subsample number of cells to subsample
#' @param seed seed for subsampling of cells
#' @return a SingleCellExperiment object
#' @export
read_10x_data <- function(path, subsample = NULL, seed = 1) {
  browser()
  barcode.path <- list.files(path, pattern = "barcodes.tsv", full.names = T)
  features.path <- list.files(path, pattern = "genes.tsv", full.names = T)
  matrix.path <- list.files(path, pattern = "matrix.mtx", full.names = T)
  mat <- Matrix::readMM(file = matrix.path)
  feature.names <- utils::read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(feature.names) <- c("gene_id", "gene_name")
  rownames(feature.names) <- feature.names$gene_id
  feature.names <- feature.names[, -1, drop = F]

  barcode.names <- utils::read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(barcode.names) <- "barcode"

  set.seed(seed)
  if (!is.null(subsample) && subsample < ncol(mat)) {
    sampled <- sample(seq(ncol(mat)), subsample)
    mat <- mat[, sampled]
    barcode.names <- barcode.names[sampled, , drop = F]
  }

  rownames(barcode.names) <- barcode.names$barcode

  sce <- SingleCellExperiment(
    assays = list(counts = mat),
    colData = barcode.names,
    rowData = feature.names
  )

  sce <- update_summary_stats(sce)

  return(sce)
}

#' Compute summary statistics on count matrix
#'
#' @param sce SingleCellExperiment object
#' @return A data frame containing the number of UMI, genes with expression,
#' mean ribosomal gene counts, mean mitochondrial gene counts, fraction of
#' ribosomal gene expression per cell, fraction of mitochrondrial gene expression per cell
#' @importFrom Matrix colSums colMeans
update_summary_stats <- function(sce) {
  mat <- assay(sce, "counts")
  feature.names <- rowData(sce)
  barcode.names <- colData(sce)
  barcode.names$nUMI <- colSums(mat)
  barcode.names$nGenes <- colSums(mat > 0)

  is_ribo <- is_ribo(feature.names$gene_name)
  barcode.names$mean_ribo <- colMeans(mat[is_ribo, ])
  barcode.names$frac_ribo <- colSums(mat[is_ribo, ])/barcode.names$nUMI

  is_mito <- is_mito(feature.names$gene_name)
  barcode.names$mean_mito <- colMeans(mat[is_mito, ])
  barcode.names$frac_mito <- colSums(mat[is_mito, ])/barcode.names$nUMI
  barcode.names <- barcode.names[, colnames(barcode.names) != "barcode"]
  colData(sce) <- barcode.names

  return(sce)
}

plot_summary_statistics <- function(data) {
  for (ind1 in seq(length(colnames(data)) - 1)) {
    for (ind2 in (ind1 + 1):(length(colnames(data)))) {
      v1 <- colnames(data)[ind1]
      v2 <- colnames(data)[ind2]
      print(ggplot(data, aes_string(v1, v2)) + geom_point())
    }
  }
}

filter_sce <- function(sce, min.genes.per.cell = 500, max.genes.per.cell.quantile = 0.98,
                       max.UMIs.per.cell.quantile = 0.98, min.cells.per.gene = 1, max.mitochondrial.frac = 0.1,
                       max.ribosomal.frac = NULL, cell.filters = NULL) {

  print.message("Original dimensions")
  print(dim(sce))

  genes.per.cell <- colSums(counts > 0)
  print.message("Original genes.per.cell")
  print(summary(genes.per.cell))
  print(stats::quantile(genes.per.cell, seq(0.01, 1, 0.01)))
  UMIs.per.cell <- colSums(counts)
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

  final.gene.filter.passed.cells <- rep(T, ncol(counts))
  if (!is.null(cell.filters)) {
    filter.genes <- as.vector(cell.filters$gene)
    if (sum(filter.genes %in% rownames(counts)) == 0)
      print.message("All filter genes do not exist in dataset - not filtering anything")
    not.exist.genes <- sum(!filter.genes %in% rownames(counts))
    if (not.exist.genes > 0)
      print.message("Not filtering based on", not.exist.genes, "filter genes do not exist in dataset")
    filter.genes <- filter.genes[filter.genes %in% rownames(counts)]
    if (length(filter.genes) == 1)
      filter.genes <- c(filter.genes, filter.genes)
    gene <- filter.genes[1]
    matched.conditions <- rep(T, ncol(counts))
  }
  sum.not.matching <- sum(matched.conditions)
  print.message("# not qualify dataset.filters", sum.not.matching,
                "out of", ncol(counts), "[", round(sum.not.matching/ncol(counts) * 100, 1), "% ]")
  final.gene.filter.passed.cells <- matched.conditions

  criteria <- genes.per.cell >= min.genes.per.cell & genes.per.cell <=
  max.genes.per.cell & UMIs.per.cell <= max.UMIs.per.cell & data$Mitochondrial.frac <=
  max.mitochondrial.frac & final.gene.filter.passed.cells
  if (!is.null(max.ribosomal.frac))
    criteria <- criteria & data$Ribosomal.frac <= max.ribosomal.frac
  print.message("# qualifying cells", sum(criteria), "# not qualifying cells",
              sum(!criteria))
  sce <- sce[, criteria]
  genes.per.cell <- colSums(counts > 0)
  print.message("Filtered genes.per.cell")
  print(summary(genes.per.cell))
  print(stats::quantile(genes.per.cell, seq(0.01, 1, 0.01)))
  UMIs.per.cell <- colSums(counts)
  print.message("Filtered UMIs.per.cell")
  print(summary(UMIs.per.cell))
  print(stats::quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))

  data <- compute_summary_stats(counts)
  print(utils::head(data))

  plot_summary_statistics(data)


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

  return(sce)
}


