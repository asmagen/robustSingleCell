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
  barcode.path <- list.files(path, pattern = "barcodes.tsv", full.names = T)
  features.path <- list.files(path, pattern = "genes.tsv", full.names = T)
  matrix.path <- list.files(path, pattern = "matrix.mtx", full.names = T)
  mat <- Matrix::readMM(file = matrix.path)
  feature.names <- utils::read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(feature.names) <- c("gene_id", "gene_name")
  rownames(feature.names) <- feature.names$gene_id
  feature.names <- feature.names[, -1, drop = F]
  gene_name <- feature.names$gene_name
  feature.names$is_ribo <- is_ribo(gene_name)
  feature.names$is_mito <- is_mito(gene_name)

  barcode.names <- utils::read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(barcode.names) <- "barcode"
  rownames(barcode.names) <- barcode.names$barcode
  barcode.names <- barcode.names[, -1, drop = F]

  set.seed(seed)
  if (!is.null(subsample) && subsample < ncol(mat)) {
    sampled <- sample(seq(ncol(mat)), subsample)
    mat <- mat[, sampled]
    barcode.names <- barcode.names[sampled, , drop = F]
  }

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
  if (is.null(barcode.names)) {
    barcode.names <- as.data.frame(matrix(, nrow = nrow(mat), ncol = 0))
  }
  barcode.names$nUMI <- colSums(mat)
  barcode.names$nGenes <- colSums(mat > 0)

  is_ribo <- feature.names$is_ribo
  barcode.names$mean_ribo <- colMeans(mat[is_ribo, ])
  barcode.names$frac_ribo <- colSums(mat[is_ribo, ])/barcode.names$nUMI

  is_mito <- feature.names$is_mito
  barcode.names$mean_mito <- colMeans(mat[is_mito, ])
  barcode.names$frac_mito <- colSums(mat[is_mito, ])/barcode.names$nUMI

  colData(sce) <- barcode.names
  return(sce)
}

plot_summary_statistics <- function(data) {
  data <- as.data.frame(data) %>%
    na.omit()
  for (ind1 in seq(length(colnames(data)) - 1)) {
    for (ind2 in (ind1 + 1):(length(colnames(data)))) {
      v1 <- colnames(data)[ind1]
      v2 <- colnames(data)[ind2]
      print(ggplot(data, aes_string(v1, v2)) + geom_jitter() + theme_classic())
    }
  }
}

#' Filter expression matrix
#'
#' @param environment SingleCellExperiment object
#' @param genome genome annotation
#' @param min.genes.per.cell minimum required number of genes per cell
#' @param max.genes.per.cell.quantile upper quantile for number of genes per cell
#' @param max.UMIs.per.cell.quantile upper quantile for number of UMIs per cell
#' @param min.cells.per.gene minimum required number of cells per gene
#' @param max.mitochondrial.frac maximum fraction of reads mapped to mitochondrial
#' genes per cell
#' @param max.ribosomal.frac maximum fraction of reads mapped to ribosomal genes per cell
#' @return A filtered SingleCellExperiment object
#' @export
filter <- function(sce, min.genes.per.cell = 500, max.genes.per.cell.quantile = 0.98,
                       max.UMIs.per.cell.quantile = 0.98, min.cells.per.gene = 1, max.mitochondrial.frac = 0.1,
                       max.ribosomal.frac = NULL, verbose = F) {
  sce <- sce[, colData(sce)$nGenes >= min.genes.per.cell]
  sce <- update_summary_stats(sce)
  counts <- assay(sce, "counts")

  genes.per.cell <- colData(sce)$nGenes
  UMIs.per.cell <- colData(sce)$nUMI
  max.genes.per.cell <- stats::quantile(genes.per.cell, max.genes.per.cell.quantile)
  max.UMIs.per.cell <- stats::quantile(UMIs.per.cell, max.UMIs.per.cell.quantile)
  frac_mito <- colData(sce)$frac_mito
  frac_ribo <- colData(sce)$frac_ribo

  if (verbose) {
    print.message("Original genes.per.cell")
    print(summary(genes.per.cell))
    print(stats::quantile(genes.per.cell, seq(0.01, 1, 0.01)))
    print.message("Original UMIs.per.cell")
    print(summary(UMIs.per.cell))
    print(stats::quantile(UMIs.per.cell, seq(0.01, 1, 0.01)))
    print.message("max.genes.per.cell", max.genes.per.cell)
    print.message("max.UMIs.per.cell", max.UMIs.per.cell)
    print.message("max.mitochondrial.frac", max.mitochondrial.frac)

    print.message("# not qualify min.genes.per.cell", sum(genes.per.cell <
                                                            min.genes.per.cell))
    print.message("# not qualify max.genes.per.cell", sum(genes.per.cell >
                                                            max.genes.per.cell))
    print.message("# not qualify max.UMIs.per.cell", sum(UMIs.per.cell >
                                                           max.UMIs.per.cell))
    print.message("# not qualify max.mitochondrial.frac", sum(frac_mito >
                                                                max.mitochondrial.frac))
  }


  criteria <- genes.per.cell >= min.genes.per.cell & genes.per.cell <=
  max.genes.per.cell & UMIs.per.cell <= max.UMIs.per.cell & frac_mito <=
  max.mitochondrial.frac

  if (!is.null(max.ribosomal.frac))
    criteria <- criteria & frac_ribo <= max.ribosomal.frac

  if (verbose) print.message("# qualifying cells", sum(criteria),
                             "\n# not qualifying cells", sum(!criteria))

  sce <- sce[, criteria]
  sce <- update_summary_stats(sce)

  return(sce)
}

#' Log normalize the expression matrix
#'
#' @param sce A SingleCellExperiment object
#' @return Return a SingleCellExperiment object
log_normalize <- function(sce, verbose = F) {
  counts <- assay(sce, "counts")
  sizeFactors(sce) <- colData(sce)$nUMI
  normalized <- counts / sizeFactors(sce)
  normalized <- normalized * 10000
  normalized <- log(normalized + 1)
  assay(sce, "normcounts") <- normalized
  return(sce)
}

