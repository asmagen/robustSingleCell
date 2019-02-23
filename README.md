robustSingleCell
================

### Overview

robustSingleCell is a pipeline designed to identify robust cell subpopulations using scRNAseq data and compare population compositions across tissues and experimental models via similarity analysis as described in Magen et al. (2019) \[1\].

### Installation and tutorial

To install the latest version from Github repository, use

``` r
if(!require(devtools))
  install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("grimbough/biomaRt")

devtools::install_github("asmagen/robustSingleCell")
```

This pipeline currently supports [slurm](https://slurm.schedmd.com)
for parallel batch jobs. 

### Getting help

If you encounter a clear bug, please submit an
[issue](https://github.com/asmagen/robustSingleCell/issues) with
reproducible example.

### Reference

1.  Magen *et al*. “Single-cell profiling of tumor-reactive CD4<sup>+</sup> T-cells reveals unexpected transcriptomic diversity” [*bioRxiv 543199*](https://doi.org/10.1101/543199)
