MagenSingleCell
================

### Overview

MagenSingleCell is a robust clustering and similarity measurement method
designed for single-cell data analysis \[1\].

### Installation and tutorial

To install the latest version from Github repository, use

``` r
if(!require(devtools))
  install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

devtools::install("JinmiaoChenLab/Rphenograph")
BiocManager::install("limma")

devtools::install_github("asmagen/magensinglecell")
```

This pipeline currently only supports [slurm](https://slurm.schedmd.com)
for parallel batch jobs. See a tutorial at
<https://asmagen.github.io/MagenSingleCell>.

### Getting help

If you encounter a clear bug, please submit an
[issue](https://github.com/asmagen/MagenSingleCell/issues) with
reproducible example.

### Reference

1.  Magen *et al*. “Single-cell profiling of tumor-reactive CD4<sup>+</sup> T-cells reveals unexpected transcriptomic diversity” [*bioRxiv 543199*](https://doi.org/10.1101/543199)
