setClass(
  "robustSingleCell",
  contains="SingleCellExperiment",
  slots=c(origin="character",
          experiment="character",
          cache="character")
) -> robustSingleCell
