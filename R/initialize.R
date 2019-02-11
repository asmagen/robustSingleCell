#' Intialize the project object
#'
#' TODO: description
#'
#' @param datasets List of dataset code names
#' @param origins List of dataset tissue origin/condition full name
#' @param experiments List of experiment design annotations
#' @param data.path Path to where the data is located
#' @param work.path path to where the analysis should create folders and store files
#' @param marker.genes Set of genes of interest for visualization purposes
#' @param clear.history Whether you would like to remove any previous project by this name
#' @param analysis.label Whether you would like to add a specific label to the analysis folder
#' @param convert.to.mouse.gene.symbols Whether you are using human gene symbols and would like to convert them
#' to mouse gene symbols
#' @export
initialize.project <- function (
  datasets,
  origins,
  experiments,
  data.path,
  work.path,
  marker.genes,
  clear.history = F,
  analysis.label = NA,
  convert.to.mouse.gene.symbols = NULL) {

  options("width"=205)
  environment = {}
  environment$datasets = datasets
  environment$origins = origins
  environment$experiments = experiments
  environment$labels = labels
  environment$data.path = data.path
  environment$project = ifelse(length(datasets)>1,paste('merged',paste(datasets,collapse='.'),sep='.'),datasets)
  if (!is.na(analysis.label)) environment$project = paste(analysis.label,environment$project,sep='_')
  if (length(convert.to.mouse.gene.symbols) > 1 || !is.null(convert.to.mouse.gene.symbols)) environment$convert.to.mouse.gene.symbols = convert.to.mouse.gene.symbols
  environment$baseline.work.path = environment$work.path = file.path(work.path,environment$project)
  if (clear.history) unlink(environment$work.path,recursive = T, force = T)
  dir.create(environment$work.path, showWarnings = F, recursive = T, mode = "700")
  dir.create(file.path(environment$work.path, 'tracking'),
             showWarnings = F, recursive = T, mode = "700")
  environment$baseline.data.path = environment$res.data.path = file.path(environment$work.path,'data')
  dir.create(environment$res.data.path, showWarnings = F, recursive = T, mode = "700")
  environment$marker.genes = marker.genes
  t = start(file.path(environment$work.path, 'tracking'))
  print(environment)
  end(t)

  return(environment)
}