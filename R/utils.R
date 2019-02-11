capwords <- function(s, strict = FALSE) {
  s = tolower(s)
  cap <- function(s) paste(toupper(substring(s, 1, 1)), {s <- substring(s, 2); if(strict) tolower(s) else s}, sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

corner <- function (matrix,n=5,m=5) {
  print(matrix[seq(min(n,nrow(matrix))),seq(min(m,nrow(matrix)))])
}

get_slurm_out <- function (sjob) {

  Sys.sleep(1)
  queued = length(system(paste('squeue -hn', sjob$jobname),intern = T)) > 0
  while(length(system(paste('squeue -hn', sjob$jobname),intern = T)) > 0) {
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
  if (length(res_files) == 0) return(NA)

  slurm_out <- lapply(res_files, readRDS)
  slurm_out <- do.call(c, slurm_out)
  slurm_out <- as.data.frame(do.call(rbind, slurm_out))

  return(slurm_out)
}

start <- function (track_dir_path, name = NA,append = F,split = F,print = T) {
  time = get.time()
  label = as.character(sys.calls()[[sys.nframe()-1]])[1]
  if (!is.na(name)) label = paste(label,name,sep='.')
  file = paste(label,'txt',sep='.')
  sink(file.path(track_dir_path, file),type = 'output',append,split )
  if (!split && print) print(as.list(sys.calls())[seq(sys.nframe()-1)])
  return(time)
}

end <- function(time = NA) {
  if (!is.na(time)) elapsed.time(time)
  sink()
}

print.message <- function (...) {
  cat(cat(...,sep=' '),'\n',sep='')
}

get.time <- function () {
  return(Sys.time())
}

elapsed.time <- function (time) {
  print(Sys.time()-time)
  cat('\n')
}

apply.by.group <- function( groups,values,... ) {
  group = group.by( groups,values )
  for( fun in list(...) ) {
    group = sapply(group,fun)
  }
  return(group)
}

group.by <- function( groups,values ) {
  return( split(x = as.vector(values), f = as.vector(groups)) )
}
