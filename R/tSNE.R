#' @import rslurm
run.tSNE <- function (perplexity,max_iter,rerun) {

  tSNE <- function (perplexity,max_iter) {

    load(data.path)

    duplicated.indices = duplicated(t(PCA))
    tSNE = Rtsne::Rtsne(t(PCA[,!duplicated.indices]), pca = F, initial_dims = nrow(PCA), perplexity = perplexity, max_iter = max_iter, verbose=T, whiten=F)$Y

    save(tSNE,file=file.path(tSNEs.dir,paste(perplexity,max_iter,'tSNE.RData',sep='.')))
  }

  tSNEs.dir = file.path(environment$res.data.path,'tSNEs')
  list.files(tSNEs.dir)
  sjob = NA

  if (rerun || !dir.exists(tSNEs.dir)) {
    t = start(split = T)
    print(tSNEs.dir)
    # unlink(tSNEs.dir,recursive=T,force=T)
    dir.create(tSNEs.dir)

    duplicated.indices = duplicated(t(environment$PCA))
    if (sum(duplicated.indices) > 0) {
      print.message('Excluding',sum(duplicated.indices),'duplicated cells from tSNE analysis')
      print(table(environment$dataset.labels[duplicated.indices]))
    }

    data.path = environment$PCA.path
    params = data.frame(expand.grid(perplexity = perplexity,max_iter = max_iter))
    if (nrow(params) > 10) {
      print.message('Warning: possibly too many tSNE parameter combinations [',nrow(params),'combinations ]')
      terminate = readline(prompt="Terminate? (y/n) ")
      if (terminate != 'n') {
        return()
      }
    }

    sopt <- list(mem = '8GB', time = '1:00:00', share = TRUE)
    sjob <- slurm_apply(tSNE, params, nodes = nrow(params), cpus_per_node = 1, add_objects = c('data.path','tSNEs.dir'), submit = TRUE, slurm_options = sopt)
    end()
  }

  return(sjob)
}