cluster <- function (knn.ratio,label,path,data.path,nPCs) {

  cat(paste('Params','\nknn.ratio = ',knn.ratio,'\nlabel = ',label,'\npath = ',path,'\ndata.path = ',data.path,'\n',sep=''))

  if (length(grep('shuffled',label))>0) {
    # PCA = apply(PCA,2,sample,nrow(PCA))
    rep = strsplit(as.character(label),split='[.]')[[1]][2]
    file = file.path(path,'shuffled.PCA',paste('shuffled.PCA.rep',rep,'RData',sep='.'))
    if (!file.exists(file)) print.message('Shuffled PCA file',rep,'does not exist [Change PCA or clustering parameters]')
    load(file)
    PCA = t(pca.perm$x)[seq(as.numeric(nPCs)),]
  } else {
    load(file.path(data.path))
  }

  library(Rphenograph,quietly=T);library(cccd,quietly=T)
  # Rphenograph_out = Rphenograph(t(PCA), k = floor(knn.ratio*ncol(PCA)))
  neighborMatrix <- find_neighbors(t(PCA), k=floor(knn.ratio*ncol(PCA)))[,-1]
  jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
  }
  links <- jaccard_coeff(neighborMatrix)
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  community <- cluster_louvain(graph.data.frame(relations, directed=FALSE))

  modularity = modularity(community)
  membership = membership(community)
  memberships = community$memberships
  nclusters = length(unique(membership))

  clustering = list(modularity=modularity,memberships=memberships,membership=membership,knn.ratio=knn.ratio,nclusters=nclusters)

  save(clustering,file=file.path(path,'clustering',paste(label,knn.ratio,'RData',sep='.')))
}

get.clustering.results <- function (clustering.dir,knn.ratio,shuffledKNN) {

  shuffled.result.files = list.files(clustering.dir,pattern = paste('shuffled.*.*.RData',sep=''));shuffled.result.files
  nshuffled = length(shuffled.result.files)
  shuffled = {}
  rep = 1
  shuffled.membership = list()
  for (rep in seq(nshuffled)) {
    load(file.path(clustering.dir,shuffled.result.files[rep]))
    shuffled = rbind(shuffled,unlist(clustering[c(1,4,5)]))
    shuffled.membership[[rep]] = as.vector(clustering[2])
  }
  shuffled = as.data.frame(shuffled)

  real.clustering.file = file.path(clustering.dir,paste('real',knn.ratio,'RData',sep='.'))
  if(!file.exists(real.clustering.file)) cat('\nERROR: Couldn\'t find real clustering file knn.ratio =',knn.ratio,'\nCHECK FOR FAILED JOBS\n\n\n')
  load(real.clustering.file);#str(clustering)
  membership = clustering$membership
  nFinalClusters = length(unique(membership))
  memberships = clustering$memberships
  modularity = clustering$modularity
  nclusters = apply(memberships,1,function(v) length(unique(v)))#length(unique(clustering$membership))
  shuffled.indices = order(abs(shuffled$nclusters-nFinalClusters))[seq(shuffledKNN)]
  shuffled.matches = shuffled[shuffled.indices,]
  shuffled.membership = shuffled.membership[shuffled.indices]

  cat('knn.ratio = ',knn.ratio,' nclusters = ',nFinalClusters,' (nClustShuffled = ',ifelse(var(shuffled.matches$nclusters)>0,paste(min(shuffled.matches$nclusters),'~',max(shuffled.matches$nclusters),sep=''),shuffled.matches$nclusters[1]),')\n',sep='')

  mean.shuffled = mean(shuffled.matches$modularity)
  max.shuffled = max(shuffled.matches$modularity)
  cat(
    'mdlrty=',round(modularity,2),
    ' mean.shfl=',round(mean.shuffled,2),
    ' max.shfl=',round(max.shuffled,2),
    ' mdlrty/mean.shfl=',round(modularity/mean.shuffled,2),
    ' mdlrty/max.shfl=',round(modularity/max.shuffled,2),
    '\n',
    sep='')

  return(list(membership=membership,memberships=memberships,modularity=modularity,shuffled=shuffled.matches,shuffled.membership=shuffled.membership,shuffled.modularity=shuffled.matches$modularity))
}

cluster.analysis <- function (
  knn.ratios = c(0.01,0.05,0.1), # Range of KNN parameters to scan (corresponding to different resolutions)
  nShuffleRuns = 10, # Number of shuffled clustering analyses to perform per KNN threshold
  shuffledKNN = 10, # Number of closest KNN shuffled analyses to include in background clustering quality computation
  loadPreviousKnn = T, # Whether to load previous analysis results
  rerun = F, # Whether to rerun
  deleteCache = F, # Whether to delete cahce files
  mem = '4GB', # HPC memory
  time = '0:15:00', # HPC time
  plot = T) { # Whether to plot the clustering qualities compared to shuffled

  clustering.dir = file.path(environment$res.data.path,'clustering')
  shuffled.result.files = list.files(clustering.dir,pattern = paste('shuffled.*.*.RData',sep=''));shuffled.result.files

  extended.knn.ratios = unique(c(knn.ratios,seq(max(knn.ratios),max(knn.ratios)*1.25,max(knn.ratios)-knn.ratios[order(knn.ratios)][length(knn.ratios)-1])))
  path = environment$res.data.path
  data.path = environment$PCA.path
  params = rbind(
    data.frame(knn.ratio=knn.ratios,label='real',path,data.path,nPCs = nrow(environment$PCA)),
    data.frame(expand.grid(knn.ratio=extended.knn.ratios,label=paste('shuffled',seq(nShuffleRuns),sep='.')),path,data.path,nPCs = nrow(environment$PCA))
  )
  cache = file.path(environment$res.data.path,'clustering.RData')
  sjob = NA

  if( rerun || !dir.exists(clustering.dir) ) {
    print.message('Computing')
    t = start()

    if (deleteCache) unlink(clustering.dir,recursive=T,force=T)
    dir.create(clustering.dir)

    print.message('knn.ratios:');print(knn.ratios)
    print.message('# params:',nrow(params))
    print.message('Head params:');print(head(params))
    print.message('Tail params:');print(tail(params))

    suppressWarnings(library(rslurm,quietly=T))
    sopt <- list(mem = mem, time = time, share = TRUE)
    sjob <- slurm_apply(cluster, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)#, add_objects = c('path','data.path')
    tryCatch({get_slurm_out(sjob);get_slurm_out(sjob);get_slurm_out(sjob)},error=function(v) v)

    end(t)
  }

  if( loadPreviousKnn && file.exists(cache) ) {
    print.message('Loading precomputed')
    load(cache)
  } else {
    t = start(name = 'KNN.stats',split = T)

    nResults = length(list.files(clustering.dir,pattern = paste('*.RData',sep='')))
    if (nResults < nrow(params)) {
      cat('\nERROR: Found just',nResults,'shuffled clusterings instead of',nrow(params),'\nCHECK FOR FAILED JOBS\n\n\n')
      terminate = readline(prompt="Terminate? (y/n) ")
      if (terminate != 'n') {
        tryCatch({cleanup_files(sjob);cleanup_files(sjob);cleanup_files(sjob)},error=function(v) v)
        stop()
      }
    }
    tryCatch({cleanup_files(sjob);cleanup_files(sjob);cleanup_files(sjob)},error=function(v) v)

    plot.stats = list()
    clusters.aggregate = {}

    for (knn.ratio in knn.ratios[length(knn.ratios):1]) {
      tryCatch({
        result = get.clustering.results (clustering.dir,knn.ratio,shuffledKNN)
        nClusters = paste(length(unique(result$membership)))
        stats = plot.stats[[nClusters]]
        if (is.null(stats)) stats = {}
        stats = rbind(stats,data.frame(knn.ratio = knn.ratio, Type = 'Original data', Modularity = result$modularity))
        stats = rbind(stats,data.frame(knn.ratio = knn.ratio, Type = 'Shuffled data', Modularity = result$shuffled.modularity))
        plot.stats[[nClusters]] = stats
        clusters.aggregate = cbind(clusters.aggregate,result$membership)
      }, error = function(v) v )
    }

    if (plot) {
      ncluster = names(plot.stats)[1]
      new.plot.stats = {}
      fold.data = {}
      for (ncluster in as.character(sort(as.numeric(names(plot.stats))))) {
        summary(plot.stats[[ncluster]])
        data = plot.stats[[ncluster]]
        original = data$Modularity[data$Type == 'Original data']
        shuffled = data$Modularity[data$Type == 'Shuffled data']
        new.plot.stats = rbind(new.plot.stats,data.frame(nClusters=ncluster, Type = 'Original Data', Modularity = mean(original),sd = sd(original)))
        new.plot.stats = rbind(new.plot.stats,data.frame(nClusters=ncluster, Type = 'Shuffled Data', Modularity = mean(shuffled),sd = sd(shuffled)))
        fold.data = rbind(fold.data,data.frame(nClusters=ncluster,Type='Original Data', Modularity = mean(original), Fold = mean(original)/mean(shuffled)))
      }
      library(ggplot2)
      pdf(file.path(environment$work.path,'Clustering.modularity.pdf'),width=10)
      print(ggplot(new.plot.stats, aes(x=nClusters, y=Modularity, fill=Type)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=Modularity-sd, ymax=Modularity+sd), width=.2, position=position_dodge(.9)) + labs(title="Clustering modularity analysis", x="Number of Clusters", y = "Modularity") + theme_classic(base_size=25) + scale_fill_manual(values=c('#999999','#E69F00')) + geom_text(data = fold.data,aes(nClusters, max(Modularity)*1.05, label = sprintf("%2.1f", Fold)),size=7))
      print(ggplot(new.plot.stats, aes(x=nClusters, y=Modularity, fill=Type)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=Modularity-sd, ymax=Modularity+sd), width=.2, position=position_dodge(.9)) + labs(title="Clustering modularity analysis", x="Number of Clusters", y = "Modularity") + theme_classic(base_size=25) + scale_fill_manual(values=c('#999999','#E69F00')) + geom_text(data = fold.data,aes(nClusters, max(Modularity)*1.05, label = sprintf("%2.2f", Fold)),size=7))
      dev.off()
    }

    knn.choose = as.numeric(readline(prompt="Select KNN ratio: "))
    print.message('knn.choose =',knn.choose)
    clustering.results = get.clustering.results (clustering.dir,knn.choose,nShuffleRuns)
    nclusters = as.numeric(readline(prompt="Select # clusters: "))
    print.message('nclusters =',nclusters)

    clustering = {}
    clustering$knn.choose = knn.choose
    clustering$membership=as.vector(clustering.results$memberships[which(apply(clustering.results$memberships,1,function(v) length(unique(v)))==nclusters),])
    clustering$nclusters = length(unique(clustering$membership))
    print.message('# clusters:',clustering$nclusters)
    print.message('Length:',length(clustering$membership))
    clustering$shuffled=clustering.results$shuffled
    clustering$shuffled.membership=clustering.results$shuffled.membership

    files = c(
      list.files(environment$work.path,pattern = '*.pdf',full.names=T),list.files(environment$work.path,pattern = '*.csv',full.names=T),
      list.files(file.path(environment$work.path,'diff.exp','main'),pattern = '*.csv',full.names=T),
      list.files(file.path(environment$work.path,'Seurat'),pattern = '*.pdf',full.names=T),
      list.files(file.path(environment$work.path,'tracking'),pattern = '*.txt',full.names=T))
    dir = file.path(environment$work.path,format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"));dir.create(dir)
    if (length(files) > 0) {
      file.copy(files, dir);file.remove(files)
      file.remove(c(list.files(environment$res.data.path,pattern = '*.diff.exp.RData',full.names=T),file.path(environment$res.data.path,c('main.diff.exp.RData','clustering.RData','homogeneity.RData','Seurat.RData'))))
    }
    unlink(file.path(environment$work.path,'diff.exp'),recursive=T,force=T)
    save(clustering,file=cache)
    end(t)
  }

  environment$clustering = clustering
  environment$cluster.names = environment$clustering$membership#apply(cbind(environment$datasets,'clust',environment$clustering$membership),1,function(v) paste(v,collapse=' '))

  print.message('Membership table:');print(table(clustering$membership))

  return(environment)
}

pearson.correlation <- function (diff1,diff2) {

  cor = cor.test(diff1,diff2,method='pearson')
  return (data.frame(similarity = cor$estimate,significance = cor$p.value))
}