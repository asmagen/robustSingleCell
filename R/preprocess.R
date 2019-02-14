#' Get highly variable genes by Heteroscedasticity controlled binning
#'
#'
#' @param min.mean Minimum mean expression per gene
#' @param min.frac.cells Minimum fraction of cells expressing each gene
#' @param min.dispersion.scaled Minimum dispersion value
#' @param rerun Whether to rerun
#' @return \code{environment} parameter containing highly variable genes selection
#' @export
get.variable.genes <- function (
  environment,
  min.mean = 0.05,
  min.frac.cells = 0,
  min.dispersion.scaled = 1,
  rerun = F) {

  cache = file.path(environment$baseline.data.path,'HVG.RData')

  if( !rerun & file.exists(cache) ) {
    print.message('Loading precomputed')
    load(cache)
  } else {
    print.message('Computing')
    t = start(file.path(environment$work.path, 'tracking'))

    normalized = environment$normalized

    datasets = environment$datasets;table(datasets)
    dataset = sort(unique(datasets))[1];dataset
    merged = {}
    for(dataset in sort(unique(datasets))) {
      print.message(dataset)
      filtered.normalized = normalized[,datasets==dataset]
      means = apply(filtered.normalized,1,function(v) log(mean(exp(v) - 1) + 1))
      frac.cells = rowSums(filtered.normalized>0)/ncol(filtered.normalized)
      vars = apply(filtered.normalized,1,function(v) log(var(exp(v) - 1) + 1))
      dispersion = apply(filtered.normalized,1,function(v) log(var(exp(v) - 1) / mean(exp(v) - 1)))
      dispersion[is.na(x = dispersion)] = 0
      means[is.na(x = means)] = 0

      pdf(file.path(environment$baseline.work.path,paste(dataset,'VariableGenes.pdf',sep='.')))
      plot.data = data.frame(gene = names(means),means,dispersion)#head(plot.data)
      print(ggplot(plot.data, aes(x=means, y=dispersion, label = gene)) + geom_text(check_overlap = TRUE,size=2))
      smoothScatter(means,dispersion)
      dev.off()

      num.bin = 20
      bins = cut(x = means, breaks = num.bin)
      names(x = bins) = names(x = means)
      mean_y = tapply(dispersion,bins,mean)
      sd_y = tapply(dispersion,bins,sd)
      dispersion.scaled = (dispersion - mean_y[as.numeric(x = bins)]) / sd_y[as.numeric(x = bins)]
      dispersion.scaled[is.na(x = dispersion.scaled)] = 0
      names(x = dispersion.scaled) = names(x = means)

      criterion = means >= min.mean & frac.cells >= min.frac.cells & dispersion.scaled >= min.dispersion.scaled
      HVG = names(means)[criterion]
      print.message('# qualifying genes',length(HVG))
      print.message('Qualifying markers');print(environment$marker.genes[environment$marker.genes%in%HVG])
      print.message('NOT qualifying markers');print(environment$marker.genes[!environment$marker.genes%in%HVG])
      genes = c(get.ribo.genes(rownames(normalized)),get.mito.genes(rownames(normalized)))
      print.message('Qualifying Ribosomal & Mitochondrial');print(genes[genes%in%HVG])
      write.csv(HVG,file=file.path(environment$baseline.work.path,paste(dataset,'VariableGenes.csv',sep='.')))
      merged = unique(c(merged,HVG))

    }
    HVG = merged

    print.message('Overall # qualifying genes',length(HVG))
    print.message('Overall qualifying markers')
    print(environment$marker.genes[environment$marker.genes%in%HVG])
    print.message('Overall NOT qualifying markers')
    print(environment$marker.genes[!environment$marker.genes%in%HVG])
    genes = c(get.ribo.genes(rownames(normalized)),get.mito.genes(rownames(normalized)))
    print.message('Overall qualifying Ribosomal & Mitochondrial');print(genes[genes%in%HVG])

    save(HVG,file=cache)

    end(t)
  }

  environment$HVG = HVG

  cat('# highly variable genes = ', length(environment$HVG), '\n', sep='')

  return(environment)
}

nUMIs <- function (environment) { return(colSums(environment$counts)) }

nGenes <- function () { return(colSums(environment$counts>0)) }

#' Add confounder variable to environment object
#' 
#'
#' @param environment The environment object
#' @param ... confounding variables
#' @return \code{environment} parameter containing added confounder variable
#' @export
add.confounder.variables <- function (environment, ...) {
  environment$confounders = data.frame(environment$confounders, data.frame(...))
  print(head(environment$confounders))
  return(environment)
}

#' Compute ribosomal score
#' 
#'
#' @param environment The environment object
#' @param control Whether to subtract the score defined by technically similar genes
#' @param knn Number of nearest neighbor
#' @return A vector of ribosomal genes activation score 
#' @export
ribosomal.score <- function (environment, control = T,knn=10) {
  t = start(file.path(environment$work.path, 'tracking'))
  genes = get.ribo.genes(environment$genes)
  print.message('Using genes:');print(genes)
  if (control) {
    score = controlled.mean.score(environment, genes,knn)
  } else {
    score = colMeans(environment$normalized[genes,])
  }
  end(t)
  return(score)
}

get.ribo.genes <- function (genes) {
  return (genes[c(grep('^Rpl',genes,ignore.case = T),grep('^Rps',genes,ignore.case = T))])
}

#' Compute mitochondrial score
#' 
#'
#' @param environment The environment object
#' @param control Whether to subtract the score defined by technically similar genes
#' @param knn Number of nearest neighbor
#' @return A vector of mitochondrial genes activation score 
#' @export
mitochondrial.score <- function (environment, control = F, knn=10) {
  #browser()
  t = start(file.path(environment$work.path, 'tracking'))
  genes = get.mito.genes(environment$genes)
  print.message('Using genes:');print(genes)
  if (control) {
    score = controlled.mean.score(environment, genes, knn)
  } else {
    score = Matrix::colMeans(environment$normalized[genes %in% rownames(environment$normalized),])
  }
  end(t)
  return(score)
}

get.mito.genes <- function (genes) {
  return (genes[grep('^Mt-',genes,ignore.case = T)])
}

#' Compute cell cycle score
#' 
#' 
#' @param environment The environment object
#' @param knn The number of nearest neighbor
#' @param cc.genes.path Optional; path to defined cell cycle genes. Default uses gene sets defined in Kowalczyk, M. S. et al. Single-cell RNA-seq reveals changes in cell cycle and differentiation programs upon aging of hematopoietic stem cells. Genome Res 25, 1860-1872, doi:10.1101/gr.192237.115 (2015).
#' @return A matrix of cell cycle genes activation scores (S, G2M and aggregated S/G2M scores, separately)
#' @export 
cell.cycle.score <- function (environment, knn = 10, cc.genes.path = NA) {
  t = start(file.path(environment$work.path, 'tracking'))

  if(is.na(cc.genes.path)) {
    cc.genes <- capwords(cell_cycle_genes)
  } else {
    cc.genes <- capwords(readLines(cc.genes.path))
  }

  s.genes = cc.genes[1:43];s.genesi = s.genes[s.genes %in% capwords(environment$genes)];print(s.genes)
  g2m.genes = cc.genes[44:98];g2m.genes = g2m.genes[g2m.genes %in% capwords(environment$genes)];print(g2m.genes)

  s.score = controlled.mean.score (environment, s.genes,knn)
  g2m.score = controlled.mean.score (environment, g2m.genes,knn)
  cell.cycle.score = controlled.mean.score (environment, c(s.genes,g2m.genes), knn)

  print.message('# s.score > 0:',sum(s.score>0),'fraction',sum(s.score>0)/length(s.score))
  print.message('# g2m.score > 0:',sum(g2m.score>0),'fraction',sum(g2m.score>0)/length(g2m.score))
  end(t)

  return(data.frame(s.score=s.score,g2m.score=g2m.score,cell.cycle.score=cell.cycle.score))
}

#' Compute controlled mean gene signatures activation score
#' 
#' 
#' @param environment The environment object
#' @param genes The gene list upon which to calculate gene signature activate
#' @param knn Number of nearest neighbors
#' @param exclude.missing.genes Whether to exclude genes with missing values
#' @param constrain.cell.universe A binary vector indicating in which subset of cells to calculate the gene signature activation. Default is all cells.
#' @return Gene signature activation scores per cell
#' @export
controlled.mean.score <- function (environment, genes,knn = 10,exclude.missing.genes = T,constrain.cell.universe = NA) {
  # similarly to http://science.sciencemag.org/content/sci/suppl/2016/04/07/352.6282.189.DC1/Tirosh.SM.pdf/Seurat to reduce association with library size or other technical

  if (is.na(constrain.cell.universe)) constrain.cell.universe = rep(T,ncol(environment$normalized))
  if (knn > 0) {
    genes = rownames(environment$normalized)[match(capwords(genes),capwords(rownames(environment$normalized)))]
    nExclude = sum(is.na(genes))
    if (nExclude > 0) {
      if (exclude.missing.genes) {
        print.message('Excluding',nExclude,'out of',length(genes),'genes not found in the dataset')
        genes = genes[genes%in%capwords(rownames(environment$normalized))]
      } else {
        print.message('Some signature genes are missing in dataset')
        return(NA)
      }
    }

    background.genes = background.genes (environment, foreground.genes = genes,knn)

    return(Matrix::colMeans(environment$normalized[genes,constrain.cell.universe])-Matrix::colMeans(environment$normalized[background.genes,constrain.cell.universe]))
  } else {
    return(Matrix::colMeans(environment$normalized[genes,constrain.cell.universe]))
  }
}

get.technically.similar.genes <- function (environment, knn = 10) {

  t = start(file.path(environment$work.path, 'tracking'))
  cache = file.path(environment$baseline.data.path,paste(knn,'technical.background.genes.distances.RData',sep='.'))

  if( file.exists(cache) ) {
    print.message('Loading precomputed')
    load(cache)
  } else {
    print.message('Computing')

    technical.variables = data.frame(means=rowMeans(environment$normalized),vars=apply(environment$normalized,1,var))
    HVG.technical.variables = technical.variables[environment$HVG,]
    scaled.technical.variables = apply(technical.variables,2,scale)
    rownames(scaled.technical.variables) = rownames(technical.variables)
    HVG.scaled.technical.variables = apply(HVG.technical.variables,2,scale)
    rownames(HVG.scaled.technical.variables) = rownames(HVG.technical.variables)

    distances = as.matrix(dist(scaled.technical.variables))
    knns = array('',c(length(environment$genes),knn))
    rownames(knns) = environment$genes
    for(index in seq(length(environment$genes))) {
      gene.dist = distances[environment$genes[index],]
      knns[index,] = names(gene.dist[order(gene.dist)[2:(knn+1)]])
    }
    save(knns,technical.variables,file=cache)
  }
  end(t)

  return(list(knns=knns,technical.variables=technical.variables))
}

background.genes <- function (environment, foreground.genes,knn) {

  t = start(file.path(environment$work.path, 'tracking'))
  foreground.genes = foreground.genes[foreground.genes%in%environment$genes]
  technically.similar.genes = get.technically.similar.genes (environment, knn)
  knns = technically.similar.genes$knns
  technical.variables = technically.similar.genes$technical.variables

  background.genes = unique(setdiff(as.vector(knns[foreground.genes,]),foreground.genes))
  print.message('Head foreground.genes technical.variables');print(head(technical.variables[foreground.genes,]))
  print.message('Head background.genes technical.variables');print(head(technical.variables[background.genes,]))
  print.message('background.genes');print(background.genes)
  end(t)
  return(background.genes)
}

regress.covariates <- function (environment, regress,data,groups,rerun = F,save = F) {

  cache = file.path(environment$res.data.path,paste(paste(colnames(regress),collapse='+'),'HVG.regressed.covariates.RData',sep='_'))

  if( !rerun && file.exists(cache) ) {
    print.message('Loading precomputed')
    load(cache)
  } else {
    print.message('Computing')
    t = start(file.path(environment$work.path, 'tracking'))

    formula.str = paste('gene',paste(colnames(regress),collapse=' + '),sep=' ~ ')
    formula = as.formula(formula.str)
    print.message('Regressing:',formula.str)
    print.message('Not regressed matrix');corner(data)
    corrected = data

    for (group in unique(groups)) {
      group.indices = groups==group
      for (gene in rownames(data)) {
        lm.data = data.frame(gene=data[gene,group.indices],regress[group.indices,])#,raw.gene=raw.data[gene,group.indices]
        colnames(lm.data)[2:ncol(lm.data)] = colnames(regress)
        model = lm(formula,lm.data)
        corrected[gene,group.indices] = model$residuals
      }
    }
    # dev.off()
    if (save) save(corrected,file=cache)
    end(t)
  }

  return (corrected)
}

#' Summarize the clustering results into figures
#' 
#'
#' @param environment The environment object
#' @param perplexity The perplexity parameter for tSNE anlaysis
#' @param max_iter Maximum iterations for tSNE
#' @param rerun Whether to rerun 
#' @param order Order in which to plot the clusters 
#' @param contrast Either 'all' indicating differential expression between one cluster against all others or 'datasets' indicating differential expression analysis comparing one cluster to all other within each dataset separately ('datasets' should be used in pooled analysis for optimal results)
#' @param min.fold Minimum fold change for filtering final differentially expressed gene lists
#' @param quantile q-value cutoff
#' @export
summarize <- function (environment, perplexity = seq(10,30,10),max_iter = 10000,rerun = F, order = NA,contrast = 'all',min.fold = 1.5,quantile = 0.95) {

  cluster.size = table(environment$cluster.names)
  if (length(order) == 1 && is.na(order)) order = names(cluster.size)[order(cluster.size,decreasing=T)]
  tSNE.job = run.tSNE (environment, perplexity, max_iter, rerun)
  plot.PCA (environment, quantile = 0.05,order)
  plot.cluster.stats (environment, membership = environment$cluster.names,order = order)
  if (length(environment$seurat.cluster.association)>1) tryCatch({ plot.cluster.stats (environment, membership = environment$seurat.cluster.association,label = 'Seurat',order = order) },error=function(v) v)

  final.diff = run.diff.expression (environment, clustering = environment$clustering,min.fold,quantile,label='main', rerun = rerun,contrast = contrast)

  order = sort(unique(environment$cluster.names))
  plot.heatmaps (environment, diff.exp = final.diff,membership = environment$cluster.names,order = order)
  if (length(environment$seurat.cluster.association)>1) tryCatch({ plot.heatmaps (environment, diff.exp = final.diff,membership = environment$seurat.cluster.association,label = 'Seurat') },error=function(v) v)
  plot.tSNE (environment, tSNE.job,perplexity,max_iter)
}
