#' Plot PCA
#'
#' TODO: description
#'
#' @param environment The environment object
#' @param quantile quantile
#' @param order order
#' @import GGally
#' @import ggrepel
plot.PCA <- function (environment, quantile,order) {

  work.path = environment$work.path
  PCA = environment$PCA
  Rotation = environment$Rotation
  cluster.names = environment$cluster.names
  dataset = environment$dataset.labels
  marker.genes = environment$marker.genes

  drivers = apply(Rotation,1,function(v) { q = quantile(v,c(quantile,1-quantile)); names(v[v<=q[1] | v>=q[2]])})

  nPCs = 3

  drivers = apply(Rotation,1,function(v) { q = quantile(v,c(quantile,1-quantile)); names(v[v<=q[1] | v>=q[2]])})
  pdf(file.path(work.path,'Rotation.PCA.pdf'))
  data = data.frame(gene = colnames(Rotation),Rot = t(Rotation))
  for (row in seq(nrow(Rotation)-1)) {
    if(typeof(drivers)=='list') {
      gene.set = as.vector(unlist(drivers[c(row,row+1)]))
    } else {
      gene.set = as.vector(drivers[,c(row,row+1)])
    }
    plot.data = data[data$gene %in% gene.set,c(1,row+1,row+2)]
    print(ggplot(plot.data, aes_string(x=colnames(plot.data)[2],y=colnames(plot.data)[3],label = 'gene')) + geom_text(check_overlap = TRUE,size = 2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")))
    print(ggplot(plot.data, aes_string(x=colnames(plot.data)[2],y=colnames(plot.data)[3],label = 'gene')) + geom_text(check_overlap = TRUE,size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")))
    plot.data = plot.data[plot.data$gene%in%marker.genes,]
    print(ggplot(plot.data, aes_string(x=colnames(plot.data)[2],y=colnames(plot.data)[3],label = 'gene')) + geom_point() + geom_text_repel() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")))
  }
  dev.off()

  pdf(file.path(work.path,'all.PCA.pdf'))
  data = data.frame(t(PCA),Cluster = factor(cluster.names),Dataset = factor(dataset));head(data)
  for (row in seq(nrow(PCA)-1)) {
    print(ggplot(data, aes_string(x=rownames(PCA)[row], y=rownames(PCA)[row+1], color='Cluster')) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")))
  }
  if (length(unique(dataset))>1) {
    for (row in seq(nrow(PCA)-1)) {
      print(ggplot(data, aes_string(x=rownames(PCA)[row], y=rownames(PCA)[row+1], color='Dataset')) + geom_point() + scale_colour_brewer(palette = 'Set3') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")))
    }
  }
  dev.off()

  data = data.frame(t(PCA),Cluster = factor(cluster.names),Dataset = factor(dataset));head(data)
  pdf(file.path(work.path,'PC.scores.histogram.pdf'),width = 10)
  for (row in seq(nrow(PCA)-1)) {
    print(ggplot(data, aes_string(x=rownames(PCA)[row], fill = 'Cluster')) + geom_density(alpha = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Density') + theme_classic(base_size=25))
  }
  if (length(unique(dataset))>1) {
    for (row in seq(nrow(PCA)-1)) {
      print(ggplot(data, aes_string(x=rownames(PCA)[row], fill = 'Dataset')) + geom_density(alpha = 0.5) + scale_fill_brewer(palette = 'Paired') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Density') + theme_classic(base_size=25))
    }
  }
  dev.off()

  pdf(file.path(work.path,'PC.scores.heatmap.pdf'),width = 5,height = 8)
  for (row in seq(nrow(PCA)-1)) {
    pc = environment$Rotation[row,]
    high = names(pc[pc>quantile(pc,0.99)]);high
    low = names(pc[pc<quantile(pc,0.01)]);low
    gene.list = list(markers=c(low,high))
    names(gene.list) = paste('PC',row)
    plot.expression.heatmap.based.on.FC.marker( environment$normalized,environment$cluster.names,gene.list=gene.list,scale='row',counts = F,order=order,filter.diff.exp=F,cellnote=F )
  }
  if (length(unique(dataset))>1) {
    for (row in seq(nrow(PCA)-1)) {
      pc = environment$Rotation[row,]
      high = names(pc[pc>quantile(pc,0.99)]);high
      low = names(pc[pc<quantile(pc,0.01)]);low
      gene.list = list(markers=c(low,high))
      names(gene.list) = paste('PC',row)
      plot.expression.heatmap.based.on.FC.marker( environment$normalized,environment$dataset.labels,gene.list=gene.list,scale='row',counts = F,order=NA,filter.diff.exp=F,cellnote=F )
    }
  }
  dev.off()
}

plot.cluster.stats <- function (environment, membership,label = NA,order = NA) {
  file.name = 'cluster.size.pdf'
  work.path = environment$work.path
  if (!is.na(label)) {
    work.path = file.path(environment$work.path,label)
    dir.create(work.path,showWarnings = F)
    file.name = paste(label,file.name,sep='.')
  }
  if (is.na(order)) order = order(table(membership),decreasing=T)
  pdf(file.path(work.path,file.name),width=8,height=5)
  data = data.frame(clustering=factor(membership,levels = order),Dataset=factor(environment$dataset.labels))
  if (length(unique(environment$dataset.labels))>1) {
    print(ggplot(data, aes(clustering,fill=Dataset)) + geom_bar() + scale_fill_brewer(palette = 'Set3') + xlab('Cluster ID') + ylab('Number of cells') + theme_classic(base_size=15) + theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + guides(fill=guide_legend(nrow=max(1,floor(length(unique(environment$dataset.labels))/2)),byrow=F)))
    print(ggplot(data, aes(clustering,fill=Dataset)) + geom_bar(aes(y = (..count..)/sum(..count..))) + scale_y_continuous(labels=scales::percent) + ylab("relative frequencies") + scale_fill_brewer(palette = 'Set3') + xlab('Cluster ID') + ylab('Number of cells') + theme_classic(base_size=15) + theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + guides(fill=guide_legend(nrow=max(1,floor(length(unique(environment$dataset.labels))/2)),byrow=F)))
    print(ggplot(data, aes(Dataset,fill=clustering)) + geom_bar() + xlab('Cluster ID') + ylab('Number of cells') + theme_classic(base_size=15) + theme(axis.text.x = element_text(angle = 25, hjust = 1)))
  }
  for (dataset in unique(environment$dataset.labels)) {
    print(ggplot(data[data$Dataset==dataset,], aes(clustering)) + geom_bar(aes(y = (..count..)/sum(..count..))) + scale_y_continuous(labels=scales::percent) + ylab("relative frequencies") + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + xlab('Cluster ID') + ylab('Number of cells') + ggtitle(dataset) + theme_classic(base_size=15))
  }
  print(ggplot(data, aes(clustering)) + geom_bar() + theme(axis.text.x = element_text(angle = 25, hjust = 1)))
  print(ggplot(data, aes(clustering)) + geom_bar(aes(y = (..count..)/sum(..count..))) + scale_y_continuous(labels=scales::percent) + ylab("relative frequencies") + theme(axis.text.x = element_text(angle = 25, hjust = 1)))

  dev.off()

  cell.confusion = table(membership,environment$dataset.labels)
  write.csv(cell.confusion,file=file.path(environment$work.path,'cell.confusion.csv'))

  confounders = environment$confounders
  data = data.frame(clustering=factor(as.vector(membership)),dataset=factor(environment$dataset.labels),confounders);head(data)
  file.name = 'confounder.stats.violin.pdf'
  if (!is.na(label)) file.name = paste(label,file.name,sep='.')
  pdf(file.path(work.path,file.name))
  for (variable in colnames(confounders)) {
    print(ggplot(data, aes_string('clustering',variable)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "width") + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ggtitle(variable))
  }
  # if (length(unique(environment$dataset.labels))>1) {
  for (variable in colnames(confounders)) {
    print(ggplot(data, aes_string('dataset',variable)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "width") + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ggtitle(variable))
  }
  # }
  dev.off()
}

plot.tSNE <- function (environment, tSNE.job,perplexity,max_iter,membership = NA) {

  if (!is.na(tSNE.job)) {
    tSNE = get_slurm_out(tSNE.job)
    tryCatch({cleanup_files(tSNE.job);cleanup_files(tSNE.job);cleanup_files(tSNE.job)},error=function(v) v)
  }

  if (is.na(membership)) membership = environment$cluster.names

  configs = expand.grid(perplexity,max_iter)

  for (row in seq(nrow(configs))) {
    perplexity = configs[row,1]
    max_iter = configs[row,2]
    tryCatch({
      load(file.path(environment$res.data.path,'tSNEs',paste(perplexity,max_iter,'tSNE.RData',sep='.')))

      duplicated.indices = duplicated(t(environment$PCA))
      data = data.frame(Cluster = factor(membership[!duplicated.indices]),Origin = factor(environment$origins[!duplicated.indices]),Experiment = factor(environment$experiments[!duplicated.indices]),tSNE = tSNE);head(data)

      pdf(file.path(environment$work.path,paste('tSNE_perplexity',perplexity,'max_iter',max_iter,'pdf',sep='.')),width=10,height=10)
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Cluster,shape=Origin)) + geom_point(data = data,size = 4,alpha = 0.6) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Cluster, shap = Origin'))# + scale_color_brewer(palette = 'Set2')
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Cluster,shape=Origin)) + geom_point(data = data,size = 5,alpha = 0.6) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Cluster, shap = Origin'))# + scale_color_brewer(palette = 'Set2')
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Cluster,shape=Origin)) + geom_point(data = data,size = 6,alpha = 0.6) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Cluster, shap = Origin'))# + scale_color_brewer(palette = 'Set2')
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Cluster,shape=Origin)) + geom_point(data = data,size = 7,alpha = 0.4) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Cluster, shap = Origin'))# + scale_color_brewer(palette = 'Set2')
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Cluster,shape=Origin)) + geom_point(data = data,size = 7,alpha = 0.6) + scale_shape(solid = F) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Cluster, shap = Origin'))# + scale_color_brewer(palette = 'Set2')
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Cluster,shape=Experiment)) + geom_point(data = data,size = 7,alpha = 0.4) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Cluster, shape = Experiment'))# + scale_color_brewer(palette = 'Set2')
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Experiment,shape=Origin)) + geom_point(data = data,size = 7,alpha = 0.4) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Experiment, shape = Origin') + scale_color_brewer(palette = 'Set2'))
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Origin,shape=Experiment)) + geom_point(data = data,size = 7,alpha = 0.4) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Origin, shape = Experiment') + scale_color_brewer(palette = 'Set2'))
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Origin,shape=Origin)) + geom_point(data = data,size = 7,alpha = 0.4) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Origin, shape = Origin') + scale_color_brewer(palette = 'Set2'))
      print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=Experiment,shape=Experiment)) + geom_point(data = data,size = 7,alpha = 0.4) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + ggtitle('color = Experiment, shape = Experiment') + scale_color_brewer(palette = 'Set2'))

      if (ncol(environment$confounders)>0) {
        data = data.frame(Cluster = factor(membership[!duplicated.indices]),Origin = factor(environment$origins[!duplicated.indices]),Experiment = factor(environment$experiments[!duplicated.indices]),environment$confounders,tSNE = tSNE);head(data)
        # name = colnames(environment$confounders)[2]
        for (name in colnames(environment$confounders)) {
          signature.activation = data[[name]]
          if (is.numeric(signature.activation)) {
            print(ggplot(data, aes(x=tSNE.1, y=tSNE.2, color=signature.activation,shape=Origin)) + geom_point(data = data,size = 4,alpha = 0.6) + scale_shape(solid = T) + xlab('tSNE 1') + ylab('tSNE 2') + theme_classic(base_size=20) + theme(legend.position="bottom") + scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + guides(color=FALSE) + ggtitle(name))
          }
        }
      }
    },error=function(e) e)
    dev.off()
  }
}

plot.expression.heatmap.based.on.FC.marker <- function( measurements,clustering,gene.list,counts = F,order = NA,RowSideColors=NA,scale='row',save=NA,filter.diff.exp=T,cellnote=T,exponent=F,doMeans=T,srtCol=45,multiplication=100,rounding=0,breaks=50,key=T,sort.rows=T,sort.cols=T,Rowv=F,Colv=F,dendrogram='none' ) {

  colors = c('#d60c0c','#ffb1ad','#f4f4f4','#aed4fc','#0050ba')#c('#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4')
  colors = colors[length(colors):1]#red should correspond to high
  color.palette = colorRampPalette(colors)
  list.row = 1
  titles = names(gene.list)
  plots = list()
  nclusters = length(unique(clustering))

  list.row=1

  for( list.row in seq(length(gene.list)) ) {
    tryCatch({
      genes = gene.list[[list.row]]

      if( length(genes) == 0 ) next
      if( length(genes) == 1 ) genes = rep(genes,2)
      # genes = c(markers,genes)
      genes = genes[genes%in%rownames(measurements)]

      measurements.diff.exp = measurements[match(genes,rownames(measurements)),]
      if (doMeans) {
        names = rownames(measurements.diff.exp)
        # if (exponent) measurements.diff.exp = exp(measurements.diff.exp)
        mat = {}
        row=1
        if (counts) measurements.diff.exp = measurements.diff.exp>0
        for( row in seq(length(genes)) ) {
          averages = apply.by.group(clustering,measurements.diff.exp[row,],mean)
          mat = rbind(mat,averages)
        }
        rownames(mat) = names
        if (sort.cols) mat = mat[,order(as.numeric(colnames(mat)))]
      } else {
        mat = measurements.diff.exp
      }

      if (exponent) mat = exp(mat)

      if (is.na(order)) order = seq(nclusters)
      others = mat
      others = others[,order]
      if (sort.rows) {
        row.order = order(unlist(apply(others,1,which.max)),decreasing=T)
      } else {
        row.order = seq(nrow(others))
      }
      others = others[row.order,]
      RowSideColors = RowSideColors[row.order]
      if( filter.diff.exp ){
        significant.diff = getDE.limma( Y = measurements[rownames(measurements)%in%rownames(others),], group = factor(clustering) )
        keep = rownames(others)%in%rownames(significant.diff)
        others = others[keep,]
        RowSideColors = RowSideColors[keep]
      }

      if(length(unique(RowSideColors))>1 & sort.rows){
        blocks = as.vector(unlist(apply(others,1,which.max)))
        index=1
        for( index in unique(blocks) ) {
          indices = which(blocks == index)
          if( length(indices) < 2 ) next
          order = order(RowSideColors[indices],as.vector(others[indices,index]))
          RowSideColors[indices] = RowSideColors[indices][order]
          others[indices,] = others[indices,][order,]
          rownames(others)[indices] = rownames(others)[indices][order]
        }
      }

      RowSideColors[is.na(RowSideColors)] = 'white'
      RowSideColors[RowSideColors=='-1'] = 'blue'
      RowSideColors[RowSideColors=='1'] = 'red'

      if( cellnote == T ) {
        cellnote = round(multiplication*others,rounding)
        plots[[titles[list.row]]] = heatmap.2(as.matrix(others), col=color.palette, scale=scale, key=key, cexRow=1, cexCol=1, density.info="none", trace="none",srtCol=srtCol,adjCol = c(1,1),Rowv = Rowv, Colv = Colv,margins = c(8,5),dendrogram = dendrogram,main = titles[list.row],RowSideColors=RowSideColors,breaks=breaks,cellnote=cellnote,notecol='white')#100*round((2^others)-1,2)
      } else {
        plots[[titles[list.row]]] = heatmap.2(as.matrix(others), col=color.palette, scale=scale, key=key, cexRow=1, cexCol=1, density.info="none", trace="none",srtCol=srtCol,adjCol = c(1,1),Rowv = Rowv, Colv = Colv,margins = c(8,5),dendrogram = dendrogram,main = titles[list.row],RowSideColors=RowSideColors,breaks=breaks)
      }

      if( !is.na(save) ) {
        write.csv(as.matrix(others),file=save)
      }
    },error = function(e) print(e))
  }

  rownames(others)
}

plot.simple.heatmap <- function (environment, name,path = NA,markers,membership = NA,normalized = NA,order = NA,width = 5, height = 5,scale='row',RowSideColors = NA,counts = F,filter.diff.exp=F,cellnote=F,key=F,save=NA,sort.rows=T,sort.cols=T,Colv=F,Rowv=F,dendrogram = 'none') {

  if (is.na(path)) path = environment$work.path
  pdf(file=file.path(path,paste(name,'heatmap.pdf',sep='.')),width = width, height = height)

  if (is.na(membership)) membership = environment$cluster.names
  if (is.na(normalized)) normalized = environment$normalized

  cluster.size = table(membership)
  if (is.na(order)) order = names(cluster.size)[order(cluster.size,decreasing=T)]

  plot.expression.heatmap.based.on.FC.marker(
    normalized,
    membership,
    gene.list = list(markers = markers),
    scale = scale,
    RowSideColors = RowSideColors,
    counts = counts,
    order = order,
    filter.diff.exp = filter.diff.exp,
    cellnote = cellnote,
    doMeans = T,
    exponent = ifelse(counts,F,T),
    multiplication = ifelse(counts,100,10),
    rounding = 0,
    key = key,
    save = save,
    sort.rows = sort.rows,
    sort.cols = sort.cols,
    Rowv = Rowv, Colv = Colv,dendrogram = dendrogram)

  if (key & length(RowSideColors)>1) {
    legend("topright",
           legend = unique(RowSideColors),
           col = unique(as.numeric(RowSideColors)),
           lty= 1,
           lwd = 5,
           cex=.7
    )
  }

  dev.off()
}

plot.violin <- function (environment, genes,types,fore1exp1,fore2exp1,fore1exp2,fore2exp2,back1exp1,back2exp1,back1exp2,back2exp2,path,height = 5,width = 5,scale = T,palette = "Greys",separate.background = F) { #Pastel2

  violin.plot.data = {}
  cluster.indices = environment$cluster.names %in% c(fore1exp1,fore2exp1,fore1exp2,fore2exp2,back1exp1,back2exp1,back1exp2,back2exp2)
  for (gene in genes) {
    if (scale) {
      expression = environment$normalized[gene,cluster.indices]
      expression = (expression-min(expression))/(max(expression) - min(expression))
    } else {
      expression = environment$normalized[gene,cluster.indices]
    }
    violin.plot.data = rbind(
      violin.plot.data,
      data.frame(cluster = environment$cluster.names[cluster.indices],gene = gene,expression = expression)
    )
  }

  violin.plot.data$cell_type = NA
  violin.plot.data$cell_type[violin.plot.data$cluster%in%c(fore1exp1,fore1exp2)] = types[1]
  violin.plot.data$cell_type[violin.plot.data$cluster%in%c(fore2exp1,fore2exp2)] = types[2]

  if (separate.background) {
    violin.plot.data$cell_type[violin.plot.data$cluster%in%c(back1exp1,back1exp2)] = 'Background 1'
    violin.plot.data$cell_type[violin.plot.data$cluster%in%c(back2exp1,back2exp2)] = 'Background 2'
    violin.plot.data$cell_type = factor(violin.plot.data$cell_type,levels = c('Background 1',types,'Background 2'))
  } else {
    violin.plot.data$cell_type[violin.plot.data$cluster%in%c(back1exp1,back2exp1,back1exp2,back2exp2)] = 'Others'
    violin.plot.data$cell_type = factor(violin.plot.data$cell_type,levels = c(types,'Others'))
  }

  tryCatch({
    pdf(file.path(path,paste(paste(types,collapse='_'),'violin.pdf',sep='.')),height = height,width = width)
    print(ggplot(violin.plot.data, aes(x = gene, y = expression, fill = cell_type)) + geom_violin(scale = "width",draw_quantiles = c(0.25, 0.5, 0.75)) + theme_classic(base_size=15) + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + scale_fill_brewer(palette = palette))
    dev.off()
  },error=function(v) {
    pdf(file.path(path,paste(paste(types,collapse='_'),'violin.pdf',sep='.')),height = height,width = width)
    print(ggplot(violin.plot.data, aes(x = gene, y = expression, fill = cell_type)) + geom_violin(scale = "width",draw_quantiles = c(0.5)) + theme_classic(base_size=15) + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + scale_fill_brewer(palette = palette))
    dev.off()
  })
}

plot.heatmaps <- function (environment, diff.exp,membership,order = NA,nTopRanked = 10,label = NA) {

  symbols = {}
  for (cluster in unique(diff.exp$cluster)) {
    top.ranked = diff.exp[diff.exp$cluster==cluster,]
    top.ranked = top.ranked[order(top.ranked$fold,decreasing=T),]
    symbols = c(symbols,head(as.vector(top.ranked$gene),nTopRanked))
  }
  symbols = unique(symbols);symbols

  clustering = as.vector(membership)

  file.name = 'diff.genes.pdf'
  work.path = environment$work.path
  if (!is.na(label)) {
    work.path = file.path(environment$work.path,label)
    dir.create(work.path,showWarnings = F)
    file.name = paste(label,file.name,sep='.')
  }
  pdf(file=file.path(work.path,file.name),width = length(unique(clustering))/2, height = length(symbols)^(1/1.5))
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=symbols),order=order,filter.diff.exp=T,cellnote = T,doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=symbols),order=order,filter.diff.exp=F,cellnote = T,doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=symbols),counts = T,order=order,filter.diff.exp=F,cellnote = T,doMeans = T,exponent = F,multiplication = 100,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=symbols),counts = T,order=order,filter.diff.exp=T,cellnote = T,doMeans = T,exponent = F,multiplication = 100,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=symbols),order=order,filter.diff.exp=F,scale='col',cellnote = T,doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=symbols),order=order,filter.diff.exp=F,scale='none',cellnote = T,doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  dev.off()

  extended.genes = environment$marker.genes
  file.name = 'all.markers.pdf'
  if (!is.na(label)) file.name = paste(label,file.name,sep='.')
  pdf(file=file.path(work.path,file.name),width = length(unique(clustering))/2, height = length(extended.genes)^(1/1.5))
  plot.expression.heatmap.based.on.FC.marker( measurements=environment$normalized,clustering,gene.list=list(markers=extended.genes),order=order,filter.diff.exp=T,doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( measurements=environment$normalized,clustering,gene.list=list(markers=extended.genes),counts = T,order=order,filter.diff.exp=T,doMeans = T,exponent = F,multiplication = 100,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=extended.genes),order=order,filter.diff.exp=T,scale='col',doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=extended.genes),order=order,filter.diff.exp=F,scale='col',doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=extended.genes),order=order,filter.diff.exp=T,scale='none',doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  plot.expression.heatmap.based.on.FC.marker( environment$normalized,clustering,gene.list=list(markers=extended.genes),order=order,filter.diff.exp=F,scale='none',doMeans = T,exponent = T,multiplication = 10,rounding = 0)
  dev.off()
}


visualize.cluster.cors.heatmaps <- function (environment, work.path,similarity) {

  name = 'cross.sample.similarities'
  work.path = file.path(environment$work.path,name)
  if (file.exists(work.path)) {
    new.dir = file.path(environment$work.path,paste(name,format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),sep='---'))
    file.rename(work.path,new.dir)
  }
  dir.create(work.path,showWarnings = T,recursive = T)


  samples.names = unique(c(as.vector(similarity$sample1),as.vector(similarity$sample2)))
  samples.names = sort(samples.names)
  sample.name1 = samples.names[3]
  sample.name2 = samples.names[4]
  for (sample1 in seq(length(samples.names)-1)) {
    sample.name1 = samples.names[sample1]
    for (sample2 in (sample1+1):length(samples.names)) {
      sample.name2 = samples.names[sample2]
      sub.work.path = file.path(work.path,sample.name1,sample.name2)
      dir.create(sub.work.path,showWarnings = T,recursive = T)
      file.name = paste(paste(strsplit(sample.name1,split=' ')[[1]],collapse='.'),'to',paste(strsplit(sample.name2,split=' ')[[1]],collapse='.'),sep='_')
      samples = c(sample.name1,sample.name2)#unique(c(as.vector(similarity$sample1[similarity$sample1 == sample.name1]),as.vector(similarity$sample2[similarity$sample2 == sample.name2])))
      sample = similarity[(similarity$sample1 %in% samples[1] & similarity$sample2 %in% samples[2]) | (similarity$sample1 %in% samples[2] & similarity$sample2 %in% samples[1]),]

      similarity.summary.df = data.frame(cluster1 = sample$name1,cluster2 = sample$name2,coef = sample$similarity)
      similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
      write.csv(sample,file=file.path(sub.work.path,paste(file.name,'similarity.csv',sep='_')))
      similarity.matrix = t(similarity.matrix)

      ocl.similarity.summary.df = data.frame(cluster1 = sample$name1,cluster2 = sample$name2,coef = -scale(sample$ocldist))
      ocl.similarity.matrix = acast(ocl.similarity.summary.df, cluster1~cluster2, value.var="coef")
      ocl.similarity.matrix = t(ocl.similarity.matrix)

      pdf(file.path(sub.work.path,paste(file.name,'similarity.heatmap.pdf',sep='_')),width=10,height=10)
      colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
      print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, srtCol = 45, scale = 'none', density.info="none", trace="none",Rowv = T, Colv = T,dendrogram = 'both',margins=c(15,15),cellnote=round(similarity.matrix,1),notecol='white',main='Cluster FC Correlations'))
      print(heatmap.2(ocl.similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, srtCol = 45, scale = 'none', density.info="none", trace="none",Rowv = F, Colv = F,dendrogram = 'none',margins=c(15,15),cellnote=round(ocl.similarity.matrix,1),notecol='white',main='Cluster mean Euclidean similarity'))
      dev.off()
    }
  }
}
