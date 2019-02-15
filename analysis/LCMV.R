source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
dir.create("~/Downloads/LCMV")
rep1_path <- GEOquery::getGEOSuppFiles("GSM3423794", baseDir = "~/Downloads/LCMV")
rep2_path <- GEOquery::getGEOSuppFiles("GSM3423795", baseDir = "~/Downloads/LCMV")
file.rename("~/Downloads/LCMV/GSM3423794", "~/Downloads/LCMV/LCMV1")
file.rename("~/Downloads/LCMV/GSM3423795", "~/Downloads/LCMV/LCMV2")


library(MagenSingleCell)
exhaustion_markers <- c('Pdcd1','Cd244','Havcr2','Ctla4','Cd160','Lag3','Tigit','Cd96')
env <- initialize.project(datasets = "LCMV1",
                          origins = "CD44+ cells",
                          experiments = "Rep1",
                          data.path = "~/Downloads/LCMV",
                          work.path = "~/Downloads/LCMV_analysis")
env <- read.data(env, subsample = 1000)
env <- get.variable.genes(env) 

ribosomal.score = ribosomal.score(env)	
mitochondrial.score = mitochondrial.score(env)
cell.cycle.score = cell.cycle.score(env)
Exhaustion = controlled.mean.score(env, exhaustion_markers)

env <- add.confounder.variables (env, ribosomal.score = ribosomal.score, 
  mitochondrial.score = mitochondrial.score,
  cell.cycle.score = cell.cycle.score,
  Exhaustion = Exhaustion)
env <- PCA(env)
env <- cluster.analysis(env)
summarize(env)
types = rbind(
                data.frame(type='Tfh',gene=c('Tcf7','Cxcr5','Bcl6')),
                data.frame(type='Th1',gene=c('Cxcr6','Ifng','Tbx21')),
                data.frame(type='Tcmp',gene=c('Ccr7','Bcl2','Tcf7')),
                data.frame(type='Treg',gene=c('Foxp3','Il2ra')),
                data.frame(type='Tmem',gene=c('Il7r','Ccr7')),
                data.frame(type='CD8',gene=c('Cd8a')),
                data.frame(type='Cycle',gene=c('Mki67','Top2a','Birc5'))
)

env_cluster_names = get.cluster.names(env, types)
env <- set.cluster.names(env, names = env_cluster_names)
summarize(env)

######

env2 <- initialize.project(datasets = "LCMV2",
                          origins = "CD44+ cells",
                          experiments = "Rep2",
                          data.path = "~/Downloads/LCMV",
                          work.path = "~/Downloads/LCMV_analysis")
env2 <- read.data(env2, subsample = 1000)
env2 <- get.variable.genes(env2) 

ribosomal.score = ribosomal.score(env2)	
mitochondrial.score = mitochondrial.score(env2)
cell.cycle.score = cell.cycle.score(env2)
Exhaustion = controlled.mean.score(env2, exhaustion_markers)

env2 <- add.confounder.variables (env2, ribosomal.score = ribosomal.score,
  mitochondrial.score = mitochondrial.score,
  cell.cycle.score = cell.cycle.score,
  Exhaustion = Exhaustion)

env2 <- PCA(env2)
env2 <- cluster.analysis(env2)
summarize(env2)
env2_cluster_names <- get.cluster.names(env2, types)
env2 <- set.cluster.names(env2, names = env2_cluster_names)
summarize(env2)

# pooled analysis
pooled_env <- initialize.project(datasets = c("LCMV1", "LCMV2"), 
                          origins = rep("CD44+ cells", 2),
                          experiments = c("Rep1", "Rep2"),
                          data.path = "~/Downloads/LCMV",
                          work.path = "~/Downloads/LCMV_analysis")

pooled_env <- read.preclustered.datasets(pooled_env)

ribosomal.score = ribosomal.score(pooled_env)	
mitochondrial.score = mitochondrial.score(pooled_env)
cell.cycle.score = cell.cycle.score(pooled_env)
Exhaustion = controlled.mean.score(pooled_env, exhaustion_markers)

pooled_env <- add.confounder.variables (pooled_env, ribosomal.score = ribosomal.score,
  mitochondrial.score = mitochondrial.score,
  cell.cycle.score = cell.cycle.score,
  Exhaustion = Exhaustion)

pooled_env <- PCA(pooled_env, clear.previously.calculated.clustering = F)

summarize(pooled_env, contrast = "datasets")

# cluster similarity analysis
cluster.similarity <- compare.cluster.similarity(pooled_env)
similarity <- cluster.similarity$similarity
map <- cluster.similarity$map
visualize.cluster.cors.heatmaps(pooled_env, pooled_env$work.path,similarity)

filtered.similarity <- get.robust.cluster.similarity(pooled_env, similarity,min.sd = qnorm(.9),max.q.val = 0.01,rerun = F)
robust.clusters <- sort(unique(c(filtered.similarity$cluster1,filtered.similarity$cluster2)));length(robust.clusters)
head(filtered.similarity);dim(filtered.similarity)
summary(filtered.similarity)
visualize.cluster.cors.heatmaps (pooled_env, pooled_env$work.path,filtered.similarity)
similarity = filtered.similarity

# visualize cluster similarity metrics
net <- igraph::graph_from_data_frame(d=similarity[similarity$similarity>0.1,c('name1','name2')], directed=F) 
deg <- igraph::degree(net, mode="all")
igraph::V(net)$size <- deg
l <- igraph::layout_with_kk(net)

similarity.summary.df = data.frame(cluster1 = similarity$name1,cluster2 = similarity$name2,coef = similarity$similarity)
mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
similarity.summary.df = rbind(similarity.summary.df,mirror)
similarity.matrix = reshape2::acast(similarity.summary.df, cluster1~cluster2, value.var="coef")

pdf(file.path(environment$work.path,paste('hclust.dist.Cor.FC.pdf',sep='_')),width=20,height=20)
similarity.matrix[1:5,1:5]
dist(similarity.matrix)
hc.dist = hclust(as.dist(1-similarity.matrix))
plot(hc.dist)
colors = rev(RColorBrewer::brewer.pal(5, "PuOr"))
color.palette = colorRampPalette(colors)
hc.dist = hclust(as.dist(1-similarity.matrix))
clusters = cutree(hc.dist, k = 8)
sort(clusters)
clusters.ordered = clusters[match(names(igraph::V(net)),names(clusters))]
mark.groups = lapply(unique(clusters.ordered),function(c) as.vector(which(clusters.ordered==c)))
lapply(unique(clusters.ordered),function(c) names(clusters.ordered[clusters.ordered==c]))
plot(net, mark.groups=mark.groups, layout=l)
print(gplots::heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Correlation',margins = c(10,10)))
print(gplots::heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = T, Colv = T,dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
dev.off()

