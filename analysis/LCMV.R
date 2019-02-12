source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
dir.create("~/Downloads/LCMV")
rep1_path <- GEOquery::getGEOSuppFiles("GSM3423794", baseDir = "~/Downloads/LCMV")
rep2_path <- GEOquery::getGEOSuppFiles("GSM3423795", baseDir = "~/Downloads/LCMV")
file.rename("~/Downloads/LCMV/GSM3423794", "~/Downloads/LCMV/LCMV1")
file.rename("~/Downloads/LCMV/GSM3423795", "~/Downloads/LCMV/LCMV2")


library(MagenSingleCell)
Exhaustion <- c('Pdcd1','Cd244','Havcr2','Ctla4','Cd160','Lag3','Tigit','Cd96')
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
Exhaustion = controlled.mean.score(env, Exhaustion)

env <- add.confounder.variables (env, ribosomal.score = ribosomal.score, 
  mitochondrial.score = mitochondrial.score,
  cell.cycle.score = cell.cycle.score,
  Exhaustion = Exhaustion)
env <- PCA(env)
env <- cluster.analysis(env)
env <- summarize(env)


# pooled analysis
env <- initialize.project(datasets = c("LCMV1", "LCMV2"), 
                          origins = rep("CD44+ cells", 2),
                          experiments = c("Rep1", "Rep2"),
                          data.path = "~/Downloads/LCMV",
                          work.path = "~/Downloads/LCMV_analysis",
                          marker.genes = unique(marker_genes$symbol))
env <- read.data(env, subsample = 1000)


