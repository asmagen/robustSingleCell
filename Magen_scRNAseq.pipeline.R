
# scRNAseq analysis pipeline methods ______________________________________________________________________________________
# 
# This pipeline was designed, written and produced by Assaf Magen as part of CD4+ T tumor-infiltrating lymphocytes scRNAseq study
# Single-cell resolution profiling of tumor-reactive CD4+ T-cells reveals immune landscape alterations
# 
# Laboratory of Immune Cell Biology, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, Maryland
# Center for Bioinformatics and Computational Biology, University of Maryland, College Park, Maryland
# 
# Using any piece of this code or its modifications require the appropriate citations and agreement to code usage terms
# Default parameters should be adjusted to account for dataset-specific data features as needed and should not be taken as is
# 
# Contact: Assaf Magen, assaf.magen@nih.gov

# for(i in seq(10)) sink() Use to shift back command line output to user interface rather than documentation tracking files

initialize.project <- function (
	datasets, # List of dataset code names
	origins, # List of dataset tissue origin/condition full name
	experiments, # List of experiment design annotations
	data.path, # Path to where the data is located
	work.path, # path to where the analysis should create folders and store files
	marker.genes, # Set of genes of interest for visualization purposes
	clear.history = F, # Whether you would like to remove any previous project by this name
	analysis.label = NA, # Whether you would like to add a specific label to the analysis folder
	convert.to.mouse.gene.symbols = NULL) { # Whether you are using human gene symbols and would like to convert them to mouse gene symbols
	
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
	setwd(environment$work.path)
	dir.create('tracking', showWarnings = F, recursive = T, mode = "700")
	environment$baseline.data.path = environment$res.data.path = file.path(environment$work.path,'data')
	dir.create(environment$res.data.path, showWarnings = F, recursive = T, mode = "700")
	environment$marker.genes = marker.genes
	t = start()
	print(environment)
	end(t)

	return(environment)
}

read.10x.data <- function( path ) {

	library(Matrix)
	barcode.path <- file.path(path,"barcodes.tsv")
	features.path <- file.path(path,"genes.tsv")
	matrix.path <- file.path(path,"matrix.mtx")
	mat <- readMM(file = matrix.path)
	feature.names = read.delim(features.path, 
	                           header = FALSE,
	                           stringsAsFactors = FALSE)
	barcode.names = read.delim(barcode.path, 
	                           header = FALSE,
	                           stringsAsFactors = FALSE)
	colnames(mat) = barcode.names$V1
	rownames(mat) = feature.names$V2

	return(mat)
}

read.data <- function(
	genome = 'mm10', # Genome annotation
	min.genes.per.cell = 500, # Minimum required number of genes per cell
	max.genes.per.cell.quantile = 0.98, # Upper quantile for number of genes per cell
	max.UMIs.per.cell.quantile = 0.98, # Upper quantile for number of UMIs per cell
	min.cells.per.gene = 1, # Minimum required number of cells per gene
	max.mitochondrial.frac = 0.1, # Maximum fraction of reads mapped to mitochondrial genes per cell
	max.ribosomal.frac = NA, # Maximum fraction of reads mapped to ribosomal genes per cell
	cell.filters = NA, # Filtering option for cells based on marker genes
	raw.data.matrices = NA, # Are you supplying the actual data matrices instead of 10X dataset format?
	rerun = F) { # Whether to rerun loading the dataset

	cache = file.path(environment$baseline.data.path,'data.RData')

  if( !rerun & file.exists(cache) ) {
  	print.message('Loading precomputed')
  	load(cache)
  } else {

  	print.message('Computing')
		t = start()

		merged = NA
		dataset = environment$datasets[1]
		merged.dataset.labels = merged.origins = merged.experiments = merged.criteria = {}
		for(dataset in environment$datasets) {
			if (is.na(raw.data.matrices)) {
				print.message('Loading',dataset)
				measurements = read.10x.data (path = file.path(environment$data.path,dataset))
				measurements = as.matrix(measurements[,colSums(measurements>0)>=min.genes.per.cell])
			} else {
				print.message('Using input',dataset)
				measurements = raw.data.matrices[[dataset]]
			}
			colnames(measurements) = rep(environment$datasets[environment$datasets==dataset],ncol(measurements))
			dataset.labels = rep(paste(environment$origins[environment$datasets==dataset],' (',environment$experiments[environment$datasets==dataset],')',sep=''),ncol(measurements))
			origins = rep(environment$origins[environment$datasets==dataset],ncol(measurements))
			experiments = rep(environment$experiments[environment$datasets==dataset],ncol(measurements))
			# corner(measurements)
			# measurements = measurements[,measurements['Cd4',]>0]

			data = data.frame(
				nUMI=colSums(measurements),
				nGenes=colSums(measurements>0),
				Ribosomal=colMeans(measurements[get.ribo.genes(rownames(measurements)),]),
				Mitochondrial=colMeans(measurements[get.mito.genes(rownames(measurements)),]),
				Ribosomal.frac=colSums(measurements[get.ribo.genes(rownames(measurements)),])/colSums(measurements),
				Mitochondrial.frac=colSums(measurements[get.mito.genes(rownames(measurements)),])/colSums(measurements)
				)
			print(head(data))
			library(ggplot2,quietly=T)
			pdf(file.path(environment$baseline.work.path,paste(dataset,'pre.filter.dataset.stats.pdf',sep='.')))
			for (ind1 in seq(length(colnames(data))-1)) {
				for (ind2 in (ind1+1):(length(colnames(data)))) {
					v1=colnames(data)[ind1]
					v2=colnames(data)[ind2]
					print(ggplot(data, aes_string(v1,v2)) + geom_point())
				}
			}
			dev.off()

			print.message('Original dimensions');print(dim(measurements))

			genes.per.cell = colSums(measurements>0)
			print.message('Original genes.per.cell');print(summary(genes.per.cell))
			print(quantile(genes.per.cell,seq(0.01,1,0.01)))
			UMIs.per.cell = colSums(measurements)
			print.message('Original UMIs.per.cell');print(summary(UMIs.per.cell))
			print(quantile(UMIs.per.cell,seq(0.01,1,0.01)))
			max.genes.per.cell = quantile(genes.per.cell,max.genes.per.cell.quantile)
			print.message('max.genes.per.cell',max.genes.per.cell)
			max.UMIs.per.cell = quantile(UMIs.per.cell,max.UMIs.per.cell.quantile)
			print.message('max.UMIs.per.cell',max.UMIs.per.cell)
			print.message('max.mitochondrial.frac',max.mitochondrial.frac)

			print.message('# not qualify min.genes.per.cell',sum(genes.per.cell<min.genes.per.cell))
			print.message('# not qualify max.genes.per.cell',sum(genes.per.cell>max.genes.per.cell))
			print.message('# not qualify max.UMIs.per.cell',sum(UMIs.per.cell>max.UMIs.per.cell))
			print.message('# not qualify max.mitochondrial.frac',sum(data$Mitochondrial.frac>max.mitochondrial.frac))

			final.gene.filter.passed.cells = rep(T,ncol(measurements))
			if (!is.na(cell.filters) && dataset%in%as.vector(cell.filters$dataset)) {
				dataset.filters = cell.filters[as.vector(cell.filters$dataset)==dataset,]
				filter.genes = as.vector(dataset.filters$gene)
				if (sum(filter.genes%in%rownames(measurements))==0) print.message('All filter genes do not exist in dataset - not filtering anything')
				not.exist.genes = sum(!filter.genes%in%rownames(measurements))
				if (not.exist.genes>0) print.message('Not filtering based on',not.exist.genes,'filter genes do not exist in dataset')
				filter.genes = filter.genes[filter.genes%in%rownames(measurements)]
				if (length(filter.genes)==1) filter.genes = c(filter.genes,filter.genes)
				gene = filter.genes[1]
				matched.conditions = rep(T,ncol(measurements))
				for (gene in filter.genes) {
					matched.conditions = matched.conditions & ((measurements[gene,] > 0) == dataset.filters$expressed[dataset.filters$gene==gene]) #check if matching condition
				}
				sum.not.matching = sum(matched.conditions)
				print.message('# not qualify dataset.filters',sum.not.matching,'out of',ncol(measurements),'[',round(sum.not.matching/ncol(measurements)*100,1),'% ]')
				final.gene.filter.passed.cells = matched.conditions
			}

			criteria = genes.per.cell>=min.genes.per.cell & genes.per.cell<=max.genes.per.cell & UMIs.per.cell<=max.UMIs.per.cell & data$Mitochondrial.frac<=max.mitochondrial.frac & final.gene.filter.passed.cells
			if (!is.na(max.ribosomal.frac)) criteria = criteria & data$Ribosomal.frac<=max.ribosomal.frac
			print.message('# qualifying cells',sum(criteria),'# not qualifying cells',sum(!criteria))
			measurements = measurements[,criteria]
			dataset.labels = dataset.labels[criteria]
			origins = origins[criteria]
			experiments = experiments[criteria]
			genes.per.cell = colSums(measurements>0)
			print.message('Filtered genes.per.cell');print(summary(genes.per.cell))
			print(quantile(genes.per.cell,seq(0.01,1,0.01)))
			UMIs.per.cell = colSums(measurements)
			print.message('Filtered UMIs.per.cell');print(summary(UMIs.per.cell))
			print(quantile(UMIs.per.cell,seq(0.01,1,0.01)))

			data = data.frame(
				nUMI=colSums(measurements),
				nGenes=colSums(measurements>0),
				Ribosomal=colMeans(measurements[get.ribo.genes(rownames(measurements)),]),
				Mitochondrial=colMeans(measurements[get.mito.genes(rownames(measurements)),]),
				Ribosomal.frac=colSums(measurements[get.ribo.genes(rownames(measurements)),])/colSums(measurements),
				Mitochondrial.frac=colSums(measurements[get.mito.genes(rownames(measurements)),])/colSums(measurements)
				)

			print(head(data))
			library(ggplot2,quietly=T)
			pdf(file.path(environment$baseline.work.path,paste(dataset,'post.filter.dataset.stats.pdf',sep='.')))
			for (ind1 in seq(length(colnames(data))-1)) {
				for (ind2 in (ind1+1):(length(colnames(data)))) {
					v1=colnames(data)[ind1]
					v2=colnames(data)[ind2]
					print(ggplot(data, aes_string(v1,v2)) + geom_point())
				}
			}
			dev.off()

			if(length(merged)==1 && is.na(merged)) {
				merged = measurements
				merged.dataset.labels = dataset.labels
				merged.origins = origins
				merged.experiments = experiments
				merged.criteria = criteria
			} else {
				if(sum(rownames(merged) != rownames(measurements))>0) stop('Feature genes mismatch - need to correct dataset binding matching')
				merged = cbind(merged,measurements)
				merged.dataset.labels = c(merged.dataset.labels,dataset.labels)
				merged.origins = c(merged.origins,origins)
				merged.experiments = c(merged.experiments,experiments)
				merged.criteria = c(merged.criteria,criteria)
			}
			print(table(colnames(merged)))
		}

		counts = merged
		dataset.labels = merged.dataset.labels
		origins = merged.origins
		experiments = merged.experiments
		criteria = merged.criteria
		rm(gbm,merged,measurements,merged.dataset.labels,merged.origins,merged.experiments,merged.criteria)

		cells.per.gene = rowSums(counts>0)
		print.message('Filtered cells.per.gene');print(summary(cells.per.gene))
		genes.filter = cells.per.gene>=min.cells.per.gene
		print.message('# qualifying genes',sum(genes.filter),'# not qualifying genes',sum(!genes.filter))

		print.message('Aggregated dataset dim');print(dim(counts))
		genes.per.cell = colSums(counts[genes.filter,]>0)
		print.message('Aggregated dataset genes.per.cell');print(summary(genes.per.cell))
		print(quantile(genes.per.cell,seq(0.01,1,0.01)))
		UMIs.per.cell = colSums(counts[genes.filter,])
		print.message('Aggregated dataset UMIs.per.cell');print(summary(UMIs.per.cell))
		print(quantile(UMIs.per.cell,seq(0.01,1,0.01)))

		rownames = data.frame(old=rownames(counts))
		if(sum(duplicated(rownames(counts)))>0) {
		  gene = unique(rownames(counts)[duplicated(rownames(counts))])[2]
		  for( gene in unique(rownames(counts)[duplicated(rownames(counts))]) ) {
		    print(gene)
		    indices = which(rownames(counts)==gene)
		    print(rownames(counts)[indices])
		    rownames(counts)[indices] = paste(rownames(counts)[indices],seq(length(indices)),sep='.')
		    print(rownames(counts)[indices])
		  }
		}
		rownames = data.frame(rownames,new=rownames(counts))

		normalized = sweep(counts,MARGIN=2,FUN="/",STATS=colSums(counts))
		# sum((counts[,1]/sum(counts[,1]))!=normalized[,1])
		normalized = normalized*10000
		normalized = log(normalized+1)
		print.message('Normalized');corner(normalized)

		t = start(append=T,split=T)
		duplicated.indices = duplicated(t(counts[genes.filter,])) | duplicated(t(counts[genes.filter,]),fromLast=T)
		if (sum(duplicated.indices) > 0) {
			print.message('\nRemoving all',sum(duplicated.indices),'duplicated cells\n')
			print(table(colnames(counts[,duplicated.indices])))
			counts = counts[,!duplicated.indices]
			normalized = normalized[,!duplicated.indices]
			dataset.labels = dataset.labels[!duplicated.indices]
			origins = origins[!duplicated.indices]
			experiments = experiments[!duplicated.indices]
			criteria[criteria==T][duplicated.indices] = FALSE
		}
		end(t)

		save(genes.filter,counts,normalized,dataset.labels,origins,experiments,criteria,file=cache)

		end(t)
	}

	counts = counts[genes.filter,]
	normalized = normalized[genes.filter,]

	environment$genes.filter = genes.filter
	environment$counts = counts
	environment$normalized = normalized
	environment$genes = rownames(normalized)
	environment$datasets = colnames(normalized)
	environment$dataset.labels = dataset.labels
	environment$origins = origins
	environment$experiments = experiments
	if (!exists('criteria')) criteria = NA
	environment$criteria = criteria
	environment$confounders = data.frame(nUMI=colSums(counts),nGenes=colSums(counts>0))
	environment$nsamples = ncol(counts)

	return(environment)
}

read.preclustered.datasets <- function (
	path = NA, # Search path for previous projects
	recursive = T, # Recursive path search
	rerun = F) { # Whether to rerun the reading process
	
	cache = file.path(environment$baseline.data.path,'preclustered.datasets.RData')

  if( !rerun & file.exists(cache) ) {
  	print.message('Loading precomputed')
  	load(cache)
  } else {

  	print.message('Computing')
		t = start(split=T)

		if (is.na(path)) path = dirname(environment$baseline.work.path)
		dataset = environment$datasets[1]
		merged.clustering = {}
		merged.original.clustering = {}
		merged.diff.exp = {}
		merged.HVG = {}
		merged.counts = {}
		merged.normalized = {}
		merged.dataset.labels = {}
		merged.origins = {}
		merged.experiments = {}
		merged.cluster.names = {}
		union.genes.filter = {}
		dataset.genes = NA
		sample.index = 1
		for (sample.index in seq(length(environment$datasets))) {
			dataset = environment$datasets[sample.index]
			origin = environment$origins[sample.index]
			experiment = environment$experiments[sample.index]
			data.files = list.files(path = file.path(path,dataset),pattern = 'clustering.RData',full.names = T,recursive = recursive)
			file.index = 1
			if (length(data.files) > 1) {
				for (row in seq(length(data.files))) print.message(row,dirname(dirname(data.files[row])))
				file.index = as.numeric(readline(prompt="Select clustering: "))
			}
			print.message('Loading',dirname(dirname(data.files[file.index])),'\n')
			load(data.files[file.index])
			if (length(merged.clustering) == 0) {
				min = 0
			} else {
				min = max(merged.clustering)
			}
			
			# names(clustering)
			# str(clustering$shuffled.membership)
			# str(clustering$shuffled.membership[[1]])
			# apply(clustering$shuffled.membership[[1]]$memberships,1,table)
			# membership = clustering$shuffled.membership[[2]]$memberships[2,]
			# limma.diff = {}
			# for( cluster in seq(length(unique(membership))) ) {
			# 	print.message('cluster =',cluster)
				
			# 	group = membership == cluster
			# 	group = factor(group)
				
			# 	diff.exp = getDE.limma( Y = normalized, group = group, filter = F )
			# 	diff.exp = diff.exp[order(diff.exp$logFC,decreasing=T),]
			# 	diff.exp = data.frame(gene=rownames(diff.exp),logFC=diff.exp$logFC,fold=exp(diff.exp$logFC),QValue=diff.exp$QValue,PValue=diff.exp$PValue,AveExpr=diff.exp$AveExpr)
			# 	limma.diff = rbind(limma.diff,data.frame(cluster=cluster,diff.exp[,c(1,3,4)]))
			# }
			# limma.diff = limma.diff[order(limma.diff$fold,decreasing=T),]
			# head(limma.diff)
			# lcdif[lcdif$gene=='Foxp3',]
			# lndif[lndif$gene=='Foxp3',]
			# res={}
			# for(i in seq(length(unique(lcdif$cluster)))){
			# 	for(j in seq(length(unique(lndif$cluster)))){
			# 		lc=lcdif[lcdif$cluster == i,]
			# 		ln=lndif[lndif$cluster == j,]
			# 		match = match(ln$gene,lc$gene)
			# 		head(lc[match,])
			# 		head(ln)
			# 		ccc=cor.test(lc$fold[match],ln$fold,method='spearman')
			# 		res = rbind(res,data.frame(ccc$estimate,ccc$p.value))
			# 	}
			# }
			# head(res)
			# res[order(res$ccc.estimate),]

			merged.clustering = c(merged.clustering,clustering$membership + min)
			merged.original.clustering = c(merged.original.clustering,clustering$membership)

			# index = 1
			# merged.shuffled.clustering[[]]
			# for (index in seq(length(clustering$shuffled.membership))) {
			# 	v = clustering$shuffled.membership[[index]]
			# 	shuffled.membership = v$memberships[2,]
			# 	if (length(shuffled.membership) == length(clustering$membership))
			# }

			load(file.path(dirname(data.files[file.index]),'main.all.diff.exp.RData'))
			if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] == T) {
				print.message('Converting from human to mouse genes')
				human.genes = toupper(limma.all$gene)
				mouse.gene.names = convertHumanGeneList(human.genes)
				if (length(dataset.genes) > 1 || !is.na(dataset.genes)) mouse.gene.names = mouse.gene.names[mouse.gene.names[,2] %in% dataset.genes,]
				limma.all$gene = mouse.gene.names[match(human.genes,mouse.gene.names[,1]),2]
				limma.all = limma.all[!is.na(limma.all$gene),]
				intersect.genes = intersect(unique(limma.all$gene),unique(merged.diff.exp$gene))
				limma.all = limma.all[limma.all$gene %in% intersect.genes,]
				merged.diff.exp = merged.diff.exp[merged.diff.exp$gene %in% intersect.genes,]
			}
			limma.all = cbind(dataset,origin,experiment,limma.all)
			merged.diff.exp = rbind(merged.diff.exp,limma.all)
			base.path = dirname(dirname(dirname(data.files[file.index])))
			load(file.path(base.path,'data/HVG.RData'))
			if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] == T) {
				human.genes = toupper(HVG)
				mouse.gene.names = convertHumanGeneList(human.genes)
				if (length(dataset.genes) > 1 || !is.na(dataset.genes)) mouse.gene.names = mouse.gene.names[mouse.gene.names[,2] %in% dataset.genes,]
				HVG = unique(mouse.gene.names[,2])
			}
			merged.HVG = unique(c(merged.HVG,HVG))
			load(file.path(base.path,'data/data.RData'))
			colnames(counts) = colnames(normalized) = rep(dataset,ncol(normalized))
			if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] == T) {
				human.genes = toupper(rownames(counts))
				mouse.gene.names = convertHumanGeneList(human.genes)
				if (length(dataset.genes) > 1 || !is.na(dataset.genes)) mouse.gene.names = mouse.gene.names[mouse.gene.names[,2] %in% dataset.genes,]
				
				rownames(counts) = mouse.gene.names[match(human.genes,mouse.gene.names[,1]),2]
				counts = counts[!is.na(rownames(counts)),]
				intersect.genes = intersect(rownames(counts),rownames(merged.counts))
				counts = counts[rownames(counts) %in% intersect.genes,]
				merged.counts = merged.counts[rownames(merged.counts) %in% intersect.genes,]
				counts = counts[match(rownames(merged.counts),rownames(counts)),]
				
				rownames(normalized) = mouse.gene.names[match(human.genes,mouse.gene.names[,1]),2]
				normalized = normalized[!is.na(rownames(normalized)),]
				intersect.genes = intersect(rownames(normalized),rownames(merged.normalized))
				normalized = normalized[rownames(normalized) %in% intersect.genes,]
				merged.normalized = merged.normalized[rownames(merged.normalized) %in% intersect.genes,]
				normalized = normalized[match(rownames(merged.normalized),rownames(normalized)),]
			}
			if (length(merged.counts)>1 && sum(rownames(merged.counts) != rownames(counts))>0) stop('Feature genes mismatch - need to correct dataset binding matching')
			merged.counts = cbind(merged.counts,counts)
			merged.normalized = cbind(merged.normalized,normalized)
			merged.dataset.labels = c(merged.dataset.labels,rep(dataset,ncol(normalized)))
			merged.origins = c(merged.origins,rep(origin,ncol(normalized)))
			merged.experiments = c(merged.experiments,rep(experiment,ncol(normalized)))
			load(file.path(dirname(data.files[file.index]),'cluster.names.RData'))
			merged.cluster.names = c(merged.cluster.names,cluster.names)
			if (!is.null(environment$convert.to.mouse.gene.symbols) && environment$convert.to.mouse.gene.symbols[sample.index] == T) {
				human.genes = names(genes.filter)
				mouse.gene.names = convertHumanGeneList(human.genes)
				if (length(dataset.genes) > 1 || !is.na(dataset.genes)) mouse.gene.names = mouse.gene.names[mouse.gene.names[,2] %in% dataset.genes,]
				names(genes.filter) = mouse.gene.names[match(human.genes,mouse.gene.names[,1]),2]
				genes.filter = genes.filter[!is.na(names(genes.filter))]
				intersect.genes = intersect(names(genes.filter),names(union.genes.filter))
				genes.filter = genes.filter[names(genes.filter) %in% intersect.genes]
				union.genes.filter = union.genes.filter[names(union.genes.filter) %in% intersect.genes]
				genes.filter = genes.filter[match(names(union.genes.filter),names(genes.filter))]
			}
			if (length(union.genes.filter) == 0) {
				union.genes.filter = genes.filter
			} else {
				union.genes.filter = union.genes.filter | genes.filter
			}

			dataset.genes = rownames(merged.counts)[apply(merged.counts,1,sd)>0]
		}

		# merged.HVG[!merged.HVG%in%rownames(normalized[genes.filter,])]
		# sum(genes.filter)
		merged.HVG = merged.HVG[merged.HVG%in%rownames(counts)]
		union.genes.filter = union.genes.filter[names(union.genes.filter)%in%rownames(counts)]

		print.message('dim(merged.normalized)');print(dim(merged.normalized))
		print(table(merged.dataset.labels))
		print(table(merged.origins))
		print(table(merged.experiments))
		print(table(merged.clustering))

		genes.filter = union.genes.filter
		counts = merged.counts
		normalized = merged.normalized
		dataset.labels = merged.dataset.labels
		origins = merged.origins
		experiments = merged.experiments
		HVG = merged.HVG
		clustering = merged.clustering
		cluster.names = merged.cluster.names
		print.message('Saving')
		save(genes.filter,counts,normalized,dataset.labels,origins,experiments,HVG,clustering,merged.diff.exp,merged.original.clustering,cluster.names,file=cache)

		end(t)
	}
	counts = counts[genes.filter,]
	normalized = normalized[genes.filter,]
	HVG = HVG[HVG %in% rownames(normalized)]

	environment$counts = counts
	environment$normalized = normalized
	environment$genes = rownames(normalized)
	environment$datasets = colnames(environment$normalized)
	environment$dataset.labels = dataset.labels
	environment$origins = origins
	environment$experiments = experiments
	environment$confounders = data.frame(nUMI=colSums(counts),nGenes=colSums(counts>0))
	environment$nsamples = ncol(counts)
	environment$HVG = HVG
	environment$clustering = {}
	environment$clustering$membership = clustering
	environment$original.clustering = merged.original.clustering
	environment$diff.exp = merged.diff.exp
	environment$cluster.names = apply(cbind(environment$datasets,cluster.names),1,function(v) paste(v,collapse='_'))

	return(environment)
}

add.confounder.variables <- function (...) {
	
	environment$confounders = data.frame(environment$confounders,data.frame(...));print(head(environment$confounders))
	print(head(environment$confounders))
	return(environment)
}

get.variable.genes <- function (
	min.mean = 0.05, # Minimum mean expression per gene
	min.frac.cells = 0, # Minimum fraction of cells expressing each gene
	min.dispersion.scaled = 1, # Minimum dispersion value
	rerun = F) { # Whether to rerun

	cache = file.path(environment$baseline.data.path,'HVG.RData')

  if( !rerun & file.exists(cache) ) {
  	print.message('Loading precomputed')
  	load(cache)
  } else {
  	print.message('Computing')
  	t = start()

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

			library(ggplot2,quietly=T)
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

	cat('# highly variable genes = ',length(environment$HVG),'\n',sep='')

	return(environment)
}

nUMIs <- function () { return(colSums(environment$counts)) }

nGenes <- function () { return(colSums(environment$counts>0)) }

ribosomal.score <- function (control = T,knn=10) {
	t = start()
	genes = get.ribo.genes(environment$genes)
	print.message('Using genes:');print(genes)
	if (control) {
		score = controlled.mean.score(genes,knn)
	} else {
		score = colMeans(environment$normalized[genes,])
	}	
	end(t)
	return(score)
}

get.ribo.genes <- function (genes) {
	return (genes[c(grep('^Rpl',genes,ignore.case = T),grep('^Rps',genes,ignore.case = T))])
}

mitochondrial.score <- function (control = F,knn=10) {
	t = start()
	genes = get.mito.genes(environment$genes)
	print.message('Using genes:');print(genes)
	if (control) {
		score = controlled.mean.score(genes,knn)
	} else {
		score = colMeans(environment$normalized[genes,])
	}	
	end(t)
	return(score)
}

get.mito.genes <- function (genes) {
	return (genes[grep('^Mt-',genes,ignore.case = T)])
}

capwords <- function(s, strict = FALSE) {
	s = tolower(s)
  cap <- function(s) paste(toupper(substring(s, 1, 1)), {s <- substring(s, 2); if(strict) tolower(s) else s}, sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

cell.cycle.score <- function (knn = 10,cc.genes.path = NA) {

	if (is.na(cc.genes.path)) cc.genes.path = file.path(environment$data.path,"regev_lab_cell_cycle_genes.txt")

	t = start()
	cc.genes = capwords(readLines(cc.genes.path))
	s.genes = cc.genes[1:43];s.genes = s.genes[s.genes%in%environment$genes];print(s.genes)
	g2m.genes = cc.genes[44:98];g2m.genes = g2m.genes[g2m.genes%in%environment$genes];print(g2m.genes)

	s.score = controlled.mean.score (s.genes,knn)
	g2m.score = controlled.mean.score (g2m.genes,knn)
	cell.cycle.score = controlled.mean.score (c(s.genes,g2m.genes),knn)

	print.message('# s.score > 0:',sum(s.score>0),'fraction',sum(s.score>0)/length(s.score))
	print.message('# g2m.score > 0:',sum(g2m.score>0),'fraction',sum(g2m.score>0)/length(g2m.score))
	end(t)

	return(data.frame(s.score=s.score,g2m.score=g2m.score,cell.cycle.score=cell.cycle.score))
}

controlled.mean.score <- function (genes,knn = 10,exclude.missing.genes = T,constrain.cell.universe = NA) {
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

		background.genes = background.genes (foreground.genes = genes,knn)

		return(colMeans(environment$normalized[genes,constrain.cell.universe])-colMeans(environment$normalized[background.genes,constrain.cell.universe]))
	} else {
		return(colMeans(environment$normalized[genes,constrain.cell.universe]))
	}
}

get.technically.similar.genes <- function (knn = 10) {
	
	t = start()
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

background.genes <- function (foreground.genes,knn) {
	
	t = start()
	foreground.genes = foreground.genes[foreground.genes%in%environment$genes]
	technically.similar.genes = get.technically.similar.genes (knn)
	knns = technically.similar.genes$knns
	technical.variables = technically.similar.genes$technical.variables

	background.genes = unique(setdiff(as.vector(knns[foreground.genes,]),foreground.genes))
	print.message('Head foreground.genes technical.variables');print(head(technical.variables[foreground.genes,]))
	print.message('Head background.genes technical.variables');print(head(technical.variables[background.genes,]))
	print.message('background.genes');print(background.genes)
	end(t)
	return(background.genes)
}

PCA <- function (
	regress = NA, # Gene signature activation scores to regress
	groups = NA, # Experimental design annotation to guide dataset-specific regression
	nShuffleRuns = 10, # Number of shuffled analyses
	threshold = 0.1, # FDR threshold
	maxPCs = 100, # Maximum number of possible PCs
	label = NA, # Optional analyses label folder
	mem = '2GB', # HPC memory
	time = '0:10:00', # HPC time
	rerun = F, # Whether to rerun
	clear.previously.calculated.clustering = T) { # Whether to clear previous clustering analysis

	if (length(regress)>1 || !is.na(regress)) {
		config = paste(colnames(regress),collapse='+')
	} else {
		config = 'not.regressed'
	}

	if(length(groups)==1 && is.na(groups)) groups = rep(1,environment$nsamples)
	if(nShuffleRuns != 10) config = paste(config,nShuffleRuns,sep='.')
	if(threshold != 0.1) config = paste(config,threshold,sep='.')
	if(maxPCs != 100) config = paste(config,maxPCs,sep='.')
	if(!is.na(label)) config = paste(config,label,sep='.')

	environment$work.path = file.path(environment$baseline.work.path,config)
	environment$PCA = environment$Rotation = environment$PCA.path = NA
	if (clear.previously.calculated.clustering) environment$clustering = environment$seurat.cluster.association = NA
	dir.create(file.path(environment$work.path,'tracking'), showWarnings = F, recursive = T, mode = "700")
	print.message('Transitioning to',config,'folder')
	setwd(environment$work.path)
	environment$res.data.path = file.path(environment$work.path,'data')
	shuffled.PCA.data.path = file.path(environment$res.data.path,'shuffled.PCA')

	dir.create(environment$res.data.path, showWarnings = F, recursive = T, mode = "700")
	
	cache = file.path(environment$res.data.path,paste(config,'PCA.RData',sep='.'))

  if( !rerun && file.exists(cache) ) {
  	print.message('Loading precomputed')
  	load(cache)
  } else {
  	print.message('Computing')
  	t = start()

  	if (clear.previously.calculated.clustering) unlink(file.path(environment$res.data.path,'clustering'),recursive=T,force=T)

  	unlink(shuffled.PCA.data.path,recursive=T,force = T)
		dir.create(shuffled.PCA.data.path, showWarnings = F, recursive = T, mode = "700")
		
		raw.data = environment$counts
		data = environment$normalized[environment$HVG,]
		print.message('Dim');print(dim(data))

		if (length(regress)>1 || !is.na(regress)) {
			corrected = regress.covariates (regress,data,groups,rerun,save=T)
			print.message('Regressed matrix');corner(corrected)
			data = corrected
		}

		n <- min(maxPCs, ncol(data))
		m <- nrow(data)
		ndf <- n - 1

		get.shuffled.var <- function (rep) {

			data.perm <- apply(data, 2, sample, replace = FALSE)
			pca.perm = prcomp(t(data.perm), retx = TRUE, center = T, scale. = T)
			var.perm = pca.perm$sdev[1:ndf]^2/sum(pca.perm$sdev[1:ndf]^2)
			save (pca.perm,var.perm,file=file.path(shuffled.PCA.data.path,paste('shuffled.PCA.rep',rep,'RData',sep='.')))
			return (var.perm)
		}

		suppressWarnings(library(rslurm,quietly=T))
		sopt <- list(mem = mem, time = time, share = TRUE)
		sjob <- slurm_apply(get.shuffled.var, data.frame(rep=seq(nShuffleRuns)), add_objects = c('shuffled.PCA.data.path','data','ndf'), pkgs = NULL, nodes = nShuffleRuns, cpus_per_node = 1, submit = TRUE, slurm_options = sopt)

		pc.time = Sys.time()
		pca = prcomp(t(data), retx = TRUE, center = T, scale. = T)
		print.message('Single PCA run time:')
		print(Sys.time()-pc.time)
		var <- pca$sdev[1:ndf]^2/sum(pca$sdev[1:ndf]^2)
		print.message('Real PCA Var');print(head(var,10))

		var.perm = get_slurm_out(sjob);dim(var.perm)

		end(t)
		start(append=T,split=T)

		if (length(var.perm) == 0 || nrow(var.perm) < nShuffleRuns) {
			print.message('JOB ERROR: Not enough shuffled results:',nrow(var.perm),'<',nShuffleRuns,'\nCHECK FOR FAILED JOBS\n\n\n')
			terminate = readline(prompt="Terminate? (y/n) ")
			if (terminate != 'n') {
				end(t)
				tryCatch({cleanup_files(sjob);cleanup_files(sjob);cleanup_files(sjob)},error=function(v) v)
				return(environment)
			}
		}
		end(t)
		start(append=T)
		tryCatch({cleanup_files(sjob);cleanup_files(sjob);cleanup_files(sjob)},error=function(v) v)
		print.message('Shuffled PCA Var');corner(var.perm,10)

		p <- rep(1, n)
		for (i in 1:ndf) {
		  p[i] <- mean(var.perm[,i] >= var[i])
		}
		p
		for (i in 2:ndf) {
		  p[i] <- max(p[(i - 1)], p[i])
		}
		nPCs <- sum(p <= threshold,na.rm=T)
		print.message('nPCs',nPCs)

		PCA = t(pca$x)[seq(nPCs),]
		print.message('PCA');corner(PCA)		
		Rotation = t(pca$rotation[,seq(nPCs)])
		print.message('Rotation');corner(Rotation)

		save(PCA,Rotation,file=cache)

		end(t)
	}

	environment$PCA = PCA
	environment$Rotation = Rotation
	environment$PCA.path = cache

	cat('# PCs = ',nrow(environment$PCA),'\n',sep='')

	return(environment)	
}

regress.covariates <- function (regress,data,groups,rerun = F,save = F) {

	cache = file.path(environment$res.data.path,paste(paste(colnames(regress),collapse='+'),'HVG.regressed.covariates.RData',sep='_'))

  if( !rerun && file.exists(cache) ) {
  	print.message('Loading precomputed')
  	load(cache)
  } else {
  	print.message('Computing')
  	t = start()

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

get.robust.cluster.similarity <- function (similarity,min.sd = qnorm(.95),max.q.val = 0.01,rerun = F) {

	cache = file.path(environment$res.data.path,'filtered.cluster.similarity.RData')
	if( !rerun && file.exists(cache) ) {
  	print.message('Loading precomputed similarity')
  	load(cache)
  } else {
  	print.message('Computing')

		library("dplyr");library("ggpubr");library(ggplot2)

		match.significance.stats = {}
		origin = unique(similarity$origin1)[2]
		origins = sort(unique(c(as.vector(similarity$origin1),as.vector(similarity$origin2))))
		for (origin in origins) {
			filter = similarity$origin1 == origin & similarity$origin2 == origin & similarity$experiment1 != similarity$experiment2
			filtered.similarity = similarity[filter,]
			cor.val = filtered.similarity$similarity
			shapiro.test.p.value = tryCatch({ shapiro.test(cor.val)$p.value }, error = function(v) return(NA))
			sd.dist = (cor.val - mean(cor.val)) / sd(cor.val)
			match.significance.stats = rbind(match.significance.stats,data.frame(origin,cor.val,sd.dist,shapiro.test.p.value,filtered.similarity))
		}

		name = paste('get.robust.cluster.similarity',sep='_')
		work.path = file.path(environment$work.path,name)
		if (file.exists(work.path)) {
			new.dir = file.path(environment$work.path,paste(name,format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),sep='---'))
			file.rename(work.path,new.dir)
		}
		dir.create(work.path,showWarnings = T,recursive = T)
		
		plot.data = data.frame(correlation = match.significance.stats$cor.val,SD = match.significance.stats$sd.dist,origin = match.significance.stats$origin)
		pdf(file.path(work.path,'cluster.matching.similarity.histogram.pdf'))
		print(ggplot(plot.data, aes(x=correlation, fill = origin)) + geom_density(alpha = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Density') + theme_classic(base_size=25))
		print(ggplot(plot.data, aes(x=SD, fill = origin)) + geom_density(alpha = 0.5) + geom_vline(xintercept = min.sd) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Density') + theme_classic(base_size=25))
		print(ggplot(plot.data, aes(x=correlation)) + geom_density(alpha = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Density') + theme_classic(base_size=25))
		print(ggplot(plot.data, aes(x=SD)) + geom_density(alpha = 0.5) + geom_vline(xintercept = 2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Density') + theme_classic(base_size=25))
		dev.off()

		sum(match.significance.stats$sd.dist >= min.sd)/length(match.significance.stats$sd.dist)
		criteria = match.significance.stats$sd.dist >= min.sd & match.significance.stats$Q.val.1 <= max.q.val
		print.message('Found',sum(criteria),'robust populations')
		robust.cluster.similarity = match.significance.stats[criteria,]
		robust.cluster.similarity = robust.cluster.similarity[order(robust.cluster.similarity$origin,robust.cluster.similarity$similarity.1,decreasing=T),]
		summary(robust.cluster.similarity)
		print(robust.cluster.similarity[1:min(nrow(robust.cluster.similarity),5),])
		write.csv(robust.cluster.similarity,file = file.path(work.path,'robust.cluster.similarity.csv'))

		robust.clusters = unique(c(robust.cluster.similarity$cluster1,robust.cluster.similarity$cluster2))
		nclusters.overall = length(unique(c(similarity$cluster1,similarity$cluster2)))
		cat('found',length(robust.clusters),'/',nclusters.overall,'robust.clusters (',round(length(robust.clusters)/nclusters.overall,2),'%)\n')
		filtered.cluster.similarity = similarity[similarity$cluster1 %in% robust.clusters & similarity$cluster2 %in% robust.clusters,]
		dim(filtered.cluster.similarity)
		summary(filtered.cluster.similarity)

		print(filtered.cluster.similarity[1:min(nrow(filtered.cluster.similarity),5),])
		write.csv(filtered.cluster.similarity,file = file.path(work.path,'filtered.cluster.similarity.csv'))

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = filtered.cluster.similarity$name1,cluster2 = filtered.cluster.similarity$name2,coef = filtered.cluster.similarity$similarity)
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.Cor.FC.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = filtered.cluster.similarity$name1,cluster2 = filtered.cluster.similarity$name2,coef = filtered.cluster.similarity$similarity.1)
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.Cor.means.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = filtered.cluster.similarity$name1,cluster2 = filtered.cluster.similarity$name2,coef = -scale(filtered.cluster.similarity$ocldist.FC))
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.ocldist.FC.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = filtered.cluster.similarity$name1,cluster2 = filtered.cluster.similarity$name2,coef = -scale(filtered.cluster.similarity$ocldist))
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.ocldist.means.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		save(filtered.cluster.similarity,file = cache)
	}

	return(filtered.cluster.similarity)
}

compare.cluster.similarity <- function (diff.exp.file = 'main.datasets.diff.exp.RData',cluster.similarity.function = pearson.correlation,label = 'pearson',rerun = F) {

	cache = file.path(environment$res.data.path,paste(label,'cluster.similarity.RData',sep='.'))
	if( !rerun && file.exists(cache) ) {
  	print.message('Loading precomputed similarity')
  	load(cache)
  } else {
  	print.message('Computing')
  	t = start(split = F)

  	name = paste('compare.cluster.similarity',label,sep='_')
		work.path = file.path(environment$work.path,name)
		if (file.exists(work.path)) {
			new.dir = file.path(environment$work.path,paste(name,format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),sep='---'))
			file.rename(work.path,new.dir)
		}
		dir.create(work.path,showWarnings = T,recursive = T)

  	load(file.path(environment$work.path,'data',diff.exp.file))
		final.diff = limma.all
		head(final.diff)
		configs = apply(final.diff,1,function(v) paste(c(v[1],v[2]),collapse = '_'))
		membership = environment$clustering$membership
		original.membership = environment$original.clustering
		nclusters = length(unique(membership))
		clusters = sort(unique(membership))
		map = data.frame(
			t(
				sapply(clusters,function(c) c(
					unique(environment$cluster.names[membership==c][1]),
					unique(environment$origins[membership==c][1]),
					unique(environment$experiments[membership==c][1]),
					unique(environment$datasets[membership==c][1]),
					unique(original.membership[membership==c][1]),
					c))))
		colnames(map) = c('name','origin','experiments','sample','original.membership','membership')
		head(map)

		cluster.descriptors = {}
		for (cluster.index in seq(nclusters)) {
			cluster = clusters[cluster.index]
			indices = membership==cluster
			diff.exp.indices = configs == paste(environment$datasets[indices][1],membership[indices][1],sep = '_ ') | configs == paste(environment$datasets[indices][1],membership[indices][1],sep = '_')
			diff = final.diff[diff.exp.indices,];rownames(diff) = diff$gene
			if (is.null(cluster.descriptors)) {
				cluster.descriptors = data.frame(diff$fold)
				rownames(cluster.descriptors) = diff$gene
			} else {
				cluster.descriptors = cbind(cluster.descriptors,diff$fold[match(rownames(cluster.descriptors),diff$gene)])
			}
		}
		colnames(cluster.descriptors) = map[match(seq(nclusters),map$membership),1]
		head(cluster.descriptors)

		pdf(file.path(environment$work.path,paste('hclust.cluster.descriptors.pdf',sep='_')),width=10,height=10)
		hc.dist = hclust(dist(t(cluster.descriptors)))
		plot(hc.dist)
		dev.off()

		# Cluster similarity
		similarity = {}
		for (cluster.index1 in seq(nclusters-1)) {
			cluster1 = clusters[cluster.index1]
			indices = membership==cluster1
			diff.exp.indices1 = configs == paste(environment$datasets[indices][1],membership[indices][1],sep = '_ ') | configs == paste(environment$datasets[indices][1],membership[indices][1],sep = '_')
			diff1 = final.diff[diff.exp.indices1,];rownames(diff1) = diff1$gene
			measurements1 = environment$normalized[,membership==cluster1]
			measurements1.means = rowMeans(measurements1)
			for (cluster.index2 in (cluster.index1+1):nclusters) {
				cluster2 = clusters[cluster.index2]
				indices = membership==cluster2
				diff.exp.indices2 = configs == paste(environment$datasets[indices][1],membership[indices][1],sep = '_ ') | configs == paste(environment$datasets[indices][1],membership[indices][1],sep = '_')
				diff2 = final.diff[diff.exp.indices2,];rownames(diff2) = diff2$gene
				matches = match(diff1$gene,diff2$gene)
				diff2 = diff2[matches,]
				if (any(rownames(diff1)!=rownames(diff2))) stop ('\n\n\n\nGene mismatch\n\n\n\n')
				measurements2 = environment$normalized[,membership==cluster2]
				measurements2.means = rowMeans(measurements2)
				ocldist = as.vector(dist(rbind(measurements1.means,measurements2.means)))
				ocldist.FC = as.vector(dist(rbind(diff1$fold,diff2$fold)))
				over1 = diff1[diff1$QValue<=0.05 & diff1$fold>=1,]
				over2 = diff2[diff2$QValue<=0.05 & diff2$fold>=1,]
				intersect = intersect(over1$gene,over2$gene);sort(intersect)
				union = union(over1$gene,over2$gene)
				jaccard = length(intersect)/length(union)
				over1 = diff1[diff1$QValue<=0.05 & diff1$fold>=1.5,]
				over2 = diff2[diff2$QValue<=0.05 & diff2$fold>=1.5,]
				intersect = intersect(over1$gene,over2$gene);sort(intersect)
				union = union(over1$gene,over2$gene)
				jaccard2 = length(intersect)/length(union)
				diff1.compact = diff1[diff1$fold>1.05 | diff1$fold<(1/1.05),]
				diff2.compact = diff2[diff2$fold>1.05 | diff2$fold<(1/1.05),]
				diff.compact.genes = union(diff1.compact$gene,diff2.compact$gene)
				diff1.compact.strict = diff1[diff1$fold>1.15 | diff1$fold<(1/1.15),]
				diff2.compact.strict = diff2[diff2$fold>1.15 | diff2$fold<(1/1.15),]
				diff.compact.strict.genes = union(diff1.compact.strict$gene,diff2.compact.strict$gene)
				similarity = rbind(similarity,data.frame(cluster1,cluster2,cluster.similarity.function (diff1$fold,diff2$fold),cluster.similarity.function (measurements1.means,measurements2.means),ocldist,ocldist.FC,jaccard,jaccard2))
			}
		}
		head(similarity)
		dim(similarity)

		# Cluster mapping
		similarity$Q.val = p.adjust(similarity$significance,method='BH')
		similarity$Q.val.1 = p.adjust(similarity$significance.1,method='BH')

		similarity = cbind(map[match(similarity$cluster1,map$membership),1:5],map[match(similarity$cluster2,map$membership),1:5],similarity)
		head(similarity)
		colnames(similarity)[1:10] = c('name1','origin1','experiment1','sample1','original.membership1','name2','origin2','experiment2','sample2','original.membership2')
		similarity = similarity[order(similarity$similarity,decreasing=T),]
		write.csv(similarity,row.names=F,file=file.path(work.path,paste(label,'cluster.similarity.csv',sep='_')))

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = similarity$name1,cluster2 = similarity$name2,coef = similarity$similarity)
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.Cor.FC.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = similarity$name1,cluster2 = similarity$name2,coef = similarity$similarity.1)
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.Cor.means.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = similarity$name1,cluster2 = similarity$name2,coef = -scale(similarity$ocldist.FC))
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.ocldist.FC.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

		library(reshape2)
		similarity.summary.df = data.frame(cluster1 = similarity$name1,cluster2 = similarity$name2,coef = -scale(similarity$ocldist))
		mirror = similarity.summary.df; cluster1 = mirror$cluster1; cluster2 = mirror$cluster2; mirror$cluster1 = cluster2; mirror$cluster2 = cluster1
		similarity.summary.df = rbind(similarity.summary.df,mirror)
		head(similarity.summary.df)
		similarity.matrix = acast(similarity.summary.df, cluster1~cluster2, value.var="coef")
		similarity.matrix[1:10,1:10]
		hc.dist = hclust(as.dist(1-similarity.matrix))
		library(RColorBrewer);library(gplots)
		pdf(file.path(work.path,paste('filtered.cluster.similarity.heatmap.ocldist.means.pdf',sep='_')),width=20,height=20)
		colors = rev(brewer.pal(5, "PuOr"));color.palette = colorRampPalette(colors)
		print(heatmap.2(similarity.matrix, col=color.palette, key=T, cexRow=1, cexCol=1, scale = 'none', density.info="none", trace="none",Rowv = as.dendrogram(hc.dist), Colv = as.dendrogram(hc.dist),dendrogram = 'both',cellnote=round(similarity.matrix,1),notecol='white',main='Cor',margins = c(10,10)))
		dev.off()

  	save(similarity,map,file=cache)
		end(t)
  }

  return(list(map=map,similarity=similarity))
}

visualize.cluster.cors.heatmaps <- function (work.path,similarity) {

	name = 'cross.sample.similarities'
	work.path = file.path(environment$work.path,name)
	if (file.exists(work.path)) {
		new.dir = file.path(environment$work.path,paste(name,format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),sep='---'))
		file.rename(work.path,new.dir)
	}
	dir.create(work.path,showWarnings = T,recursive = T)

	library(reshape2)

	samples.names = unique(c(as.vector(similarity$sample1),as.vector(similarity$sample2)))
	samples.names = sort(samples.names)
	sample.name1 = samples.names[3]
	sample.name2 = samples.names[4]
	library(gplots,quietly=T);library(RColorBrewer)
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

		suppressWarnings(library(rslurm,quietly=T))
		sopt <- list(mem = '8GB', time = '1:00:00', share = TRUE)
		sjob <- slurm_apply(tSNE, params, nodes = nrow(params), cpus_per_node = 1, add_objects = c('data.path','tSNEs.dir'), submit = TRUE, slurm_options = sopt)
		end()
	}

	return(sjob)
}

getDE.limma <- function( Y, group,filter = T ){
  require(limma)
  design = model.matrix( ~group )
  my.lm <- limma::lmFit(Y, design = design )
  my.lm <- eBayes(my.lm)
  topTable = topTable(my.lm,number=Inf, coef = seq(2,ncol(design)))
  if( filter ) topTable = topTable[topTable$"adj.P.Val" < 0.1,]
  topTable = topTable[order(topTable$"adj.P.Val",decreasing = F),]
  colnames(topTable)[colnames(topTable)=='P.Value'] = 'PValue'
  colnames(topTable)[colnames(topTable)=='adj.P.Val'] = 'QValue'

  return(topTable)
}

run.diff.expression <- function (clustering,min.fold,quantile,label,robust = F,rerun = F,contrast = 'all',contrast.groups = NA) {

	get.diff.exp.stats <- function (id) {
		t=Sys.time()
		if (is.na(id)) {
			label = 'real'
			clustering = as.vector(membership) 
		} else {
			label = 'shuffled'
			clustering = as.vector(unlist(shuffled.membership[[id]]))
		}

		fold <- function (v,group) { exp(mean(v[group==T])-mean(v[group==F])) }

		cat(label,id,'\n')

		stats = {}
		for( cluster in sort(unique(clustering)) ) {

			cat('cluster',cluster,'\n')
			tt=Sys.time()
			group = factor(clustering == cluster)

			stats = rbind(stats,data.frame(cluster=cluster,gene=rownames(matrix),fold=apply(matrix,1,fold,group)))

			print(Sys.time()-tt)
		}

		stats = cbind(label,stats)

		print(Sys.time()-t)
		return(stats)
	}

	summarize.diff.exp.stats <- function (job.portion,min.fold,quantile) {

		real.stats = stats[stats$label=='real',];print(head(real.stats))
		shuffled.stats = stats[stats$label=='shuffled',];print(head(shuffled.stats))
		genes = genes[job.portions==job.portion]
		results = {}
		empirical.diff = {}
		for (gene in genes) {
			real = real.stats[real.stats$gene==gene,]
			real = real[real$fold>=min.fold,] 
			if (nrow(real)==0) next
			shuffled = shuffled.stats$fold[shuffled.stats$gene==gene]
			# quantile(shuffled,seq(0.01,1,0.01));quantile(shuffled,0.95);real
			significant = real[real$fold>=quantile(shuffled,quantile),]
			empirical.diff = rbind(empirical.diff,significant)
		}
		rownames(empirical.diff) = NULL
		empirical.diff = empirical.diff[order(empirical.diff$fold,decreasing=T),]
		print(head(empirical.diff))
		return(empirical.diff)
	}

	cache = file.path(environment$res.data.path,paste(label,contrast,'diff.exp.RData',sep='.'))

  if( !rerun && file.exists(cache) ) {
  	print.message('Loading precomputed')
  	load(cache)
  } else {
  	print.message('Computing')
  	t = start(split = T)
  	membership = as.vector(clustering$membership)

  	empirical.diff = NA

  	diff.exp.dir = file.path(environment$work.path,'diff.exp',label,contrast)
		unlink(diff.exp.dir,recursive=T,force=T)
		dir.create(diff.exp.dir,recursive = T)
		limma.diff = {}

		if (contrast == 'all') {
			contrast.groups = rep(1,length(membership))
		} else if (contrast == 'datasets') {
			contrast.groups = environment$datasets
		}

		print.message('contrast =',contrast)

		print(table(contrast.groups))
		contrast.group = unique(contrast.groups)[1]
		table(membership,contrast.groups)
		for (contrast.group in unique(contrast.groups)) {
			print.message('contrast.group =',contrast.group)
			contrast.group.indices = contrast.groups == contrast.group
			sum(contrast.group.indices)
			contrast.group.clusters = sort(unique(membership[contrast.group.indices]))
			cluster = contrast.group.clusters[1]
			for( cluster in contrast.group.clusters ) {
				print.message('cluster =',cluster)
				
				group = membership[contrast.group.indices] == cluster
				group = factor(group)
				
				diff.exp = getDE.limma( Y = environment$normalized[,contrast.group.indices], group = group, filter = F )
				diff.exp = diff.exp[order(diff.exp$logFC,decreasing=T),]
				diff.exp = data.frame(gene=rownames(diff.exp),logFC=diff.exp$logFC,fold=exp(diff.exp$logFC),QValue=diff.exp$QValue,PValue=diff.exp$PValue,AveExpr=diff.exp$AveExpr)
				write.csv(diff.exp,file=file.path(diff.exp.dir,paste('cluster',cluster,'csv',sep='.')))
				markers.diff = diff.exp[diff.exp$gene%in%environment$marker.genes,]
				write.csv(markers.diff,file=file.path(diff.exp.dir,paste('markers.cluster',cluster,'csv',sep='.')))
				limma.diff = rbind(limma.diff,data.frame(contrast.group=contrast.group,cluster=cluster,diff.exp[,c(1,3,4)]))
			}
		}
		rownames(limma.diff) = NULL
		limma.all = limma.diff
		limma.diff = limma.diff[limma.diff$QValue<=(1-quantile) & limma.diff$fold>=min.fold,]
		final.diff = limma.diff = limma.diff[order(limma.diff$fold,decreasing=T),]
		print.message('head(limma.diff):');print(head(limma.diff))
		write.csv(limma.diff,file=file.path(diff.exp.dir,'limma.diff.csv'))
		write.csv(final.diff,file=file.path(diff.exp.dir,'final.diff.csv'))

		save(final.diff,limma.diff,empirical.diff,limma.all,file=cache)
		end(t)
	}

	return(final.diff)
}

plot.PCA <- function (quantile,order) {

	work.path = environment$work.path
	PCA = environment$PCA
	Rotation = environment$Rotation
	cluster.names = environment$cluster.names
	dataset = environment$dataset.labels
	marker.genes = environment$marker.genes

	library(ggplot2,quietly=T);library(GGally,quietly=T);library(ggrepel,quietly=T)

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

plot.cluster.stats <- function (membership,label = NA,order = NA) {

	library(ggplot2,quietly=T)
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
	library(ggplot2,quietly=T)
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

plot.tSNE <- function (tSNE.job,perplexity,max_iter,membership = NA) {
	
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
		
		library(ggplot2,quietly=T);library(GGally,quietly=T);library(ggrepel,quietly=T)
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

	library(gplots,quietly=T)
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
			library(xlsx)
			write.xlsx(as.matrix(others),file=save)
		}
		},error = function(e) print(e))
	}

	rownames(others)
}

plot.simple.heatmap <- function (name,path = NA,markers,membership = NA,normalized = NA,order = NA,width = 5, height = 5,scale='row',RowSideColors = NA,counts = F,filter.diff.exp=F,cellnote=F,key=F,save=NA,sort.rows=T,sort.cols=T,Colv=F,Rowv=F,dendrogram = 'none') {
	
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

plot.violin <- function (genes,types,fore1exp1,fore2exp1,fore1exp2,fore2exp2,back1exp1,back2exp1,back1exp2,back2exp2,path,height = 5,width = 5,scale = T,palette = "Greys",separate.background = F) { #Pastel2

	library(ggplot2)
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

plot.heatmaps <- function (diff.exp,membership,order = NA,nTopRanked = 10,label = NA) {

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

summarize <- function (perplexity = seq(10,30,10),max_iter = 10000,rerun = F,robust = F,order = NA,contrast = 'all',min.fold = 1.5,quantile = 0.95) {

	cluster.size = table(environment$cluster.names)
	if (is.na(order)) order = names(cluster.size)[order(cluster.size,decreasing=T)]
	tSNE.job = run.tSNE (perplexity,max_iter,rerun)
	plot.PCA (quantile = 0.05,order)
	plot.cluster.stats (membership = environment$cluster.names,order = order)
	if (length(environment$seurat.cluster.association)>1) tryCatch({ plot.cluster.stats (membership = environment$seurat.cluster.association,label = 'Seurat',order = order) },error=function(v) v)

  final.diff = run.diff.expression (clustering = environment$clustering,min.fold,quantile,label='main',robust = robust,rerun = rerun,contrast = contrast)
	
	order = sort(unique(environment$cluster.names))
	plot.heatmaps (diff.exp = final.diff,membership = environment$cluster.names,order = order)
	if (length(environment$seurat.cluster.association)>1) tryCatch({ plot.heatmaps (diff.exp = final.diff,membership = environment$seurat.cluster.association,label = 'Seurat') },error=function(v) v)
	plot.tSNE (tSNE.job,perplexity,max_iter)
}

get.cluster.names <- function (types,min.fold = 1.25,max.Qval = 0.1,print = T) {

	load(file.path(environment$res.data.path,paste('main','all','diff.exp.RData',sep='.')))
	if(print) print(summary(limma.all))
	diff.exp = limma.all[limma.all$fold > min.fold & limma.all$QValue < max.Qval,]
	if(print) print(summary(diff.exp))

	cluster = 1
	cluster.names = array('Unknown',environment$clustering$nclusters)
	for (cluster in seq(environment$clustering$nclusters)) {
		cluster.diff = diff.exp[diff.exp$cluster == cluster,]
		cluster.name = get.cluster.names.with.diff (cluster.diff,types,print)
		if (print) print(cluster.name)
		if (!is.na(cluster.name)) cluster.names[cluster] = paste(cluster.name,collapse = '_')
	}
	cluster.names

	for (name in unique(cluster.names)) {
		match = cluster.names == name
		if (sum(match) > 1)
		cluster.names[match] = paste(name,seq(sum(match)),sep='_')
	}

	return(cluster.names)
}

get.cluster.names.with.diff <- function (cluster.diff,types,print) {

	types$gene = as.vector(types$gene)
	minimum.genes.to.qualify = table(types$type)/2

	expression = cbind(types,cluster.diff[match(types$gene,cluster.diff$gene),])
	if(print) print(expression[!is.na(expression$fold),])
	table = sort(table(expression$type[!is.na(expression$fold)])-minimum.genes.to.qualify,decreasing=T)
	if(print) print(table)
	if (sum(table>0) == 0) return(NA)
	table = table[table>0]
	cluster.name = names(table)
	
	return (cluster.name)
}

set.cluster.names <- function (names) {
	
	cluster.name.map = data.frame(id = seq(length(names)),name = names)
	environment$cluster.names = cluster.names = names[environment$clustering$membership]
	save(cluster.names,cluster.name.map,file = file.path(environment$res.data.path,'cluster.names.RData'))
	write.csv(cluster.name.map,file = file.path(environment$work.path,'cluster.name.map.csv'))
	print(table(environment$cluster.names))
	return(environment)
}

load.cluster.names <- function () {

	load(file.path(environment$res.data.path,'cluster.names.RData'))
	environment$cluster.names = cluster.names
	print(table(environment$cluster.names))
	return(environment)
}

remove.cluster.names <- function () {

	environment$cluster.names = environment$clustering$membership

	return(environment)
}

filter.cluster.data <- function (remove.clusters) {
	membership  = as.vector(environment$clustering$membership)
	keep = !membership %in% remove.clusters
	filter.data(keep)
}

filter.data <- function (keep) {
	data.file = file.path(environment$baseline.data.path,'data.RData')
	load(data.file)
	file.rename(data.file,paste(data.file,format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),sep='---'))
	counts = counts[,keep]
	genes.filter = genes.filter & apply(counts,1,var)>0
	normalized = normalized[,keep]
	dataset.labels = dataset.labels[keep]
	origins = origins[keep]
	experiments = experiments[keep]
	unlink(environment$baseline.data.path,recursive=T,force=T)
	dir.create(environment$baseline.data.path)
	cache = file.path(environment$baseline.data.path,'data.RData')
	save(genes.filter,counts,normalized,dataset.labels,origins,experiments,file=cache)
	file.rename(environment$work.path,paste(environment$work.path,'pre.filter',format(Sys.time(), "%a_%b_%e_%Y__%H_%M_%S"),sep='_'))
}

filter.robust.clusters <- function (robust.clusters) {

	load(file.path(environment$baseline.data.path,'preclustered.datasets.RData'))

	membership  = as.vector(environment$clustering$membership)
	keep = membership %in% robust.clusters

	counts = counts[,keep]
	normalized = normalized[,keep]
	genes.filter = genes.filter & apply(counts,1,var)>0
	dataset.labels = dataset.labels[keep]
	origins = origins[keep]
	experiments = experiments[keep]
	HVG = NA
	clustering = clustering[keep]
	merged.original.clustering = merged.original.clustering[keep]
	merged.diff.exp = NA

	dir = dirname(environment$work.path)
	new.dir = file.path(dirname(dir),paste('filtered',basename(dir),sep='_'),'data')
	dir.create(new.dir,recursive=T)

	save(genes.filter,counts,normalized,dataset.labels,origins,experiments,HVG,clustering,merged.diff.exp,merged.original.clustering,file=file.path(new.dir,'preclustered.datasets.RData'))
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

start <- function (name = NA,append = F,split = F,print = T) {
	time = get.time()
	label = as.character(sys.calls()[[sys.nframe()-1]])[1]
	if (!is.na(name)) label = paste(label,name,sep='.')
	file = paste(label,'txt',sep='.')
	sink( file.path('tracking',file),type = 'output',append,split )
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
