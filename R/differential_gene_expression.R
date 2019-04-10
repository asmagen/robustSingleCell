getDE.limma <- function(Y, group, filter = T) {
    design <- stats::model.matrix(~group)
    my.lm <- limma::lmFit(Y, design = design)
    my.lm <- limma::eBayes(my.lm)
    topTable <- limma::topTable(my.lm, number = Inf, coef = seq(2, ncol(design)))
    if (filter)
        topTable <- topTable[topTable$adj.P.Val < 0.1, ]
    topTable <- topTable[order(topTable$adj.P.Val, decreasing = F), ]
    colnames(topTable)[colnames(topTable) == "P.Value"] <- "PValue"
    colnames(topTable)[colnames(topTable) == "adj.P.Val"] <- "QValue"

    return(topTable)
}

run.diff.expression <- function(environment, clustering, min.fold, quantile, label,
    rerun = F, contrast = "all", contrast.groups = NA) {

    cache <- file.path(environment$res.data.path, paste(label, contrast, "diff.exp.rds",
        sep = "."))

    if (!rerun && file.exists(cache)) {
        print.message("Loading precomputed")
        precomputed <- readRDS(cache)
        final.diff <- precomputed$final.diff
        limma.diff <- precomputed$limma.diff
        empirical.diff <- precomputed$empirical.diff
        limma.all <- precomputed$limma.all
        rm(precomputed)
    } else {
        print.message("Computing")
        t <- start(file.path(environment$work.path, "tracking"), split = T)
        on.exit(end(t))
        membership <- as.vector(clustering$membership)

        empirical.diff <- NA

        diff.exp.dir <- file.path(environment$work.path, "diff.exp", label, contrast)
        unlink(diff.exp.dir, recursive = T, force = T)
        dir.create(diff.exp.dir, recursive = T)
        limma.diff <- {
        }

        if (contrast == "all") {
            contrast.groups <- rep(1, length(membership))
        } else if (contrast == "datasets") {
            contrast.groups <- environment$dataset.labels
        }

        print.message("contrast =", contrast)

        print(table(contrast.groups))
        contrast.group <- unique(contrast.groups)[1]
        table(membership, contrast.groups)
        for (contrast.group in unique(contrast.groups)) {
            print.message("contrast.group =", contrast.group)
            contrast.group.indices <- contrast.groups == contrast.group
            sum(contrast.group.indices)
            contrast.group.clusters <- sort(unique(membership[contrast.group.indices]))
            cluster <- contrast.group.clusters[1]
            for (cluster in contrast.group.clusters) {
                print.message("cluster =", cluster)

                group <- membership[contrast.group.indices] == cluster
                group <- factor(group)

                diff.exp <- getDE.limma(Y = environment$normalized[, contrast.group.indices],
                  group = group, filter = F)
                diff.exp <- diff.exp[order(diff.exp$logFC, decreasing = T), ]
                diff.exp <- data.frame(gene = rownames(diff.exp), logFC = diff.exp$logFC,
                  fold = exp(diff.exp$logFC), QValue = diff.exp$QValue, PValue = diff.exp$PValue,
                  AveExpr = diff.exp$AveExpr)
                utils::write.csv(diff.exp, file = file.path(diff.exp.dir, paste("cluster",
                  cluster, "csv", sep = ".")))
                markers.diff <- diff.exp[diff.exp$gene %in% environment$marker.genes,
                  ]
                utils::write.csv(markers.diff, file = file.path(diff.exp.dir, paste("markers.cluster",
                  cluster, "csv", sep = ".")))
                limma.diff <- rbind(limma.diff, data.frame(contrast.group = contrast.group,
                  cluster = cluster, diff.exp[, c(1, 3, 4)]))
            }
        }
        rownames(limma.diff) <- NULL
        limma.all <- limma.diff
        limma.diff <- limma.diff[limma.diff$QValue <= (1 - quantile) & limma.diff$fold >=
            min.fold, ]
        final.diff <- limma.diff <- limma.diff[order(limma.diff$fold, decreasing = T),
            ]
        print.message("head(limma.diff):")
        print(utils::head(limma.diff))
        utils::write.csv(limma.diff, file = file.path(diff.exp.dir, "limma.diff.csv"))
        utils::write.csv(final.diff, file = file.path(diff.exp.dir, "final.diff.csv"))

        saveRDS(list(final.diff = final.diff, limma.diff = limma.diff, empirical.diff = empirical.diff,
            limma.all = limma.all), file = cache)
    }

    return(final.diff)
}


#' Get Robust Marker
#'
#' Analysis of robust subpopulation marker prioritization
#'
#' @param environment \code{environment} object
#' @param cluster_group1 cluster group 1 to be used as a foreground
#' @param cluster_group2 cluster group 2 to be used as a background
#' @param group1_label label for group 1
#' @param group2_label label for group 2
#' @param annotate.genes specific gene names to annotate in figure in addition to novel markers
#' @param min.fold.diff average expression fold change cutoff
#' @param min.ratio.diff detection ratio fold change cutoff
#' @param QValue Qvalue cutoff
#' @import ggrepel
#' @importFrom graphics plot
#' @export
#' @examples
#' \donttest{
#' LCMV1 <- setup_LCMV_example()
#' LCMV1 <- get.variable.genes(LCMV1, min.mean = 0.1, min.frac.cells = 0,
#' min.dispersion.scaled = 0.1)
#' LCMV1 <- PCA(LCMV1)
#' LCMV1 <- cluster.analysis(LCMV1)
#' diff_exp <- get.robust.markers(LCMV1,
#' cluster_group1 = c('1','2'),
#' cluster_group2 = c('3','4'),
#' group1_label = 'CD4 T Cells',
#' group2_label = 'CD8 T Cells')
#' }
get.robust.markers <- function (
    environment,
    cluster_group1,
    cluster_group2,
    group1_label,
    group2_label,
    annotate.genes = NA,
    min.fold.diff = 1.5,
    min.ratio.diff = 3,
    QValue = 0.05) {

    indices = environment$cluster.names %in% c(cluster_group1,cluster_group2)
    measurements = environment$normalized[,indices]
    groups = environment$cluster.names[indices] %in% cluster_group1
    diff.exp = getDE.limma( Y = measurements, group = groups, filter = F )
    diff.exp = diff.exp[order(diff.exp$logFC,decreasing=T),]
    diff.exp = data.frame(gene=rownames(diff.exp),logFC=diff.exp$logFC,fold=exp(diff.exp$logFC),QValue=diff.exp$QValue,PValue=diff.exp$PValue,AveExpr=diff.exp$AveExpr,stringsAsFactors = F)
    diff.exp = diff.exp[((diff.exp$fold >= min.fold.diff | diff.exp$fold <= 1/min.fold.diff) & diff.exp$QValue <= QValue) | diff.exp$gene %in% annotate.genes,]
    rownames(diff.exp) = diff.exp$gene
    diff.exp$detected.cluster_group1 = rowMeans(measurements[match(diff.exp$gene,rownames(measurements)),groups]>0)
    diff.exp$detected.cluster_group2 = rowMeans(measurements[match(diff.exp$gene,rownames(measurements)),!groups]>0)
    diff.exp$detected.fold = diff.exp$detected.cluster_group1/diff.exp$detected.cluster_group2
    diff.exp = diff.exp[order(diff.exp$detected.cluster_group1/diff.exp$detected.cluster_group2,decreasing=T),]
    diff.exp = diff.exp[order(diff.exp$fold,decreasing=T),]

    diff.exp = diff.exp[,c('gene','fold','QValue','detected.cluster_group1','detected.cluster_group2','detected.fold')]
    diff.exp = diff.exp[rowSums(is.na(diff.exp))==0,]
    if (!is.na(annotate.genes)) print(diff.exp[annotate.genes,])

    figure.data = data.frame(Gene = diff.exp$gene,Fold = diff.exp$fold,Foreground.ratio = diff.exp$detected.cluster_group1,Background.ratio = diff.exp$detected.cluster_group2,Detected.fold = diff.exp$detected.fold)
    label.indices = (figure.data$Fold >= min.fold.diff | figure.data$Fold <= 1/min.fold.diff) & (figure.data$Detected.fold >= min.ratio.diff | figure.data$Detected.fold < 1/min.ratio.diff)
    label.indices.annotated = (figure.data$Fold >= min.fold.diff | figure.data$Fold <= 1/min.fold.diff) & (figure.data$Detected.fold >= min.ratio.diff | figure.data$Detected.fold < 1/min.ratio.diff) | figure.data$Gene %in% annotate.genes

    path = file.path(environment$work.path,paste(group1_label,group2_label,'Robust.diff.detection',sep='_'))
    dir.create(path)

    grDevices::pdf(file = file.path(path,paste(group1_label,group2_label,'robust.diff.exp.pdf',sep='_')),width=20,height=10)
    print(ggplot(data = figure.data, aes(Foreground.ratio, Background.ratio,label=Gene)) + geom_point(data = figure.data,aes(alpha=0.25,size=Fold)) + theme_classic(base_size=30) + xlab(paste('Detection ratio in',group1_label)) + ylab(paste('Detection ratio in',group2_label)) + geom_text_repel(data = figure.data[label.indices,],size=10,box.padding = 0.25,point.padding = 0.25) + theme(legend.position="right"))
    grDevices::dev.off()

    diff.exp = diff.exp[order(diff.exp$detected.fold,decreasing=T),]
    rownames(diff.exp) = NULL
    colnames(diff.exp) = c('Gene','Fold Change Mean Expression','QValue','Detection ratio in group1','Detection ratio in group2','Fold Change detection ratio')

    write.csv(diff.exp,file = file.path(path,paste(group1_label,group2_label,'robust.diff.exp.csv',sep='_')))

    return(diff.exp)
}
