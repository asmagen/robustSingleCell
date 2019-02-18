extract_pages () {
  source=$1
  dest=$2
  page=$3
  pdftk $source cat $page output $dest
}


extract_pages ~/LCMV/LCMV_analysis/LCMV2/LCMV2.pre.filter.dataset.stats.pdf ~/Documents/repos/magensinglecell/vignettes/figs/pre_stats.pdf 3
extract_pages ~/LCMV/LCMV_analysis/LCMV2/LCMV2.post.filter.dataset.stats.pdf ~/Documents/repos/magensinglecell/vignettes/figs/post_stats.pdf 3

extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/all.PCA.pdf figs/PCA.pdf 1
extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/Rotation.PCA.pdf figs/rotation.pdf 2
extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/PC.scores.heatmap.pdf figs/PC1_heatmap.pdf 1
extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/PC.scores.heatmap.pdf figs/PC2_heatmap.pdf 2

extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/cluster.size.pdf figs/cluster.size.pdf 3

extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/confounder.stats.violin.pdf figs/violin.pdf 7

extract_pages ~/LCMV/LCMV_analysis/LCMV2/not.regressed/tSNE_perplexity.20.max_iter.10000.pdf figs/tsne.pdf 1

# heatmap todo extract_pages

extract_pages ~/LCMV/LCMV_analysis/merged.LCMV1.LCMV2/not.regressed/cross.sample.similarities/LCMV1/LCMV2/LCMV1_to_LCMV2_similarity.heatmap.pdf figs/LCMV1_LCMV2_similarity.pdf 1

extract_pages ~/LCMV/LCMV_analysis/merged.LCMV1.LCMV2/not.regressed/compare.cluster.similarity_pearson/filtered.cluster.similarity.heatmap.Cor.FC.pdf figs/cor_FC.pdf 1
extract_pages ~/LCMV/LCMV_analysis/merged.LCMV1.LCMV2/not.regressed/compare.cluster.similarity_pearson/filtered.cluster.similarity.heatmap.Cor.means.pdf figs/cor_means.pdf 1

