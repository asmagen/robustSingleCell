extract_pages () {
  source=$1
  dest=$2
  page=$3
  pdftk $source cat $page output $dest
}


extract_pages ~/LCMV/LCMV_analysis/LCMV1/LCMV1.pre.filter.dataset.stats.pdf figs/pre_stats.pdf 3
extract_pages ~/LCMV/LCMV_analysis/LCMV1/LCMV1.post.filter.dataset.stats.pdf figs/post_stats.pdf 3

extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/all.PCA.pdf figs/PCA.pdf 1
extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/Rotation.PCA.pdf figs/rotation.pdf 2
extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/PC.scores.heatmap.pdf figs/PC1_heatmap.pdf 1
extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/PC.scores.heatmap.pdf figs/PC2_heatmap.pdf 2

extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/cluster.size.pdf figs/cluster.size.pdf 3

extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/confounder.stats.violin.pdf figs/violin.pdf 7

extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/tSNE_perplexity.20.max_iter.10000.pdf figs/tsne.pdf 1

extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/canonical.heatmap.pdf figs/diff.genes.pdf 1

extract_pages ~/LCMV/LCMV_analysis/merged.LCMV1.LCMV2/not.regressed/cross.sample.similarities/LCMV1/LCMV2/LCMV1_to_LCMV2_similarity.heatmap.pdf figs/LCMV1_LCMV2_similarity.pdf 1

extract_pages ~/LCMV/LCMV_analysis/merged.LCMV1.LCMV2/not.regressed/hclust.dist.Cor.FC.pdf figs/cor_FC.pdf 3

extract_pages ~/LCMV/LCMV_analysis/LCMV1/not.regressed/Thu_Feb_21_2019__16_50_54/Clustering.modularity.pdf figs/Clustering.modularity.pdf 1

cd figs

find . -type f -name '*.pdf' -print0 |
  while IFS= read -r -d '' file
    do convert -verbose -density 500 -resize 800 "${file}" "${file%.*}.png"
  done
