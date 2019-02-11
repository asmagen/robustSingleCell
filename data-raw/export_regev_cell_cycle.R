cell_cycle_genes <- readr::read_csv('data/regev_lab_cell_cycle_genes.txt', col_names = F)
cell_cycle_genes <- unlist(cell_cycle_genes$X1)
marker_genes <- readr::read_csv("data/markers.csv")
usethis::use_data(marker_genes, cell_cycle_genes, internal = TRUE, overwrite = T)
