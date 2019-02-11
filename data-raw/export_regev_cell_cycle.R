cell_cycle_genes <- readr::read_csv('data/regev_lab_cell_cycle_genes.txt', col_names = F)
cell_cycle_genes <- unlist(cell_cycle_genes$X1)
usethis::use_data(cell_cycle_genes, internal = TRUE, overwrite = T)
