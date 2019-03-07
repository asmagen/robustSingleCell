context("Get pairwise distance")
library(robustSingleCell)

names <- sample(LETTERS, 5, replace = F)
dist_mat <- dist(structure(c(1, 1, 2, 3, 5), names = names))

test_that("the distance for row i is", {
  expect_equal(get_dist(dist_mat, 1, 5, names), as.matrix(dist_mat)[1,])
  expect_equal(get_dist(dist_mat, 3, 5, names), as.matrix(dist_mat)[3,])
  expect_equal(get_dist(dist_mat, 5, 5, names), as.matrix(dist_mat)[5,])
})
