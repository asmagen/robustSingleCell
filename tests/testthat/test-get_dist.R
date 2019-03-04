context("Get pairwise distance")
library(robustSingleCell)

dist_mat <- dist(c(1, 1, 2, 3, 5))

test_that("the distance for row i is", {
  expect_equal(get_dist(dist_mat, 1, 5), unname(as.matrix(dist_mat)[1,]))
  expect_equal(get_dist(dist_mat, 3, 5), unname(as.matrix(dist_mat)[3,]))
  expect_equal(get_dist(dist_mat, 5, 5), unname(as.matrix(dist_mat)[5,]))
})
