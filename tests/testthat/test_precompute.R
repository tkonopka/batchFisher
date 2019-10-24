# Tests for pre-calculation of tables with p-values and odds ratios


count.columns = paste0("count_", c("11", "10", "01", "00"))


test_that("precompute produces a small table with p-values and odds.ratios", {
  result = precompute_fisher(50, c(0,1,2))
  # general structure should be a tbale with fixed columns
  expect_true("data.table" %in% class(result))
  expect_equal(colnames(result), c(count.columns, "p.value", "odds.ratio"))
  # all rows should add up to the universe size
  result.sums = apply(as.matrix(result)[, count.columns], 1, sum)
  expect_true(all(result.sums==50))
})


test_that("precompute makes all combinations", {
  result = precompute_fisher(50, c(0,1))
  expect_equal(sort(unique(result$count_11)), 0:1)
  expect_equal(sort(unique(result$count_10)), 0:1)
  expect_equal(sort(unique(result$count_01)), 0:1)
  expect_equal(sort(unique(result$count_00)), 48:50)
})


test_that("precompute takes different ranges for a and b values", {
  result = precompute_fisher(50, a_vals=c(0,1,2,3), b_vals=c(0,1))
  expect_equal(sort(unique(result$count_11)), c(0,1))
  # a values should be between 0 and 3
  expect_equal(sort(unique(result$count_10)), c(0,1,2,3))
  # b values were specified as 0,1 - but the output is symmetrized
  expect_equal(sort(unique(result$count_01)), c(0,1,2,3))
  # some combinations will not occur because of the a_vals and b_vals ranges
  expect_equal(nrow(result[count_01==3 & count_10==3]), 0)
})

