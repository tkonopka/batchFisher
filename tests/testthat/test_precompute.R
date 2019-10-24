# Tests for pre-calculation of tables with p-values and odds ratios


count.columns = paste0("count_", c("11", "10", "01", "00"))


test_that("precompute produces a small table with p-values and odds.ratios", {
  result = precompute_fisher2x2(50, vals=c(0,1,2))
  # general structure should be a tbale with fixed columns
  expect_true("data.table" %in% class(result))
  expect_equal(colnames(result), c(count.columns, "p.value", "odds.ratio"))
  # all rows should add up to the universe size
  result.sums = apply(as.matrix(result)[, count.columns], 1, sum)
  expect_true(all(result.sums==50))
})


test_that("precompute makes all combinations, without rotation", {
  result = precompute_fisher2x2(50, vals=c(0,1))
  expect_equal(sort(unique(result$count_11)), 0:1)
  expect_equal(sort(unique(result$count_10)), 0:1)
  expect_equal(sort(unique(result$count_01)), 0:1)
  expect_equal(sort(unique(result$count_00)), 47:50)
})


test_that("precompute makes all combinations, with rotation", {
  result = precompute_fisher2x2(50, vals=c(0,1), rotate=TRUE)
  counts = c(0, 1, 47:50)
  expect_equal(sort(unique(result$count_11)), counts)
  expect_equal(sort(unique(result$count_10)), counts)
  expect_equal(sort(unique(result$count_01)), counts)
  expect_equal(sort(unique(result$count_00)), counts)
})


