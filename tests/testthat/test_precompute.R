# Tests for pre-calculation of tables with p-values and odds ratios


count.columns = paste0("count_", c("11", "10", "01", "00"))




# calculation of fisher values from set descriptions

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


test_that("precompute can limit configuration by odds-ratio", {
  result.all = precompute_fisher(800, 0:6)
  result.or = precompute_fisher(800, 0:6, max_or=30)
  expect_true(nrow(result.or)<nrow(result.all))
})


test_that("precompute takes different ranges for a and b values", {
  result = precompute_fisher(50, a_sizes=c(0,1,2,3), b_sizes=c(0,1))
  expect_equal(sort(unique(result$count_11)), c(0,1))
  # a values should be between 0 and 3
  expect_equal(sort(unique(result$count_10)), c(0,1,2,3))
  # b values were specified as 0,1 - but the output is symmetrized
  expect_equal(sort(unique(result$count_01)), c(0,1,2,3))
  # some combinations will not occur because of the a_vals and b_vals ranges
  expect_equal(nrow(result[count_01==3 & count_10==3]), 0)
})




# calculation of fisher values from configurations

test_that("precomputing from contingency table checks structure", {
  conf0 = data.table(count_11=4, count_10=2, count_01=3)
  expect_error(precompute_fisher_contingency(conf0), "missing count")
})


test_that("precomputing from contingency table works with single configuration", {
  conf1 = data.table(count_11=4, count_10=2, count_01=3, count_00=30)
  result = precompute_fisher_contingency(conf1)
  expect_equal(dim(result), c(1, 6))
  expected = fisher.test(matrix(c(4,2,3,30), nrow=2))
  expect_equal(result$p.value, expected$p.value)
  expect_equal(result$odds.ratio, as.numeric(expected$estimate))  
})


test_that("precomputing from a table works with a couple of configuration", {
  conf2 = data.table(count_11=c(4,5), count_10=2, count_01=3, count_00=30)
  result = precompute_fisher_contingency(conf2)
  expect_equal(dim(result), c(2, 6))
  expected1 = fisher.test(matrix(c(4,2,3,30), nrow=2))
  expected2 = fisher.test(matrix(c(5,2,3,30), nrow=2))
  expect_equal(result$p.value, c(expected1$p.value, expected2$p.value))
})



