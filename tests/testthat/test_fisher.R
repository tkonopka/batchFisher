# Tests for calculation of fisher values using pre-calculated tables 




# tests for calculation of fisher values from configurations

test_that("fisher2x2 can process a single configuration", {
  conf1 = data.table(count_11=4, count_10=2, count_01=3, count_00=30)
  result = fisher2x2(conf1)
  expect_equal(dim(result), c(1, 6))
  expected = fisher.test(matrix(c(4,2,3,30), nrow=2))
  expect_equal(result$p.value, expected$p.value)
  expect_equal(result$odds.ratio, as.numeric(expected$estimate))  
})


test_that("fisher2x2 can process a couple of configuration", {
  conf2 = data.table(count_11=c(4,5), count_10=2, count_01=3, count_00=30)
  result = fisher2x2(conf2)
  expect_equal(dim(result), c(2, 6))
  expected1 = fisher.test(matrix(c(4,2,3,30), nrow=2))
  expected2 = fisher.test(matrix(c(5,2,3,30), nrow=2))
  expect_equal(result$p.value, c(expected1$p.value, expected2$p.value))
})




# tests for batch assignment of p-values

test_that("batch calculation checks that all necessary columns exist", {
  missing.cols = data.table(count_11=0, count_01=c(1,2,3,4))
  expect_error(batch_fisher(missing.cols), "count")
})


test_that("batch calculation checks p.values and odds.ratio are not already present", {
  with.p = data.table(count_11=0, count_01=c(1,2,3,4), count_10=3, count_00=20,
                      p.value=1)
  expect_error(batch_fisher(with.p), "p.value")
  with.ratio = data.table(count_11=0, count_01=c(1,2,3,4), count_10=3, count_00=20,
                          odds.ratio=1)
  expect_error(batch_fisher(with.ratio), "ratio")
})


test_that("batch calculation should return p.values and odds ratios", {
  counts = data.table(count_11=0, count_10=1:4, count_01=1:4, count_00=50)
  result = batch_fisher(counts)
  expect_true(all(c("p.value", "odds.ratio") %in% colnames(result)))
})


test_that("batch calculation should use relevant precomputed values", {
  # fudge a precomputed matrix - this is a signal that those values will be used!
  precomp = precompute_fisher(50, seq(0, 8))
  precomp$p.value = 2+seq_len(nrow(precomp))

  # define a matrix with a mixture of configurations that are in precomp
  counts.1 = data.table(count_11=2, count_10=1:4, count_01=1:4, count_00 = 48-seq(2, 8, by=2))
  counts.0 = data.table(count_11=2, count_10=1:4, count_01=1:4, count_00 = 100)
  counts = rbind(counts.0, counts.1)

  result = batch_fisher(counts, precomputed=precomp)

  ## all p-values that have 00==100 should have realistic de-novo computed p-values
  good.p = result[count_00==100]$p.value
  expect_true(all(good.p > 0 & good.p < 1))
  ## all p-values that have 00==50 should come from precomp, which was corrupted
  bad.p = result[count_00<100]$p.value
  expect_true(all(bad.p > 1))
})
  

test_that("batch calculation with and without precomputed values should be equivalnet", {
  counts = data.table(count_11=5, count_10=1:4, count_01=1:4, count_00 = 45-seq(2,8, by=2))
  result.0 = batch_fisher(counts)
  precomp = precompute_fisher(50, seq(0,6))
  result.1 = batch_fisher(counts, precomp)
  expect_equal(result.0$p.value, result.1$p.value)
})


test_that("batch calculation should preserve labels", {
  data = data.table(a=letters[1:4],
                    count_11=5, count_10=1:4, count_01=1:4, count_00 = 45-seq(2,8, by=2))
  result.0 = batch_fisher(data)
  expect_equal(colnames(result.0)[1], "a")
  precomp = precompute_fisher(50, seq(0,6))
  result.1 = batch_fisher(data, precomp)
  expect_equal(colnames(result.1)[1], "a")
})
