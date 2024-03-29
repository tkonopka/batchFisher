# Tests for obtaining contingency vectors from sets


# Some data to use within the tests

if (!exists("letter_sets")) {
  letter_sets = lapply(1:20, function(x) { letters[1:x] })
  names(letter_sets) = LETTERS[1:length(letter_sets)]
}



# calculation of contingency vectors

test_that("contingency vec works with an empty set (character)", {
  result.a = contingency_vector(letters[1:3], c(), 40)
  result.b = contingency_vector(c(), letters[1:3], 40)
  expect_equal(result.a, c(0, 3, 0, 37))
  expect_equal(result.b, c(0, 0, 3, 37))
})


test_that("contingency vec works with an empty set (int)", {
  result.a = contingency_vector(1:3, c(), 40)
  result.b = contingency_vector(c(), 1:3, 40)
  expect_equal(result.a, c(0, 3, 0, 37))
  expect_equal(result.b, c(0, 0, 3, 37))
})


test_that("contingency vec counts with characters", {
  alonger = contingency_vector(letters[3:6], letters[1:5], 40)
  expect_equal(alonger, c(3, 1, 2, 34))
  blonger = contingency_vector(letters[1:5], letters[3:6], 40)
  expect_equal(blonger, c(3, 2, 1, 34))
})


test_that("contingency vec counts with integers", {
  alonger = contingency_vector(3:6, 1:5, 40)
  expect_equal(alonger, c(3, 1, 2, 34))             
  blonger = contingency_vector(1:5, 3:6, 40)
  expect_equal(blonger, c(3, 2, 1, 34))
})


test_that("contingency vec counts with complete overlap (character)", {
  alonger = contingency_vector(letters[1:10], letters[3:6], 40)
  expect_equal(alonger, c(4, 6, 0, 30))
  blonger = contingency_vector(letters[3:6], letters[1:10], 40)
  expect_equal(blonger, c(4, 0, 6, 30))
})


test_that("contingency vec counts with off-by-one (character)", {
  alonger = contingency_vector(letters[3:7], letters[2:6], 40)
  expect_equal(alonger, c(4, 1, 1, 34))
  blonger = contingency_vector(letters[2:6], letters[3:7], 40)
  expect_equal(blonger, c(4, 1, 1, 34))
})




# batch contingency tables

test_that("batch_contingency requires equal length a, b", {
  expect_error(batch_contingency(1:4, 1:6, 10, sets=letter_sets),
               "incompatible")
})


test_that("batch_contingency requires same types of inputs for a and b", {
  expect_error(batch_contingency(1:4, LETTERS[1:4], 10, sets=letter_sets),
               "incompatible")
})


test_that("batch_contingency requires lists or characters", {
  bad.input = matrix(0, ncol=2, nrow=2)
  expect_error(batch_contingency(bad.input, bad.input, 10, sets=letter_sets),
               "invalid inputs")
})


test_that("batch_contingency detects missing sets", {
  # when a/b are vectors, sets must be supplied
  expect_error(batch_contingency(1:4, 1:4, 10),
               "insufficient")
})


test_that("batch_contingency does not accept empty ", {
  expect_error(batch_contingency(c(), c(), 100, sets=letter_sets), "invalid")
})


test_that("batch_contingency can process batch of size 1 (from lists)", {
  result = batch_contingency(letter_sets[1], letter_sets[2], 100)
  expect_equal(nrow(result), 1)
})


test_that("batch_contingency can process batch of size 1 (from named sets)", {
  result = batch_contingency("A", "B", 100, sets=letter_sets)
  expect_equal(nrow(result), 1)
})


test_that("batch_contingency gives counts as integers, not numerics", {
  result = batch_contingency(letter_sets[1:4], letter_sets[2:5], 100)
  expect_equal(class(result$count_00), "integer")
  expect_equal(class(result$count_01), "integer")
  expect_equal(class(result$count_10), "integer")
  expect_equal(class(result$count_11), "integer")
})


test_that("batch_contingency can process (from lists)", {
  result = batch_contingency(letter_sets[1:4], letter_sets[2:5], 100)
  expect_equal(result$count_11, 1:4)
  expect_equal(result$count_01, rep(1, length=4))
  expect_equal(result$count_10, rep(0, length=4))  
})

test_that("batch_contingency can process (from named sets)", {
  setnames = names(letter_sets)
  result = batch_contingency(setnames[1:4], setnames[2:5], 100, sets=letter_sets)
  expect_equal(result$count_11, 1:4)
  expect_equal(result$count_01, rep(1, length=4))
  expect_equal(result$count_10, rep(0, length=4))
})


test_that("batch_contingency tracks labels (from lists)", {
  result = batch_contingency(letter_sets[1:4], letter_sets[2:5], 100)
  expect_equal(result$a, names(letter_sets)[1:4])
  expect_equal(result$b, names(letter_sets)[2:5])
})


test_that("batch_contingency tracks labels (from named sets)", {
  anames = names(letter_sets[1:3])
  bnames = names(letter_sets[5:7])
  result = batch_contingency(anames, bnames, 100, sets=letter_sets)
  expect_equal(result$a, anames)
  expect_equal(result$b, bnames)
})


test_that("batch_contingency tracks labels (from set indexes)", {
  anames = names(letter_sets[1:3])
  bnames = names(letter_sets[5:7])
  result.ints = batch_contingency(1:3, 5:7, 100, sets=letter_sets)
  result.names = batch_contingency(anames, bnames, 100, sets=letter_sets)
  expect_equal(result.ints$count_01, result.names$count_01)
  expect_equal(result.ints$count_10, result.names$count_10)
  expect_equal(result.ints$b, 5:7)
  expect_equal(result.names$b, bnames)
})


test_that("batch contingency provides answers in same order as inputs", {
  # define a series where some indexes in ** a ** are repeated
  a.indexes = as.integer(c(2,4,2,2,5))
  b.indexes = as.integer(c(1,2,5,8,4))
  asets = letter_sets[a.indexes]
  bsets = letter_sets[b.indexes]
  u = c(100, 120, 100, 120, 100)
  result.ints = batch_contingency(a.indexes, b.indexes, u, sets=letter_sets)
  result.sets = batch_contingency(asets, bsets, u)
  # the counts must be the same
  # (this is nontrivial because the batch_contingency with indexes does not
  # perform the calculations in order. So it is important that it outputs
  # the results as they are intended)
  expect_equal(result.ints$count_11, result.sets$count_11)
  expect_equal(result.ints$count_00, result.sets$count_00)
  expect_equal(result.ints$count_01, result.sets$count_01)
  expect_equal(result.ints$count_10, result.sets$count_10)
})


test_that("batch contingency when all sets are the same", {
  # define a series where some indexes in ** a ** are repeated
  a.indexes = as.integer(rep(3, 5))
  b.indexes = as.integer(c(1,2,8,5,4))
  asets = letter_sets[a.indexes]
  bsets = letter_sets[b.indexes]
  result.ints = batch_contingency(a.indexes, b.indexes, 200, sets=letter_sets)
  result.sets = batch_contingency(asets, bsets, 200)
  expect_equal(result.ints$count_11, result.sets$count_11)
  expect_equal(result.ints$count_00, result.sets$count_00)
  expect_equal(result.ints$count_01, result.sets$count_01)
  expect_equal(result.ints$count_10, result.sets$count_10)
})

