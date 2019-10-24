# pre-calculation of fisher.test outputs
#

#' precompute a data table with fisher p-values and odds ratios 
#'
#' @param universe_size integer, large number
#' @param vals vector of integer, values to consider in first three cells in the confusion matrix
#' @param rotate logical, set TRUE to obtain full table of combinations, leave FALSE
#' to consider column "d" as the special column that holds constrained counts
#'
#' @return data.table with columns (a,b,c,d) that describe entries in a contingency
#' table and columns (p.value, odds.ratio) that capture output from fisher.test
precompute_fisher2x2 = function(universe_size, vals=seq(0, 20), rotate=FALSE) {

  # prepare a table of all configurations
  vals = unique(pmin(universe_size, pmax(vals, 0)))
  configs = list(count_11=vals, count_10=vals, count_01=vals)
  configs = data.table(expand.grid(configs))
  configs$count_00 = universe_size - (configs$count_11 + configs$count_01 + configs$count_10)
  # consider only positive entries
  # also consider only when 10>01 (will symmetrize afterward)
  configs = configs[count_00 >=0 & count_10 >= count_01]
  
  # brute-force compute fisher outputs, then symmetrize
  result = fisher2x2(configs)
  other = copy(result[count_10 > count_01])
  temp = other$count_01
  other$count_01 = other$count_10
  other$count_10 = temp
  
  result = rbind(result, other)
  if (rotate) {
    result_1 = rotate_abcd(result)
    result_2 = rotate_abcd(result_1)
    result_3 = rotate_abcd(result_2)
    result = rbind(result, result_1, result_2, result_3)
  }
  
  unique(result)[order(count_11, count_10, count_01, count_00)]
}



#' produce a new table with columns rotated.
#'
#' This creates configurations that are expected to have the same properties under
#' fisher.test.
#'
#' @param x data table with columns a,b,c,d
#'
#' @return data table with same structure as x, but columns a,b,c,d rotated
rotate_abcd = function(x) {
  temp = x$count_11
  x$count_11 = x$count_10
  x$count_10 = x$count_01
  x$count_01 = x$count_00
  x$count_00= temp
  x
}

