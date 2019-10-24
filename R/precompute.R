# pre-calculation of fisher.test outputs
#

#' precompute a data table with fisher p-values and odds ratios 
#'
#' Nomenclature: set a is compare against set b, the contingency matrix is
#' composed of count_11, count_10, count_01, and count_00
#'
#' @param universe_size integer, large number
#' @param a_vals vector of integers, values to consider for sizes of set a
#' @param b_vals vector of integers, values to consider for szies of set b
#' @export
#'
#' @return data.table with columns (a,b,c,d) that describe entries in a contingency
#' table and columns (p.value, odds.ratio) that capture output from fisher.test
precompute_fisher = function(universe_size, a_vals=seq(0, universe_size), b_vals=a_vals) {

  # prepare a table of all configurations
  a_vals = unique(pmin(universe_size, pmax(a_vals, 0)))
  b_vals = unique(pmin(universe_size, pmax(b_vals, 0)))
  configs = list(count_11=seq(0, min(max(a_vals), max(b_vals))),
                 count_10=a_vals,
                 count_01=b_vals)
  configs = data.table(expand.grid(configs))
  configs = configs[(count_11 + count_10) %in% a_vals, ]
  configs = configs[(count_11 + count_01) %in% b_vals, ]
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
  unique(result)[order(count_11, count_10, count_01, count_00)]
}

