# pre-calculation of fisher.test outputs
#

#' precompute fisher p-values and odds ratios
#'
#' Imagine that two sets a and b are composed of items from some superset.
#' These sets are compared, which leads to a contingency table with components
#' count_11, count_10, count_01, and count_00. The contingency table can be
#' processed using a fisher test to obtain a p-value and an odds-ratio.
#'
#' @param universe_size integer, large number
#' @param a_sizes vector of integers, values to consider for sizes of set a
#' @param b_sizes vector of integers, values to consider for szies of set b
#' @param max_or numeric, maximal population odds-ratio;
#' When max_or is NULL, the precomputed configurations contain all scenarios
#' that have sets a and b compatible with the previous arguments.
#' When max_or is a number, some configurations are omitted if they are unlikely.
#' @import data.table
#' @export
#'
#' @return data.table with columns (a,b,c,d) that describe entries in a contingency
#' table and columns (p.value, odds.ratio) that capture output from fisher.test
precompute_fisher = function(universe_size,
                             a_sizes=seq(0, universe_size), b_sizes=a_sizes,
                             max_or=NULL) {

  # prepare a table of all configurations
  a_sizes = unique(pmin(universe_size, pmax(a_sizes, 0)))
  b_sizes = unique(pmin(universe_size, pmax(b_sizes, 0)))
  ab_sizes = unique(c(a_sizes, b_sizes))
  pairs = expand.grid(list(x=ab_sizes, y=ab_sizes))
  pairs = split(pairs, seq_len(nrow(pairs)))
  configs = rbindlist(lapply(pairs, function(z) {
    x = z$x
    y = z$y
    # find the maximal overlap between the sets
    if (is.null(max_or)) {
      # consider all possible overlaps
      solution = min(x, y)
    } else {
      # capping by odds-ratio
      # This scenario gives a quadratic equation with the following coefficients
      aa =  max_or-1
      bb = -universe_size - (x + y)*(max_or-1)
      cc = max_or * x * y
      solutions = (-bb + c(1, -1)* sqrt(bb*bb - 4*aa*cc))/(2*aa)
      solution = ceiling(solutions[solutions <= min(x, y)])
    }
    data.table(count_11=seq(0, solution),
               count_10=x-seq(0, solution),
               count_01=y-seq(0, solution))
  }))
  configs$count_00 = universe_size - (configs$count_11 + configs$count_01 + configs$count_10)
  # consider only positive entries,
  # Also pick to compute only when count_10 > count_01 (will symmetrize later)
  configs = configs[count_00 >= 0 & count_10 >= count_01]
  
  # brute-force compute fisher outputs, then symmetrize
  result = precompute_fisher_contingency(configs)
  other = copy(result[count_10 > count_01])
  other[, c("count_01", "count_10") := list(count_10, count_01)]
  result = rbind(result, other)
  result = result[(count_11 + count_10) %in% a_sizes]
  result = result[(count_11 + count_01) %in% b_sizes]
  result[order(count_11, count_10, count_01, count_00)]
}


#' precompute fisher p-values and odds-ratios from a table of contingency values
#'
#' @param x data table with count columns (count_11, count_10, count_01, count_00)
#' @import data.table
#' @importFrom stats fisher.test
#' @export
#'
#' @return a table with the contingency values, and new columns with p.value and odds.ratio
precompute_fisher_contingency = function(x) {
  count_columns = paste0("count_", c("11", "10", "01", "00"))
  if (!all(count_columns %in% colnames(x))) {
    stop("input data object has missing count columns\n")
  }  
  po = function(v) {
    f = fisher.test(matrix(v, nrow=2))
    list(f$p.value, as.numeric(f$estimate))
  }
  result = data.table(x)
  result[, c("p.value", "odds.ratio") := po(c(count_11, count_10, count_01, count_00)),
                  by=count_columns]
}

