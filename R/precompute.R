# pre-calculation of fisher.test outputs
#

#' precompute fisher p-values and odds ratios 
#'
#' Nomenclature: set a is compared against set b, the contingency matrix is
#' composed of count_11, count_10, count_01, and count_00
#'
#' @param universe_size integer, large number
#' @param a_vals vector of integers, values to consider for sizes of set a
#' @param b_vals vector of integers, values to consider for szies of set b
#' @param max_or numeric, maximal population odds-ratio;
#' When max_or is NULL, the precomputed configurations contain all scenarios
#' that have sets a and b compatible with the previous arguments.
#' When max_or is a number, some configurations are omitted if they are unlikely.
#' @export
#'
#' @return data.table with columns (a,b,c,d) that describe entries in a contingency
#' table and columns (p.value, odds.ratio) that capture output from fisher.test
precompute_fisher = function(universe_size,
                             a_vals=seq(0, universe_size), b_vals=a_vals,
                             max_or=NULL) {

  # prepare a table of all configurations
  a_vals = unique(pmin(universe_size, pmax(a_vals, 0)))
  b_vals = unique(pmin(universe_size, pmax(b_vals, 0)))
  if (is.null(max_or)) {
    # construct exhaustive set of configurations
    configs = list(count_11=seq(0, min(max(a_vals), max(b_vals))),
                   count_10=a_vals,
                   count_01=b_vals)
    configs = data.table(expand.grid(configs))
    configs = configs[(count_11 + count_10) %in% a_vals &
                      (count_11 + count_01) %in% b_vals, ]
  } else {
    # consider pairings of sizes a and b, for each construct a small set
    temp = expand.grid(list(count_10=a_vals, count_01=b_vals))
    temp = split(temp, seq_len(nrow(temp)))
    configs = rbindlist(lapply(temp, function(x) {
      # for a given size of set a and set b, compute the maximal overlap
      # that provides an odds-ratio of max_or
      # This calculation gives a quadratic formula with the following coefficients
      aa = max_or-1
      bb = -universe_size - (x$count_10 + x$count_01)*(max_or-1)
      cc = max_or * x$count_10 * x$count_01
      solutions = (-bb + c(1, -1)* sqrt(bb*bb - 4*aa*cc))/(2*aa)
      solution = ceiling(solutions[solutions <= x$count_01 & solutions <= x$count_10])
      x.configs = list(count_11=seq(0, solution),
                       count_10=x$count_10-seq(0, solution),
                       count_01=x$count_01-seq(0, solution))
      x.out = data.table(expand.grid(x.configs))
      x.out[(count_11+count_10) == x$count_10 & (count_11+count_01)==x$count_01]
    }))
  }
  # both above branches produce too many configurations, so trim
  configs$count_00 = universe_size - (configs$count_11 + configs$count_01 + configs$count_10)
  # consider only positive entries
  # also consider only when 10>01 (will symmetrize afterward)
  configs = unique(configs[count_00 >= 0 & count_10 >= count_01])
  
  # brute-force compute fisher outputs, then symmetrize
  result = precompute_fisher_contingency(configs)
  other = copy(result[count_10 > count_01])
  temp = other$count_01
  other$count_01 = other$count_10
  other$count_10 = temp
  
  result = rbind(result, other)
  unique(result)[order(count_11, count_10, count_01, count_00)]
}


#' precompute fisher p-values and odds-ratios from a table of contingency values
#'
#' @param x data table with columns count_11, count_10, count_01, count_00
#' @export
#'
#' @return a table with the contingency values, and new columns with p.value and odds.ratio
precompute_fisher_contingency = function(x) {  
  p.values = rep(1.0, nrow(x))
  odds.ratios = rep(1.0, nrow(x))
  xt = t(as.matrix(x[, c("count_11", "count_10", "count_01", "count_00")]))
  for (i in 1:ncol(xt)) {
    ff = fisher.test(matrix(xt[,i], nrow=2))
    p.values[i] = ff$p.value
    odds.ratios[i] = as.numeric(ff$estimate)
  }
  result = data.table(t(xt))
  result$p.value = p.values
  result$odds.ratio = odds.ratios
  result
}

