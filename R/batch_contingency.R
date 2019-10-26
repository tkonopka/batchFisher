# Calculation of contingency vectors
#
# Contingency vectors are length-4 vectors of integers that can be put
# into a 2x2 matrix to give a traditional contingency table.


#' summarize overlap of two sets into a vector that is ready for a contingency matrix
#'
#' @param a vector, a set of elements, can be 'numeric' or 'character'
#' @param b vector, a set of elements, can be 'numeric' or 'character'
#' @param universe_size integer, total number of elements 
#' @export
#'
#' @return vector with four components (count_11, count_10, count_01, count_00)
contingency_vector = function(a, b, universe_size) {
  common = sum(a %in% b)
  a_only = length(a) - common
  b_only = length(b) - common
  c(common, a_only, b_only, universe_size - common - a_only - b_only)
}


#' summarize overlap of two sets into a vector (partial contingency vector)
#'
#' @keywords internal
#' @param a vector, a set of elements, can be 'numeric' or 'character'
#' @param b vector, a set of elements, can be 'numeric' or 'character'
#'
#' @return vector with three components (count_11, count_10, count_01)
contingency_triple = function(a, b) {
  common = sum(a %in% b)
  c(common, length(a) - common, length(b) - common)
}


#' construct a table with contingency values comparing many configuration of type (a vs. b)
#'
#' @param a sets a. if of type list, each component is interpreted as a set.
#' If a simple vector of characters or integers, interpreted as indexes into
#' the list of sets. It is advantageous when many a sets are the same.
#' @param b sets b to match a;
#' if type list, each component is interpreted as a set;
#' if a simple vector, each component is interpreted as an index into 'sets'
#' @param sets list with vectors, used when 'a' or 'b' are not lists
#' @param universe_size integer, the total number of elements that sets a, b were drawn from
#' @import data.table
#' @export
#'
#' @return data table with one row per item, column id and a,b,c,d
batch_contingency = function(a, b, universe_size, sets=NULL) {
  
  n = length(a)
  if (n != length(b)) {
    stop("incompatible components 'a' and 'b' - must be of equal length\n")
  }
  if (!identical(class(a), class(b))) {
    stop("incompatible components 'a' and 'b' - must of the same class\n")
  }
  if (!(class(a) %in% c("list", "character", "integer"))) {
    stop("invalid inputs - must be in form of list or character/integer vectors\n")
  }
  if (!is.list(a) & !is.list(sets)) {
    stop("insufficient information - sets must be specified if a, b are vectors\n")
  }
  
  universe_size = rep(universe_size, length.out=n)
  
  # compute the contingency tables
  result = matrix(as.integer(0), ncol=n, nrow=3)
  if (is.list(a) & is.list(b)) {
    # all configurations are provided explicitly - no shortcuts
    for (i in seq_len(n)) {
      result[, i] = contingency_triple(a[[i]], b[[i]])
    }
  } else {
    # configurations are provided through sets
    # perform the calculations in an order that avoids sets lookup for "a" when possible
    if (length(unique(a))==1) {
      # when all the a sets are the same, no need to look them up within the loop
      aset = sets[[a[1]]]
      for (i in seq_len(n)) {
        bset = sets[[b[i]]]
        result[, i] = contingency_triple(aset, bset)
      }
    } else {
      # when a sets come in different types, group them to avoid some lookups
      fwd.order = order(a)
      rev.order = rank(a, ties="first")
      a = a[fwd.order]
      b = b[fwd.order]
      anew = !duplicated(a)
      for (i in seq_len(n)) {
        if (anew[i]) {
          aset = sets[[a[i]]]
        }
        bset = sets[[b[i]]]
        result[, i] = contingency_triple(aset, bset)
      }
      result = result[, rev.order, drop=FALSE]
      a = a[rev.order]
      b = b[rev.order]
    }
  }
  count_00 = as.integer(universe_size - colSums(result))
  result = data.table(t(result))
  count.cols = c("count_11", "count_10", "count_01", "count_00")
  colnames(result) = count.cols[1:3]
  result$count_00 = count_00
  
  # add annotations to keep track of which contingency tables belong to which inputs
  if (is.list(a) & is.list(b)) {
    result$a = names(a)
    result$b = names(b)
  } else {
    result$a = a
    result$b = b
  }
  setcolorder(result, intersect(c("a", "b", count.cols), colnames(result)))
  
  result
}

