# Calculation of contingency vectors
#
# Contingency vectors are length-4 vectors of integers that can be put
# into a 2x2 matrix to give a traditional contingency table.


#' summarize overlap of two sets
#'
#' @param a vector, a set of elemnts, can be 'numeric' or 'character'
#' @param b vector, a set of elements, can be 'numeric' or 'character'
#' @param universe_size integer, total number of elements 
#' @export
#'
#' @return vector with four components (num_common, num_unique_a, num_unique_b, num_others)
contingency_vector = function(a, b, universe_size) {
  common = sum(a %in% b)
  a_only = sum(!(a %in% b))
  b_only = sum(!(b %in% a))
  c(common, a_only, b_only, universe_size-common-a_only-b_only)
}


#' construct a table with contingency values comparing many configuration of type (a vs. b)
#'
#' @param a sets a; if of type list, each component is interpreted as a set
#' @param b sets b to match a;
#' if type list, each component is interpreted as a set;
#' if a simple vector, each component is interpreted as an index into 'sets'
#' @param sets list with vectors, used when 'a' or 'b' are not lists
#' @param universe_size integer, the total number of elements that sets a, b were drawn from
#' @export
#'
#' @return data table with one row per item, column id and a,b,c,d
batch_contingency = function(a, b, universe_size, sets=NULL) {
  
  n = length(a)
  if (n != length(b)) {
    stop("incompatible components 'a' and 'b' - must be of equal length\n")
  }
  if (class(a) != class(b)) {
    stop("incompatible components 'a' and 'b' - must of the same class\n")
  }
  if (!class(a) %in% c("list", "character", "integer")) {
    stop("invalid inputs - must be in form of list or character/integer vectors\n")
  }
  
  universe_size = rep(universe_size, length.out=n)
  
  # compute the contingency tables
  result = base::vector("list", n)
  if (is.list(a) & is.list(b)) {
    # all configurations are provided explicitly
    for (i in seq_len(n)) {
      result[[i]] = contingency_vector(a[[i]], b[[i]], universe_size[i])
    }
  } else {
    # configuration are provided through sets
    for (i in seq_len(n)) {
      aset = sets[[a[i]]]
      bset = sets[[b[i]]]
      result[[i]] = contingency_vector(aset, bset, universe_size[i])
    }
  }
  # put into one large table - this relies on rbind preserving the order of a and b
  result = data.table(do.call(rbind, result))
  count.cols = c("count_11", "count_10", "count_01", "count_00")
  colnames(result) = count.cols

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


