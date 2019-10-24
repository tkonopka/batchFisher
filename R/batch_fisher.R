# obtain fisher.test outputs for a batch of configurations
#


#' compute fisher comparisons of paired
#'
#' @param data data.table with contingency values as output by batchContingency
#' @param precomputed data table with precomputed p-values and odds-ratios
#' @export
#'
#' @return data table with contingency data and fisher output for each item in data
batch_fisher = function(data, precomputed=NULL) {

  count_columns = paste0("count_", c("11", "10", "01", "00"))

  if (!all(count_columns %in% colnames(data))) {
    stop("input data object has missing count columns\n")
  }
  if ("p.value" %in% colnames(data)) {
    stop("input already has a p.value column\n")
  }
  if ("odds.ratio" %in% colnames(data)) {
    stop("input already has an odds-ratio column\n")
  }

  # identify the types of configuration that have not already been computed
  if (is.null(precomputed)) {
    relevant = NULL
    missing = unique(data[, count_columns, with=FALSE])
  } else {
    # identify the configurations in the data that are not available as precomputed
    temp = merge(unique(data[, count_columns, with=FALSE]), precomputed,
                 by=count_columns, all.x=TRUE)
    missing = temp[is.na(p.value)]
    relevant = temp[!is.na(p.value)]
  }
  
  if (nrow(missing)>0) {
    newcomputed = fisher2x2(missing)
    relevant = rbind(relevant, newcomputed)
  }

  result = data
  result$.index = seq_along(result$count_11)
  result = merge(result, relevant, by=count_columns, all.x=TRUE)
  result = result[order(.index)]
  result$.index = NULL
  result
}




#' compute fisher p-values and odds-ratios from a table of contingency values
#'
#' @keywords internal
#' @param x data table with column $a, $b, $c $d
#'
#' @return an augmented data table with $p.value and $odds.ratio
fisher2x2 = function(x) {  
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
