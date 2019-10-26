# batchFisher

The Fisher exact test is routinely used to compare pairs of sets. In situations when many such tests must be performed in a systematic fashion, the calculation of Fisher p-values can be a bottleneck. The package `batchFisher` streamlines such batch calculations. It provides substantial performance gains in situations when many tests in the batch rely on the same contingency matrix.




## Strategy

The package uses the following optimization strategies:

 - Large batches of set comparisons can contain several instances described by the same contingency matrices. The package can identify the unique configurations and perform the p-value and odds-ratios calculation fewer times.
 - Calculations that handle several batches can encounter the same contingency matrices in distinct batches. The package provides a mechanism to precompute outputs for certain configuration and then re-use the results when appropriate.
 - Batches can include comparisons of one set against many others. The repeated set can be re-used for many comparisons instead of fetching it anew every time.
 - Calculations that begin with sets must compute contingency matrices before evaluating the Fisher test. Out of several possible implementations (including Rcpp), the package uses an implementation in base R that works well for sets of moderate size.


The above techniques exploit properties of the batches to avoid repetition. The results are themselves computed using `fisher.test`. It is in principle possible to optimize the calculations that take place *within* `fisher.test`. Below is a list of possibilities, but these are *not* used in the package.

 - The calculation of Fisher p-values involves evaluating many combinations - in the mathematical sense, e.g. `choose(5,2) = 10`. Dynamic programming is likely used within `fisher.test`, but cached values are not re-used on subsequent invokations of `fisher.test`. A re-implementation from scratch might provide some speedup.
 - Function `fisher.test` provides other information in addition to p-values and odds ratios. While convenient during interactive use, these additional data are typically not relevant to batch calculations.
 - P-values for contingency matrices with large entries might be well-approximated by alternative statistical tests, i.e. chi-squared. 




## Performance

An example calculation is presented in the vignette. The summary graph shows performance as a function of batch size (number of tests performed). Performance is measured by total running time as well as by the number of tests performed per second.

<img src="https://github.com/tkonopka/batchFisher/blob/master/images/readme_simulation.png?raw=true" alt="Running times and operations per second"></img>

The two methods provided by the `batchFisher` package are denoted as 'batch' and 'batch with precomp.`. They differ by whether or not they use a table or pre-calculated results. These methods are compared against a third approach denoted as 'simple', which computes one test at a time. 



## Limitations

- While the Fisher test is a general procedure that works on contingency tables of any size, the package only supports 2x2 matrices.


