#' LogSumExp function
#'
#' @param xs a vector of values to be summed up in exp scale
#'
#' @return a value
#' @export
#'
#' @examples
#' xs <- c(1000, 1001, 999)
#'
#' ## go infinity because of numerical overflow
#' is.infinite(log(sum(exp(x))))
#'
#' ## Use lse function
#' lse(xs)
lse <- function(xs, na.rm=FALSE) {
  m <- max(xs)
  return(log(sum(exp(xs - m), na.rm = na.rm)) + m)
}


weighted.sd <- function(x, wt) {
  wt <- wt / sum(wt)
  mu <- weighted.mean(x, wt)
  tau <- sqrt(sum(wt * (x - mu) ^ 2))
  tau
}
