#' Fit simulation model to data with Likelihood-free methods
#'
#' @param lf
#' @param n_posterior
#' @param method
#' @param ... arguments for a specific method
#'
#' @return
#' @export
#'
#' @examples
fit <- function(lf, n_posterior, method = c("abc", "abcsmc", "abcpmc"), ...) {
  method <- match.arg(method)

  alg <- switch(method,
    abc = fit_abc, abcsmc = fit_abc_smc, abcpmc = fit_abc_pmc
  )
  post <- alg(lf, n_posterior = n_posterior, ...)
  return(post)
}
