#' Calculate the distance between data and simulated dynamics
#'
#' @param sim
#' @param pars
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calc_dist <- function(sim, pars, ...) {
  UseMethod("calc_dist", sim)
}


#' @rdname calc_dist
#' @export
calc_dist.likefree_model <- function(sim, pars) {
  res <- simulate(sim, pars=pars)
  return(calc_dist(res))
}


#' @rdname calc_dist
#' @export
calc_dist.sim_results <- function(sim) {
  simulated <- sim$ys
  rownames(simulated) <- simulated[, 1]
  simulated <- simulated[res$model$ts_to_fit, res$model$cols_to_fit]

  dist <- - sum(dnorm(simulated, sim$model$data[, -1], log=T), na.rm = T)
  return(dist)
}


#' Effective Sample Size (ESS)
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
ess <- function(obj) {
  UseMethod("ess", obj)
}


#' @rdname ess
#' @export
ess.numeric <- function(obj) {
  wts <- exp(obj - lse(obj, na.rm = T))
  return(sum(wts) ^ 2 /sum(wts^2))
}


#' @rdname ess
#' @export
ess.fitted_abc <- function(obj) {
  return (obj$meta$ess)
}
