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
calc_dist.sim_model <- function(sim, pars, dat) {
  lf <- compile_model_likefree(dat, sim)
  return(calc_dist(lf, pars))
}


#' @rdname calc_dist
#' @export
calc_dist.sim_model_likefree <- function(lf, pars) {
  rs <- simulate(lf, pars = pars, warmup = lf$Model$WarmupStage == "Yes")
  return(calc_dist(rs, lf))
}


#' @rdname calc_dist
#' @export
calc_dist.sim_results <- function(rs, lf) {
  simulated <- rs$Ys
  rownames(simulated) <- simulated[, 1]
  simulated <- simulated[lf$Ts2fit, lf$Cols2fit]

  dist <- - sum(dnorm(simulated, lf$Data[, -1], log=T), na.rm = T)
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
