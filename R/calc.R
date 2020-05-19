#' Calculate the distance between data and simulated dynamics
#'
#' @param sim a sim_model object
#' @param lf a sim_model_likefree object
#' @param lik a sim_model_likelihood object
#' @param dat a data for composing likelihood-free distance
#' @param pars parameter set
#'
#' @return
#' @export
#'
#' @examples
calc_dist <- function(sim, pars) {
  UseMethod("calc_dist", sim)
}


#' @rdname calc_dist
#' @export
calc_dist.sim_model <- function(sim, pars, dat) {
  if (is.infinite(sim$d_prior(pars))) return(Inf)

  if (is.function(dat)) {
    lf <- compile_model_likelihood(dat, sim)
  } else {
    lf <- compile_model_likefree(dat, sim)
  }
  return(calc_dist(lf, pars))
}


#' @rdname calc_dist
#' @export
calc_dist.sim_model_likefree <- function(lf, pars) {
  rs <- tryCatch({
    simulate(lf, pars = pars, warmup = lf$Model$WarmupStage == "Yes")
  }, error = function(e) e$message)
  if (class(rs) != "sim_results") return(Inf)
  return(calc_dist(rs, lf))
}


#' @rdname calc_dist
#' @export
calc_dist.sim_model_likelihood <- function(lik, pars) {
  rs <- tryCatch({
    simulate(lik, pars = pars, warmup = lik$Model$WarmupStage == "Yes")
  }, error = function(e) e$message)
  if (class(rs) != "sim_results") return(Inf)
  return(-lik$FnLike(rs$Ys))
}


#' @rdname calc_dist
#' @export
calc_dist.sim_results <- function(rs, lf) {
  simulated <- rs$Ys
  if (class(lf)[1] == "sim_model_likefree") {
    sel <- simulated[, 1] %in% lf$Ts2fit
    rows <- simulated[sel, 1]
    simulated <- simulated[sel, lf$Cols2fit]
    dat <- lf$Data[lf$Data[, 1] %in% rows, -1]
    dist <- sum((simulated - dat) ^ 2, na.rm = T)
  } else {
    dist <- -lf$FnLike(simulated)
  }
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
