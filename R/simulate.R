#' Simulate a system dynamic
#'
#' @param sim
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
simulate.likefree_model <- function(sim, y0=NA, pars=NA, times=NA) {
  if(all(is.na(pars))) pars <- sim$r_prior()

  y_eq <- warmup(sim, y0, pars)

  ys <- project(sim, y_eq, pars, times)

  res <- list(
    model = sim,
    y_eq = y_eq,
    pars = pars,
    ys = ys
  )
  class(res) <- "sim_results"
  return(res)
}
