
#' Projecting from a simulation model
#'
#' @param sim
#' @param y_eq
#' @param pars
#' @param times
#' @param ...
#'
#' @return a data.frame of simulated data
#' @export
#'
#' @examples
project <- function(sim, y_eq, pars, ...) {
  UseMethod("project", sim)
}


#' @rdname project
#' @export
project.likefree_model <- function(sim, y_eq, pars, times=NA) {
  if (any(is.na(times))) times <- sim$ts_sim

  env <- pars
  env$y0 <- y_eq
  m <- sim$cm_sim
  m$set_user(user=env)
  ys <- m$run(times)

  return (ys)
}
