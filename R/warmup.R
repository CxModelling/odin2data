#' Warming-up a model to a steady state
#'
#' @param sim
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
warmup <- function(sim, ...) {
  UseMethod("warmup", sim)
}


#' @rdname warmup
#' @export
warmup.likefree_model <- function(sim, y0=NA, pars, t_warmup=NA) {
  if (all(is.na(y0))) y0 <- sim$y0

  if (is.na(sim$cm_warmup)) {
    return(y0)
  }

  if (!is.na(t_warmup)) {
    ts_warmup <- seq(min(sim$ts_sim) - t_warmup, min(sim$ts_sim), by=1)
  } else {
    ts_warmup <- sim$ts_warmup
  }

  env <- pars
  env$y0 <- y0
  m <- sim$cm_warmup
  m$set_user(user=env)
  ys <- m$run(ts_warmup)
  y_eq <- ys[nrow(ys), 1:length(y0) + 1]

  return(y_eq)
}
