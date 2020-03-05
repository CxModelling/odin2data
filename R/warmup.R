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
warmup.sim_model <- function(sim, y0, pars, times = sim$TS_wp) {
  stopifnot(sim$WarmupStage == "Yes")

  if (missing(pars)) {
    pars <- sim$r_prior()
  }

  if (missing(y0)) {
    y0 <- sim$Y0_wp
  }

  inp <- pars
  inp$Y0 <- y0

  cm_wp <- sim$CM_wp
  cm_wp$set_user(user = inp)

  ys <- cm_wp$run(times, method = sim$Method)

  if ("Checker" %in% names(sim)) {
    stopifnot(sim$Checker(ys))
  }

  if (is.array(y0)) {
    dim0 <- dim(y0)
    y1 <- array(ys[nrow(ys), 1 + 1:prod(dim0)], dim0)
  } else if (is.vector(y0)) {
    y1 <- ys[nrow(ys), 1 + 1:length(y0)]
  }


  y0 <- sim$Linker(ys)

  res <- list(
    Yend = y1,
    Y0 = y0,
    Parameters = pars
  )

  class(res) <- "Y_eq"

  return(res)
}


#' @rdname warmup
#' @export
update.Y_eq <- function(yeq, sim, nforward = length(sim$TS_wp)) {
  times <- sim$Time_wp[2] - nforward:0
  return(warmup(sim, y0 = yeq$Yend, pars = yeq$Parameters, times = times))
}
