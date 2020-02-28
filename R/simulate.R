#' Simulate a single run
#'
#' @param sim
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
simulate.sim_model <- function(sim, y0, pars, warmup = F) {
  if (missing(pars)) {
    pars <- sim$r_prior()
  }

  if (warmup & sim$WarmupStage == "Yes") {
    if (missing(y0)) {
      y0 <- sim$Y0_wp
    }
    temp <- warmup(sim, y0 = y0, pars = pars)
    if (class(temp) != "Y_eq") return("Eq check failed")

    pars <- temp$Parameters
    y0 <- temp$Y0
  } else {
    if (missing(y0)) {
      y0 <- sim$Y0_sim
    }
  }

  inp <- pars
  inp$Y0 <- y0

  cm_sim <- sim$CM_sim
  cm_sim$set_user(user = inp)

  st <- system.time({ ys <- cm_sim$run(sim$TS_sim) })

  res <- list(
    Y0 = y0,
    Parameters = pars,
    Ys = ys,
    ProcTime = st
  )
  class(res) <- "sim_results"
  return(res)
}


#' @rdname simulate.sim_model
#' @export
simulate.sim_model_likefree <- function(sim, y0, pars, warmup = F) {
  return(simulate(sim$Model, y0, pars, warmup))
}
