#' Simulate a single run
#'
#' @param sim simulation model, sim_model
#' @param y0 initial values
#' @param pars parameters in a list
#'
#' @return
#' @export
#'
#' @examples
simulate.sim_model <- function(sim, y0, pars, warmup = T) {
  if (missing(pars)) {
    pars <- sim$r_prior()
  }

  stopifnot(is.finite(sim$d_prior(pars)))

  if (warmup & sim$WarmupStage == "Yes") {
    if (missing(y0)) {
      y0 <- sim$Y0_wp
    }
    temp <- tryCatch({ warmup(sim, y0 = y0, pars = pars) }, error = function(e) e$message)
    stopifnot(is.list(temp))

    pars <- temp$Parameters
    y0 <- temp$Y0
  } else if (missing(y0)) {
    y0 <- sim$Y0_sim
  }

  inp <- sim$r_sto_sim(pars, sim$Input_sim)
  inp$Y0 <- y0

  cm_sim <- sim$CM_sim
  cm_sim$set_user(user = inp)

  st <- system.time({ ys <- cm_sim$run(sim$TS_sim, method = sim$Method) })

  res <- list(
    Y0 = y0,
    Parameters = pars,
    Ys = ys,
    ProcTime = st
  )
  if (warmup & sim$WarmupStage == "Yes") {
    res$Y1wp <- temp$Yend
  }

  class(res) <- "sim_results"
  return(res)
}


#' @rdname simulate.sim_model
#' @export
simulate.sim_model_likefree <- function(sim, y0, pars, warmup = T) {
  return(simulate(sim$Model, y0, pars, warmup))
}


#' @rdname simulate.sim_model
#' @export
simulate.sim_model_likelihood <- function(sim, y0, pars, warmup = T) {
  return(simulate(sim$Model, y0, pars, warmup))
}


a_sample <- function(sim, max_attempt = 100) {
  n_attempt <- 1
  while (T) {
    rs <- tryCatch({
      simulate(sim)
    }, error = function(e) e$message)

    if ("sim_results" %in% class(rs)) {
      rs$N_attempt <- n_attempt
      return (rs)
    } else {
      n_attempt <- n_attempt + 1
    }
    stopifnot(n_attempt <= max_attempt)
  }
}


mutate_sample <- function(sim, ptc, tau, max_attempt = 100, max_move = 1000) {
  n_attempt <- 0
  n_move <- 0

  pars <- ptc$Parameters
  if (sim$WarmupStage == "Yes") {
    y0 <- ptc$Y1wp
  } else {
    y0 <- sim$Y0_sim
  }

  while (T) {
    n_move <- n_move + 1
    prop <- lapply(1:length(pars), function(i) {
      p <- pars[[i]]
      sig <- tau[i]
      p + rnorm(length(p), 0, sig)
    })
    names(prop) <- names(pars)

    if (is.infinite(sim$d_prior(prop))) {
      next
    }

    n_attempt <- n_attempt + 1
    rs <- tryCatch({
      simulate(sim, y0 = y0, pars = prop)
    }, error = function(e) e$message)

    if ("sim_results" %in% class(rs)) {
      rs$N_attempt <- n_attempt
      return (rs)
    }
    stopifnot(n_attempt <= max_attempt)
    stopifnot(n_move <= max_move)
  }
}
