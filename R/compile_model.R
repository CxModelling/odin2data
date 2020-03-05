#' @rdname compile_model
#' @export
test_model <- function(d_prior, r_prior, y0, inp_sim, ts_sim, m_sim, inp_wp, t_wp, m_wp, fn_pass_y0, fn_check, method) {
  res <- list()

  fp <- r_prior()
  res$FreeParameters <- fp

  stopifnot(is.finite(d_prior(fp)))

  res$d_prior <- d_prior
  res$r_prior <- r_prior


  if (missing(t_wp) | missing(m_wp) | missing(fn_pass_y0)) {
    res$WarmupStage <- "No"
  } else {
    res$WarmupStage <- "Yes"

    pars <- fp
    if (!missing(inp_wp) & !any(is.null(inp_wp))) {
      pars <- c(inp_wp, pars)
    } else {
      inp_wp <- NULL
    }

    pars$Y0 <- y0

    cm_wp <- m_wp(user = pars)
    ts_wp <- ts_sim[1] - (t_wp:0)

    st <- system.time({ ys0 <- cm_wp$run(ts_wp) })
    cat("Warm-up time:\n")
    print(st)

    ys0 <- ys0[ts_wp == round(ts_wp), ]

    res <- c(res, list(
      Input_wp = inp_wp,
      Y0_wp = y0,
      Ys_wp = ys0,
      CM_wp = cm_wp,
      Time_wp = range(ts_wp),
      TS_wp = ts_wp
    ))

    if(!missing(fn_check)) {
      stopifnot(fn_check(ys0))
      res$Checker <- fn_check
    }
    y0 <- fn_pass_y0(ys0)
    res$Linker = fn_pass_y0
  }

  pars <- fp
  if (!missing(inp_sim) & !any(is.null(inp_sim))) {
    pars <- c(inp_sim, pars)
  } else {
    inp_sim <- NULL
  }
  pars$Y0 <- y0

  cm_sim <- m_sim(user = pars)

  if (is.null(method)) method = "rk4"
  st <- system.time({ ys1 <- cm_sim$run(ts_sim, method = method) })
  cat("Simulation time:\n")
  print(st)

  ys1 <- ys1[ts_sim == round(ts_sim),]

  res <- c(res, list(
    Input_sim = inp_sim,
    Y0_sim = y0,
    Ys_sim = ys1,
    CM_sim = cm_sim,
    Time_sim = range(ts_sim),
    TS_sim = ts_sim,
    Method = method
  ))
  return(res)
}


#' Compile a simulation model given all components
#'
#' @param d_prior a probability density function supporting prior distributions
#' @param r_prior a function for generating parameters from their prior distributions
#' @param y0 initial values
#' @param inp_sim input data for the simulation stage
#' @param ts_sim timespan for simulation
#' @param m_sim an odin model for simulation
#' @param inp_wp input data for the warm-up stage
#' @param t_wp length of warm-up stage
#' @param m_wp an odin model for waru-up
#' @param fn_pass_y0 a function for bringing the states at the end of warm-up to simulation initials
#' @param fn_check a function for checking if a parameter set can generate validated output
#'
#' @return
#' @export
#'
#' @examples
compile_model <- function(d_prior, r_prior, y0, inp_sim = NULL, ts_sim, m_sim, inp_wp = NULL, t_wp, m_wp = m_sim,
                          fn_pass_y0, fn_check, method = NULL, max_attempt = 10) {
  n_attempt <- 0

  if (missing(fn_check)) {
    fn_check <- function(x) { T }
  }

  while(T) {
    tested <- tryCatch({
      test_model(d_prior, r_prior, y0, inp_sim, ts_sim, m_sim, inp_wp, t_wp, m_wp, fn_pass_y0, fn_check, method)
    }, error = function(e) e$message)


    n_attempt <- n_attempt + 1
    if (is.list(tested)) break

    stopifnot(n_attempt < max_attempt)
  }
  class(tested) <- "sim_model"
  return(tested)
}


#' Compile a simulation model with a likelihood-free link to data
#'
#' @param dat a dataframe of data to be fitted, "t" as the indicator of time
#' @param sim a compile model, see compile_model
#' @param y0 initial values of the model
#'
#' @return
#' @export
#'
#' @examples
compile_model_likefree <- function(dat, sim) {

  vars <- intersect(colnames(sim$Ys_sim), colnames(dat))
  vars <- vars[vars != "t"]

  res <- list(
    Data = dat,
    Model = sim,
    Cols2fit = vars,
    Ts2fit = dat[, "t"]
  )

  class(res) <- "sim_model_likefree"
  return(res)
}
