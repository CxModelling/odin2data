#' Approximate Bayesian Computation (ABC) for System dynamic models
#'
#' @param lf a "model_likefree" object with a compile model and data
#' @param n_posterior number of posteriors to collect
#' @param epsilon threshold for collecting
#' @param keep simulation results to keep; "Y0": initial values, "Ys": trajectories, "both": both, "none": none
#' @param target_acc target acceptance rate
#' @param n_test number of simulations for testing threshold
#'
#' @return
#' @export
#'
#' @examples
#' library(odin2data)
#'
#' ## Example 1: model without warmup stage
#' test_data <- data.frame(
#'   t = 1:5,
#'   incidence = c(20, 49, 109, 184, 206) / 1000
#' )
#'
#' ### Load the model file
#' f <- system.file("models/SIR.txt", package = "odin2data")
#' test_m <- odin::odin(f, verbose=F)
#'
#'
#' ### Generate a prior and set up it as a list
#' r_prior <- function() {
#'   list(beta = runif(1, 1, 10), gamma = runif(1, .1, 1))
#' }
#'
#' d_prior <- function(pars) {
#'   dunif(pars$beta, 1, 10, log = T) + dunif(pars$gamma, .1, 1, log = T)
#' }
#'
#' times = seq(0, 10, 0.2)
#' y0 <- c(995, 5, 0)
#'
#' ### Compile all elements as a simulation model
#' sim <- odin2data::compile_model(d_prior, r_prior, y0, ts_sim = times, m_sim = test_m)
#'
#' ### Compile the model with data
#' lf <- odin2data::compile_model_likefree(test_data, sim)
#'
#' ### Model fitting
#' post <- odin2data::fit_abc(lf, 100, target_acc = 0.2)
#' summary(post)
#'
#'
#' ## Example 2: model with warmup stage
#' fn_pass_y0 <- function(ys) {
#' ys[nrow(ys), 2:4]
#' }
#'
#' test_data <- data.frame(
#'   t = 1:5,
#'   incidence = rep(0.009580023, 5)
#' )
#'
#' sim <- odin2data::compile_model(d_prior, r_prior, y0, ts_sim = times, m_sim = test_m,
#'                                 m_wp = test_m, t_wp = 100, fn_pass_y0 = fn_pass_y0)
#'
#' lf <- odin2data::compile_model_likefree(test_data, sim)
#'
#' ### Model fitting
#' post <- odin2data::fit_abc(lf, 100, target_acc = 0.2)
#' summary(post)
#'
fit_abc <- function(lf, n_posterior = 300, epsilon = NA, keep = c("Y0", "Ys", "both", "none"),
                    target_acc = 0.1, n_test = 100) {
  sim <- lf$Model

  keep <- match.arg(keep)
  keep_y0 <- keep %in% c("Y0", "both")
  keep_ys <- keep %in% c("Ys", "both")

  # Finding epsilon --------------------------
  if (is.na(epsilon)) {
    cat("Test threshold\n")
    stopifnot(target_acc < 1 & target_acc > 0)

    ds <- rep(0, n_test)
    for (i in 1:n_test) {
      pars <- sim$r_prior()
      ds[i] <- tryCatch({ calc_dist(lf, pars) }, error = function(e) Inf )
    }

    epsilon <- as.numeric(quantile(ds, target_acc))
    stopifnot(is.finite(epsilon))
  }

  cat("Collect posteriors\n")
  # Collect posterior ------------------------
  n_collected <- 0
  n_run <- 0
  posteriors <- list()

  while(n_collected < n_posterior) {
    pars <- sim$r_prior()

    res <- tryCatch({
      simulate(lf, pars = pars)
    }, error = function(e) "Fail")

    if (!is.list(res)) {
      n_run <- n_run + 1
      next
    }
    dist <- calc_dist(res, lf)
    n_run <- n_run + 1

    if (dist < epsilon) {
      n_collected <- n_collected + 1

      posterior <- list(
        parameters = pars,
        distance = dist
      )

      if (keep_ys) {
        posterior$Ys <- res$Ys
      }
      if (keep_y0) {
        posterior$Y0 <- res$Y0
      }
      posteriors[[n_collected]] <- posterior
    }
  }

  lis <- - sapply(posteriors, function(x) x$distance)

  # Collect meta data ------------------------
  meta <- list(
    n_iter = n_posterior,
    n_run = n_run,
    epsilon = epsilon,
    p_acc = n_posterior / n_run,
    ess = ess(lis)
  )

  pss <- sapply(posteriors, function(x) unlist(x$parameters))
  pss <- t(pss)

  res <- list(
    Posteriors = posteriors,
    Parameters = pss,
    Meta = meta,
    Data = lf$Data
  )
  class(res) <- "posterior_likefree"

  return(res)
}
