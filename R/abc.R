#' Approximate Bayesian Computation (ABC) for System dynamic models
#'
#' @param sim
#' @param n_posterior
#' @param epsilon
#' @param keep_all
#' @param target_acc
#'
#' @return
#' @export
#'
#' @examples
#' library(odin2data)
#'
#' test_data <- data.frame(
#'   t = 1:5,
#'   incidence = c(20, 49, 109, 184, 206) / 1000
#' )
#'
#'
#' test_m <- odin::odin({
#'   deriv(sus) <- - foi
#'   deriv(inf) <- foi - gamma * inf
#'   deriv(rec) <- gamma * inf
#'
#'   initial(sus) <- y0[1]
#'   initial(inf) <- y0[2]
#'   initial(rec) <- y0[3]
#'
#'   output(incidence) <- foi / n
#'
#'   # parameters
#'   n <- sus + inf + rec
#'   beta <- user(1.5)
#'   gamma <- user(0.5)
#'
#'   foi <- beta * sus * inf / n
#'
#'   # data
#'   y0[] <- user()
#'   dim(y0) <- 3
#' }, verbose=F)
#'
#'
#' r_prior <- function() {
#'   list(
#'     beta = runif(1, 1, 10),
#'     gamma = runif(1, .1, 1)
#'   )
#' }
#'
#' d_prior <- function(pars) {
#'   dunif(pars$beta, 1, 10, log = T) + dunif(pars$gamma, .1, 1, log = T)
#' }
#'
#'
#' times = seq(0, 10, 0.2)
#' y0 <- c(995, 5, 0)
#'
#' sim <- compile_likefree_model(test_data, test_m, y0, rprior = r_prior, dprior = d_prior, times = times)
#'
#' fitted <- fit_abc(sim)
#'
#' summary(fitted)
fit_abc <- function(sim, n_posterior = 300, epsilon = NA, keep_all = FALSE, target_acc = 0.05) {
  # Finding epsilon --------------------------
  if (is.na(epsilon)) {
    ds <- rep(0, 300)
    for (i in 1:300) {
      pars <- sim$r_prior()
      ds[i] <- calc_dist(sim, pars)
    }

    stopifnot(target_acc < 1 & target_acc > 0)
    epsilon <- as.numeric(quantile(ds, target_acc))
  }

  # Collect posterior ------------------------
  n_collected <- 0
  n_run <- 0
  posteriors <- list()

  while(n_collected < n_posterior) {
    pars <- sim$r_prior()
    res <- simulate(sim, pars = pars)
    dist <- calc_dist(res)
    n_run <- n_run + 1

    if (dist < epsilon) {
      n_collected <- n_collected + 1

      posterior <- list(
        parameters = pars,
        distance = dist
      )

      if (keep_all) {
        posterior$ys <- res$ys
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
    posteriors = posteriors,
    parameters = pss,
    meta = meta
  )
  class(res) <- "fitted_abc"

  return(res)
}


#' @rdname fit_abc
#' @export
summary.fitted_abc <- function(fitted) {
  ans <- fitted$meta

  su <- apply(fitted$parameters, 2, function(ps) {
    c(
      mean(ps), sd(ps), min(ps),
      quantile(ps, c(0.025, 0.25, 0.5, 0.75, 0.975)),
      max(ps)
    )
  })

  su <- t(su)

  colnames(su)[c(1:3, 9)] <- c("Mean", "Std", "Min", "Max")

  ans$Post <- su
  class(ans) <- "summary.abc"
  return(ans)
}
