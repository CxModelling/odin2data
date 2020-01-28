#' Title
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
