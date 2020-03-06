find_eps <- function(ds, eps0, alpha) {
  e0 <- mean(ds < eps0)
  et <- alpha * e0
  for (eps1 in rev(sort(ds))) {
    e1 <- mean(ds < eps1)
    if (e1 <= et) {
      return (eps1)
    }
  }
  return (eps0)
}

#' Approximate Bayesian Computation Sequential Monte Carlo
#'
#' @param lf a sim_model_likefree object
#' @param n_posterior number of posterior particles to collect
#' @param alpha percentage particles surviving
#' @param p_thres proportion of threshold population
#' @param keep keep simulation (Ys), intial values (Y0), or both
#' @param max_round max round, reach to stop
#' @param max_stay max staying round of epsilon, reach to stop
#' @param verbose show log or not
#'
#' @references Del Moral, P., Doucet, A., & Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation. Statistics and Computing, 22(5), 1009-1020.
#' @export
#'
#' @examples
fit_abc_smc <- function(lf, n_posterior, alpha = 0.9, p_thres = 0.6, max_round = 20, max_stay = 3, verbose = T) {
  ess_thres <- n_posterior * p_thres

  # initialization
  if (verbose) {
    cat("Initialise\n")
  }
  posteriors <- list()
  theta_0 <- list()
  d_0 <- rep(0, n_posterior)

  for (i in 1:n_posterior) {
    d <- Inf
    while(is.infinite(d)) {
      rs <- a_sample(sim)

      if (class(rs) != "sim_results") next
      d <- calc_dist(rs, lf)
    }

    posteriors[[i]] <- rs
    theta_0[[i]] <- rs$Parameters
    d_0[i] <- d
  }

  wts <- rep(1, n_posterior) / n_posterior

  eps_0 <- Inf
  n_round <- 0
  n_stay <- 0
  n_eval <- n_posterior

  traj <- c(Round = n_round, Eval = n_eval, Eps = eps_0, ESS = 1 / sum(wts ^ 2), Acc = 1)

  if (verbose) {
    cat("Round ", n_round, ", ESS ", round(1 / sum(wts ^ 2), 1), ", Epsilon ", eps_0, "\n")
  }

  while (TRUE) {
    n_round <- n_round + 1
    eps_1 <- find_eps(d_0, eps_0, alpha)
    if (eps_1 > eps_0) {
      n_stay <- n_stay + 1
      eps_1 <- eps_0
    } else {
      n_stay <- 0
    }

    # Step 1 updating weight
    act_np_0 <- (d_0 < eps_0) + 0
    act_np_1 <- (d_0 < eps_1) + 0

    a <- act_np_0 > 0
    wts[a] <- wts[a] * act_np_1[a] / act_np_0[a]
    wts[act_np_0 == 0] <- 0
    wts <- wts / sum(wts)

    # Step 2 resampling
    if (ess_thres * sum(wts ^2) > 1) {
      stopifnot(sum(wts > 0) > 2)
      alive <- wts > 0
      ind <- (1:n_posterior)[alive]
      re_index <- sample(ind, n_posterior, prob = wts[alive], replace = T)
      theta_1 <- theta_0[re_index]

      posteriors <- posteriors[re_index]

      d_1 <- d_0[re_index]
      wts <- rep(1, n_posterior) / n_posterior
    } else {
      theta_1 <- theta_0
      d_1 <- d_0
    }

    # Step 3 MCMC proposal
    act_np_0 <- (d_0 < eps_1) + 0
    alive <- act_np_0 > 0

    ## MH step
    tau <- apply(sapply(theta_1, unlist), 1, weighted.sd, wt = wts)

    theta_p <- list()
    y_p <- posteriors
    d_p <- rep(0, n_posterior)

    for (i in 1:n_posterior) {
      rs <- mutate_sample(sim, posteriors[[i]], tau = tau)

      theta_p[[i]] <- rs$Parameters
      d_p[i] <- d <- calc_dist(rs, lf)

      rs$Distance <- d

      y_p[[i]] <- rs
      n_eval <- n_eval + rs$N_attempt
    }

    act_np_p <- (d_p < eps_1) + 0

    ## MH acceptance ratio
    acc <- rep(0, n_posterior)
    acc[alive] <- act_np_p[alive] / act_np_0[alive]
    acc[is.infinite(d_p)] <- 0

    ## Update accepted proposals
    a <- runif(n_posterior) < acc
    if (sum(a) > 2) {
      theta_1[a] <- theta_p[a]
      d_1[a] <- d_p[a]
      posteriors[a] <- y_p[a]
    }

    # Collect history
    traj <- rbind(traj,
                  c(Round = n_round, Eval = n_eval, Eps = eps_1,
                    ESS = 1 / sum(wts ^ 2), Acc = sum(a) / sum(alive)))

    if (verbose) {
      cat("Round ", n_round, ", ESS ", round(1 / sum(wts ^ 2), 1), ", Epsilon ", eps_1,
          ", Acceptance ", round(100 * sum(a) / sum(alive), 2), "%\n")
    }

    if (n_round >= max_round | (n_stay >= max_stay & n_round > 3)) {
      break
    }

    theta_0 <- theta_1
    d_0 <- d_1
    eps_0 <- eps_1
  }

  alive <- wts > 0
  ind <- (1:n_posterior)[alive]
  re_index <- sample(ind, n_posterior, prob = wts[alive], replace = T)


  # Collect meta data ------------------------
  meta <- list(
    n_iter = n_posterior,
    epsilon = eps_1,
    ess = 1 / sum(wts ^ 2),
    trajectory = traj
  )

  pss <- sapply(posteriors, function(x) unlist(x$Parameters))
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
