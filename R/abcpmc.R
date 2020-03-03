weighted.sd <- function(x, wt) {
  wt <- wt / sum(wt)
  mu <- weighted.mean(x, wt)
  tau <- sqrt(sum(wt * (x - mu) ^ 2))
  tau
}

#' Title
#'
#' @param lf
#' @param n_posterior
#' @param k
#' @param keep
#' @param max_round
#' @param q_stop
#' @param verbose
#'
#' @references Simola, U., Cisewski-Kehe, J., Gutmann, M. U., & Corander, J. (2019). Adaptive Approximate Bayesian Computation Tolerance Selection. arXiv preprint arXiv:1907.01505.
#' @return
#' @export
#'
#' @examples
fit_abc_pmc <- function(lf, n_posterior, k = 5, keep = c("Y0", "Ys", "both", "none"), max_round = 20, q_stop = 0.99, verbose = T) {
  sim <- lf$Model

  keep <- match.arg(keep)
  keep_y0 <- keep %in% c("Y0", "both")
  keep_ys <- keep %in% c("Ys", "both")

  ##### Initialisation -----
  n_init <- n_posterior * k
  posteriors <- list()

  ##### Step 0 -----
  if (verbose) {
    cat("Initialising\n")
  }
  theta0 <- lapply(1:n_init, function(x) sim$r_prior())
  ds0 <- sapply(theta0, function(x) calc_dist(lf, pars = x))
  wts0 <- rep(1, n_posterior) / n_posterior

  eps <- sort(ds0)[n_posterior + 1]
  ess <- n_posterior

  theta0 <- theta0[order(ds0) <= n_posterior]
  ds0 <- ds0[order(ds0) <= n_posterior]

  traj <- c(round = 0, ess = ess, eps = eps, qt = 1/k, acc = 1/k, n = n_init)

  if (verbose) {
    cat("Round ", 0, ", ESS ", round(ess, 1), ", Epsilon ", eps, ", Acceptance ", round(100 * 1/k, 2), "%\n")
  }

  n_round <- 1
  qt = 1/k

  while(n_round <= max_round) {
    ##### Step 1+ -----
    tau <- apply(sapply(theta0, unlist), 1, weighted.sd, wt = wts0)

    theta1 <- list()
    ds1 <- rep(0, n_posterior)
    wts1 <- rep(0, n_posterior)

    posteriors <- list()
    n_try <- 0

    for (j in 1:n_posterior) {
      dj <- eps + 1
      dpr <- -Inf
      while (dj > eps) {
        n_try <- n_try + 1
        the <- theta0[[sample.int(n_posterior, 1, prob = wts0)]]
        the_prop <- as.list(rnorm(length(the), unlist(the), tau))
        names(the_prop) <- names(the)
        dpr <- sim$d_prior(the_prop)
        if (is.finite(dpr)) {
          y_prop <- simulate(sim, pars = the_prop)
          dj <- calc_dist(y_prop, lf)
        }
      }

      post <- list(Parameters = y_prop$Parameters, Distance = dj)
      if (keep_ys) {
        post$Ys <- y_prop$Ys
      } else {
        post$Y0 <- y_prop$Y0
      }
      posteriors[[j]] <- post
      theta1[[j]] <- the_prop

      phi <- sapply(theta0, function(x) sum(dnorm(unlist(x), mean = unlist(the_prop), sd = tau, log = T)))
      ds1[j] <- dj
      wts1[j] <- dpr - lse(log(wts0) + phi)
    }
    wts1 <- exp(wts1 - lse(wts1))

    ess <- sum(wts1) ^ 2 / sum(wts1^2)
    acc <- n_posterior / n_try
    traj <- rbind(traj, c(round = n_round, ess = ess, eps = eps, qt = qt, acc = acc, n = n_try))

    if (verbose) {
      cat("Round ", n_round, ", ESS ", round(ess, 1), ", Epsilon ", eps,
          ", Acceptance ", round(100 * acc, 2), "%\n")
    }
    n_round <- n_round + 1

    x1 <- t(sapply(theta1, unlist))
    dens <- densratio::densratio(x1, t(sapply(theta0, unlist)),
                                 method = "KLIEP", kernel_num = 100, verbose = F)
    ct <- max(max(dens$compute_density_ratio(x1)), 1)
    qt <- 1/ct
    eps <- quantile(ds1, qt)

    ##### Go to next round -----
    theta0 <- theta1
    ds0 <- ds1
    wts0 <- wts1 / sum(wts1)

    if (qt > q_stop & n_round >= 3) {
      break
    }
  }

  # Collect meta data ------------------------
  meta <- list(
    n_iter = n_posterior,
    epsilon = eps,
    ess = ess,
    trajectory = traj
  )

  pss <- sapply(posteriors, function(x) unlist(x$Parameters))
  pss <- t(pss)

  res <- list(
    posteriors = posteriors,
    parameters = pss,
    meta = meta
  )
  class(res) <- "fitted_abcpmc"

  return(res)
}


#' @rdname fit_abc_pmc
#' @export
summary.fitted_abcpmc <- function(fitted) {
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
  class(ans) <- "summary.abcpmc"
  return(ans)
}
