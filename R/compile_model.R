#' Compiling a likelihood-free simulation model
#'
#' @param dat a dataframe of data to be fitted, "t" as the indicator of time
#' @param m_sim an odin model
#' @param y0 initial values of the model
#' @param rprior a function for simulating prior parameters
#' @param dprior a function for calulating prior probability in log scale
#' @param times a vectors of times for projecting resutls
#' @param m_warmup a odin model for warm-up period
#' @param t_warmup the length of warm-up period
#'
#' @return
#' @export
#'
#' @examples
compile_likefree_model <- function(dat, m_sim, y0, rprior, dprior, times=dat[, "t"],
                                   m_warmup=NA, t_warmup=100) {
  
  sim <- list(data=dat, y0=y0,
              r_prior=rprior, d_prior=dprior,
              m_sim=m_sim, m_warmup=m_warmup)
  
  
  stopifnot(is.finite(dprior(rprior())))
  
  if (is.na(m_warmup)) {
    sim$cm_warmup <- NA
    sim$ts_warmup <- NA
  } else {
    env <- r_prior()
    env$y0 <- y0
    
    sim$cm_warmup <- m_warmup(user=env)
    sim$ts_warmup <- seq(min(times) - t_warmup, min(times), by=1)
  }
  
  env <- r_prior()
  env$y0 <- y0
  
  sim$cm_sim <- m_sim(user=env)
  sim$ts_sim <- sort(times)
  
  class(sim) <- "likefree_model"
  return(sim)
}