#' @rdname fit.likefree
#' @export
summary.posterior_likefree <- function(fitted) {
  ans <- fitted$Meta

  su <- apply(fitted$Parameters, 2, function(ps) {
    c(
      mean(ps), sd(ps), min(ps),
      quantile(ps, c(0.025, 0.25, 0.5, 0.75, 0.975)),
      max(ps)
    )
  })

  su <- t(su)

  colnames(su)[c(1:3, 9)] <- c("Mean", "Std", "Min", "Max")

  ans$Post <- su
  class(ans) <- "summary_posterior_likefree"
  return(ans)
}


#' @rdname fit.likefree
#' @export
print.summary_posterior_likefree <- function(su) {
  print(su$Post)
}

