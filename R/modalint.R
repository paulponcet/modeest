
# Author: P. Poncet

modalint <-
function(x,                          # sample (the data)
         bw = seq(0.01, 0.99, 0.01), # sequence of fractions of the observations defining the modal interval
         plot.it = FALSE,
         tie.action = "mean",
         tie.limit = 0.05)
{
  ny <- length(x)

  if (!is.null(bw)) {
    if (any(bw <= 0) | any(bw > 1)) stop("argument 'bw' must belong to (0, 1]")
    K <- ceiling(bw*ny) - 1
  }

  y <- sort(x)

  modal <- c()
  sh <- c()
  for (k in K) {
    inf <- y[1:(ny-k)]
    sup <- y[(k+1):ny]
    diffs <- sup - inf
    i <- which(diffs==min(diffs))

    ## Ties?
    if (length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit)

    ## Modal interval and its mid-point
    modal <- c(modal, y[i], y[i+k])
    sh <- c(sh, mean(y[i:(i+k)]))
  }
  if (plot.it) {
    plot(x = modal, y = rep(bw, each = 2), type = "l", ylim = c(0,1), ylab = "Bandwidth", main = "Modal plot", xlab = "x")
    lines(x = sh, y = bw, col = 2, lwd = 2)
  }
  return(modal)
}
