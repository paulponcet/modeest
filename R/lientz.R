#' @title 
#' The empirical Lientz function and the Lientz mode estimator
#' 
#' @description 
#' The Lientz mode estimator is nothing but the value minimizing the empirical 
#' Lientz function. A 'plot' and a 'print' methods are provided. 
#' 
#' @references 
#' \itemize{
#'   \item Lientz B.P. (1969).
#'   On estimating points of local maxima and minima of density functions.
#'   \emph{Nonparametric Techniques in Statistical Inference (ed. M.L. Puri, Cambridge University Press}, p.275-282.
#'   
#'   \item Lientz B.P. (1970).
#'   Results on nonparametric modal intervals.
#'   \emph{SIAM J. Appl. Math.}, \bold{19}:356-366.
#'   
#'   \item Lientz B.P. (1972).
#'   Properties of modal intervals.
#'   \emph{SIAM J. Appl. Math.}, \bold{23}:1-5.
#' }
#' 
#' @details 
#' The Lientz function is the smallest non-negative quantity \eqn{S(x,\beta)}{S(x,b)}, 
#' where \eqn{\beta}{b} = \code{bw}, such that 
#' \deqn{F(x+S(x,\beta)) - F(x-S(x,\beta)) \geq \beta.}{F(x+S(x,b)) - F(x-S(x,b)) >= b.} 
#' Lientz (1970) provided a way to estimate \eqn{S(x,\beta)}{S(x,b)}; this estimate 
#' is what we call the empirical Lientz function. 
#' 
#' @note 
#' The user may call \code{mlv.lientz} through 
#' \code{mlv(x, method = "lientz", ...)}. 
#' 
#' @param x 
#' numeric (vector of observations) or an object of class \code{"lientz"}.
#' 
#' @param bw
#' numeric. The smoothing bandwidth to be used. 
#' Should belong to (0, 1). Parameter 'beta' in Lientz (1970) function.
#' 
#' @param abc
#' logical. If \code{FALSE} (the default), the Lientz empirical function 
#' is minimised using \code{\link[stats]{optim}}.
#' 
#' @param par
#' numeric. The initial value used in \code{\link[stats]{optim}}.
#' 
#' @param optim.method
#' character. If \code{abc = FALSE}, the method used in 
#' \code{\link[stats]{optim}}.
#' 
#' @param zoom
#' logical. If \code{TRUE}, one can zoom on the graph created.
#' 
#' @param digits
#' numeric. Number of digits to be printed.
#' 
#' @param ...
#' if \code{abc = FALSE}, further arguments to be passed to 
#' \code{\link[stats]{optim}}, or further arguments to be passed to 
#' \code{\link[graphics]{plot}}.
#' 
#' @return 
#' \code{lientz} returns an object of class \code{c("lientz", "function")}; 
#' this is a function with additional attributes:
#' \itemize{
#'   \item{x}{ the \code{x} argument}
#'   \item{bw}{ the \code{bw} argument }
#'   \item{call}{ the call which produced the result }
#' }
#' 
#' \code{mlv.lientz} returns a numeric value, the mode estimate. 
#' If \code{abc = TRUE}, the \code{x} value minimizing the Lientz empirical 
#' function is returned. Otherwise, the \code{\link[stats]{optim}} method is 
#' used to perform minimization, and the attributes: 'value', 'counts', 
#' 'convergence' and 'message', coming from the \code{\link[stats]{optim}} 
#' method, are added to the result.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation; 
#' \code{\link[modeest]{shorth}} for the shorth estimate of the mode
#' 
#' @export
#' @aliases Lientz
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rbeta(1000,23,4)
#' 
#' ## True mode
#' betaMode(23, 4)
#' 
#' ## Lientz object
#' f <- lientz(x, 0.2)
#' print(f)
#' plot(f)
#' 
#' ## Estimate of the mode
#' mlv(f)              # optim(shorth(x), fn = f)
#' mlv(f, abc = TRUE)  # x[which.min(f(x))]
#' mlv(x, method = "lientz", bw = 0.2)
#' 
#' # Bimodal distribution
#' x <- c(rnorm(1000,5,1), rnorm(1500, 22, 3))
#' f <- lientz(x, 0.1)
#' plot(f)
#' 
lientz <-
function(x,
         bw = NULL)
{
  if (bw <= 0 || bw >= 1) {
    stop("argument 'bw' must belong to (0, 1)", call. = FALSE)
  }
  
  y <- sort(x)
  ny <- length(y)
  k <- ceiling(bw*ny) - 1
  
  if (k==0) {
    f <-
    function(z)
    {
      yy <- (y + c(y[-1],Inf))/2
      i <- sapply(z, FUN = function(zz) min(which(zz <= yy)))      
      return(abs(y[i] - z))    
    }
  
  } else if (k>0) {
    f <-
    function(z)
    {
      yy1 <- c(y[(k+1):ny],rep(Inf,k))   # k+1 lags
      yy2 <- c(y[(k+2):ny],rep(Inf,k+1))   # k+2 lags
      yy <- sort(c((y+yy1)/2,(y+yy2)/2))
      
      i <- sapply(z, FUN = function(zz) min(which(zz <= yy)))
      j <- i%%2
      yy <- y[j*k + (i+j)/2]
      return(ifelse(j==1, yy-z, z-yy))
    }
  }
  
  class(f) <- c("lientz", class(f))
  attr(f, "call") <- sys.call()
  attr(f, "x") <- x
  attr(f, "bw") <- bw
  attr(f, "source") <- NULL
  f
}


#' @importFrom graphics plot points legend locator
#' @export
#' @rdname lientz
#' 
plot.lientz <-
function(x,            # an object of class 'lientz'
         zoom = FALSE, # if TRUE, one can zoom on the graph created
         ...)
{
  if (!inherits(x, "lientz")) {
    stop("argument 'x' must inherit from class 'lientz'")
  }

  arg <- list(...)
  ylim <- arg$ylim
  main <- arg$main
  xlab <- arg$xlab
  ylab <- arg$ylab

  xx <- attr(x, "x")
  #bw <- attr(x, "bw") 
  
  inf <- min(xx)
  sup <- max(xx)
  z <- seq(inf, sup, (sup - inf)/1024)
  lz <- x(z)

  if (is.null(ylim)) ylim <- range(lz)    
  if (is.null(main)) main <- "Empirical Lientz's function"
  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "Sn(x)"
    
  graphics::plot(z, lz, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
  graphics::points(xx, rep(ylim[1],length(xx)), pch = "'", col = 4)
  graphics::legend("topleft",legend = c("Regular grid", "x"), col = c(1,4), pch = 19, bg = "white")
  
  if (zoom) {
    cat("you can zoom on the graph (press 'Esc' to escape)\n")  
    lc <- graphics::locator(2)
    while (!is.null(lc)) {
      xlim <- sort(c(lc$x[1], lc$x[2]))
      ylim <- sort(c(lc$y[1], lc$y[2]))
      plot.lientz(x, zoom = FALSE, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
      lc <- graphics::locator(2)
    }  
  }
  invisible(NULL)
}


#' @export
#' @rdname lientz
#' @method print lientz
#' 
print.lientz <-
function(x,             # an object of class 'lientz'
         digits = NULL,
         ...)
{
  if (!inherits(x, "lientz")) {
    stop("argument 'x' must inherit from class 'lientz'", 
         call. = FALSE)
  }
  bw <- attr(x, "bw")
  call <- attr(x, "call")
  #xx <- attr(x, "x")
  cat("Empirical Lientz function\n")
  cat("Call:",deparse(call),"\n")
  cat("bw =", format(bw, digits = digits), "\n")
}


#' @importFrom stats optim
#' @export
#' @rdname lientz
#' 
mlv.lientz <-
function(x,                       # sample (the data) or object of class 'lientz'
         bw = NULL,               # bandwidth
         abc = FALSE,            # if FALSE, 'optim' is used
         par = shorth(x),         # initial value used in 'optim'
         optim.method = "BFGS",   # method used in 'optim'
         ...)
{
  ## Initialization
  if (!inherits(x, "lientz")) {
    Sn <- lientz(x, bw)
  } else {
    Sn <- x
    x <- attr(Sn, "x")
  }
    
  if (!abc) {
    mini <- stats::optim(par, fn = Sn, method = optim.method, control=list(fnscale=1),...)
    M <- mini$par
    attr(M, "value") <- mini$value
    attr(M, "counts") <- mini$counts
    attr(M, "convergence") <- mini$convergence
    attr(M, "message") <- mini$message
  } else {
    Sn <- Sn(x)
    M <- mean(x[Sn == min(Sn)])
  }
  M
}
