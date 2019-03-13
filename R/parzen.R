#' @title 
#' Parzen's Kernel mode estimator
#' 
#' @description 
#' Parzen's kernel mode estimator is the value 
#' maximizing the kernel density estimate. 
#' 
#' @details 
#' If \code{kernel = "uniform"}, the \code{\link[modeest]{naive}} mode estimate is returned.
#' 
#' @note 
#' The user may call \code{parzen} through 
#' \code{mlv(x, method = "kernel", ...)} or \code{mlv(x, method = "parzen", ...)}. 
#' 
#' Presently, \code{parzen} is quite slow.
#' 
#' @references 
#' \itemize{ 
#' \item Parzen E. (1962). 
#' On estimation of a probability density function and mode. 
#' \emph{Ann. Math. Stat.}, \bold{33}(3):1065--1076. 
#' 
#' \item Konakov V.D. (1973). 
#' On the asymptotic normality of the mode of multidimensional distributions. 
#' \emph{Theory Probab. Appl.}, \bold{18}:794-803. 
#' 
#' \item Eddy W.F. (1980). 
#' Optimum kernel estimators of the mode. 
#' \emph{Ann. Statist.}, \bold{8}(4):870-882.
#' 
#' \item Eddy W.F. (1982). 
#' The Asymptotic Distributions of Kernel Estimators of the Mode. 
#' \emph{Z. Wahrsch. Verw. Gebiete}, \bold{59}:279-290. 
#' 
#' \item Romano J.P. (1988). 
#' On weak convergence and optimality of kernel density estimates of the mode. 
#' \emph{Ann. Statist.}, \bold{16}(2):629-647. 
#' 
#' \item Abraham C., Biau G. and Cadre B. (2003). 
#' Simple Estimation of the Mode of a Multivariate Density. 
#' \emph{Canad. J. Statist.}, \bold{31}(1):23-34. 
#' 
#' \item Abraham C., Biau G. and Cadre B. (2004). 
#' On the Asymptotic Properties of a Simple Estimate of the Mode. 
#' \emph{ESAIM Probab. Stat.}, \bold{8}:1-11. 
#' }
#' 
#' @param x 
#' numeric. Vector of observations.
#' 
#' @param bw
#' numeric. The smoothing bandwidth to be used. 
#' 
#' @param kernel
#' character. The kernel to be used. For available kernels see 
#' \code{\link[statip]{densityfun}} in package \pkg{statip}.
#' 
#' @param abc
#' logical. If \code{FALSE} (the default), the kernel density estimate 
#' is maximised using \code{\link[stats]{optim}}. 
#' 
#' @param tolerance
#' numeric. Desired accuracy in the \code{\link[stats]{optimize}} function.
#' 
#' @param ... 
#' If \code{abc = FALSE}, further arguments to be passed to \code{\link[stats]{optim}}. 
#' 
#' @return 
#' \code{parzen} returns a numeric value, the mode estimate. 
#' If \code{abc = TRUE}, the \code{x} value maximizing the density 
#' estimate is returned. Otherwise, the \code{\link[stats]{optim}} 
#' method is used to perform maximization, and the attributes: 
#' 'value', 'counts', 'convergence' and 'message', coming from 
#' the \code{\link[stats]{optim}} method, are added to the result.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}}, \code{\link[modeest]{naive}}
#' 
#' @importFrom stats optimize bw.SJ bw.nrd0 bw.nrd bw.ucv bw.bcv
#' @importFrom statip kernelfun
#' @importFrom statip .kernelsList
#' @export
#' @aliases Parzen
#' 
#' @examples 
#' # Unimodal distribution 
#' x <- rlnorm(10000, meanlog = 3.4, sdlog = 0.2) 
#' 
#' ## True mode 
#' lnormMode(meanlog = 3.4, sdlog = 0.2) 
#' 
#' ## Estimate of the mode 
#' mlv(x, method = "kernel", kernel = "gaussian", bw = 0.3, par = shorth(x)) 
#' 
parzen <-
function(x,
         bw = NULL,
         kernel = "gaussian",
         abc = FALSE,
         tolerance = .Machine$double.eps^0.25,
         ...)
{

  if (pmatch(tolower(kernel), "normal", nomatch = 0)) {
    kernel <- "gaussian"
  }

  if (is.null(bw)) bw <- "nrd0"
  f <- statip::densityfun(x, bw = bw, kernel = kernel, ...)
  
  if (!abc) {
    M <- stats::optimize(f, range(x), maximum = TRUE, tol = tolerance)$maximum
  } else {
    f <- f(x)
    M <- x[f == max(f)]
  }
  M
}
