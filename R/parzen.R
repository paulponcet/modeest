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
#' The user should preferentially call \code{parzen} through 
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
#' character. The kernel to be used. Available kernels are 
#' \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"}, 
#' \code{"gaussian"}, \code{"optcosine"}, \code{"rectangular"}, 
#' \code{"triangular"}, \code{"uniform"}. 
#' See \code{\link[stats]{density}} for more details on some of these kernels.
#' 
#' @param abc
#' logical. If \code{FALSE} (the default), the kernel density estimate 
#' is maximised using \code{\link[stats]{optim}}. 
#' 
#' @param par
#' numeric. The initial value used in \code{\link[stats]{optim}}.
#' 
#' @param optim.method
#' character. If \code{abc = FALSE}, the method used in \code{\link[stats]{optim}}.
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
#' @importFrom stats ecdf optim bw.SJ bw.nrd0 bw.nrd bw.ucv bw.bcv
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
function(x,                       # sample (the data)
         bw = NULL,               # bandwidth
         kernel = "gaussian",     # kernel used
         abc = FALSE,             # if FALSE, 'optim' is used
         par = shorth(x),         # initial value used in 'optim'
         optim.method = "BFGS",   # method used in 'optim'
         ...)
{

  if (pmatch(tolower(kernel), "normal", nomatch = 0)) {
    kernel <- "gaussian"
  } else {
    kernel <- match.arg(tolower(kernel), c(statip::.kernelsList(), "uniform"))
  }

  if (kernel == "uniform") {
    if (is.null(bw)) bw <- 1/2
    if (bw <= 0 | bw >= 1) stop("argument 'bw' must belong to (0, 1)")    
    Fn <- stats::ecdf(x)
    fn <- Fn(x+bw) - Fn(x-bw) # the estimate of the density is (Fn(x+bw) - Fn(x-bw))/(2*bw)
    ## Remark : optimization fails since fn is not regular enough
    ## (any point is viewed as a local maxima). The following is the only solution:
    M <- x[fn == max(fn)]
    #! problem: possibly length(M) > 1 !!
  
  } else {  
    
    ## Initialization
    nx <- length(x)
    if (is.null(bw)) bw <- "nrd0"
    if (is.character(bw)) {
      if (nx < 2)
        stop("need at least 2 points to select a bandwidth automatically")
      bw <- switch(tolower(bw), 
                   nrd0 = stats::bw.nrd0(x), nrd = bw.nrd(x), 
                   ucv = bw.ucv(x), 
                   bcv = stats::bw.bcv(x), 
                   sj = , 
                   `sj-ste` = stats::bw.SJ(x, method = "ste"), 
                   `sj-dpi` = stats::bw.SJ(x, method = "dpi"), 
                   stop("unknown bandwidth rule"))
    }
   
    fn <- 
    function(z)
    {
      mat <- z/bw - x/bw
      k <- statip::kernelfun(kernel)(mat)
      return(sum(k))
    }
    
    #FN <- 
    #function(z)
    #{
    #  mat <- kronecker(z/bw, t(-x/bw), FUN = "+")
    #  k <- do.call(paste(".kernel.", kernel, sep = ""), list(mat))$k
    #  return(rowSums(k))
    #}
  
    if (!abc) {
      maxi <- stats::optim(par, fn, method = optim.method, control=list(fnscale=-1), ...)
      M <- maxi$par
      attr(M, "value") <- maxi$value
      attr(M, "counts") <- maxi$counts
      attr(M, "convergence") <- maxi$convergence
      attr(M, "message") <- maxi$message
    } else {
      f <- Vectorize(fn)
      f <- f(x)
      M <- x[f == max(f)]
    }
  }
  M
}
