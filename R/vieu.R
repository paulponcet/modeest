#' @title 
#' Vieu's mode estimator
#' 
#' @description 
#' Vieu's mode estimator is the value at which the kernel density derivative 
#' estimate is null. 
#' 
#' @note 
#' The user should preferentially call \code{vieu} through 
#' \code{mlv(x, method = "vieu", ...)}. 
#' 
#' Presently, \code{vieu} is quite slow.
#' 
#' @references 
#' \itemize{ 
#'   \item Vieu P. (1996). A note on density mode estimation. 
#'   \emph{Statistics \& Probability Letters}, \bold{26}:297--307.
#' }
#' 
#' @param x 
#' numeric. Vector of observations. 
#' 
#' @param bw
#' numeric. The smoothing bandwidth to be used. 
#' 
#' @param kernel
#' character. The kernel to be used. Available kernels are \code{"biweight"}, 
#' \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"}, \code{"gaussian"}, 
#' \code{"optcosine"}, \code{"rectangular"}, \code{"triangular"}, 
#' \code{"uniform"}. See \code{\link[stats]{density}} for more details on some 
#' of these kernels. 
#' 
#' @param abc
#' logical. If \code{FALSE} (the default), the root of the density derivate 
#' estimate is searched with \code{\link[stats]{uniroot}}. 
#' 
#' @param ... 
#' If \code{abc = FALSE}, further arguments to be passed to 
#' \code{\link[stats]{uniroot}}. 
#' 
#' @return 
#' \code{vieu} returns a numeric value, the mode estimate. If \code{abc = TRUE}, 
#' the \code{x} value at which the density derivative estimate is null is 
#' returned. Otherwise, the \code{\link[stats]{uniroot}} method is used.
#' 
#' @importFrom stats uniroot
#' @importFrom statip .kernelsList kernelfun
#' @export
#' @aliases Vieu
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}}, \code{\link[modeest]{parzen}}. 
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rlnorm(10000, meanlog = 3.4, sdlog = 0.2)
#' 
#' ## True mode
#' lnormMode(meanlog = 3.4, sdlog = 0.2)
#' 
#' ## Estimate of the mode
#' mlv(x, method = "vieu", kernel = "gaussian")
#' 
vieu <-
function(x,
         bw = NULL,
         kernel = "gaussian",
         abc = FALSE,
         ...)
{
  if (pmatch(tolower(kernel), "normal", nomatch = 0)) {
    kernel <- "gaussian"
  } else {
    kernel <- match.arg(tolower(kernel), c(statip::.kernelsList(), "uniform"))
  }
    
  #nx <- length(x)
  if (is.null(bw)) bw <- bw.SJ(x)
  
  fn <- 
  function(z)
  {
    k <- statip::kernelfun(kernel, derivative = TRUE)((z-x)/bw)
    sum(k)
  }
  
  #FN <- 
  #function(z)
  #{
  #  mat <- kronecker(z/bw, t(-x/bw), FUN = "+")
  #  k <- do.call(paste(".kernel.d", kernel, sep = ""), list(mat))$k
  #  return(rowSums(k))
  #}
  
  if (!abc) {
    r <- stats::uniroot(f = fn, interval = c(min(x), max(x)), ...)
    M <- r$root
  } else {
    FN <- Vectorize(fn)
    f <- abs(FN(x)) 
    M <- x[f == min(f)]
  }
  M
}
