#' @title 
#' The Grenander mode estimator
#' 
#' @description 
#' This function computes the Grenander mode estimator.
#' 
#' @details 
#' The Grenander estimate is defined by 
#' \deqn{ \frac{ \sum_{j=1}^{n-k} \frac{(x_{j+k} + x_{j})}{2(x_{j+k} - x_{j})^p} }
#' { \sum_{j=1}^{n-k} \frac{1}{(x_{j+k} - x_{j})^p} } }{ ( sum_{j=1}^{n-k} (x_{j+k} + x_{j})/(2(x_{j+k} - x_{j})^p) ) / ( sum_{j=1}^{n-k} 1/((x_{j+k} - x_{j})^p) ) } 
#' 
#' If \eqn{p}{p} tends to infinity, this estimate tends to the Venter mode estimate; 
#' this justifies to call \code{\link[modeest]{venter}} if \code{p = Inf}. 
#' 
#' The user should either give the bandwidth \code{bw} or the argument \code{k}, 
#' \code{k} being taken equal to \code{ceiling(bw*n) - 1} if missing.
#' 
#' @references 
#' \itemize{
#'   \item Grenander U. (1965). 
#'   Some direct estimates of the mode. 
#'   \emph{Ann. Math. Statist.}, \bold{36}:131-138.
#'   
#'   \item Dalenius T. (1965). 
#'   The Mode - A Negleted Statistical Parameter. 
#'   \emph{J. Royal Statist. Soc. A}, \emph{128}:110-117. 
#'   
#'   \item Adriano K.N., Gentle J.E. and Sposito V.A. (1977). 
#'   On the asymptotic bias of Grenander's mode estimator. 
#'   \emph{Commun. Statist.-Theor. Meth. A}, \bold{6}:773-776. 
#'   
#'   \item Hall P. (1982). 
#'   Asymptotic Theory of Grenander's Mode Estimator. 
#'   \emph{Z. Wahrsch. Verw. Gebiete}, \bold{60}:315-334.
#' }
#' 
#' @note 
#' The user should preferentially call \code{grenander} through 
#' \code{mlv(x, method = "grenander", bw, k, p, ...)}. 
#'
#' @param x 
#' numeric. Vector of observations.
#' 
#' @param bw
#' numeric. The bandwidth to be used. Should belong to (0, 1].
#' 
#' @param k
#' numeric. Paramater 'k' in Grenander's mode estimate, see below. 
#' 
#' @param p
#' numeric. Paramater 'p' in Grenander's mode estimate, see below. 
#' If \code{p = Inf}, the function \code{\link[modeest]{venter}} is used. 
#' 
#' @param ...
#' Additional arguments to be passed to \code{\link[modeest]{venter}}. 
#' 
#' @return 
#' A numeric value is returned, the mode estimate. 
#' If \code{p = Inf}, the \code{\link[modeest]{venter}} mode estimator is returned. 
#' 
#' @author D.R. Bickel for the original code, 
#' P. Poncet for the slight modifications introduced.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation; 
#' \code{\link[modeest]{venter}} for the Venter mode estimate. 
#' 
#' @export
#' @aliases Grenander
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rnorm(1000, mean = 23, sd = 0.5) 
#' 
#' ## True mode
#' normMode(mean = 23, sd = 0.5) # (!)
#' 
#' ## Parameter 'k'
#' k <- 5
#' 
#' ## Many values of parameter 'p'
#' ps <- seq(0.1, 4, 0.01)
#' 
#' ## Estimate of the mode with these parameters
#' M <- sapply(ps, function(p) grenander(x, p = p, k = k))
#' 
#' ## Distribution obtained
#' plot(density(M), xlim = c(22.5, 23.5))
#'
grenander <-
function(x, 
         bw = NULL,
         k,
         p,
         ...)
{
  if (p == Inf) {
    cat("argument 'p' is infinite. Venter's mode estimator is used")
    return(venter(x = x, bw = bw, k = k, ...)) 
  }

  ny <- length(x)    

  if (missing(k) & !is.null(bw)) {
    if (bw <= 0 | bw > 1) stop("argument 'bw' must belong to (0, 1]")
    k <- ceiling(bw*ny) - 1
  } else if (missing(k) & is.null(bw)) {
    k <- ceiling(ny/2) - 1
  }
  
  if (k < 0 | k >= ny) stop("argument 'k' must belong to [0, length('x'))") 

  y <- sort(x)

  inf <- y[1:(ny-k)]
  sup <- y[(k+1):ny]
  diff <- sup - inf
  tot <- inf + sup
  if (any(diff==0)) {
    warning("limiting value of Grenander mode used") #! ??
    M <- mean(ifelse(diff==0, tot, NA), na.rm = TRUE)/2
  } else {
    b <- sum(tot/diff^p)/2
    a <- sum(1/diff^p)
    if (is.finite(b/a)) {
      M <- b/a
    } else {
      stop("function 'grenander' failed. Argument 'p' may be too large")
    }
  }
  M
}
