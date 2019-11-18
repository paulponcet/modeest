#' @title 
#' The Meanshift mode estimator
#' 
#' @description 
#' The Meanshift mode estimator. 
#' 
#' @note 
#' The user should preferentially call \code{meanshift} through 
#' \code{mlv(x, method = "meanshift", ...)}. 
#' 
#' @references 
#' \itemize{ 
#'   \item Fukunaga, K. and Hostetler, L. (1975).  
#'   The estimation of the gradient of a density function, 
#'   with applications in pattern recognition. 
#'   \emph{IEEE Transactions on Information Theory}, \bold{21}(1):32--40. 
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
#' \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, 
#' \code{"epanechnikov"}, \code{"gaussian"}, \code{"optcosine"}, 
#' \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}. 
#' See \code{\link[stats]{density}} for more details on some of these kernels. 
#' 
#' @param par
#' numeric. The initial value used in the meanshift algorithm. 
#' 
#' @param iter
#' numeric. Maximal number of iterations. 
#' 
#' @param tolerance
#' numeric. Stopping criteria.
#' 
#' @return 
#' \code{meanshift} returns a numeric value, the mode estimate, 
#' with an attribute \code{"iterations"}. 
#' The number of iterations can be less than \code{iter} 
#' if the stopping criteria specified by \code{eps} is reached. 
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}}, \code{\link[modeest]{tsybakov}}. 
#' 
#' @importFrom stats bw.SJ bw.nrd0 bw.nrd bw.ucv bw.bcv
#' @importFrom statip kernelfun
#' @export
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rweibull(100, shape = 12, scale = 0.8)
#' 
#' ## True mode
#' weibullMode(shape = 12, scale = 0.8)
#' 
#' ## Estimate of the mode
#' mlv(x, method = "meanshift", par = mean(x))
#' 
meanshift <-
function(x,
         bw = NULL,
         kernel = "gaussian",
         par = shorth(x),
         iter = 1000,
         tolerance = sqrt(.Machine$double.eps))
{
  if (is.null(bw)) bw <- "nrd0"
  if (is.character(bw)) {
    if (length(x) < 2L) 
      stop("need at least 2 points to select a bandwidth automatically", 
           call. = FALSE)
    bw <- switch(tolower(bw), 
                 nrd0 = stats::bw.nrd0(x), 
                 nrd = stats::bw.nrd(x), 
                 ucv = stats::bw.ucv(x), 
                 bcv = stats::bw.bcv(x), 
                 sj = , 
                 `sj-ste` = stats::bw.SJ(x, method = "ste"), 
                 `sj-dpi` = stats::bw.SJ(x, method = "dpi"), 
                 stop("unknown bandwidth rule", call. = FALSE))
  }
  s <- 0
  for (j in seq_len(iter)) {
    z <- (x - par)/bw
    k <- statip::kernelfun(kernel)(z)
    M <- crossprod(x, k)/sum(k)
    if (is.nan(M)) stop("sum(k) is zero in the meanshift function. Change the 
                         bandwidth 'bw' or the initial value 'par'.", 
                        call. = FALSE)
    th <- abs(M/par-1)
    if (th < tolerance) {
      s <- j
      break
    }
    par <- as.vector(M)
  }
  attr(par, "iterations") <- s
  par
}
