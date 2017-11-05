#' @title 
#' The Chernoff or 'naive' mode estimator
#' 
#' @description 
#' This estimator, also called the *naive* mode estimator, is defined as the 
#' center of the interval of given length containing the most observations. 
#' It is identical to Parzen's kernel mode estimator, when the kernel is chosen 
#' to be the uniform kernel.
#' 
#' @note 
#' The user should preferentially call \code{naive} through 
#' \code{mlv(x, method = "naive", bw)}. 
#' 
#' @references 
#' \itemize{ 
#'   \item Chernoff H. (1964). 
#'   Estimation of the mode. 
#'   \emph{Ann. Inst. Statist. Math.}, \bold{16}:31-41.
#'   
#'   \item Leclerc J. (1997). 
#'   Comportement limite fort de deux estimateurs du mode : 
#'   le shorth et l'estimateur naif. 
#'   \emph{C. R. Acad. Sci. Paris, Serie I}, \bold{325}(11):1207-1210.
#' }
#' 
#' @param x
#' numeric. Vector of observations.
#' 
#' @param bw 
#' numeric. The smoothing bandwidth to be used. Should belong to (0, 1). See below.
#' 
#' @return 
#' A numeric vector is returned, the mode estimate, 
#' which is the center of the interval of length \code{2*bw} 
#' containing the most observations.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation; 
#' \code{\link[modeest]{parzen}} for Parzen's kernel mode estimation. 
#' 
#' @export
#' @aliases Chernoff chernoff
# #' @rdname parzen 
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rf(10000, df1 = 40, df2 = 30)
#' 
#' ## True mode
#' fMode(df1 = 40, df2 = 30)
#' 
#' ## Estimate of the mode
#' mean(naive(x, bw = 1/4))
#' mlv(x, method = "naive", bw = 1/4)
#' 
naive <-
function(x,
         bw = 1/2)
{
  parzen(x = x, bw = bw, kernel = "uniform", abc = TRUE)
}
