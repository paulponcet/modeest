#' @title 
#' Half sample mode estimator
#' 
#' @description 
#' This function computes the Robertson-Cryer mode estimator 
#' described in Robertson and Cryer (1974), 
#' also called half sample mode (if \code{bw = 1/2}) 
#' or fraction sample mode (for some other \code{bw}) by Bickel (2006). 
#' 
#' @details 
#' The modal interval, i.e. the shortest interval among 
#' intervals containing \code{k+1} observations, is computed 
#' iteratively, until only one value is found, the mode estimate. 
#' At each step \eqn{i}{i}, one takes \code{k = ceiling(bw*n) - 1}, 
#' where \code{n} is the length of the modal interval computed 
#' at step \eqn{i-}{i-}\code{1}. 
#' If \code{bw} is of class \code{"function"}, 
#' then \code{k = ceiling(bw(n)) - 1} instead.
#' 
#' @note 
#' The user may call \code{hsm} through 
#' \code{mlv(x, method = "hsm", ...)}. 
#' 
#' @references 
#' \itemize{ 
#'   \item Robertson T. and Cryer J.D. (1974). 
#'   An iterative procedure for estimating the mode. 
#'   \emph{J. Amer. Statist. Assoc.}, \bold{69}(348):1012-1016.
#'   
#'   \item Bickel D.R. and Fruehwirth R. (2006). 
#'   On a Fast, Robust Estimator of the Mode: Comparisons to 
#'   Other Robust Estimators with Applications. 
#'   \emph{Computational Statistics and Data Analysis}, \bold{50}(12):3500-3530.
#' }
#' 
#' @param x 
#' numeric. Vector of observations. 
#' 
#' @param bw
#' numeric or function. 
#' The bandwidth to be used. Should belong to (0, 1].
#' 
#' @param k
#' numeric. See 'Details'.
#' 
#' @param tie.action
#' character. The action to take if a tie is encountered. 
#' 
#' @param tie.limit
#' numeric. A limit deciding whether or not a warning 
#' is given when a tie is encountered. 
#' 
#' @param ...
#' Additional arguments.
#' 
#' @return 
#' A numeric value is returned, the mode estimate.
#' 
#' @author 
#' D.R. Bickel for the original code, 
#' P. Poncet for the slight modifications introduced.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation; 
#' \code{\link[modeest]{venter}} for the Venter mode estimate. 
#' 
#' @export
#' @aliases HSM
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rweibull(10000, shape = 3, scale = 0.9)
#' 
#' ## True mode
#' weibullMode(shape = 3, scale = 0.9)
#' 
#' ## Estimate of the mode
#' bandwidth <- function(n, alpha) {1/n^alpha}
#' hsm(x, bw = bandwidth, alpha = 2)
#' mlv(x, method = "hsm", bw = bandwidth, alpha = 2)
#' 
hsm <-
function(x,
         bw = NULL,
         k,
         tie.action = "mean",
         tie.limit = 0.05,
         ...)
{
  if (!missing(k) && is.null(bw)) {
    bw <- (k+1)/length(x)
  } else if (missing(k) && is.null(bw)) {
    bw <- 1/2
  }

  if (is.numeric(bw)) {
    if (bw <= 0 || bw > 1) stop("argument 'bw' must belong to (0, 1]")
  }
  
  y <- sort(x)
   
  while (length(y) >= 4) {
    ny <- length(y)
    if (is.function(bw)) {
      k <- ceiling(bw(ny, ...)*ny) - 1
    } else {
      k <- ceiling(bw*ny) - 1
    }

    inf <- y[1:(ny-k)]
    sup <- y[(k+1):ny]
    diffs <- sup - inf
    i <- which(diffs==min(diffs))
    
    ## Ties?
    if(length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit) 
   
    if (diffs[i]==0) {
      y <- y[i]
    } else {
      y <- y[i:(i+k)]
    }
    #y <- ifelse(diffs[i]==0, y[i], y[i:(i+k)])
  }
  if (length(y) == 3) {
    z <- 2*y[2] - y[1] - y[3]
    M <- switch(as.character(sign(z)), 
                "-1" =  mean(y[1:2]), 
                "1" = mean(y[2:3]), 
                "0" = y[2])
  } else {
    M <- mean(y)
  }
  M
}


#! Recursive estimator
#mlv.hsm.rec <-
#function(x,  # sample
#         bw, # fraction of the observations to consider
#         tie.action = "mean",
#         tie.limit = 0.05)
#{
#  if (missing(bw)) bw <- 1/2
#  if (bw <= 0 || bw > 1) stop("Argument 'bw' must belong to (0, 1].")

#  aux <-
#  function(y, ny, k)
#  {
#    if (ny == 3) {
#      z <- 2*y[2] - y[1] - y[3]
#      return(switch(as.character(sign(z)), "-1" =  mean(y[1:2]), "1" = mean(y[2:3]), "0" = y[2]))
#    } else {
#      if (ny < 3) {
#        return(mean(y))
#      } else {
#        diffs <- y[(k+1):ny] - y[1:(ny-k)]
#        i <- which(diffs==min(diffs))
#        if(length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit)        
#        y <- ifelse(diffs[i]==0, y[i], y[i:(i+k)])
#        ny <- length(y)
#        Recall(y, ny, ceiling(ny*bw)-1)
      
#      }
#    }
#  }
  
  ## Output
#  nx <- length(x)
#  return(aux(sort(x), nx, ceiling(bw*nx)-1))
   
#}
