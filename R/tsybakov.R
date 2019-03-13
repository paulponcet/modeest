#' @title 
#' The Tsybakov mode estimator
#' 
#' @description 
#' This mode estimator is based on a gradient-like recursive algorithm, 
#' more adapted for online estimation. 
#' It includes the Mizoguchi-Shimura (1976) mode estimator, 
#' based on the window training procedure. 
#' 
#' @details 
#' If \code{bw} or \code{a} is missing, a default 
#' value advised by Djeddour et al (2003) is used: 
#' \code{bw = (1:length(x))^(-1/7)} and \code{a = (1:length(x))^(-alpha)}. 
#' (with \code{alpha = 0.9} if \code{alpha} is missing).
#' 
#' @note 
#' The user may call \code{tsybakov} through 
#' \code{mlv(x, method = "tsybakov", ...)}. 
#' 
#' @section
#' Warning: The Tsybakov mode estimate as it is presently 
#' computed does not work very well. 
#' The reasons of this inefficiency should be further investigated. 
#' 
#' @references 
#' \itemize{ 
#'   \item Mizoguchi R. and Shimura M. (1976).
#'   Nonparametric Learning Without a Teacher Based on Mode Estimation.
#'   \emph{IEEE Transactions on Computers}, \bold{C25}(11):1109-1117.
#'   
#'   \item Tsybakov A. (1990).
#'   Recursive estimation of the mode of a multivariate distribution.
#'   \emph{Probl. Inf. Transm.}, \bold{26}:31-37.
#'   
#'   \item Djeddour K., Mokkadem A. et Pelletier M. (2003).
#'   Sur l'estimation recursive du mode et de la valeur modale d'une densite de 
#'   probabilite.
#'   \emph{Technical report 105}.
#'   
#'   \item Djeddour K., Mokkadem A. et Pelletier M. (2003).
#'   Application du principe de moyennisation a l'estimation recursive du mode 
#'   et de la valeur modale d'une densite de probabilite.
#'   \emph{Technical report 106}.
#' }
#' 
#' @param x
#' numeric. Vector of observations. 
#' 
#' @param bw
#' numeric. Vector of length \code{length(x)} 
#' giving the sequence of smoothing bandwidths to be used.
#' 
#' @param a
#' numeric. Vector of length \code{length(x)} used in the 
#' gradient algorithm
#' 
#' @param alpha
#' numeric. An alternative way of specifying \code{a}. See 'Details'. 
#' 
#' @param kernel
#' character. The kernel to be used. Available kernels are 
#' \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, 
#' \code{"epanechnikov"}, \code{"gaussian"}, \code{"optcosine"}, 
#' \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}. 
#' See \code{\link[stats]{density}} for more details on some 
#' of these kernels. 
#' 
#' @param dmp
#' logical. If \code{TRUE}, Djeddour et al. 
#' version of the estimate is used. 
#' 
#' @param par
#' numeric. Initial value in the gradient algorithm. 
#' Default value is \code{\link[modeest]{shorth}(x)}. 
#' 
#' @return 
#' A numeric value is returned, the mode estimate.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation. 
#' 
#' @importFrom statip .kernelsList kernelfun
#' @export
#' @aliases Tsybakov
#' 
#' @examples 
#' x <- rbeta(1000, shape1 = 2, shape2 = 5)
#' 
#' ## True mode:
#' betaMode(shape1 = 2, shape2 = 5)
#' 
#' ## Estimation:
#' tsybakov(x, kernel = "triangular")
#' tsybakov(x, kernel = "gaussian", alpha = 0.99)
#' mlv(x, method = "tsybakov", kernel = "gaussian", alpha = 0.99)
#' 
tsybakov <-
function(x,
         bw = NULL,
         a,
         alpha = 0.9,
         kernel = "triangular", 
         dmp = TRUE, 
         par = shorth(x))
{
  if (pmatch(tolower(kernel), "normal", nomatch = 0)) {
    kernel <- "gaussian"
  } else {
    kernel <- match.arg(tolower(kernel), statip::.kernelsList())
  }
  
  K <- statip::kernelfun(kernel, derivative = TRUE)
  nx <- length(x)
    
  if (missing(a)) {
    a <- (1:nx)^(-alpha)
  }
  
  if (is.null(bw)) bw <- (1:nx)^(-1/7)
    
  ## Initialization
  M.dmp <- M <- par
  b <- a/(bw^2)
  p <- bw^3
  p <- p/cumsum(p)
  
  for (n in 1:nx) {
    M <- M + b[n]*K((M-x[n])/bw[n])
    M.dmp <- M.dmp + p[n]*(M - M.dmp)
  }
  
  ifelse(dmp, M.dmp, M)
}
