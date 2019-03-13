#' @title 
#' The Venter / Dalenius / LMS mode estimator
#' 
#' @description
#' This function computes the Venter mode estimator, also called the Dalenius, 
#' or LMS (Least Median Square) mode estimator. 
#' 
#' @details 
#' The modal interval, i.e. the shortest interval among intervals containing 
#' \code{k+1} observations, is first computed. (In dimension > 1, this question 
#' is known as a 'k-enclosing problem'.)
#' The user should either give the bandwidth \code{bw} or the argument \code{k}, 
#' \code{k} being taken equal to \code{ceiling(bw*n) - 1} if missing, so 
#' \code{bw} can be seen as the fraction of the observations to be considered 
#' for the shortest interval. 
#' 
#' If \code{type = 1}, the midpoint of the modal interval is returned.
#' If \code{type = 2}, the \code{floor((k+1)/2)}th element of the modal 
#' interval is returned.
#' If \code{type = 3} or \code{type = "dalenius"}, the median of the modal 
#' interval is returned.
#' If \code{type = 4} or \code{type = "shorth"}, the mean of the modal interval 
#' is returned.
#' If \code{type = 5} or \code{type = "ekblom"}, Ekblom's 
#' \eqn{L_{-\infty}}{L_{-infinity}} estimate is returned, see Ekblom (1972). 
#' If \code{type = 6} or \code{type = "hsm"}, the half sample mode (hsm) is 
#' computed, see \code{\link{hsm}}.
#' 
#' @note 
#' The user may call \code{venter} through 
#' \code{mlv(x, method = "venter", ...)}. 
#' 
#' @references 
#' \itemize{
#'   \item Dalenius T. (1965). 
#'   The Mode - A Negleted Statistical Parameter. 
#'   \emph{J. Royal Statist. Soc. A}, \emph{128}:110-117.
#'   
#'   \item Venter J.H. (1967). 
#'   On estimation of the mode. 
#'   \emph{Ann. Math. Statist.}, \bold{38}(5):1446-1455. 
#'   
#'   \item Ekblom H. (1972). 
#'   A Monte Carlo investigation of mode estimators in small samples. 
#'   \emph{Applied Statistics}, \bold{21}:177-184.
#
#\item Rousseeuw and Leroy, 1987  #(ou bien Andrews ?)
#' 
#'   \item Leclerc J. (1997). 
#'   Comportement limite fort de deux estimateurs du mode : le shorth et l'estimateur naif. 
#'   \emph{C. R. Acad. Sci. Paris, Serie I}, \bold{325}(11):1207-1210.
#' }
#' 
#' @param x 
#' numeric. Vector of observations. 
#' 
#' @param bw
#' numeric. The bandwidth to be used. Should belong to (0, 1]. See 'Details'.
#' 
#' @param k
#' numeric. See 'Details'.
#' 
#' @param iter
#' numeric. Number of iterations.
#' 
#' @param type 
#' numeric or character. The type of Venter estimate to be computed. See 'Details'.
#' 
#' @param tie.action
#' character. The action to take if a tie is encountered.
#' 
#' @param tie.limit
#' numeric. A limit deciding whether or not a warning is given when a tie is 
#' encountered.
#' 
#' @param warn
#' logical. If \code{TRUE}, a warning is thrown when a tie is encountered. 
#' 
#' @param ...
#' Further arguments.
#' 
#' @return 
#' A numeric value is returned, the mode estimate.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation, 
#' \code{\link[modeest]{hsm}} for the half sample mode. 
#' 
#' @importFrom stats median
#' @export 
#' @aliases Venter
#' 
#' @examples 
#' library(evd)
#' 
#' # Unimodal distribution
#' x <- rgev(1000, loc = 23, scale = 1.5, shape = 0)
#' 
#' ## True mode
#' gevMode(loc = 23, scale = 1.5, shape = 0)
#' 
#' ## Estimate of the mode
#' venter(x, bw = 1/3)
#' mlv(x, method = "venter", bw = 1/3)
#' 
venter <-
function(x,
         bw = NULL, 
         k,
         iter = 1,
         type = 1,
         tie.action = "mean",
         tie.limit = 0.05, 
         warn = FALSE)
{

  ny <- length(x)
    
  ## Initialization
  type <- match.arg(tolower(as.character(type)), 
                    c("-inf", "1", "2", "3", "dalenius", "4", "shorth", "5", "ekblom", "6", "hsm"))
  if (type == "3") type <- "dalenius"
  if (type == "4") type <- "shorth"
  if (type == "-Inf" || type == "5") type <- "ekblom"
  if (type == "6") type <- "hsm"
  
  if (type == "hsm") return(hsm(x = x, bw = bw, k = k, tie.action = tie.action, 
                                tie.limit = tie.limit))
  
  if (missing(k) && !is.null(bw)) {
    if (bw <= 0 || bw > 1) stop("argument 'bw' must belong to (0, 1]")
    k <- ceiling(bw*ny) - 1
  } else if (missing(k) & is.null(bw)) {
    if (type == "ekblom") {
      k <- 1
    } else {
      k <- ceiling(ny/2) - 1
    }
  }
    
  if (k < 0 || k >= ny) stop("argument 'k' must belong to [0, length('x'))") 
  
  y <- sort(x)

  inf <- y[1:(ny-k)]
  sup <- y[(k+1):ny]
  diffs <- sup - inf
  i <- which(diffs==min(diffs))
  
  ## Ties?
  if (length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit, warn = warn)
  
  ## Output
  M <- switch(type,
              "1" = (y[i] + y[i+k])/2,
              "2" = y[i+floor((k+1)/2)],
              "dalenius" = stats::median(y[i:(i+k)]),
              "shorth" = mean(y[i:(i+k)]),
              "ekblom" = ifelse(y[i+2]-y[i+1]>y[i]-y[i-1], y[i], y[i+1]))
  
  if (iter > 1) {
    M <- Recall(x=y[i:(i+k)], bw=(k+1)/ny, iter = iter-1, type=type, 
                tie.action=tie.action, tie.limit=tie.limit)
  }
    
  #attr(M, "inf") <- y[i]
  #attr(M, "sup") <- y[i+k]
  M
}


#' @export 
#' @rdname venter
#' 
shorth <- function(x, ...) venter(x, type = "shorth", ...)
