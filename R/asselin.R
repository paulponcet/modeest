#' @title 
#' The Asselin de Beauville mode estimator
#' 
#' @description 
#' This mode estimator is based on the algorithm 
#' described in Asselin de Beauville (1978).
#' 
#' @note 
#' The user should preferentially call \code{asselin} through 
#' \code{mlv(x, method = "asselin", ...)}. 
#' 
#' @references 
#' \itemize{
#'   \item Asselin de Beauville J.-P. (1978). 
#'   Estimation non parametrique de la densite et du mode, 
#'   exemple de la distribution Gamma. 
#'   \emph{Revue de Statistique Appliquee}, \bold{26}(3):47-70.
#' }
#' 
#' @param x 
#' numeric. Vector of observations. 
#' 
#' @param bw
#' numeric. A number in \code{(0, 1]}. 
#' If \code{bw = 1}, the selected 'modal chain' may be too long. 
#' 
#' @param ... 
#' further arguments to be passed to the \code{\link[stats]{quantile}} function.
#' 
#' @return 
#' A numeric value is returned, the mode estimate.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation. 
#' 
#' @importFrom stats median quantile
#' @export
#' @aliases Asselin
#' 
#' @examples 
#' x <- rbeta(1000, shape1 = 2, shape2 = 5)
#' 
#' ## True mode:
#' betaMode(shape1 = 2, shape2 = 5)
#' 
#' ## Estimation:
#' asselin(x, bw = 1)
#' asselin(x, bw = 1/2)
#' mlv(x, method = "asselin")
#' 
asselin <-
function(x,
         bw = NULL, # bw = 1 donne une chaine modale longue, bw < 1 est plus severe 
         ...)
{
# TODO: look at 'na.contiguous' 

  if (is.null(bw)) bw <- 1
  
  nx <- length(x)
  kmax <- floor(ifelse(nx < 30, 10, 15)*log(nx))
  
  y <- sort(x)

  ok1 <- FALSE
  
  while (!ok1) {

    ny <- length(y)
    if (ny==1) return(y)

    qy <- stats::quantile(y, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), 
                          names = FALSE, ...)
    delta <- min(qy[5] - qy[4], qy[2] - qy[1])

    a <- qy[1] - 3*delta
    b <- qy[5] + 3*delta
    yab <- y[y>=a & y <= b]

    k <- kmax
    ok2 <- FALSE

    while (!ok2) {

      #hy <- hist(yab, breaks = k, plot = FALSE);b <- hy$breaks;n <- c(hy$counts, 0)
      b <- seq(from = min(yab), to = max(yab), length = k+1)
      n <- c(tabulate(findInterval(yab, b[-(k+1)])), 0)

      N <- sum(n)
      v <- as.numeric(n >= N/k)

      ## Beginning of the first chain
      w <- which.max(v)
      v2 <- v[w:(k+1)]

      ## End of the first chain
      w2 <- which.min(v2) + w - 1
      v3 <- v[w2:(k+1)]

      ## Length of the first chain
      nc <- sum(n[w:(w2-1)])
      
      ## There exists another chain, and the first chain has only one element
      if (any(v3==1) & nc==1) {
        if (k > 3) {
          k <- k-1
        } else if (k==3) {
          if (n[3] > 1) {
            w <- 3
            w2 <- 4
          } 
          ok2 <- TRUE
        } else {
          stop("k < 3")
        }
      
      ## There exists another chain, and the first chain has more than one element
      } else if (any(v3==1) & nc > 1) {
        if (k > 3) {
          k <- k-1
        
        ### In this case, w = 1 necessarily
        } else if (k==3) {
          if (n[3] > 1) {           
            p1 <- (1/n[1])*prod(diff(yab[yab >= b[w] & yab <= b[w2]])) # here, n[1] = length(first chain)
            p2 <- (1/n[3])*prod(diff(yab[yab >= b[3] & yab <= b[4]])) # and n[3] = length(second chain)
            if (p1 > p2) {
              w <- 3
              w2 <- 4
            }
          }
          ok2 <- TRUE
        } else {
          stop("k < 3")
        }
      
      ## There is no other chain: the modal chain is found!
      } else if (!any(v3==1)) {
        ok2 <- TRUE
      }

    }

    ## Update 'nc'
    nc <- sum(n[w:(w2-1)])
    #cat("Modal chain length = ", nc, "\n")
    
    d <- abs((qy[4] + qy[2] - 2*qy[3])/(qy[4] - qy[2]))
    nc2 <-  ny*(1-d)
    #cat("d = ", d, "\n")

    y <- yab[yab >= b[w] & yab <= b[w2]]
    if (nc == ny) {
      ok1 <- TRUE
    } else if (nc <= ifelse(nx < 30, nx/3, bw*nc2)) {
      ok1 <- TRUE
    } else {
      ok1 <- FALSE 
    }
  
  }
  stats::median(y)
}
