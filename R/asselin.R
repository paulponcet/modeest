#' @title
#' The Asselin de Beauville mode estimator
#'
#' @description
#' This mode estimator is based on the algorithm
#' described in Asselin de Beauville (1978).
#'
#' @note
#' The user may call \code{asselin} through
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
asselin <- function(x,
                    bw = NULL, # bw = 1 donne une chaine modale longue, bw < 1 est plus severe
                    ...) {

  # TODO: look at 'na.contiguous'

  if (is.null(bw)) bw <- 1.

  nx <- length(x)
  kmax <- floor(ifelse(nx < 30L, 10L, 15L) * log(nx))

  y <- sort(x)

  ok1 <- FALSE

  while (!ok1) {

    ny <- length(y)
    if (ny == 1L) return(y)

    qy <- stats::quantile(y, probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
                          names = FALSE, ...)
    delta <- min(qy[5L] - qy[4L], qy[2L] - qy[1L])

    a <- qy[1L] - 3. * delta
    b <- qy[5L] + 3. * delta
    yab <- y[y >= a & y <= b]

    k <- kmax
    ok2 <- FALSE

    while (!ok2) {

      #hy <- hist(yab, breaks = k, plot = FALSE);b <- hy$breaks;n <- c(hy$counts, 0)
      b <- seq(from = min(yab), to = max(yab), length = k + 1L)
      n <- c(tabulate(findInterval(yab, b[-(k + 1L)])), 0L)

      N <- sum(n)
      v <- as.numeric(n >= N / k)

      ## Beginning of the first chain
      w <- which.max(v)
      v2 <- v[w:(k + 1L)]

      ## End of the first chain
      w2 <- which.min(v2) + w - 1L
      v3 <- v[w2:(k + 1L)]

      ## Length of the first chain
      nc <- sum(n[w:(w2 - 1L)])

      ## There exists another chain, and the first chain has only one element
      if (any(v3 == 1L) && nc == 1L) {
        if (k > 3L) {
          k <- k - 1L
        } else if (k == 3L) {
          if (n[3L] > 1L) {
            w <- 3L
            w2 <- 4L
          }
          ok2 <- TRUE
        } else {
          stop("k < 3", call. = FALSE)
        }

      ## There exists another chain, and the first chain has more than one element
      } else if (any(v3 == 1L) && nc > 1L) {
        if (k > 3L) {
          k <- k - 1L

        ### In this case, w = 1 necessarily
        } else if (k == 3L) {
          if (n[3L] > 1L) {
            p1 <- (1L / n[1L]) * prod(diff(yab[yab >= b[w] & yab <= b[w2]])) # here, n[1] = length(first chain)
            p2 <- (1L / n[3L]) * prod(diff(yab[yab >= b[3L] & yab <= b[4L]])) # and n[3] = length(second chain)
            if (p1 > p2) {
              w <- 3L
              w2 <- 4L
            }
          }
          ok2 <- TRUE
        } else {
          stop("k < 3", call. = FALSE)
        }

      ## There is no other chain: the modal chain is found!
      } else if (!any(v3 == 1L)) {
        ok2 <- TRUE
      }

    }

    ## Update 'nc'
    nc <- sum(n[w:(w2 - 1L)])
    #cat("Modal chain length = ", nc, "\n")

    d <- abs((qy[4L] + qy[2L] - 2L * qy[3L]) / (qy[4L] - qy[2L]))
    nc2 <-  ny * (1L - d)
    #cat("d = ", d, "\n")

    y <- yab[yab >= b[w] & yab <= b[w2]]
    if (nc == ny) {
      ok1 <- TRUE
    } else if (nc <= ifelse(nx < 30L, nx / 3L, bw * nc2)) {
      ok1 <- TRUE
    } else {
      ok1 <- FALSE
    }

  }
  stats::median(y)
}
