#' @title 
#' Skewness
#' 
#' @description 
#' This function encodes different methods 
#' to calculate the skewness from a vector of 
#' observations. 
#' 
#' @references 
#' \itemize{
#'   \item Bickel D.R. (2002).
#'   Robust estimators of the mode and skewness of continuous data.
#'   \emph{Computational Statistics and Data Analysis}, \bold{39}:153-163.
#'   
#'   \item Bickel D.R. et Fruehwirth R. (2006).
#'   On a Fast, Robust Estimator of the Mode: Comparisons to Other Robust Estimators with Applications.
#'   \emph{Computational Statistics and Data Analysis}, \bold{50}(12):3500-3530.
#' }
#' 
#' @param x
#' numeric. Vector of observations.
#' 
#' @param na.rm
#' logical. Should missing values be removed? 
#' 
#' @param method
#' character. Specifies the method of computation. 
#' These are either \code{"moment"}, \code{"fisher"} or \code{"bickel"}. 
#' The \code{"moment"} method is based on the definition of 
#' skewness for distributions; this form should 
#' be used when resampling (bootstrap or jackknife). The 
#' \code{"fisher"} method corresponds to the usual "unbiased" 
#' definition of sample variance, although in the case of skewness 
#' exact unbiasedness is not possible. 
#' 
#' @param M
#' numeric. (An estimate of) the mode of the observations \code{x}. 
#' Default value is \code{\link[modeest]{shorth}(x)}.
#' 
#' @param ...
#' Additional arguments. 
#' 
#' @return 
#' \code{skewness} returns a numeric value. 
#' An attribute reports the method used.
#' 
#' @export
#' @importFrom stats sd
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for general mode estimation; 
#' \code{\link[modeest]{shorth}} for the shorth estimate of the mode
#' 
#' @author 
#' Diethelm Wuertz and contributors for the original \code{skewness} function 
#' from package \pkg{fBasics}. 
#' 
#' @examples 
#' ## Skewness = 0
#' x <- rnorm(1000)
#' skewness(x, method = "bickel", M = shorth(x))
#' 
#' ## Skewness > 0 (left skewed case)
#' x <- rbeta(1000, 2, 5)
#' skewness(x, method = "bickel", M = betaMode(2, 5))
#' 
#' ## Skewness < 0 (right skewed case)
#' x <- rbeta(1000, 7, 2)
#' skewness(x, method = "bickel", M = hsm(x, bw = 1/3))
#' 
skewness <-
function(x,
         na.rm = FALSE,
         method = c("moment", "fisher", "bickel"),
         M,
         ...)
{
  method <- match.arg(tolower(method), c("moment", "fisher", "bickel"))
  if (!is.numeric(x)) {
    stop("argument 'x' is must be numeric", call. = FALSE)
  }
  if (na.rm) x <- x[!is.na(x)]
  nx <- length(x)
  if (missing(M)) M <- shorth(x)
  if (method == "moment") {
    sk <- sum((x - mean(x))^3/stats::sd(x)^3)/nx
  }
  if (method == "fisher") {
    if (nx < 3) sk <- NA
    else sk <- ((sqrt(nx * (nx - 1))/(nx - 2)) * (sum(x^3)/nx))/((sum(x^2)/nx)^(3/2))
  }
  if (method == "bickel") {
    cdf.value <- (length(x[x < M]) + length(x[x == M])/2)/nx
    sk <- 1-2*cdf.value
  }
  sk
}
