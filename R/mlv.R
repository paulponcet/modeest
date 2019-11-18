
# TODO: packages 'ks' et 'kedd', notamment pour les kernel density derivative estimates

#' @title 
#' Estimation of the Mode(s) or Most Likely Value(s)
#' 
#' @description 
#' \code{mlv} is a generic function for estimating the mode of a univariate distribution. 
#' Different estimates (or methods) are provided: 
#' \itemize{
#'   \item \code{\link{mfv}}, which returns the most frequent value(s) in a given numerical vector, 
#'   \item the \code{\link{Lientz}} mode estimator, which is the value minimizing the Lientz function estimate, 
#'   \item the Chernoff mode estimator, also called \code{\link{naive}} mode estimator, 
#'   which is defined as the center of the interval of given length containing the most observations, 
#'   \item the \code{\link{Venter}} mode estimator, including the \code{\link{shorth}}, i.e. the midpoint of the modal interval, 
#'   \item the \code{\link{Grenander}} mode estimator, 
#'   \item the half sample mode (\code{\link{HSM}}) and the half range mode (\code{\link{HRM}}), which are iterative versions of the Venter mode estimator, 
#'   \item \code{\link{Parzen}}'s kernel mode estimator, which is the value maximizing the kernel density estimate, 
#'   \item the \code{\link{Tsybakov}} mode estimator, based on a gradient-like recursive algorithm, 
#'   \item the \code{\link{Asselin}} de Beauville mode estimator, based on a algorithm detecting chains and holes in the sample, 
#'   \item the \code{\link{Vieu}} mode estimator, 
#'   \item the \code{\link{meanshift}} mode estimator. 
#' }
#' 
#' \code{mlv} can also be used to compute the mode of a given distribution, with \code{mlv.character}.
#' 
#' @details 
#' For the default method of \code{mlv}, available methods are \code{"lientz"}, 
#' \code{"naive"}, \code{"venter"}, 
#' \code{"grenander"}, \code{"hsm"}, \code{"parzen"}, 
#' \code{"tsybakov"}, \code{"asselin"}, and \code{"meanshift"}. 
#' See the description above and the associated links. 
#' 
#' If \code{x} is of class \code{"character"} (with length > 1), 
#' \code{"factor"}, or \code{"integer"}, then the most frequent value found in 
#' \code{x} is returned using \code{\link[statip]{mfv}} from package 
#' \pkg{statip}. 
#' 
#' If \code{x} is of class \code{"character"} (with length 1), 
#' \code{x} should be one of \code{"beta"}, \code{"cauchy"}, \code{"gev"}, etc. 
#' i.e. a character for which a function \code{*Mode} exists 
#' (for instance \code{betaMode}, \code{cauchyMode}, etc.). 
#' See \code{\link[modeest]{distrMode}} for the available functions. 
#' The mode of the corresponding distribution is returned. 
#' 
# #' If \code{x} is of class \code{"density"}, the value where the density is 
# #' maximised is returned. 
# #' 
#' If \code{x} is of class \code{mlv.lientz}, see \code{\link[modeest]{Lientz}} 
#' for more details. 
#' 
#' @references 
#' See the references on mode estimation on the \code{\link[modeest]{modeest-package}}'s page. 
#' 
#' @param x
#' numeric (vector of observations), or an object of class \code{"factor"}, \code{"integer"}, etc. 
#' 
#' @param bw
#' numeric. The bandwidth to be used. 
#' This may have different meanings regarding the \code{method} used. 
#' 
#' @param method
#' character. One of the methods available for computing the mode estimate. See 'Details'. 
#' 
#' @param na.rm
#' logical. Should missing values be removed? 
#' 
# #' @param all
# #' logical. 
# #' 
# #' @param abc
# #' logical. If \code{FALSE} (the default), the estimate of the density function 
# #' is maximised using \code{\link{optim}}. 
# #' 
#' @param ...
#' Further arguments to be passed to the function called for computation. 
#' 
#' @return 
#' A vector of the same type as \code{x}. 
#' Be aware that the length of this vector can be \code{> 1}. 
#' 
#' @seealso 
#' \code{\link[statip]{mfv}}, 
#' \code{\link[modeest]{parzen}}, 
#' \code{\link[modeest]{venter}}, 
#' \code{\link[modeest]{meanshift}},
#' \code{\link[modeest]{grenander}}, 
# #' \code{\link[modeest]{hrm}}, 
#' \code{\link[modeest]{hsm}}, 
#' \code{\link[modeest]{lientz}}, 
#' \code{\link[modeest]{naive}}, 
#' \code{\link[modeest]{tsybakov}}, 
#' \code{\link[modeest]{skewness}}
#' 
#' @export
#' 
#' @examples 
#' # Unimodal distribution
#' x <- rbeta(1000,23,4)
#' 
#' ## True mode
#' betaMode(23, 4)
#' # or
#' mlv("beta", shape1 = 23, shape2 = 4)
#' 
#' ## Be aware of this behaviour: 
#' mlv("norm") # returns 0, the mode of the standard normal distribution
#' mlv("normal") # returns 0 again, since "normal" is matched with "norm"
#' mlv("abnormal") # returns "abnormal", since the input vector "abrnormal" 
#' # is not recognized as a distribution name, hence is taken as a character 
#' # vector from which the most frequent value is requested. 
#' 
#' ## Estimate of the mode
#' mlv(x, method = "lientz", bw = 0.2)
#' mlv(x, method = "naive", bw = 1/3)
#' mlv(x, method = "venter", type = "shorth")
#' mlv(x, method = "grenander", p = 4)
# #' mlv(x, method = "hrm", bw = 0.3)
#' mlv(x, method = "hsm")
#' mlv(x, method = "parzen", kernel = "gaussian")
#' mlv(x, method = "tsybakov", kernel = "gaussian")
#' mlv(x, method = "asselin", bw = 2/3)
#' mlv(x, method = "vieu")
#' mlv(x, method = "meanshift")
#' 
mlv <-
function(x,
         ...)
{
  UseMethod("mlv")
}


#' @importFrom statip name2distr
#' @export
#' @rdname mlv
#' 
mlv.character <- 
function(x, 
         na.rm = FALSE,
         ...)
{
  stopifnot(is.character(x))
  if (length(x)==1L && statip::name2distr(x) %in% .distributionsList()) {
    distrMode(x, ...)
  } else {
    mfv(x, na_rm = na.rm)
  }
}


#' @export
#' @rdname mlv
#' 
mlv.factor <-
function(x,
         na.rm = FALSE,
         ...)
{
  stopifnot(is.factor(x))
  mfv(x, na_rm = na.rm)
}


#' @export
#' @rdname mlv
#' 
mlv.logical <-
function(x,
         na.rm = FALSE,
         ...)
{
  stopifnot(is.logical(x))
  mfv(x, na_rm = na.rm)
}


#' @export
#' @rdname mlv
#' 
mlv.integer <-
function(x,
         na.rm = FALSE,
         ...)
{
  stopifnot(is.integer(x))
  mfv(x, na_rm = na.rm)
}


# #' @importFrom bazar is_wholenumber
#' @export
#' @rdname mlv
#' 
mlv.default <-
function(x,
         bw = NULL,
         method,
         na.rm = FALSE,
         ...)
{

  stopifnot(is.numeric(x))
  x <- as.vector(x)

  #test <- bazar::is_wholenumber(x)
  #if (is.na(test) || test) {
  #  return(mfv(as.integer(round(x)), na_rm = na.rm))
  #}
  
  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm) {
      x <- x[!x.na]
    } else {
      stop("argument 'x' contains missing values", 
           call. = FALSE)
    }
  }

  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
  }

  if (missing(method)) {
    warning("argument 'method' is missing. Data are supposed to be continuous. 
            Default method 'shorth' is used", 
            call. = FALSE)
    method <- "shorth"
  #} else if (tolower(method) == "mfv") {
  #  stop("incorrect 'method' argument")
  } else if (pmatch(tolower(method), c("density", "kernel"), nomatch = 0)) {
    method <- "parzen"
  } else method <- match.arg(tolower(method), .methodsList())
  
  if (method == "lientz") method <- "mlv.lientz"
  
  do.call(method, list(x = x, bw = bw, ...)) # possibly length > 1
  #mean(theta)
}


# #' @export
# #' @rdname mlv
# #' 
# mlv.density <-
# function(x,
#          all = TRUE, 
#          abc = FALSE,
#          ...)
# {
#   # TODO: A MODIFIER EN MEME TEMPS QUE 'parzen'
#   
#   if (!inherits(x, "density")) stop("argument 'x' must inherit from class 'density'")
#   
#   y <- x$y
#   x <- x$x
# 
#   den.s <- stats::smooth.spline(x, y, all.knots=TRUE, spar=spar)
#   s.1 <- stats::predict(den.s, den.s$x, deriv = 1)
#   s.0 <- stats::predict(den.s, den.s$x, deriv = 0)
#   
#   den.sign <- sign(s.1$y)
#   b <- rle(den.sign)$values
#   nmodes <- length(b)/2
#   #if (nmodes > 10) { nmodes <- 10 }
#   if (is.na(nmodes)) { nmodes <- 0 }
#   
#   a <- c(1,1+which(diff(den.sign)!=0))
#   df <- data.frame(a,b)
#   df <- df[which(df$b %in% -1),]
#   modes <- s.1$x[df$a]
#   density <- s.0$y[df$a]
#   df2 <- data.frame(modes,density)
#   df2 <- df2[with(df2, order(-density)), ] # ordered by density
#   df2
#   
#   
# #-------------------  
#   
#   idx <- y == max(y)
#   M <- x[idx]
# 
#   if (all) {
#     yy <- c(0, y, 0)
#     ny <- length(yy)
#     idx <- (yy[2:(ny - 1)] > yy[1:(ny - 2)]) & (yy[2:(ny - 1)] > yy[3:ny])
#     M <- unique(c(x[idx], M))
#   }
#    
#   M
# }


#' @export
#' @rdname mlv
#'
mlv1 <-
function(x, 
         ...)
{
  mlv(x, ...)[[1L]]
}
