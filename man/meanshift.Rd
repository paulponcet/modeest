\name{meanshift}
\alias{meanshift}
\alias{Meanshift}

\title{ Meanshift Mode Estimator }
\description{
Meanshift mode estimator
}
\usage{
  meanshift(x, 
            bw = NULL, 
            kernel = "gaussian", 
            par = shorth(x), 
            iter = 1000, 
            eps = 1e-08)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The smoothing bandwidth to be used. }
  \item{kernel}{ character. The kernel to be used. Available kernels are 
    \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"}, 
    \code{"gaussian"}, \code{"optcosine"}, \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}. 
    See \code{\link{density.default}} for more details on some of these kernels. }
  \item{par}{ numeric. The initial value used in the meanshift algorithm. }
  \item{iter}{ numeric. Maximal number of iterations. }
  \item{eps}{ numeric. Stopping criteria. } 
}
\value{
\code{meanshift} returns a numeric value, the mode estimate, 
with an attribute \code{"iterations"}. 
The number of iterations can be less than \code{iter} 
if the stopping criteria specified by \code{eps} is reached. 
}
\note{
The user should preferentially call \code{meanshift} through 
 \code{mlv(x, method = "meanshift", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{ 
    \item Fukunaga, K. and Hostetler, L. (1975). 
    The estimation of the gradient of a density function, with applications in pattern recognition. 
    \emph{IEEE Transactions on Information Theory}, \bold{21}(1):32--40. 
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}}, \code{\link{tsybakov}} }
\examples{
# Unimodal distribution
x <- rweibull(100, shape = 12, scale = 0.8)

## True mode
weibullMode(shape = 12, scale = 0.8)

## Estimate of the mode
M <- mlv(x, method = "meanshift", kernel = "gaussian", bw = bw.SJ(x), par = shorth(x))
print(M)
plot(M)
}
\keyword{ univar }
