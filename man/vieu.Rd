\name{vieu}
%\alias{mlv.vieu}
\alias{vieu}
\alias{Vieu}

\title{ Vieu's Mode Estimator }
\description{
Vieu's mode estimator is the value at which the kernel density derivative estimate is null.
}
\usage{
  vieu(x, 
       bw = NULL, 
       kernel = "gaussian", 
       abc = FALSE, 
       ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The smoothing bandwidth to be used. }
  \item{kernel}{ character. The kernel to be used. Available kernels are 
    \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"}, 
    \code{"gaussian"}, \code{"optcosine"}, \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}. 
    See \code{\link{density.default}} for more details on some of these kernels. }
  \item{abc}{ logical. If \code{FALSE} (the default), the root of the density derivate estimate  
  is searched with \code{\link{uniroot}}. }
  \item{\dots}{ if \code{abc = FALSE}, further arguments to be passed to \code{\link{uniroot}}. }
}
\value{
\code{vieu} returns a numeric value, the mode estimate. If \code{abc = TRUE}, 
the \code{x} value at which the density derivative estimate is null is returned. Otherwise, 
the \code{\link{uniroot}} method is used.
}
\note{
The user should preferentially call \code{vieu} through 
 \code{mlv(x, method = "vieu", ...)}.
This returns an object of class \code{\link{mlv}}.

Presently, \code{vieu} is quite slow.
}
\references{
  \itemize{ 
    \item Vieu P. (1996).
    A note on density mode estimation.
    \emph{Statistics \& Probability Letters}, \bold{26}:297--307.
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}}, \code{\link{parzen}} }
\examples{
# Unimodal distribution
x <- rlnorm(10000, meanlog = 3.4, sdlog = 0.2)

## True mode
lnormMode(meanlog = 3.4, sdlog = 0.2)

## Estimate of the mode
M <- mlv(x, method = "vieu", kernel = "gaussian")
print(M)
plot(M)
}
\keyword{ univar }
