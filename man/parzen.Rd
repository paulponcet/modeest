\name{parzen}
%\alias{mlv.kernel}
\alias{parzen}
\alias{Parzen}

\title{ Parzen's Kernel Mode Estimator }
\description{
Parzen's kernel mode estimator is the value maximizing the kernel density estimate.
}
\usage{
  parzen(x, 
         bw = NULL, 
         kernel = "gaussian", 
         abc = FALSE, 
         par = shorth(x), 
         optim.method = "BFGS", 
         ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The smoothing bandwidth to be used. }
  \item{kernel}{ character. The kernel to be used. Available kernels are 
    \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"}, 
    \code{"gaussian"}, \code{"optcosine"}, \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}. 
    See \code{\link{density.default}} for more details on some of these kernels. }
  \item{abc}{ logical. If \code{FALSE} (the default), the kernel density estimate  
  is maximised using \code{\link{optim}}. }
  \item{par}{ numeric. The initial value used in \code{\link{optim}}. }
  \item{optim.method}{ character. If \code{abc = FALSE}, the method used in \code{\link{optim}}. }
  \item{\dots}{ if \code{abc = FALSE}, further arguments to be passed to \code{\link{optim}}. }%, or further arguments to be passed to \code{\link{plot.default}}. }
}
\details{
If \code{kernel = "uniform"}, the \code{\link{naive}} mode estimate is returned.
}
\value{
\code{parzen} returns a numeric value, the mode estimate. If \code{abc = TRUE}, 
the \code{x} value maximizing the density estimate is returned. Otherwise, 
the \code{\link{optim}} method is used to perform maximization, and the  
attributes: 'value', 'counts', 'convergence' and 'message', coming from 
the \code{\link{optim}} method, are added to the result.
}
\note{
The user should preferentially call \code{parzen} through 
 \code{mlv(x, method = "kernel", ...)} 
 or \code{mlv(x, method = "parzen", ...)}.
This returns an object of class \code{\link{mlv}}.

Presently, \code{parzen} is quite slow.
}
\references{
  \itemize{ 
    \item Parzen E. (1962).
    On estimation of a probability density function and mode.
    \emph{Ann. Math. Stat.}, \bold{33}(3):1065--1076.
  
    \item Konakov V.D. (1973).
    On the asymptotic normality of the mode of multidimensional distributions.
    \emph{Theory Probab. Appl.}, \bold{18}:794-803.

    \item Eddy W.F. (1980).
    Optimum kernel estimators of the mode.
    \emph{Ann. Statist.}, \bold{8}(4):870-882.
    
    \item Eddy W.F. (1982).
    The Asymptotic Distributions of Kernel Estimators of the Mode.
    \emph{Z. Wahrsch. Verw. Gebiete}, \bold{59}:279-290.

    \item Romano J.P. (1988).
    On weak convergence and optimality of kernel density estimates of the mode.
    \emph{Ann. Statist.}, \bold{16}(2):629-647.
  
    \item Abraham C., Biau G. and Cadre B. (2003).
    Simple Estimation of the Mode of a Multivariate Density.
    \emph{Canad. J. Statist.}, \bold{31}(1):23-34.

    \item Abraham C., Biau G. and Cadre B. (2004).
    On the Asymptotic Properties of a Simple Estimate of the Mode.
    \emph{ESAIM Probab. Stat.}, \bold{8}:1-11.
  
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}}, \code{\link{naive}} }
\examples{
# Unimodal distribution
x <- rlnorm(10000, meanlog = 3.4, sdlog = 0.2)

## True mode
lnormMode(meanlog = 3.4, sdlog = 0.2)

## Estimate of the mode
M <- mlv(x, method = "kernel", kernel = "gaussian", bw = 0.3, par = shorth(x))
print(M)
plot(M)
}
\keyword{ univar }
