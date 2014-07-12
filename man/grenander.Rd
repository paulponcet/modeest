\name{grenander}
%\alias{mlv.grenander}
\alias{grenander}
\alias{Grenander}

\title{ The Grenander Mode Estimator }
\description{
This function computes the Grenander mode estimator.
}
\usage{
  grenander(x, 
            bw = NULL, 
            k, 
            p, 
            ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The bandwidth to be used. Should belong to (0, 1].}
  \item{k}{ numeric. Paramater 'k' in Grenander's mode estimate. }
  \item{p}{ numeric. Paramater 'p' in Grenander's mode estimate. 
            If \code{p = Inf}, function \code{\link{venter}} is used. }
  \item{...}{ further arguments to be passed to \code{link{venter}} }
}
\details{
The Grenander estimate is defined by
  \deqn{ \frac{ \sum_{j=1}^{n-k} \frac{(x_{j+k} + x_{j})}{2(x_{j+k} - x_{j})^p} }{ \sum_{j=1}^{n-k} \frac{1}{(x_{j+k} - x_{j})^p} } }{ ( sum_{j=1}^{n-k} (x_{j+k} + x_{j})/(2(x_{j+k} - x_{j})^p) ) / ( sum_{j=1}^{n-k} 1/((x_{j+k} - x_{j})^p) ) }

If \eqn{p}{p} tends to infinity, this estimate tends to the Venter mode estimate; 
this justifies to call \code{\link{venter}} if \code{p = Inf}.

The user should either give the bandwidth \code{bw} or the argument \code{k}, 
\code{k} being taken equal to \code{ceiling(bw*n) - 1} if missing.
}
\value{
A numeric value is returned, the mode estimate.
If \code{p = Inf}, the \code{\link{Venter}} mode estimator is returned.
}
\note{
The user should preferentially call \code{grenander} through 
\code{mlv(x, method = "grenander", bw, k, p, ...)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{
    \item Grenander U. (1965).
    Some direct estimates of the mode.
    \emph{Ann. Math. Statist.}, \bold{36}:131-138.
    
    \item Dalenius T. (1965).
    The Mode - A Negleted Statistical Parameter.
    \emph{J. Royal Statist. Soc. A}, \emph{128}:110-117.
    
    \item Adriano K.N., Gentle J.E. and Sposito V.A. (1977).
    On the asymptotic bias of Grenander's mode estimator.
    \emph{Commun. Statist.-Theor. Meth. A}, \bold{6}:773-776.
    
    \item Hall P. (1982).
    Asymptotic Theory of Grenander's Mode Estimator.
    \emph{Z. Wahrsch. Verw. Gebiete}, \bold{60}:315-334.
  }
}
\author{ D.R. Bickel for the original code, \cr
P. Poncet for the slight modifications introduced.
}
\seealso{ \code{\link{mlv}} for general mode estimation;
\code{\link{venter}} for the Venter mode estimate }
\examples{
# Unimodal distribution
x <- rnorm(1000, mean = 23, sd = 0.5)

## True mode
normMode(mean = 23, sd = 0.5) # (!)

## Parameter 'k'
k <- 5

## Many values of parameter 'p'
p <- seq(0.1, 4, 0.01)

## Estimate of the mode with these parameters
M <- sapply(p, function(pp) grenander(x, p = pp, k = k))

## Distribution obtained
plot(density(M), xlim = c(22.5, 23.5))
}
\keyword{ univar }
\keyword{ nonparametric }
