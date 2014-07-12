\name{hsm}
%\alias{mlv.hsm}
\alias{hsm}
\alias{HSM}
%\alias{fsm}
%\alias{robertson}
%\alias{cryer}
%Robertson & Cryer ??  OUI (a ceci pres que Robertson et Cryer proposent d'instaurer
%un arret de la procedure iterative apres un certain nombre d'iterations)

\title{ Half Sample Mode }
\description{
This function computes the Robertson-Cryer mode estimator described in Robertson and Cryer (1974), 
also called half sample mode (if \code{bw = 1/2}) or fraction sample mode (for some other \code{bw}) by Bickel (2006). 
}
\usage{
  hsm(x, 
      bw = NULL, 
      k, 
      tie.action = "mean", 
      tie.limit = 0.05, 
      ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric or function. The bandwidth to be used. Should belong to (0, 1]. }
  \item{k}{ numeric. See 'Details'. }
  \item{tie.action}{ character. The action to take if a tie is encountered. }
  \item{tie.limit}{ numeric. A limit deciding whether or not a warning is given when a tie is encountered. }
  \item{...}{ further arguments. }
}
\details{
The modal interval, i.e. the shortest interval among intervals containing \code{k+1} observations, 
is computed iteratively, until only one value is found, the mode estimate. 
At each step \eqn{i}{i}, one takes \code{k = ceiling(bw*n) - 1}, where \code{n} is the length of the modal interval computed at step \eqn{i-}{i-}\code{1}.
If \code{bw} is of class \code{"function"}, then \code{k = ceiling(bw(n)) - 1} instead.
}
\value{
A numeric value is returned, the mode estimate.
}
\note{
The user should preferentially call \code{hsm} through 
 \code{mlv(x, method = "hsm", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{ 
    \item Robertson T. and Cryer J.D. (1974).
    An iterative procedure for estimating the mode.
    \emph{J. Amer. Statist. Assoc.}, \bold{69}(348):1012-1016.
    
    \item Bickel D.R. and Fruehwirth R. (2006).
    On a Fast, Robust Estimator of the Mode: Comparisons to Other Robust Estimators with Applications.
    \emph{Computational Statistics and Data Analysis}, \bold{50}(12):3500-3530.
  }
}
\author{ D.R. Bickel for the original code, \cr
P. Poncet for the slight modifications introduced.
}
\seealso{ \code{\link{mlv}} for general mode estimation;
\code{\link{venter}} for the Venter mode estimate}
\examples{
# Unimodal distribution
x <- rweibull(10000, shape = 3, scale = 0.9)

## True mode
weibullMode(shape = 3, scale = 0.9)

## Estimate of the mode
bandwidth <- function(n, alpha) {1/n^alpha}
hsm(x, bw = bandwidth, alpha = 2)
M <- mlv(x, method = "hsm", bw = bandwidth, alpha = 2)
print(M)
plot(M)
}
\keyword{ univar }
\keyword{ robust }
\keyword{ nonparametric }
