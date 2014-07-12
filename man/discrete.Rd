\name{mfv}
%\alias{mlv.discrete}
\alias{discrete}
%\alias{mlv.mfv}
\alias{mfv}

\title{ Estimate of the Mode of a Discrete Distribution (Most Frequent Value) }
\description{
  This function returns the most frequent value(s) in a given numerical vector.
}
\usage{
  mfv(x, ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{...}{ further arguments, which will be ignored. }
}
\details{
Argument \code{x} is to come from a discrete distribution. 
This function uses function \code{\link{tabulate}} of \R.
}
\value{
The most frequent value(s) found in \code{x} is (are) returned.
}
\note{
The user should preferentially call \code{mfv} through 
 \code{mlv(x, method = "mfv")} (or  \code{mlv(x, method = "discrete")}). 
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{ 
    \item Dutta S. and Goswami A. (2010).
    Mode estimation for discrete distributions.
    \emph{Mathematical Methods of Statistics}, \bold{19}(4):374--384.
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}} for general mode estimation;
 \code{\link{geomMode}}, \code{\link{poisMode}}, etc. for computation of the mode of the usual discrete distributions}
\examples{
# Unimodal distribution
x <- rbinom(100, size = 10, prob = 0.8)

## True mode
binomMode(size = 10, prob = 0.8)

## Most frequent value
mfv(x)
mlv(x, method = "discrete")

# Bimodal distribution
x <- rpois(100, lambda = 7)

## True mode
poisMode(lambda = 7)

## Most frequent value
mfv(x)
M <- mlv(x, method = "discrete")
print(M)
plot(M)
}
\keyword{ univar }
